#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MOS Equation Lab v3 — tiny schematic editor, multi-device DC solver, equation panels,
λ (channel-length modulation), PMOS, DC current source, capacitor (DC-open), netlist export.

Design choices
- Graph is node-centric; devices connect by node names.
- You can add/move nodes and devices. Connectivity is edited via per-device dropdowns.
- DC solver assembles KCL at all non-fixed nodes, iterates over MOS region assignments,
  solves for node voltages with SymPy, and checks region consistency.
- Equation panels are draggable Toplevel windows.
- Netlist export targets NGSpice/LTSpice Level-1 models using per-device .model cards:
    Id_sat = 0.5*KP*(W/L)*(Vgs - Vto)^2*(1 + LAMBDA*Vds)
  We set W/L=1 and KP := k (your μCox·W/L lumped parameter), so equations align.

Limitations (for now)
- DC only (capacitors act open). Transient is planned.
- One global VDD source to GND (adjustable). Add’l sources ≈ easy to generalize later.
- Gate currents are zero; body tied to source.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import sympy as sp
from math import isfinite
import itertools
import uuid

# ------------- symbols & helpers ------------------------------------------------
def sym(name): return sp.Symbol(name, real=True)
Ground0 = sp.Integer(0)

def pretty(e): return sp.pretty(e)

# Utility for numeric testing with tolerance
def within(x, lo, hi, eps=1e-12): return (x >= lo - eps) and (x <= hi + eps)

# ------------- device primitives ------------------------------------------------

class Device:
    kind = "Device"
    def __init__(self, name):
        self.name = name
        self.x, self.y = 200, 200   # canvas location
        self.selected = False
    @property
    def pins(self):
        """Return ordered list of (pin_name, default_node_name) for UI."""
        return []
    def draw(self, C, nodes_pos):
        """Override to draw on canvas; nodes_pos gives dict name->(x,y)."""
        pass
    def kcl_contribs(self, node_syms, node_to_sym):
        """
        Return dict node_name -> (symbolic current leaving node).
        Positive means current leaving node; contributes to KCL sum(node)=0.
        """
        return {}

class Resistor(Device):
    kind = "R"
    def __init__(self, name, a="N1", b="N2", R=1000.0):
        super().__init__(name)
        self.a, self.b, self.R = a, b, float(R)
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Va = node_to_sym[self.a]; Vb = node_to_sym[self.b]
        Iab = (Va - Vb)/self.R  # current from a->b
        return { self.a:  Iab, self.b: -Iab }
    def draw(self, C, nodes_pos):
        ax, ay = nodes_pos[self.a]; bx, by = nodes_pos[self.b]
        C.create_line(ax, ay, bx, by, fill="#aaa", width=2)
        # simple text in middle
        C.create_text((ax+bx)/2, (ay+by)/2 - 14, text=f"{self.name}={self.R:g}Ω",
                      fill="#aaa", font=("Consolas", 10))

class CurrentSourceDC(Device):
    kind = "I"
    def __init__(self, name, a="N1", b="N2", I=1e-3):
        super().__init__(name)
        self.a, self.b, self.I = a, b, float(I)   # +I from a to b
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def kcl_contribs(self, node_syms, node_to_sym):
        return { self.a:  sp.Float(self.I), self.b: -sp.Float(self.I) }
    def draw(self, C, nodes_pos):
        ax, ay = nodes_pos[self.a]; bx, by = nodes_pos[self.b]
        C.create_line(ax, ay, bx, by, fill="#cc8", width=2)
        C.create_text((ax+bx)/2, (ay+by)/2 - 14, text=f"{self.name}={self.I:g}A",
                      fill="#cc8", font=("Consolas", 10))

class Capacitor(Device):
    kind = "C"
    def __init__(self, name, a="N1", b="N2", Cval=1e-9):
        super().__init__(name)
        self.a, self.b, self.C = a, b, float(Cval)
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def kcl_contribs(self, node_syms, node_to_sym):
        # DC open => no current; placeholder for transient later.
        return {}
    def draw(self, C, nodes_pos):
        ax, ay = nodes_pos[self.a]; bx, by = nodes_pos[self.b]
        C.create_line(ax, ay, (ax+bx)/2-6, (ay+by)/2, fill="#8ac", width=2)
        C.create_line((ax+bx)/2-6, (ay+by)/2-12, (ax+bx)/2-6, (ay+by)/2+12, fill="#8ac", width=2)
        C.create_line((ax+bx)/2+6, (ay+by)/2-12, (ax+bx)/2+6, (ay+by)/2+12, fill="#8ac", width=2)
        C.create_line((ax+bx)/2+6, (ay+by)/2, bx, by, fill="#8ac", width=2)
        C.create_text((ax+bx)/2, (ay+by)/2 - 18, text=f"{self.name}={self.C:g}F",
                      fill="#8ac", font=("Consolas", 10))

# MOS helper: region-aware NMOS
def nmos_Id(Vgs, Vds, k, Vt, lam, region):
    Id_lin = k*((Vgs - Vt)*Vds - sp.Rational(1,2)*Vds**2)
    Id_sat = sp.Rational(1,2)*k*(Vgs - Vt)**2
    if lam and lam != 0:
        Id_sat = Id_sat*(1 + lam*Vds)
    if region == "cutoff": return sp.Integer(0)
    if region == "linear": return Id_lin
    if region == "sat":    return Id_sat
    raise ValueError("bad region")

class NMOS(Device):
    kind = "NMOS"
    def __init__(self, name, d="N1", g="N2", s="GND", k=2e-4, Vt=0.6, lam=0.0):
        super().__init__(name)
        self.d, self.g, self.s = d, g, s
        self.k, self.Vt, self.lam = float(k), float(Vt), float(lam)
        self.region = None
    @property
    def pins(self): return [("d", self.d), ("g", self.g), ("s", self.s)]
    def Id_expr(self, node_to_sym, region):
        Vd = node_to_sym[self.d]; Vg = node_to_sym[self.g]; Vs = node_to_sym[self.s]
        Vgs = Vg - Vs; Vds = Vd - Vs
        return nmos_Id(Vgs, Vds, sp.Float(self.k), sp.Float(self.Vt), sp.Float(self.lam), region), Vgs, Vds
    def kcl_contribs(self, node_syms, node_to_sym):
        # region handled by solver; here return placeholders (0); solver calls Id_expr directly
        return {}
    def draw(self, C, nodes_pos):
        # safe pin lookup with fallbacks near the device's own (x,y)
        dx, dy = nodes_pos.get(self.d, (self.x, self.y - 30))
        gx, gy = nodes_pos.get(self.g, (self.x - 40, self.y))
        sx, sy = nodes_pos.get(self.s, (self.x, self.y + 30))
        cx = (dx + gx + sx)/3; cy = (dy + gy + sy)/3
        w, h = 90, 60
        C.create_rectangle(cx-w/2, cy-h/2, cx+w/2, cy+h/2, outline="#7fd", width=2)
        C.create_text(cx, cy-h/2-12, text=f"{self.name} (NMOS)", fill="#7fd", font=("Consolas", 10))
        C.create_line(dx, dy, cx, cy-h/2, fill="#aaa", width=2); C.create_oval(dx-3,dy-3,dx+3,dy+3,fill="#fff",outline="")
        C.create_text(dx+6, dy-10, text="D", fill="#aaa", font=("Consolas", 9), anchor="w")
        C.create_line(gx, gy, cx-w/2, cy, fill="#aaa", width=2); C.create_oval(gx-3,gy-3,gx+3,gy+3,fill="#fff",outline="")
        C.create_text(gx+6, gy-10, text="G", fill="#aaa", font=("Consolas", 9), anchor="w")
        C.create_line(sx, sy, cx, cy+h/2, fill="#aaa", width=2); C.create_oval(sx-3,sy-3,sx+3,sy+3,fill="#fff",outline="")
        C.create_text(sx+6, sy-10, text="S", fill="#aaa", font=("Consolas", 9), anchor="w")


# PMOS: current flows S->D when on. We keep k (|kp|) and Vt as positive magnitude.
def pmos_Id_mag(Vsg, Vsd, k, Vt, lam, region):
    Id_lin = k*((Vsg - Vt)*Vsd - sp.Rational(1,2)*Vsd**2)         # magnitude S->D
    Id_sat = sp.Rational(1,2)*k*(Vsg - Vt)**2
    if lam and lam != 0:
        Id_sat = Id_sat*(1 + lam*Vsd)
    if region == "cutoff": return sp.Integer(0)
    if region == "linear": return Id_lin
    if region == "sat":    return Id_sat
    raise ValueError("bad region")

class PMOS(Device):
    kind = "PMOS"
    def __init__(self, name, d="N1", g="N2", s="VDD", k=2e-4, Vt=0.6, lam=0.0):
        super().__init__(name)
        self.d, self.g, self.s = d, g, s
        self.k, self.Vt, self.lam = float(k), float(Vt), float(lam)
        self.region = None
    @property
    def pins(self): return [("d", self.d), ("g", self.g), ("s", self.s)]
    def Id_d_leaving(self, node_to_sym, region):
        Vd = node_to_sym[self.d]; Vg = node_to_sym[self.g]; Vs = node_to_sym[self.s]
        Vsg = Vs - Vg; Vsd = Vs - Vd
        Idmag = pmos_Id_mag(Vsg, Vsd, sp.Float(self.k), sp.Float(self.Vt), sp.Float(self.lam), region)
        # magnitude is S->D; current leaving D equals -Idmag
        return -Idmag, Vsg, Vsd
    def draw(self, C, nodes_pos):
        dx, dy = nodes_pos.get(self.d, (self.x, self.y - 30))
        gx, gy = nodes_pos.get(self.g, (self.x - 40, self.y))
        sx, sy = nodes_pos.get(self.s, (self.x, self.y + 30))
        cx = (dx + gx + sx)/3; cy = (dy + gy + sy)/3
        w, h = 90, 60
        C.create_rectangle(cx-w/2, cy-h/2, cx+w/2, cy+h/2, outline="#f8a", width=2)
        C.create_text(cx, cy-h/2-12, text=f"{self.name} (PMOS)", fill="#f8a", font=("Consolas", 10))
        C.create_line(dx, dy, cx, cy-h/2, fill="#aaa", width=2); C.create_oval(dx-3,dy-3,dx+3,dy+3,fill="#fff",outline="")
        C.create_text(dx+6, dy-10, text="D", fill="#aaa", font=("Consolas", 9), anchor="w")
        C.create_line(gx, gy, cx-w/2, cy, fill="#aaa", width=2); C.create_oval(gx-3,gy-3,gx+3,gy+3,fill="#fff",outline="")
        C.create_text(gx+6, gy-10, text="G", fill="#aaa", font=("Consolas", 9), anchor="w")
        C.create_line(sx, sy, cx, cy+h/2, fill="#aaa", width=2); C.create_oval(sx-3,sy-3,sx+3,sy+3,fill="#fff",outline="")
        C.create_text(sx+6, sy-10, text="S", fill="#aaa", font=("Consolas", 9), anchor="w")


# ------------- circuit model ---------------------------------------------------

class Circuit:
    def __init__(self):
        # nodes: dict name -> (x,y) position + fixed flag/value
        self.nodes = {
            "VDD": {"x": 140, "y": 80,  "fixed": True,  "value": 3.3},
            "GND": {"x": 140, "y": 520, "fixed": True,  "value": 0.0},
            "X":   {"x": 360, "y": 260, "fixed": False},
            "Y":   {"x": 560, "y": 260, "fixed": False},
        }
        self.devices = []
        # seed with current mirror + Rb + Rload
        self.devices.append(Resistor("Rb", a="VDD", b="X", R=10000.0))
        self.devices.append(Resistor("Rload", a="VDD", b="Y", R=20000.0))
        self.devices.append(NMOS("M1", d="X", g="X", s="GND", k=2e-4, Vt=0.6, lam=0.0))
        self.devices.append(NMOS("M2", d="Y", g="X", s="GND", k=2e-4, Vt=0.6, lam=0.0))
        # transient state
        self.t = 0.0
        self.prev_node_voltages = {}  # dict name->float from last step (free nodes only)

    def _snapshot_voltages(self, solved_dict):
        """Save last-step voltages for free nodes."""
        self.prev_node_voltages = {n: solved_dict.get(n, 0.0)
                                for n, inf in self.nodes.items()
                                if not inf.get("fixed", False)}

    # ----- net utilities
    def node_syms(self):
        """Map nodes to SymPy symbols/constants."""
        mapping = {}
        for name, info in self.nodes.items():
            if info.get("fixed", False):
                mapping[name] = sp.Float(info.get("value", 0.0))
            else:
                mapping[name] = sym(f"V_{name}")
        return mapping

    # inside class Circuit:
    def ensure_node(self, name, x=220, y=200):
        """Create a free (unfixed) node if it doesn't exist yet."""
        if name not in self.nodes:
            self.nodes[name] = {"x": x, "y": y, "fixed": False}
    def _ensure_all_device_nodes(self):
        for d in self.devices:
            if hasattr(d, "pins"):
                for _, node in d.pins:
                    if node not in self.nodes and node not in ("GND", "VDD"):
                        self.ensure_node(node)

    def solve_transient_step(self, dt):
        """
        One Backward-Euler step:
        KCL(node) + sum_C{ I_C,leaving(node) } = 0
        I_C(a->b) = C * ((Va - Vb) - (Va_prev - Vb_prev)) / dt
        Region search for MOS as in DC, but including the capacitor currents.
        """
        self._ensure_all_device_nodes()
        node_map = self.node_syms()
        unknown_nodes = [n for n, inf in self.nodes.items() if not inf.get("fixed", False)]
        unknown_syms = [node_map[n] for n in unknown_nodes]

        # prev voltages (fixed nodes use fixed values)
        prev = {}
        for n, inf in self.nodes.items():
            if inf.get("fixed", False):
                prev[n] = float(inf.get("value", 0.0))
            else:
                prev[n] = float(self.prev_node_voltages.get(n, 0.0))

        def empty_kcl():
            return {n: sp.Integer(0) for n in unknown_nodes}

        base_kcl = empty_kcl()

        # Resistors + DC sources as before
        for dev in self.devices:
            if isinstance(dev, (Resistor, CurrentSourceDC)):
                contrib = dev.kcl_contribs(self.nodes, node_map)
                for node, cur in contrib.items():
                    if node in base_kcl:
                        base_kcl[node] += cur

        # Capacitors (Backward Euler)
        for dev in self.devices:
            if isinstance(dev, Capacitor):
                a, b, Cval = dev.a, dev.b, dev.C
                Va, Vb = node_map[a], node_map[b]
                Iab = sp.Float(Cval) * ((Va - Vb) - (sp.Float(prev[a]) - sp.Float(prev[b]))) / sp.Float(dt)
                if a in base_kcl: base_kcl[a] += Iab       # leaving a
                if b in base_kcl: base_kcl[b] -= Iab       # entering a => leaving b negative

        mos_list = [d for d in self.devices if isinstance(d, (NMOS, PMOS))]
        region_opts = ["sat", "linear", "cutoff"] if mos_list else []

        if not mos_list:
            eqs = [sp.Eq(base_kcl[n], 0) for n in unknown_nodes]
            sol = sp.solve(eqs, unknown_syms, dict=True)
            if not sol: return None
            s = sol[0]
            num = {n: float(s[node_map[n]]) for n in unknown_nodes}
            self.t += dt
            self._snapshot_voltages(num)
            return num, {"kcl_sym": {n: pretty(sp.Eq(sp.simplify(base_kcl[n]), 0)) for n in unknown_nodes}}, eqs

        # Region search with caps included
        for choice in itertools.product(region_opts, repeat=len(mos_list)):
            kcl = {k: v for k, v in base_kcl.items()}
            preds = []
            for dev, rgn in zip(mos_list, choice):
                if isinstance(dev, NMOS):
                    Id, Vgs, Vds = dev.Id_expr(node_map, rgn)
                    if dev.d in kcl: kcl[dev.d] += Id
                    if dev.s in kcl: kcl[dev.s] -= Id
                    vt = dev.Vt
                    if rgn == "cutoff": preds.append(sp.Le(Vgs, vt))
                    elif rgn == "linear": preds += [sp.Ge(Vgs, vt), sp.Le(Vds, Vgs - vt)]
                    elif rgn == "sat":    preds += [sp.Ge(Vgs, vt), sp.Ge(Vds, Vgs - vt)]
                else:
                    Id_d, Vsg, Vsd = dev.Id_d_leaving(node_map, rgn)
                    if dev.d in kcl: kcl[dev.d] += Id_d
                    if dev.s in kcl: kcl[dev.s] -= Id_d
                    vt = dev.Vt
                    if rgn == "cutoff": preds.append(sp.Le(Vsg, vt))
                    elif rgn == "linear": preds += [sp.Ge(Vsg, vt), sp.Le(Vsd, Vsg - vt)]
                    elif rgn == "sat":    preds += [sp.Ge(Vsg, vt), sp.Ge(Vsd, Vsg - vt)]

            eqs = [sp.Eq(kcl[n], 0) for n in unknown_nodes]
            try:
                sol_list = sp.solve(eqs, unknown_syms, dict=True)
            except Exception:
                sol_list = []

            for sol in sol_list:
                try:
                    num = {n: float(sol[node_map[n]]) for n in unknown_nodes}
                except Exception:
                    continue
                if not all(isfinite(v) for v in num.values()):
                    continue

                # Check predicates (with fixed + numeric)
                ok = True
                subs = {node_map[n]: num[n] for n in unknown_nodes}
                for n, inf in self.nodes.items():
                    if inf.get("fixed", False):
                        subs[node_map.get(n, sym("dummy"))] = float(inf.get("value", 0.0))
                for p in preds:
                    try:
                        if not bool(p.subs(subs)): ok = False; break
                    except Exception:
                        ok = False; break
                if ok:
                    self.t += dt
                    self._snapshot_voltages(num)
                    return num, {"regions": dict(zip([m.name for m in mos_list], choice)),
                                "kcl_sym": {n: pretty(sp.Eq(sp.simplify(kcl[n]), 0)) for n in unknown_nodes}}, eqs
        return None, None, None
    def solve_dc(self):
        self._ensure_all_device_nodes()
        node_map = self.node_syms()
        unknown_nodes = [n for n, inf in self.nodes.items() if not inf.get("fixed", False)]
        unknown_syms = [node_map[n] for n in unknown_nodes]

        # KCL assembly templates (dict node->expr)
        def empty_kcl():
            return {n: sp.Integer(0) for n in unknown_nodes}

        # Build KCL excluding MOS for now (we add them per-region below)
        base_kcl = empty_kcl()
        for dev in self.devices:
            if isinstance(dev, (Resistor, CurrentSourceDC, Capacitor)):
                contrib = dev.kcl_contribs(self.nodes, node_map)
                for node, cur in contrib.items():
                    if node in base_kcl:
                        base_kcl[node] += cur

        # Collect MOS devices
        mos_list = [d for d in self.devices if isinstance(d, (NMOS, PMOS))]
        if not mos_list:
            # Solve linear system (R + I only)
            eqs = [sp.Eq(base_kcl[n], 0) for n in unknown_nodes]
            sol = sp.solve(eqs, unknown_syms, dict=True)
            if not sol: return None
            s = sol[0]
            return {n: float(s[node_map[n]]) for n in unknown_nodes}, {}, eqs

        # Region search
        region_opts = ["sat", "linear", "cutoff"]
        for choice in itertools.product(region_opts, repeat=len(mos_list)):
            kcl = {k: v for k, v in base_kcl.items()}
            preds = []
            # Add MOS currents per region
            for dev, rgn in zip(mos_list, choice):
                if isinstance(dev, NMOS):
                    Id, Vgs, Vds = dev.Id_expr(node_map, rgn)
                    # currents leaving nodes:
                    # drain leaving: +Id ; source leaving: -Id ; gate 0
                    if dev.d in kcl: kcl[dev.d] += Id
                    if dev.s in kcl: kcl[dev.s] -= Id
                    # region inequalities to check later
                    vt = dev.Vt
                    if rgn == "cutoff": preds.append(sp.Le(Vgs, vt))
                    elif rgn == "linear": preds += [sp.Ge(Vgs, vt), sp.Le(Vds, Vgs - vt)]
                    elif rgn == "sat":    preds += [sp.Ge(Vgs, vt), sp.Ge(Vds, Vgs - vt)]
                elif isinstance(dev, PMOS):
                    Id_d_leave, Vsg, Vsd = dev.Id_d_leaving(node_map, rgn)
                    # At drain: leaving current = Id_d_leave ; at source: opposite
                    if dev.d in kcl: kcl[dev.d] += Id_d_leave
                    if dev.s in kcl: kcl[dev.s] -= Id_d_leave
                    vt = dev.Vt
                    if rgn == "cutoff": preds.append(sp.Le(Vsg, vt))
                    elif rgn == "linear": preds += [sp.Ge(Vsg, vt), sp.Le(Vsd, Vsg - vt)]
                    elif rgn == "sat":    preds += [sp.Ge(Vsg, vt), sp.Ge(Vsd, Vsg - vt)]

            # Solve this region's KCL
            eqs = [sp.Eq(kcl[n], 0) for n in unknown_nodes]
            try:
                sol_list = sp.solve(eqs, unknown_syms, dict=True)
            except Exception:
                sol_list = []

            for sol in sol_list:
                try:
                    num = {n: float(sol[node_map[n]]) for n in unknown_nodes}
                except Exception:
                    continue
                if not all(isfinite(v) for v in num.values()):
                    continue
                # Check region predicates numerically
                ok = True
                subs = {node_map[n]: num[n] for n in unknown_nodes}
                # include fixed nodes:
                for n, inf in self.nodes.items():
                    if inf.get("fixed", False):
                        subs[node_map.get(n, sym("dummy"))] = float(inf.get("value", 0.0))
                for p in preds:
                    try:
                        if not bool(p.subs(subs)): ok = False; break
                    except Exception:
                        ok = False; break
                if ok:
                    # build symbolic KCL strings for panels
                    kcl_sym = {n: pretty(sp.Eq(sp.simplify(kcl[n]), 0)) for n in unknown_nodes}
                    kcl_num = {n: pretty(sp.Eq(sp.nsimplify(kcl[n].subs(subs), rational=False), 0)) for n in unknown_nodes}
                    # device numeric info
                    dnum = {}
                    for dev, rgn in zip(mos_list, choice):
                        if isinstance(dev, NMOS):
                            Id, Vgs, Vds = dev.Id_expr(node_map, rgn)
                            dnum[dev.name] = {
                                "region": rgn,
                                "Vgs": float(Vgs.subs(subs)),
                                "Vds": float(Vds.subs(subs)),
                                "Id_expr": Id,
                                "Id_num_str": pretty(sp.nsimplify(Id.subs(subs), rational=False))
                            }
                        else:
                            Id_d, Vsg, Vsd = dev.Id_d_leaving(node_map, rgn)
                            dnum[dev.name] = {
                                "region": rgn,
                                "Vsg": float(Vsg.subs(subs)),
                                "Vsd": float(Vsd.subs(subs)),
                                "Id_d_expr": Id_d,
                                "Id_d_num_str": pretty(sp.nsimplify(Id_d.subs(subs), rational=False))
                            }
                    return num, {"regions": dict(zip([m.name for m in mos_list], choice)),
                                 "kcl_sym": kcl_sym, "kcl_num": kcl_num, "dev": dnum}, eqs
        return None, None, None

    # ----- netlist export (NGSpice/LTSpice Level-1)
    def export_netlist(self, path):
        lines = []
        lines.append("* Auto-generated by MOS Equation Lab v3")
        vdd = self.nodes["VDD"]["value"]
        lines.append(f"VDD VDD 0 {vdd}")
        # Create per-device model cards for MOS so KP=k & W/L=1
        model_cards = []
        for d in self.devices:
            if isinstance(d, NMOS):
                mname = f"MOD_{d.name}"
                model_cards.append(f".model {mname} NMOS (VTO={d.Vt} KP={d.k} LAMBDA={d.lam})")
                lines.append(f"M{d.name} {d.d} {d.g} {d.s} 0 {mname} W=1u L=1u")
            elif isinstance(d, PMOS):
                mname = f"MOD_{d.name}"
                # PMOS VTO should be negative in SPICE; we pass -Vt
                model_cards.append(f".model {mname} PMOS (VTO={-d.Vt} KP={d.k} LAMBDA={d.lam})")
                lines.append(f"M{d.name} {d.d} {d.g} {d.s} 0 {mname} W=1u L=1u")
        for d in self.devices:
            if isinstance(d, Resistor):
                lines.append(f"R{d.name} {d.a} {d.b} {d.R}")
            elif isinstance(d, CurrentSourceDC):
                # I from a->b equals +I; SPICE current source Iname N+ N- value
                lines.append(f"I{d.name} {d.a} {d.b} {d.I}")
            elif isinstance(d, Capacitor):
                lines.append(f"C{d.name} {d.a} {d.b} {d.C}")
        lines.append(".op")
        lines += model_cards
        lines.append(".end")
        with open(path, "w") as f:
            f.write("\n".join(lines))

# ------------- UI (editor + panels) -------------------------------------------

class DraggablePanel(tk.Toplevel):
    def __init__(self, master, title, fg="#eee"):
        super().__init__(master)
        self.title(title)
        self.configure(bg="#151515")
        self.resizable(True, True)
        self.text = tk.Text(self, width=70, height=16, bg="#0d0d0d", fg=fg,
                            insertbackground=fg, font=("Consolas", 11), bd=0, highlightthickness=0)
        self.text.pack(fill="both", expand=True)
        self.bind("<ButtonPress-1>", self._start_move)
        self.bind("<B1-Motion>", self._on_move)
        self._ox = self._oy = 0
    def _start_move(self, e): self._ox, self._oy = e.x, e.y
    def _on_move(self, e): self.geometry(f"+{self.winfo_x()+e.x-self._ox}+{self.winfo_y()+e.y-self._oy}")
    def set(self, s):
        self.text.config(state="normal"); self.text.delete("1.0", "end"); self.text.insert("end", s); self.text.config(state="disabled")

class Editor(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MOS Equation Lab v3 — editor")
        self.geometry("1280x800")
        self.configure(bg="#101010")
        self.circ = Circuit()

        self._build_left()
        self._build_canvas()
        self._build_panels()

        self.selected_node = None
        self.selected_dev  = None
        self.dragging_node = None
        self.bind("<Delete>", lambda e: self.delete_selected())

        self.refresh()
        self.running = False
        self.var_dt = tk.DoubleVar(value=1e-6)
        self.var_tmax = tk.DoubleVar(value=1e-3)
        self.var_csv = tk.StringVar(value="")


    # Left panel: palette & inspector
    def _build_left(self):
        p = tk.Frame(self, bg="#171717"); p.pack(side="left", fill="y", padx=6, pady=6)
        # VDD control
        tk.Label(p, text="VDD (V)", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        self.var_vdd = tk.DoubleVar(value=self.circ.nodes["VDD"]["value"])
        tk.Spinbox(p, textvariable=self.var_vdd, from_=0.0, to=20.0, increment=0.01,
                   width=12, font=("Consolas", 11), command=self.on_vdd).pack(anchor="w")
        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)

        # palette
        tk.Label(p, text="Add:", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        tk.Button(p, text="+ Node", command=self.add_node).pack(fill="x", pady=2)
        tk.Button(p, text="+ Resistor", command=lambda:self.add_dev("R")).pack(fill="x", pady=2)
        tk.Button(p, text="+ NMOS", command=lambda:self.add_dev("NMOS")).pack(fill="x", pady=2)
        tk.Button(p, text="+ PMOS", command=lambda:self.add_dev("PMOS")).pack(fill="x", pady=2)
        tk.Button(p, text="+ CurrentSrc (DC)", command=lambda:self.add_dev("I")).pack(fill="x", pady=2)
        tk.Button(p, text="+ Capacitor", command=lambda:self.add_dev("C")).pack(fill="x", pady=2)
        tk.Button(p, text="+ Depletion NMOS", command=lambda:self.add_dev("NMOS_DEP")).pack(fill="x", pady=2)

        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        tk.Button(p, text="Solve DC", command=self.refresh).pack(fill="x", pady=3)
        tk.Button(p, text="Export netlist…", command=self.export_netlist).pack(fill="x", pady=3)
        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        tk.Button(p, text="Delete Selected", command=self.delete_selected).pack(fill="x", pady=3)
        tk.Button(p, text="Cleanup Unused Nodes", command=self.cleanup_unused_nodes).pack(fill="x", pady=3)

        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        # Inspector (selected device)
        self.inspector = tk.Frame(p, bg="#171717"); self.inspector.pack(fill="x", pady=6)
        tk.Label(self.inspector, text="Inspector", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        self.ins_text = tk.Text(self.inspector, width=30, height=16, bg="#0d0d0d", fg="#ccc",
                                insertbackground="#ccc", font=("Consolas", 10), bd=0, highlightthickness=0)
        self.ins_text.pack(fill="both", expand=True)
        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        tk.Label(p, text="Transient", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        row = tk.Frame(p, bg="#171717"); row.pack(fill="x", pady=2)
        tk.Label(row, text="dt (s)", fg="#ddd", bg="#171717").pack(side="left")
        tk.Entry(row, textvariable=self.var_dt, width=10).pack(side="left", padx=4)
        row2 = tk.Frame(p, bg="#171717"); row2.pack(fill="x", pady=2)
        tk.Label(row2, text="t_max (s)", fg="#ddd", bg="#171717").pack(side="left")
        tk.Entry(row2, textvariable=self.var_tmax, width=10).pack(side="left", padx=4)
        row3 = tk.Frame(p, bg="#171717"); row3.pack(fill="x", pady=2)
        tk.Label(row3, text="CSV path", fg="#ddd", bg="#171717").pack(side="left")
        tk.Entry(row3, textvariable=self.var_csv, width=18).pack(side="left", padx=4)

        tk.Button(p, text="Step dt", command=self.transient_step).pack(fill="x", pady=2)
        tk.Button(p, text="Run", command=self.transient_run).pack(fill="x", pady=2)
        tk.Button(p, text="Stop", command=self.transient_stop).pack(fill="x", pady=2)

    def transient_step(self):
        dt = float(self.var_dt.get())
        res = self.circ.solve_transient_step(dt)
        if res is None:
            messagebox.showwarning("Transient", "No solution this step.")
            return
        sol, meta, _ = res
        # append CSV if path provided
        path = self.var_csv.get().strip()
        if path:
            headers = ["t"] + [n for n in sol.keys()]
            import os
            write_header = not os.path.exists(path)
            with open(path, "a") as f:
                if write_header:
                    f.write(",".join(headers) + "\n")
                f.write(",".join([f"{self.circ.t}"] + [f"{sol[n]}" for n in sol.keys()]) + "\n")
        # refresh panels to show latest node voltages/KCL
        self.refresh()

    def transient_run(self):
        if self.running: return
        self.running = True
        tmax = float(self.var_tmax.get())
        def _loop():
            if not self.running or self.circ.t >= tmax:
                self.running = False
                return
            self.transient_step()
            # schedule next iteration
            self.after(1, _loop)
        _loop()

    def transient_stop(self):
        self.running = False


    def _build_canvas(self):
        self.canvas = tk.Canvas(self, width=900, height=780, bg="#0d0d0d", highlightthickness=0)
        self.canvas.pack(side="right", fill="both", expand=True, padx=6, pady=6)
        self.canvas.bind("<Button-1>", self.on_canvas_click)
        self.canvas.bind("<B1-Motion>", self.on_canvas_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_canvas_release)

    def _build_panels(self):
        self.panel_kcl = DraggablePanel(self, "Node equations (KCL)", fg="#cfe"); self.panel_kcl.geometry("+160+60")
        self.panel_dev = DraggablePanel(self, "Device equations", fg="#cfc"); self.panel_dev.geometry("+780+60")

    # ---------- palette actions
    def on_vdd(self):
        self.circ.nodes["VDD"]["value"] = float(self.var_vdd.get()); self.refresh()

    def add_node(self):
        name = f"N{len(self.circ.nodes)+1}"
        self.circ.nodes[name] = {"x": 220, "y": 200, "fixed": False}
        self.refresh()

    def add_dev(self, kind):
        nid = uuid.uuid4().hex[:4].upper()
        if kind == "R":
            # make sure default endpoints exist
            self.circ.ensure_node("N1")
            self.circ.ensure_node("N2")
            self.circ.devices.append(Resistor(f"R{nid}", "N1", "N2", 10000))
        elif kind == "NMOS":
            self.circ.ensure_node("N1")
            self.circ.ensure_node("N2")
            # default: D=N2, G=N1, S=GND
            self.circ.devices.append(NMOS(f"M{nid}", "N2", "N1", "GND", 2e-4, 0.6, 0.0))
        elif kind == "PMOS":
            self.circ.ensure_node("N1")
            self.circ.ensure_node("N2")
            # default: D=N2, G=N1, S=VDD
            self.circ.devices.append(PMOS(f"MP{nid}", "N2", "N1", "VDD", 2e-4, 0.6, 0.0))
        elif kind == "I":
            self.circ.ensure_node("N1")
            self.circ.ensure_node("N2")
            # +I flows from a->b
            self.circ.devices.append(CurrentSourceDC(f"I{nid}", "N1", "N2", 1e-3))
        elif kind == "C":
            self.circ.ensure_node("N1")
            self.circ.devices.append(Capacitor(f"C{nid}", "N1", "GND", 1e-9))
        elif kind == "NMOS_DEP":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2")
            self.circ.devices.append(NMOS(f"MD{nid}", "N2", "N1", "GND", 2e-4, -0.6, 0.0))

        self.refresh()

def delete_selected(self):
    """Delete the currently selected device or (unused, non-fixed) node."""
    if self.selected_dev:
        try:
            self.circ.devices.remove(self.selected_dev)
        except ValueError:
            pass
        self.selected_dev = None
        self.refresh()
        return

    if self.selected_node:
        n = self.selected_node
        info = self.circ.nodes.get(n, {})
        if info.get("fixed", False):
            messagebox.showwarning("Delete node", f"'{n}' is fixed and cannot be deleted.")
            return
        # check references
        used = []
        for d in self.circ.devices:
            pins = []
            if hasattr(d, "pins"): pins = [node for _, node in d.pins]
            if n in pins:
                used.append(d.name)
        if used:
            messagebox.showwarning("Delete node",
                                   f"Node '{n}' is still used by: {', '.join(used)}.\n"
                                   "Rewire those devices first.")
            return
        self.circ.nodes.pop(n, None)
        self.selected_node = None
        self.refresh()

    def cleanup_unused_nodes(self):
        """Remove all non-fixed nodes that no device references."""
        referenced = set()
        for d in self.circ.devices:
            if hasattr(d, "pins"):
                for _, node in d.pins:
                    referenced.add(node)
        to_delete = [n for n, inf in self.circ.nodes.items()
                    if not inf.get("fixed", False) and n not in referenced]
        for n in to_delete:
            self.circ.nodes.pop(n, None)
        messagebox.showinfo("Cleanup", f"Removed {len(to_delete)} unused node(s).")
        self.refresh()

    # ---------- canvas interactions
    def nodes_pos(self):
        pos = {n:(inf["x"], inf["y"]) for n, inf in self.circ.nodes.items()}
        # ensure placeholders for any referenced-but-missing nodes
        for d in self.circ.devices:
            if hasattr(d, "pins"):
                for _, node in d.pins:
                    if node not in pos and node not in ("GND", "VDD"):
                        # place near canvas center
                        pos[node] = (450, 300)
                        # and add it to the circuit so it's editable later
                        self.circ.ensure_node(node, *pos[node])
        return pos


    def find_node_at(self, x, y, r=10):
        for name, inf in self.circ.nodes.items():
            nx, ny = inf["x"], inf["y"]
            if (x-nx)**2 + (y-ny)**2 <= r*r:
                return name
        return None

    def on_canvas_click(self, e):
        n = self.find_node_at(e.x, e.y)
        if n:
            self.selected_node = n; self.selected_dev = None; self.dragging_node = n
            self.update_inspector(); return
        # select nearest device by bounding box proximity (simple)
        self.selected_node = None
        self.selected_dev = self.pick_device_near(e.x, e.y)
        self.dragging_node = None
        self.update_inspector()

    def on_canvas_drag(self, e):
        if self.dragging_node:
            inf = self.circ.nodes[self.dragging_node]
            inf["x"], inf["y"] = e.x, e.y
            self.refresh(redraw_only=True)

    def on_canvas_release(self, e):
        self.dragging_node = None

    def pick_device_near(self, x, y):
        # Choose device whose centroid is closest to (x,y) via current node positions
        nodes_pos = self.nodes_pos()
        mind, best = 99999, None
        for d in self.circ.devices:
            # estimate centroid by averaging pin coords (or stored x,y later)
            ps = []
            if isinstance(d, Resistor) or isinstance(d, CurrentSourceDC) or isinstance(d, Capacitor):
                ps = [nodes_pos[d.a], nodes_pos[d.b]]
            elif isinstance(d, (NMOS, PMOS)):
                ps = [nodes_pos[d.d], nodes_pos[d.g], nodes_pos[d.s]]
            if not ps: continue
            cx = sum(p[0] for p in ps)/len(ps); cy = sum(p[1] for p in ps)/len(ps)
            dist = (cx-x)**2 + (cy-y)**2
            if dist < mind: mind, best = dist, d
        return best

    # ---------- rendering & panels
    def refresh(self, redraw_only=False):
        C = self.canvas; C.delete("all")
        nodes_pos = self.nodes_pos()

        # Wires are drawn by each device; draw VDD rail label + ground symbols
        C.create_text(self.circ.nodes["VDD"]["x"]-30, self.circ.nodes["VDD"]["y"],
                      text="VDD", fill="#58a6ff", font=("Consolas", 12))
        # draw devices
        for d in self.circ.devices:
            d.draw(C, nodes_pos)

        # draw nodes (draggable)
        for name, inf in self.circ.nodes.items():
            x, y = inf["x"], inf["y"]
            color = "#58a6ff" if name == "VDD" else ("#888" if name=="GND" else "#fff")
            C.create_oval(x-5, y-5, x+5, y+5, fill=color, outline="")
            C.create_text(x+8, y-12, text=name, fill=color, font=("Consolas", 10), anchor="w")

        if redraw_only: return

        # Solve and update panels
        sol, meta, eqs = self.circ.solve_dc()
        if sol is None:
            self.panel_kcl.set("No DC solution (for current params/regions).")
            self.panel_dev.set("No DC solution.")
        else:
            # Node KCL panel
            s1 = ["Node voltages:"]
            for n, v in sol.items(): s1.append(f"  {n} = {v:.9f} V")
            s1.append("\nKCL (symbolic):")
            for n, txt in meta["kcl_sym"].items(): s1.append(f"[{n}]\n{txt}")
            s1.append("\nKCL (numeric):")
            for n, txt in meta["kcl_num"].items(): s1.append(f"[{n}]\n{txt}")
            self.panel_kcl.set("\n".join(s1))

            # Device equations panel
            s2 = []
            for name, dinfo in meta["dev"].items():
                s2.append(f"{name}:  region={dinfo['region']}")
                if name.startswith("M") and "Vgs" in dinfo:
                    s2.append(f"  VGS={dinfo['Vgs']:.9f} V, VDS={dinfo['Vds']:.9f} V")
                    s2.append("  Id (numeric): " + dinfo["Id_num_str"] + " A\n")
                else:
                    s2.append(f"  VSG={dinfo['Vsg']:.9f} V, VSD={dinfo['Vsd']:.9f} V")
                    s2.append("  I_d_leaving (numeric): " + dinfo["Id_d_num_str"] + " A\n")
            self.panel_dev.set("\n".join(s2))

        # Inspector view (selected)
        self.update_inspector()

    def update_inspector(self):
        t = []
        if self.selected_node:
            inf = self.circ.nodes[self.selected_node]
            t.append(f"[Node] {self.selected_node}")
            t.append(f"  fixed: {inf.get('fixed', False)}")
            if inf.get("fixed", False): t.append(f"  value: {inf.get('value', 0.0)} V")
        elif self.selected_dev:
            d = self.selected_dev
            t.append(f"[{d.kind}] {d.name}")
            if isinstance(d, Resistor):
                t += [f"  a: {d.a}", f"  b: {d.b}", f"  R: {d.R:g} Ω"]
            elif isinstance(d, CurrentSourceDC):
                t += [f"  a: {d.a}", f"  b: {d.b}", f"  I: {d.I:g} A"]
            elif isinstance(d, Capacitor):
                t += [f"  a: {d.a}", f"  b: {d.b}", f"  C: {d.C:g} F (DC open)"]
            elif isinstance(d, NMOS):
                t += [f"  d: {d.d}", f"  g: {d.g}", f"  s: {d.s}",
                      f"  k: {d.k:g} A/V²", f"  Vt: {d.Vt:g} V", f"  λ: {d.lam:g} 1/V"]
            elif isinstance(d, PMOS):
                t += [f"  d: {d.d}", f"  g: {d.g}", f"  s: {d.s}",
                      f"  k: {d.k:g} A/V²", f"  |Vt|: {d.Vt:g} V", f"  λ: {d.lam:g} 1/V"]
            t.append("\nTip: edit selection via the dropdowns below.")
        else:
            t.append("Click a node to drag it.\nClick a device to inspect/edit.")

        self.ins_text.config(state="normal"); self.ins_text.delete("1.0","end"); self.ins_text.insert("end","\n".join(t)); self.ins_text.config(state="disabled")

        # build simple editing widgets under the inspector
        for w in self.inspector.pack_slaves():
            if isinstance(w, tk.Frame) and getattr(w, "_is_editor", False):
                w.destroy()
        if self.selected_dev:
            ed = tk.Frame(self.inspector, bg="#171717"); ed._is_editor=True; ed.pack(fill="x", pady=4)
            nodes = list(self.circ.nodes.keys())
            def add_dropdown(lbl, getter, setter):
                row = tk.Frame(ed, bg="#171717"); row.pack(fill="x", pady=2)
                tk.Label(row, text=lbl, fg="#ddd", bg="#171717", font=("Consolas", 10)).pack(side="left")
                var = tk.StringVar(value=getter())
                om = ttk.OptionMenu(row, var, getter(), *nodes, command=lambda v:setter(v) or self.refresh())
                om.pack(side="left")
            def add_entry(lbl, getter, setter):
                row = tk.Frame(ed, bg="#171717"); row.pack(fill="x", pady=2)
                tk.Label(row, text=lbl, fg="#ddd", bg="#171717", font=("Consolas", 10)).pack(side="left")
                var = tk.DoubleVar(value=getter())
                ent = tk.Spinbox(row, textvariable=var, from_=-1e6, to=1e6, increment=0.001,
                                 width=12, font=("Consolas", 10), command=lambda:setter(var.get()) or self.refresh())
                ent.pack(side="left")

            d = self.selected_dev
            if isinstance(d, Resistor):
                add_dropdown("a", lambda:d.a, lambda v:setattr(d,"a",v))
                add_dropdown("b", lambda:d.b, lambda v:setattr(d,"b",v))
                add_entry("R (Ω)", lambda:d.R, lambda v:setattr(d,"R",float(v)))
            elif isinstance(d, CurrentSourceDC):
                add_dropdown("a", lambda:d.a, lambda v:setattr(d,"a",v))
                add_dropdown("b", lambda:d.b, lambda v:setattr(d,"b",v))
                add_entry("I (A)", lambda:d.I, lambda v:setattr(d,"I",float(v)))
            elif isinstance(d, Capacitor):
                add_dropdown("a", lambda:d.a, lambda v:setattr(d,"a",v))
                add_dropdown("b", lambda:d.b, lambda v:setattr(d,"b",v))
                add_entry("C (F)", lambda:d.C, lambda v:setattr(d,"C",float(v)))
            elif isinstance(d, NMOS):
                add_dropdown("d", lambda:d.d, lambda v:setattr(d,"d",v))
                add_dropdown("g", lambda:d.g, lambda v:setattr(d,"g",v))
                add_dropdown("s", lambda:d.s, lambda v:setattr(d,"s",v))
                add_entry("k (A/V²)", lambda:d.k, lambda v:setattr(d,"k",float(v)))
                add_entry("Vt (V)",   lambda:d.Vt, lambda v:setattr(d,"Vt",float(v)))
                add_entry("λ (1/V)",  lambda:d.lam, lambda v:setattr(d,"lam",float(v)))
            elif isinstance(d, PMOS):
                add_dropdown("d", lambda:d.d, lambda v:setattr(d,"d",v))
                add_dropdown("g", lambda:d.g, lambda v:setattr(d,"g",v))
                add_dropdown("s", lambda:d.s, lambda v:setattr(d,"s",v))
                add_entry("k (A/V²)", lambda:d.k, lambda v:setattr(d,"k",float(v)))
                add_entry("|Vt| (V)", lambda:d.Vt, lambda v:setattr(d,"Vt",float(v)))
                add_entry("λ (1/V)",  lambda:d.lam, lambda v:setattr(d,"lam",float(v)))

    def export_netlist(self):
        path = filedialog.asksaveasfilename(defaultextension=".cir", filetypes=[("Spice netlist",".cir"),("All","*.*")])
        if not path: return
        try:
            self.circ.export_netlist(path)
        except Exception as e:
            messagebox.showerror("Export failed", str(e))
            return
        messagebox.showinfo("Export", f"Netlist written to:\n{path}")

if __name__ == "__main__":
    Editor().mainloop()
