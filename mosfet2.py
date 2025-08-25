#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MOS Equation Lab v2 — multi-node, region-aware, draggable equation windows.

Topology in this prototype:
  VDD -- Rb --> node X
  node X -> NMOS M1 (diode-connected): D=X, G=X, S=GND
  node Y -> NMOS M2 (mirrored):        D=Y, G=X, S=GND
  (Y is the output node, can be loaded by Rload to VDD)

We assemble KCL at nodes X and Y and solve for {Vx, Vy} under consistent
region assignments for each MOSFET (cutoff / linear / saturation).
Equation panels are floating (draggable).
"""

import tkinter as tk
from tkinter import ttk
from math import isfinite
import itertools
import sympy as sp

# --------------------------- Symbolic primitives ------------------------------

# Node voltages (unknowns)
Vx, Vy = sp.symbols('Vx Vy', real=True)

# Parameters
VDD, Rb, Rload = sp.symbols('VDD Rb Rload', positive=True, real=True)

# Common MOS params (per-device will have k, Vt; λ optional)
# For extensibility we keep device-local params.
def nmos_currents(Vgs, Vds, k, Vt, lam=None, region=None):
    """
    Return (Id_expr, region_predicates, pretty_name) for NMOS.
    - region can be 'cutoff', 'linear', 'sat' to force branch (skip preds),
      or None to return predicates to be used by caller.
    """
    # Square-law
    Id_lin = k * ((Vgs - Vt) * Vds - sp.Rational(1, 2) * Vds**2)
    Id_sat = sp.Rational(1, 2) * k * (Vgs - Vt) ** 2
    if lam is not None and lam != 0:
        Id_sat = Id_sat * (1 + lam * Vds)

    # Predicates
    preds = {
        "cutoff": (Vgs <= Vt),
        "linear": sp.And(Vgs >= Vt, Vds <= Vgs - Vt),
        "sat":    sp.And(Vgs >= Vt, Vds >= Vgs - Vt),
    }

    if region == "cutoff":
        return sp.Integer(0), [], "Id=0"
    elif region == "linear":
        return Id_lin, [], "Id = k*((Vgs-Vt)·Vds - 0.5·Vds²)"
    elif region == "sat":
        return Id_sat, [], "Id = 0.5·k·(Vgs-Vt)² · (1+λVds)" if lam else "Id = 0.5·k·(Vgs-Vt)²"

    # If no forced region, return piecewise with predicates (unused here)
    return sp.Piecewise(
        (0, preds["cutoff"]),
        (Id_lin, preds["linear"]),
        (Id_sat, preds["sat"])
    ), preds, "piecewise"

# --------------------------- Circuit model ------------------------------------

class NMOS:
    def __init__(self, name, d, g, s, k, Vt, lam=0.0):
        self.name = name
        self.d, self.g, self.s = d, g, s
        self.k = float(k)
        self.Vt = float(Vt)
        self.lam = float(lam)

    def Id_symbolic(self, Vmap, region):
        Vd = Vmap[self.d]; Vg = Vmap[self.g]; Vs = Vmap[self.s]
        Vgs = Vg - Vs
        Vds = Vd - Vs
        k = sp.Symbol(f"k_{self.name}")
        Vt = sp.Symbol(f"Vt_{self.name}")
        lam = sp.Symbol(f"lam_{self.name}")
        Id, _, pretty = nmos_currents(Vgs, Vds, k, Vt, lam=lam if self.lam!=0 else None, region=region)
        subs = {k: self.k, Vt: self.Vt}
        if self.lam != 0:
            subs[lam] = self.lam
        Id = sp.simplify(Id.subs(subs))
        return sp.simplify(Id), Vgs, Vds, pretty

class Resistor:
    def __init__(self, name, a, b, R):
        self.name = name
        self.a, self.b = a, b
        self.R = float(R)
    def I_from_a_to_b(self, Vmap):
        Va = Vmap[self.a]; Vb = Vmap[self.b]
        return (Va - Vb) / self.R

class VSourceDC:
    def __init__(self, name, pos, neg, value):
        self.name = name
        self.pos, self.neg = pos, neg
        self.value = float(value)

# --------------------------- The specific topology ----------------------------

class MirrorCircuit:
    """
    VDD -- Rb --> X
    M1: NMOS D=X, G=X, S=GND (diode-connected)
    M2: NMOS D=Y, G=X, S=GND
    Optional: Rload from VDD to Y
    """
    def __init__(self):
        # default params
        self.params = {
            "VDD": 3.3,
            "Rb":  10000.0,
            "Rload": 20000.0,
            "k_M1": 2e-4,
            "k_M2": 2e-4,
            "Vt_M1": 0.6,
            "Vt_M2": 0.6,
            "lam_M1": 0.0,
            "lam_M2": 0.0
        }
        # elements
        self.Rb = Resistor("Rb", "VDD", "X", self.params["Rb"])
        self.Rload = Resistor("Rload", "VDD", "Y", self.params["Rload"])
        self.M1 = NMOS("M1", d="X", g="X", s="GND", k=self.params["k_M1"], Vt=self.params["Vt_M1"], lam=self.params["lam_M1"])
        self.M2 = NMOS("M2", d="Y", g="X", s="GND", k=self.params["k_M2"], Vt=self.params["Vt_M2"], lam=self.params["lam_M2"])
        self.VDD = VSourceDC("VDD", "VDD", "GND", self.params["VDD"])

    def vmap(self):
        return {"X": Vx, "Y": Vy, "VDD": VDD, "GND": sp.Integer(0)}

    def kcl_equations(self, regions):
        """
        Build KCL at X and Y given region assignment dict:
          regions = {"M1": "cutoff|linear|sat", "M2": ...}
        Returns symbolic equations [eq_X, eq_Y], device info.
        """
        Vm = self.vmap()
        # Currents *leaving* node positive
        # Node X:
        I_Rb_X = - self.Rb.I_from_a_to_b(Vm)   # current leaving X to VDD is negative of VDD->X
        Id_M1, Vgs1, Vds1, pretty1 = self.M1.Id_symbolic(Vm, regions["M1"])
        # Node Y:
        I_Rload_Y = - self.Rload.I_from_a_to_b(Vm)
        Id_M2, Vgs2, Vds2, pretty2 = self.M2.Id_symbolic(Vm, regions["M2"])

        # KCL: sum currents leaving node = 0
        eqX = sp.Eq(I_Rb_X + Id_M1, 0)
        eqY = sp.Eq(I_Rload_Y + Id_M2, 0)

        dev_info = {
            "M1": {"Id": Id_M1, "Vgs": Vgs1, "Vds": Vds1, "pretty": pretty1, "region": regions["M1"]},
            "M2": {"Id": Id_M2, "Vgs": Vgs2, "Vds": Vds2, "pretty": pretty2, "region": regions["M2"]},
        }
        return [eqX, eqY], dev_info

    def solve(self):
        """
        Try all consistent region assignments and return the first that yields a valid solution
        satisfying the region inequalities. Returns (solution dict, dev_info, equations).
        """
        # ordered attempt preference: try sat/linear first (typical for mirrors)
        region_options = ["sat", "linear", "cutoff"]
        Vm = self.vmap()

        # substitutions for numeric values
        subs = {
            VDD: self.params["VDD"],
        }
        # refresh element params to latest values
        self.Rb.R = float(self.params["Rb"])
        self.Rload.R = float(self.params["Rload"])
        self.M1.k = float(self.params["k_M1"]); self.M1.Vt = float(self.params["Vt_M1"]); self.M1.lam = float(self.params["lam_M1"])
        self.M2.k = float(self.params["k_M2"]); self.M2.Vt = float(self.params["Vt_M2"]); self.M2.lam = float(self.params["lam_M2"])

        for r1, r2 in itertools.product(region_options, repeat=2):
            eqs, dinfo = self.kcl_equations({"M1": r1, "M2": r2})
            try:
                # Solve linear or nonlinear system for (Vx, Vy)
                sol = sp.solve([eq.subs(subs) for eq in eqs], (Vx, Vy), dict=True)
            except Exception:
                sol = []
            # sp.solve may give param families; we only accept numeric
            for s in sol:
                try:
                    vx = float(s[Vx]); vy = float(s[Vy])
                except Exception:
                    continue
                if not (isfinite(vx) and isfinite(vy)):
                    continue

                # Check region consistency
                ok = True
                for name, info in dinfo.items():
                    Vgs = float(info["Vgs"].subs({Vx: vx, Vy: vy, VDD: self.params["VDD"]}))
                    Vds = float(info["Vds"].subs({Vx: vx, Vy: vy, VDD: self.params["VDD"]}))
                    Vt = self.params[f"Vt_{name}"]
                    if info["region"] == "cutoff" and not (Vgs <= Vt + 1e-12):
                        ok = False
                    elif info["region"] == "linear" and not (Vgs >= Vt - 1e-12 and Vds <= (Vgs - Vt) + 1e-12):
                        ok = False
                    elif info["region"] == "sat" and not (Vgs >= Vt - 1e-12 and Vds >= (Vgs - Vt) - 1e-12):
                        ok = False
                    if not ok:
                        break
                if ok and 0 <= vx <= self.params["VDD"] + 1e-9 and 0 <= vy <= self.params["VDD"] + 1e-9:
                    # Build numeric KCL strings
                    eqs_sym, _ = self.kcl_equations({"M1": r1, "M2": r2})
                    kclX_sym = sp.pretty(eqs_sym[0])
                    kclY_sym = sp.pretty(eqs_sym[1])
                    kclX_num = sp.pretty(sp.Eq(sp.nsimplify(eqs_sym[0].lhs.subs({Vx: vx, Vy: vy, VDD: self.params["VDD"]}), rational=False),
                                               sp.nsimplify(eqs_sym[0].rhs, rational=False)))
                    kclY_num = sp.pretty(sp.Eq(sp.nsimplify(eqs_sym[1].lhs.subs({Vx: vx, Vy: vy, VDD: self.params["VDD"]}), rational=False),
                                               sp.nsimplify(eqs_sym[1].rhs, rational=False)))

                    # Device currents numeric
                    for name, info in dinfo.items():
                        Id_num = float(info["Id"].subs({Vx: vx, Vy: vy, VDD: self.params["VDD"]}))
                        info["Id_num"] = Id_num
                        info["Vgs_num"] = float(info["Vgs"].subs({Vx: vx, Vy: vy, VDD: self.params["VDD"]}))
                        info["Vds_num"] = float(info["Vds"].subs({Vx: vx, Vy: vy, VDD: self.params["VDD"]}))
                        info["Id_num_str"] = sp.pretty(sp.nsimplify(info["Id"].subs({Vx: vx, Vy: vy, VDD: self.params["VDD"]}), rational=False))

                    return {
                        "Vx": vx, "Vy": vy,
                        "regions": {"M1": r1, "M2": r2},
                        "kclX_sym": kclX_sym, "kclY_sym": kclY_sym,
                        "kclX_num": kclX_num, "kclY_num": kclY_num
                    }, dinfo, eqs_sym
        # If nothing worked:
        return None, None, None

# --------------------------- GUI with draggable panels -------------------------

class DraggablePanel(tk.Toplevel):
    def __init__(self, master, title, fg="#eee"):
        super().__init__(master)
        self.title(title)
        self.configure(bg="#141414")
        self.resizable(True, True)
        self.attributes("-topmost", False)
        self.text = tk.Text(self, width=60, height=14, bg="#0d0d0d", fg=fg,
                            insertbackground=fg, font=("Consolas", 11), bd=0, highlightthickness=0)
        self.text.pack(fill="both", expand=True)
        # drag by title bar substitute
        self.bind("<ButtonPress-1>", self._start_move)
        self.bind("<B1-Motion>", self._on_move)
        self._ox = self._oy = 0

    def _start_move(self, e):
        self._ox = e.x; self._oy = e.y
    def _on_move(self, e):
        self.geometry(f"+{self.winfo_x() + e.x - self._ox}+{self.winfo_y() + e.y - self._oy}")

    def set_content(self, s):
        self.text.config(state="normal")
        self.text.delete("1.0", "end")
        self.text.insert("end", s)
        self.text.config(state="disabled")

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MOS Equation Lab v2 — Current Mirror")
        self.geometry("1080x720")
        self.configure(bg="#101010")

        self.circ = MirrorCircuit()
        self._build_left_panel()
        self._build_canvas()
        self._build_floating_panels()
        self.refresh()

    def _build_left_panel(self):
        p = tk.Frame(self, bg="#171717")
        p.pack(side="left", fill="y", padx=6, pady=6)

        def add_spin(label, varname, from_, to_, step, fmt="{:.4g}", unit=""):
            row = tk.Frame(p, bg="#171717"); row.pack(fill="x", pady=5)
            tk.Label(row, text=label, fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
            v = tk.DoubleVar(value=self.circ.params[varname])
            entry = tk.Spinbox(row, textvariable=v, from_=from_, to=to_, increment=step,
                               width=14, font=("Consolas", 11), command=self.refresh)
            entry.pack(side="left")
            tk.Label(row, text=unit, fg="#888", bg="#171717", font=("Consolas", 10)).pack(side="left", padx=4)
            setattr(self, f"var_{varname}", v)

        add_spin("VDD",   "VDD",   0.0, 10.0, 0.01, unit="V")
        add_spin("Rb",    "Rb",    100.0, 1e6, 10.0, unit="Ω")
        add_spin("Rload", "Rload", 100.0, 1e6, 10.0, unit="Ω")
        add_spin("k_M1",  "k_M1",  1e-6, 5e-3, 1e-6, unit="A/V²")
        add_spin("k_M2",  "k_M2",  1e-6, 5e-3, 1e-6, unit="A/V²")
        add_spin("Vt_M1", "Vt_M1", 0.0,  2.0,  0.001, unit="V")
        add_spin("Vt_M2", "Vt_M2", 0.0,  2.0,  0.001, unit="V")
        add_spin("lam_M1","lam_M1",0.0,  0.1,  0.001, unit="1/V")
        add_spin("lam_M2","lam_M2",0.0,  0.1,  0.001, unit="1/V")

        tk.Button(p, text="Refresh", command=self.refresh).pack(pady=8)

    def _build_canvas(self):
        self.canvas = tk.Canvas(self, width=680, height=700, bg="#0d0d0d", highlightthickness=0)
        self.canvas.pack(side="right", fill="both", expand=True, padx=6, pady=6)

    def _build_floating_panels(self):
        self.panel_nodeX = DraggablePanel(self, "Node X — KCL", fg="#cfe")
        self.panel_nodeY = DraggablePanel(self, "Node Y — KCL", fg="#cfe")
        self.panel_M1    = DraggablePanel(self, "M1 Equation",  fg="#cfc")
        self.panel_M2    = DraggablePanel(self, "M2 Equation",  fg="#cfc")
        # initial placements
        self.panel_nodeX.geometry("+120+60")
        self.panel_nodeY.geometry("+120+340")
        self.panel_M1.geometry("+840+80")
        self.panel_M2.geometry("+840+360")

    def refresh(self):
        # pull params
        for k in self.circ.params.keys():
            v = getattr(self, f"var_{k}", None)
            if isinstance(v, tk.DoubleVar):
                self.circ.params[k] = float(v.get())

        sol, dinfo, eqs = self.circ.solve()
        C = self.canvas
        C.delete("all")

        # Draw simple schematic
        left_x = 90; right_x = 520; mid_x = 305
        y_top = 120; y_mid = 330; y_bot = 560

        # VDD rail
        C.create_line(left_x, y_top, right_x, y_top, fill="#58a6ff", width=2)
        C.create_text(left_x-30, y_top, text="VDD", fill="#58a6ff", font=("Consolas", 12))

        # Rb from VDD to X
        C.create_line(left_x+40, y_top, left_x+40, y_mid-40, fill="#aaa", width=2)
        # resistor squiggle
        x = left_x+40; y = (y_top + y_mid-40)/2; rz = 10
        for i in range(6):
            x1 = x - 12 + i*12
            C.create_line(x1, y - rz if i%2==0 else y + rz, x1+12, y + rz if i%2==0 else y - rz, fill="#aaa", width=2)
        C.create_line(left_x+40, y_mid-40, mid_x, y_mid-40, fill="#aaa", width=2)
        C.create_text(left_x+12, y, text=f"Rb={self.circ.params['Rb']:.0f}Ω", fill="#aaa", font=("Consolas", 10), anchor="w")

        # Node X dot
        C.create_line(mid_x, y_mid-40, mid_x, y_mid, fill="#aaa", width=2)
        C.create_oval(mid_x-3, y_mid-3, mid_x+3, y_mid+3, fill="#fff", outline="")

        # M1 diode-connected box
        bx, by, bw, bh = mid_x-60, y_mid, 120, 80
        C.create_rectangle(bx, by, bx+bw, by+bh, outline="#7fd", width=2)
        C.create_text(bx+bw/2, by-14, text="NMOS M1 (D=G=X,S=GND)", fill="#7fd", font=("Consolas", 10))
        # M1 drain/gate tied to node X
        C.create_line(mid_x, y_mid-40, mid_x, by, fill="#aaa", width=2)
        # M1 source to ground
        C.create_line(bx+bw/2, by+bh, bx+bw/2, y_bot, fill="#aaa", width=2)
        # ground symbol
        gx = bx+bw/2; gy = y_bot
        C.create_line(gx-14, gy, gx+14, gy, fill="#888", width=2)
        C.create_line(gx-10, gy+6, gx+10, gy+6, fill="#888", width=2)
        C.create_line(gx-6,  gy+12, gx+6,  gy+12, fill="#888", width=2)

        # M2 box at right
        bx2, by2, bw2, bh2 = right_x-160, y_mid, 120, 80
        C.create_rectangle(bx2, by2, bx2+bw2, by2+bh2, outline="#7fd", width=2)
        C.create_text(bx2+bw2/2, by2-14, text="NMOS M2 (D=Y,G=X,S=GND)", fill="#7fd", font=("Consolas", 10))
        # Gate from X to M2 gate midline
        C.create_line(mid_x, y_mid-0, bx2, y_mid-0, fill="#aaa", width=2)
        # M2 source to ground
        C.create_line(bx2+bw2/2, by2+bh2, bx2+bw2/2, y_bot, fill="#aaa", width=2)
        C.create_line(bx2+bw2/2-14, y_bot, bx2+bw2/2+14, y_bot, fill="#888", width=2)
        C.create_line(bx2+bw2/2-10, y_bot+6, bx2+bw2/2+10, y_bot+6, fill="#888", width=2)
        C.create_line(bx2+bw2/2-6,  y_bot+12, bx2+bw2/2+6,  y_bot+12, fill="#888", width=2)

        # Rload from VDD to Y
        C.create_line(right_x-40, y_top, right_x-40, y_mid-40, fill="#aaa", width=2)
        xR = right_x-40; yR = (y_top + y_mid-40)/2; rz = 10
        for i in range(6):
            x1 = xR - 12 + i*12
            C.create_line(x1, yR - rz if i%2==0 else yR + rz, x1+12, yR + rz if i%2==0 else yR - rz, fill="#aaa", width=2)
        C.create_line(right_x-40, y_mid-40, bx2+bw2/2, y_mid-40, fill="#aaa", width=2)
        C.create_text(right_x-120, yR, text=f"Rload={self.circ.params['Rload']:.0f}Ω", fill="#aaa", font=("Consolas", 10))

        # Node Y dot
        C.create_line(bx2+bw2/2, y_mid-40, bx2+bw2/2, by2, fill="#aaa", width=2)
        C.create_oval(bx2+bw2/2-3, y_mid-40-3, bx2+bw2/2+3, y_mid-40+3, fill="#fff", outline="")

        # Results
        if sol:
            C.create_text(120, 40, text=f"Vx={sol['Vx']:.6f} V, Vy={sol['Vy']:.6f} V",
                          fill="#cfc", font=("Consolas", 12), anchor="w")
            C.create_text(120, 62, text=f"Regions: M1={sol['regions']['M1']}  |  M2={sol['regions']['M2']}",
                          fill="#8ad", font=("Consolas", 11), anchor="w")

            # Populate draggable panels
            self.panel_nodeX.set_content(
                "KCL @ Node X (symbolic):\n" + sol["kclX_sym"] + "\n\n" +
                "KCL @ Node X (numeric):\n" + sol["kclX_num"]
            )
            self.panel_nodeY.set_content(
                "KCL @ Node Y (symbolic):\n" + sol["kclY_sym"] + "\n\n" +
                "KCL @ Node Y (numeric):\n" + sol["kclY_num"]
            )
            # Devices
            m1 = dinfo["M1"]; m2 = dinfo["M2"]
            self.panel_M1.set_content(
                f"M1 region: {m1['region']}\n" +
                f"VGS={m1['Vgs_num']:.6f} V, VDS={m1['Vds_num']:.6f} V\n" +
                f"Equation: {m1['pretty']}\n\n" +
                "Id (symbolic):\n" + sp.pretty(m1["Id"]) + "\n\n" +
                "Id (numeric):\n" + m1["Id_num_str"] + "  A"
            )
            self.panel_M2.set_content(
                f"M2 region: {m2['region']}\n" +
                f"VGS={m2['Vgs_num']:.6f} V, VDS={m2['Vds_num']:.6f} V\n" +
                f"Equation: {m2['pretty']}\n\n" +
                "Id (symbolic):\n" + sp.pretty(m2["Id"]) + "\n\n" +
                "Id (numeric):\n" + m2["Id_num_str"] + "  A"
            )
        else:
            C.create_text(120, 40, text="No consistent DC solution for current parameters.",
                          fill="#f88", font=("Consolas", 12), anchor="w")
            self.panel_nodeX.set_content("No solution.")
            self.panel_nodeY.set_content("No solution.")
            self.panel_M1.set_content("No solution.")
            self.panel_M2.set_content("No solution.")


if __name__ == "__main__":
    App().mainloop()
