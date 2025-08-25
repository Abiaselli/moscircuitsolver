#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MOS-Equation Lab: minimal, math-aware schematic prototype.
- NMOS pull-down + resistor load to VDD
- Region-aware Id equations (square-law) shown inline
- KCL at drain node solved analytically
- Live GUI controls for VDD, VG, R, Vt, k (μ·Cox·W/L)

Extendable toward PySpice / NGSpice + real netlists later.
"""
import tkinter as tk
from tkinter import ttk
from math import isfinite
import sympy as sp

# ---------- SYMBOLIC CORE -----------------------------------------------------

Vd = sp.symbols('Vd', real=True)  # unknown drain node voltage
VDD, VG, R, Vt, k = sp.symbols('VDD VG R Vt k', positive=True, real=True)

# NMOS with source at 0 V (ground), body at 0
VGS = VG
VDS = Vd  # because Vs=0

Id_lin_expr = k * ((VGS - Vt) * VDS - sp.Rational(1, 2) * VDS**2)
Id_sat_expr = sp.Rational(1, 2) * k * (VGS - Vt)**2

# Resistor current from VDD -> node (KCL sign convention: currents leaving node positive)
Ir_expr = (VDD - Vd) / R

def solve_node(vdd, vg, r_ohm, vt, k_val):
    """
    Solve KCL at the drain node:
      Ir(Vd) = Id(Vd) with correct MOSFET region selection.
    Returns: dict with Vd, Id, region, equation strings, KCL string.
    """
    # Try saturation first: Id is constant; solve Vd from resistor only.
    subs = {VDD: vdd, VG: vg, R: r_ohm, Vt: vt, k: k_val}
    Id_sat = float(Id_sat_expr.subs(subs))
    Vd_sat = float((VDD - R * Id_sat).subs(subs))

    region = None
    chosen_Id_expr = None
    chosen_eq_text = ""

    def within(x, lo, hi):
        return (x >= lo - 1e-12) and (x <= hi + 1e-12)

    # Saturation consistency: VDS >= VGS - Vt, Vd within [0, VDD]
    if within(Vd_sat, 0.0, vdd) and (Vd_sat >= (vg - vt) - 1e-12) and (vg >= vt):
        # Saturation is self-consistent
        region = "saturation"
        Vd_val = Vd_sat
        Id_val = Id_sat
        chosen_Id_expr = Id_sat_expr
        chosen_eq_text = "Id = 0.5*k*(VGS - Vt)^2"
    else:
        # Solve triode (linear) region: Ir(Vd) = Id_lin(Vd), with constraint Vd < VGS - Vt
        eq = sp.Eq(Ir_expr, Id_lin_expr)
        # Reduce to polynomial in Vd
        poly_eq = sp.expand(eq.lhs - eq.rhs)
        # Solve symbolically, then filter numerically
        sols = list(sp.solve(poly_eq.subs(subs), Vd, dict=False))
        valid = []
        for s in sols:
            try:
                sv = float(s)
            except TypeError:
                continue
            if isfinite(sv) and within(sv, 0.0, vdd) and (sv <= (vg - vt) + 1e-12) and (vg >= vt):
                valid.append(sv)
        if valid:
            # Pick the physically reasonable one (closer to VGS-Vt is often stable, but either is fine)
            Vd_val = max(valid)
            Id_val = float(Id_lin_expr.subs(subs).subs({Vd: Vd_val}))
            region = "linear (triode)"
            chosen_Id_expr = Id_lin_expr
            chosen_eq_text = "Id = k*((VGS - Vt)*VDS - 0.5*VDS^2)"
        else:
            # Cutoff or no conduction: VGS <= Vt or no valid triode solution
            region = "cutoff" if vg <= vt else "no-solution (params)"
            Vd_val = vdd  # node floats up via R (no current)
            Id_val = 0.0
            chosen_Id_expr = sp.Integer(0)
            chosen_eq_text = "Id = 0 (cutoff)"

    # Build pretty strings with numeric substitutions
    pretty = sp.pretty
    # KCL at node: (VDD - Vd)/R = Id(...)
    kcl_lhs = Ir_expr
    kcl_rhs = chosen_Id_expr
    kcl_sym = sp.Eq(kcl_lhs, kcl_rhs)
    kcl_num = kcl_sym.subs(subs).subs({Vd: Vd_val})
    kcl_sym_str = pretty(kcl_sym)
    kcl_num_str = pretty(sp.Eq(sp.nsimplify(kcl_lhs.subs(subs).subs({Vd: Vd_val}), rational=False),
                               sp.nsimplify(kcl_rhs.subs(subs).subs({Vd: Vd_val}), rational=False)))

    # Device equation (region-specific) with numeric substitution
    dev_sym_str = chosen_eq_text
    dev_num_str = pretty(sp.nsimplify(chosen_Id_expr.subs(subs).subs({Vd: Vd_val}), rational=False))

    return {
        "Vd": Vd_val,
        "Id": Id_val,
        "region": region,
        "kcl_sym": kcl_sym_str,
        "kcl_num": kcl_num_str,
        "dev_sym": dev_sym_str,
        "dev_num": dev_num_str
    }

# ---------- GUI ----------------------------------------------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MOS Equation Lab — NMOS + Resistor")
        self.geometry("980x620")
        self.configure(bg="#111")

        # Parameters (with sensible defaults)
        self.vdd = tk.DoubleVar(value=3.3)
        self.vg  = tk.DoubleVar(value=2.0)
        self.r   = tk.DoubleVar(value=10_000.0)
        self.vt  = tk.DoubleVar(value=0.6)
        self.k   = tk.DoubleVar(value=2e-4)  # A/V^2 (e.g., 200 µA/V^2)

        self.build_ui()
        self.refresh()

    def build_ui(self):
        # Controls panel
        ctrl = tk.Frame(self, bg="#1a1a1a")
        ctrl.pack(side="left", fill="y", padx=8, pady=8)

        def add_slider(label, var, from_, to_, step, unit=""):
            row = tk.Frame(ctrl, bg="#1a1a1a")
            row.pack(fill="x", pady=6)
            tk.Label(row, text=label, fg="#ddd", bg="#1a1a1a", font=("Consolas", 11)).pack(anchor="w")
            s = ttk.Scale(row, variable=var, from_=from_, to=to_, orient="horizontal", command=lambda *_: self.refresh())
            s.pack(fill="x")
            box = tk.Entry(row, textvariable=var, width=10, font=("Consolas", 11))
            box.pack(side="left", pady=2)
            tk.Label(row, text=unit, fg="#888", bg="#1a1a1a", font=("Consolas", 10)).pack(side="left", padx=4)

        style = ttk.Style(self)
        style.theme_use("default")
        style.configure("TScale", background="#1a1a1a")

        add_slider("VDD (V)", self.vdd, 0.0, 10.0, 0.01, "V")
        add_slider("VG (V)",  self.vg,  0.0, 10.0, 0.01, "V")
        add_slider("R (Ω)",   self.r,   100.0, 1e6, 1.0, "Ω")
        add_slider("Vt (V)",  self.vt,  0.0,  2.0, 0.001, "V")
        add_slider("k (A/V^2)", self.k, 1e-6, 5e-3, 1e-6, "A/V²")

        # Schematic canvas-ish area (just simple rectangles/lines)
        self.canvas = tk.Canvas(self, width=640, height=600, bg="#0d0d0d", highlightthickness=0)
        self.canvas.pack(side="right", fill="both", expand=True, padx=8, pady=8)

        # Text panels
        self.txt_wire = tk.Text(self.canvas, width=70, height=6, bg="#0d0d0d", fg="#cce", insertbackground="#cce",
                                font=("Consolas", 11), bd=0, highlightthickness=0)
        self.txt_dev  = tk.Text(self.canvas, width=70, height=6, bg="#0d0d0d", fg="#cec", insertbackground="#cec",
                                font=("Consolas", 11), bd=0, highlightthickness=0)

    def refresh(self):
        vdd = float(self.vdd.get())
        vg  = float(self.vg.get())
        r   = float(self.r.get())
        vt  = float(self.vt.get())
        kval= float(self.k.get())

        res = solve_node(vdd, vg, r, vt, kval)

        # Draw schematic elements (simple)
        C = self.canvas
        C.delete("all")

        # Coordinates
        left_x = 80
        node_x = 360
        top_y  = 120
        bot_y  = 420

        # VDD rail and resistor to node
        C.create_line(left_x, top_y, node_x, top_y, fill="#5af", width=2)  # VDD wire
        C.create_text(left_x - 30, top_y, text="VDD", fill="#5af", font=("Consolas", 12))

        # Resistor symbol
        rz = 14
        x = node_x - 80
        y = top_y
        for i in range(6):
            x1 = x + i*12
            C.create_line(x1, y - rz if i%2==0 else y + rz, x1+12, y + rz if i%2==0 else y - rz, fill="#aaa", width=2)
        C.create_line(node_x - 152, y, node_x - 80, y, fill="#aaa", width=2)  # left lead
        C.create_line(node_x, y, node_x - 8, y, fill="#aaa", width=2)        # right lead
        C.create_text(node_x - 120, y - 26, text=f"R = {r:,.1f} Ω", fill="#aaa", font=("Consolas", 10))

        # Node vertical to MOS drain
        C.create_line(node_x, top_y, node_x, bot_y - 120, fill="#aaa", width=2)
        C.create_oval(node_x-3, top_y-3, node_x+3, top_y+3, fill="#fff", outline="")

        # NMOS block
        box_w, box_h = 140, 100
        box_x, box_y = node_x - box_w/2, bot_y - box_h
        C.create_rectangle(box_x, box_y, box_x+box_w, box_y+box_h, outline="#7fd", width=2)
        C.create_text(box_x+box_w/2, box_y-16, text=f"NMOS  (region: {res['region']})",
                      fill="#7fd", font=("Consolas", 11))

        # Drain connection to node
        C.create_line(node_x, bot_y - 120, node_x, box_y, fill="#aaa", width=2)

        # Source to ground
        C.create_line(box_x+box_w/2, box_y+box_h, box_x+box_w/2, bot_y+18, fill="#aaa", width=2)
        C.create_line(box_x-40, bot_y+18, box_x+box_w/2, bot_y+18, fill="#aaa", width=2)
        # ground symbol
        gx = box_x-40; gy = bot_y+18
        C.create_line(gx-14, gy, gx+14, gy, fill="#888", width=2)
        C.create_line(gx-10, gy+6, gx+10, gy+6, fill="#888", width=2)
        C.create_line(gx-6,  gy+12, gx+6,  gy+12, fill="#888", width=2)
        C.create_text(gx-26, gy+2, text="GND", fill="#888", font=("Consolas", 10))

        # Gate drive
        gate_y = box_y + box_h/2
        C.create_line(box_x - 120, gate_y, box_x, gate_y, fill="#aaa", width=2)
        C.create_text(box_x - 140, gate_y, text=f"VG = {vg:.3f} V", fill="#aaa", font=("Consolas", 10))

        # Node readout
        C.create_text(node_x + 70, top_y + 10, text=f"Node Vd = {res['Vd']:.6f} V",
                      fill="#cfc", font=("Consolas", 11), anchor="w")
        C.create_text(node_x + 70, top_y + 28, text=f"Id = {res['Id']:.6e} A",
                      fill="#cfc", font=("Consolas", 11), anchor="w")
        C.create_text(node_x + 70, top_y + 46, text=f"VGS={vg:.3f} V, VDS={res['Vd']:.3f} V, Vt={vt:.3f} V, k={self.k.get():.3e} A/V²",
                      fill="#8ad", font=("Consolas", 10), anchor="w")

        # Wire (KCL) equations panel
        self.txt_wire.place(x=40, y=top_y+60)
        self.txt_wire.config(state="normal")
        self.txt_wire.delete("1.0", "end")
        self.txt_wire.insert("end", "KCL on the node (symbolic):\n")
        self.txt_wire.insert("end", res["kcl_sym"] + "\n\n")
        self.txt_wire.insert("end", "KCL with numeric substitution:\n")
        self.txt_wire.insert("end", res["kcl_num"] + "\n")
        self.txt_wire.config(state="disabled")

        # Device equation panel
        self.txt_dev.place(x=40, y=bot_y - 160)
        self.txt_dev.config(state="normal")
        self.txt_dev.delete("1.0", "end")
        self.txt_dev.insert("end", "MOSFET region equation:\n")
        self.txt_dev.insert("end", res["dev_sym"] + "\n\n")
        self.txt_dev.insert("end", "With numeric substitution:\n")
        self.txt_dev.insert("end", "Id = " + res["dev_num"] + "  A\n")
        self.txt_dev.config(state="disabled")


if __name__ == "__main__":
    App().mainloop()
