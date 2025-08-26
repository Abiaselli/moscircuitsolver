#!/usr/bin/env python3
# mosfet10.py
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import sympy as sp
from math import isfinite
import uuid, os

# ============================== Utils ==============================
def sym(name): return sp.Symbol(name, real=True)
def pretty(e): return sp.pretty(e)

# Smooth helpers (avoid region combinatorics)
SMOOTH_EPS = 1e-9
BETA_SWITCH = 6.0   # smoothness of triode↔sat blend

def relu_smooth(x, eps=SMOOTH_EPS):
    # 0.5*(x + sqrt(x^2 + eps)) ~ smooth ReLU
    return sp.Rational(1,2)*(x + sp.sqrt(x*x + eps))

def blend_smooth(z, beta=BETA_SWITCH):
    # 0..1 smooth step: 0.5*(1 + tanh(beta*z))
    return sp.Rational(1,2)*(1 + sp.tanh(beta*z))

# ========================== Scrollable sidebar ==========================
class ScrollableFrame(tk.Frame):
    def __init__(self, master, width=360, height=760, **kw):
        super().__init__(master, **kw)
        self.canvas = tk.Canvas(self, width=width, height=height, bg=kw.get("bg", "#171717"), highlightthickness=0)
        self.canvas.pack(side="left", fill="both", expand=False)
        self.vsb = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.vsb.pack(side="right", fill="y")
        self.inner = tk.Frame(self.canvas, bg=kw.get("bg", "#171717"))
        self.canvas_window = self.canvas.create_window((0,0), window=self.inner, anchor="nw")
        self.inner.bind("<Configure>", self._on_frame_configure)
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.canvas.bind("<Configure>", self._on_canvas_configure)
    def _on_frame_configure(self, _):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
    def _on_canvas_configure(self, event):
        self.canvas.itemconfigure(self.canvas_window, width=event.width)

# ============================== Devices ==============================
class Device:
    kind = "Device"
    def __init__(self, name): self.name=name; self.x,self.y=300,300
    @property
    def pins(self): return []
    def kcl_contribs(self, node_syms, node_to_sym): return {}
    def draw(self, C, nodes_pos): pass

class Resistor(Device):
    kind="R"
    def __init__(self, name, a="N1", b="N2", R=1e3):
        super().__init__(name); self.a,self.b=a,b; self.R=float(R)
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Va,Vb = node_to_sym[self.a], node_to_sym[self.b]
        Iab = (Va - Vb)/self.R
        return {self.a: Iab, self.b: -Iab}
    def draw(self, C, pos):
        ax,ay = pos.get(self.a,(self.x-40,self.y)); bx,by = pos.get(self.b,(self.x+40,self.y))
        C.create_line(ax,ay,bx,by, fill="#aaa", width=2)
        C.create_text((ax+bx)/2, (ay+by)/2 - 12, text=f"{self.name}={self.R:g}Ω", fill="#aaa", font=("Consolas",10))

class CurrentSourceDC(Device):
    kind="I"
    def __init__(self, name, a="N1", b="N2", I=1e-3):
        super().__init__(name); self.a,self.b=a,b; self.I=float(I)
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def kcl_contribs(self, node_syms, node_to_sym):
        return {self.a: sp.Float(self.I), self.b: -sp.Float(self.I)}
    def draw(self, C, pos):
        ax,ay = pos.get(self.a,(self.x-40,self.y)); bx,by = pos.get(self.b,(self.x+40,self.y))
        C.create_line(ax,ay,bx,by, fill="#cc8", width=2)
        C.create_text((ax+bx)/2,(ay+by)/2-12, text=f"{self.name}={self.I:g}A", fill="#cc8", font=("Consolas",10))

class Capacitor(Device):
    kind="C"
    def __init__(self, name, a="N1", b="N2", Cval=1e-12):
        super().__init__(name); self.a,self.b=a,b; self.C=float(Cval)
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def draw(self, C, pos):
        ax,ay = pos.get(self.a,(self.x-40,self.y)); bx,by = pos.get(self.b,(self.x+40,self.y))
        mx,my = (ax+bx)/2, (ay+by)/2
        C.create_line(ax,ay, mx-6,my, fill="#8ac", width=2)
        C.create_line(mx-6,my-12, mx-6,my+12, fill="#8ac", width=2)
        C.create_line(mx+6,my-12, mx+6,my+12, fill="#8ac", width=2)
        C.create_line(mx+6,my, bx,by, fill="#8ac", width=2)
        C.create_text(mx, my-18, text=f"{self.name}={self.C:g}F", fill="#8ac", font=("Consolas",10))

# Smooth SH NMOS/PMOS
def nmos_Id_smooth(Vgs, Vds, k, Vt, lam):
    Vov  = relu_smooth(Vgs - Vt)
    Id_lin = k*( Vov*Vds - sp.Rational(1,2)*Vds**2 )
    Id_lin = relu_smooth(Id_lin)
    Id_sat = sp.Rational(1,2)*k*Vov**2 * (1 + lam*Vds)
    s = blend_smooth(Vds - Vov)  # ~1 => sat, 0 => triode
    return (1 - s)*Id_lin + s*Id_sat

class NMOS(Device):
    kind="NMOS"
    def __init__(self, name, d="N1", g="N2", s="GND", k=2e-4, Vt=0.6, lam=0.0):
        super().__init__(name); self.d,self.g,self.s=d,g,s; self.k=float(k); self.Vt=float(Vt); self.lam=float(lam)
    @property
    def pins(self): return [("d", self.d), ("g", self.g), ("s", self.s)]
    def Id_expr(self, node_to_sym):
        Vd,Vg,Vs = node_to_sym[self.d], node_to_sym[self.g], node_to_sym[self.s]
        Vgs, Vds = Vg - Vs, Vd - Vs
        Id = nmos_Id_smooth(Vgs, Vds, sp.Float(self.k), sp.Float(self.Vt), sp.Float(self.lam))
        return Id       # convention: positive leaves D node
    def draw(self, C, pos):
        dx,dy = pos.get(self.d,(self.x,self.y-30)); gx,gy = pos.get(self.g,(self.x-40,self.y)); sx,sy = pos.get(self.s,(self.x,self.y+30))
        cx,cy = (dx+gx+sx)/3, (dy+gy+sy)/3; w,h=90,60
        C.create_rectangle(cx-w/2, cy-h/2, cx+w/2, cy+h/2, outline="#7fd", width=2)
        C.create_text(cx, cy-h/2-12, text=f"{self.name} (NMOS)", fill="#7fd", font=("Consolas",10))
        C.create_line(dx,dy, cx,cy-h/2, fill="#aaa", width=2); C.create_text(dx+6,dy-10,text="D",fill="#aaa",font=("Consolas",9),anchor="w")
        C.create_line(gx,gy, cx-w/2,cy,  fill="#aaa", width=2); C.create_text(gx+6,gy-10,text="G",fill="#aaa",font=("Consolas",9),anchor="w")
        C.create_line(sx,sy, cx,cy+h/2, fill="#aaa", width=2); C.create_text(sx+6,sy-10,text="S",fill="#aaa",font=("Consolas",9),anchor="w")

def pmos_Id_leavingD_smooth(Vsg, Vsd, k, Vt, lam):
    # magnitude model, current leaves drain towards source
    Vov  = relu_smooth(Vsg - Vt)
    Id_lin = k*( Vov*Vsd - sp.Rational(1,2)*Vsd**2 )
    Id_lin = relu_smooth(Id_lin)
    Id_sat = sp.Rational(1,2)*k*Vov**2 * (1 + lam*Vsd)
    s = blend_smooth(Vsd - Vov)
    Idmag = (1 - s)*Id_lin + s*Id_sat
    return Idmag

class PMOS(Device):
    kind="PMOS"
    def __init__(self, name, d="N1", g="N2", s="VDD", k=2e-4, Vt=0.6, lam=0.0):
        super().__init__(name); self.d,self.g,self.s=d,g,s; self.k=float(k); self.Vt=float(Vt); self.lam=float(lam)
    @property
    def pins(self): return [("d", self.d), ("g", self.g), ("s", self.s)]
    def Id_d_leaving(self, node_to_sym):
        Vd,Vg,Vs = node_to_sym[self.d], node_to_sym[self.g], node_to_sym[self.s]
        Vsg, Vsd = Vs - Vg, Vs - Vd
        Idmag = pmos_Id_leavingD_smooth(Vsg, Vsd, sp.Float(self.k), sp.Float(self.Vt), sp.Float(self.lam))
        return -Idmag   # negative means current leaving drain (toward node d) is -Idmag; we *add* to KCL at d
    def draw(self, C, pos):
        dx,dy = pos.get(self.d,(self.x,self.y-30)); gx,gy = pos.get(self.g,(self.x-40,self.y)); sx,sy = pos.get(self.s,(self.x,self.y+30))
        cx,cy = (dx+gx+sx)/3, (dy+gy+sy)/3; w,h=90,60
        C.create_rectangle(cx-w/2, cy-h/2, cx+w/2, cy+h/2, outline="#f8a", width=2)
        C.create_text(cx, cy-h/2-12, text=f"{self.name} (PMOS)", fill="#f8a", font=("Consolas",10))
        C.create_line(dx,dy, cx,cy-h/2, fill="#aaa", width=2); C.create_text(dx+6,dy-10,text="D",fill="#aaa",font=("Consolas",9),anchor="w")
        C.create_line(gx,gy, cx-w/2,cy,  fill="#aaa", width=2); C.create_text(gx+6,gy-10,text="G",fill="#aaa",font=("Consolas",9),anchor="w")
        C.create_line(sx,sy, cx,cy+h/2, fill="#aaa", width=2); C.create_text(sx+6,sy-10,text="S",fill="#aaa",font=("Consolas",9),anchor="w")

# Smooth comparator current source (no Piecewise)
class ComparatorToCurrent(Device):
    kind="CMP"
    def __init__(self, name, plus="N1", minus="N2", out="N3", Iout=1e-3, offset=0.0, gain=100.0):
        super().__init__(name); self.plus,self.minus,self.out=plus,minus,out
        self.Iout=float(Iout); self.offset=float(offset); self.gain=float(gain)
    @property
    def pins(self): return [("plus", self.plus), ("minus", self.minus), ("out", self.out)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Vp = node_to_sym[self.plus]; Vm = node_to_sym[self.minus]
        x  = self.gain * (Vp - Vm - sp.Float(self.offset))
        I  = sp.Float(self.Iout) * sp.tanh(x)
        return { self.out: I }
    def draw(self, C, pos):
        px,py = pos.get(self.plus,(self.x-30,self.y-10)); mx,my = pos.get(self.minus,(self.x-30,self.y+10)); ox,oy = pos.get(self.out,(self.x+50,self.y))
        C.create_polygon(self.x-24,self.y-20,self.x-24,self.y+20,self.x+12,self.y, outline="#ffa", fill="", width=2)
        C.create_text(self.x-10,self.y-28,text=f"{self.name}",fill="#ffa",font=("Consolas",10))
        C.create_line(px,py, self.x-24,self.y-10, fill="#ffa", width=2)
        C.create_line(mx,my, self.x-24,self.y+10, fill="#ffa", width=2)
        C.create_line(self.x+12,self.y, ox,oy, fill="#ffa", width=2)
        C.create_text((ox+self.x+12)/2, oy-12, text=f"±{self.Iout:g}A", fill="#ffa", font=("Consolas",9))

# (Kept for completeness—unused by default template)
class PolyCurrentToGND(Device):
    kind="POLY"
    def __init__(self, name, n="N1", g1=0.0, g2=0.0, I0=0.0):
        super().__init__(name); self.n=n; self.g1=float(g1); self.g2=float(g2); self.I0=float(I0)
    @property
    def pins(self): return [("n", self.n)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Vn = node_to_sym[self.n]
        return { self.n: sp.Float(self.g1)*Vn + sp.Float(self.g2)*Vn**2 + sp.Float(self.I0) }
    def draw(self, C, pos):
        x,y = pos.get(self.n,(self.x,self.y))
        C.create_oval(x-12,y-12,x+12,y+12, outline="#9cf", width=2)
        C.create_text(x,y-22, text=f"{self.name}\nI=g1·V+g2·V²+I0", fill="#9cf", font=("Consolas",9), anchor="s")

class CrossConductance(Device):
    kind="XG"
    def __init__(self, name, a="N1", ctrl="N2", g=0.0):
        super().__init__(name); self.a=a; self.ctrl=ctrl; self.g=float(g)
    @property
    def pins(self): return [("a", self.a), ("ctrl", self.ctrl)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Va = node_to_sym[self.a]; Vc = node_to_sym[self.ctrl]
        return { self.a: sp.Float(self.g)*Vc }
    def draw(self, C, pos):
        ax,ay = pos.get(self.a,(self.x,self.y)); cx,cy = pos.get(self.ctrl,(self.x+60,self.y))
        C.create_line(ax,ay,cx,cy, fill="#fc9", width=2, dash=(4,2))
        C.create_text((ax+cx)/2, (ay+cy)/2-12, text=f"{self.name}: I({self.a})=g·V({self.ctrl})", fill="#fc9", font=("Consolas",9))

# ============================== Circuit core ==============================
class Circuit:
    def __init__(self):
        self.nodes = {
            "VDD": {"x":140,"y":80,"fixed":True,"value":3.3},
            "GND": {"x":140,"y":520,"fixed":True,"value":0.0},
        }
        self.devices = []
        self.t = 0.0
        self.prev_node_voltages = {}

    # node helpers
    def ensure_node(self, name, x=220, y=200):
        if name not in self.nodes:
            self.nodes[name] = {"x":x,"y":y,"fixed":False}
    def _ensure_all_device_nodes(self):
        for d in self.devices:
            if hasattr(d,"pins"):
                for _,node in d.pins:
                    if node not in ("GND","VDD") and node not in self.nodes:
                        self.ensure_node(node, x=450, y=300)
    def node_syms(self):
        m={}
        for n,inf in self.nodes.items():
            m[n] = sp.Float(inf["value"]) if inf.get("fixed",False) else sym(f"V_{n}")
        return m
    def _snapshot_voltages(self, solved):
        self.prev_node_voltages = {n: float(solved.get(n,0.0))
                                   for n,inf in self.nodes.items() if not inf.get("fixed",False)}

    # DC solve (numeric, smooth models)
    def solve_dc(self):
        self._ensure_all_device_nodes()
        node_map = self.node_syms()
        unknown_nodes = [n for n,inf in self.nodes.items() if not inf.get("fixed",False)]
        if not unknown_nodes:
            return {}, {"kcl_sym":{}, "kcl_num":{}}, []

        # Build KCL with R / I / (optionally XG/POLY), exclude C and keep CMP out of DC
        kcl = {n: sp.Integer(0) for n in unknown_nodes}
        def add_dict(d):
            for node,cur in d.items():
                if node in kcl: kcl[node] += cur

        for dev in self.devices:
            if isinstance(dev, Resistor) or isinstance(dev, CurrentSourceDC):
                add_dict(dev.kcl_contribs(self.nodes, node_map))

        # (optional) include ideal math if you want DC influence:
        # for dev in self.devices:
        #     if isinstance(dev, (PolyCurrentToGND, CrossConductance)):
        #         add_dict(dev.kcl_contribs(self.nodes, node_map))

        for dev in self.devices:
            if isinstance(dev, NMOS):
                Id = dev.Id_expr(node_map)
                if dev.d in kcl: kcl[dev.d] += Id
                if dev.s in kcl: kcl[dev.s] -= Id
            elif isinstance(dev, PMOS):
                Id_d = dev.Id_d_leaving(node_map)
                if dev.d in kcl: kcl[dev.d] += Id_d
                if dev.s in kcl: kcl[dev.s] -= Id_d
            # Comparator skipped in DC

        eqs = [sp.Eq(kcl[n], 0) for n in unknown_nodes]
        unknown_syms = [node_map[n] for n in unknown_nodes]

        # initial guess from prev or 0
        guess = [self.prev_node_voltages.get(n, 0.0) for n in unknown_nodes]

        def _try_nsolve(g):
            root = sp.nsolve(eqs, unknown_syms, g, tol=1e-12, maxsteps=80, prec=60)
            nums = {}
            for i,n in enumerate(unknown_nodes):
                val = float(root[i])
                if not isfinite(val): val = 0.0
                nums[n] = val
            return nums

        try:
            num = _try_nsolve(guess)
        except Exception:
            # fallback guesses: small random-ish perturbations and mid-rail
            try:
                guess2 = [gv if isfinite(gv) else 0.0 for gv in guess]
                guess2 = [g + (0.05*(i%3-1)) for i,g in enumerate(guess2)]
                num = _try_nsolve(guess2)
            except Exception:
                try:
                    vdd = float(self.nodes["VDD"]["value"])
                    num = _try_nsolve([0.5*vdd for _ in unknown_nodes])
                except Exception:
                    # benign empty
                    return {}, {"kcl_sym":{n:pretty(eqs[i]) for i,n in enumerate(unknown_nodes)}, "kcl_num":{}}, eqs

        self._snapshot_voltages(num)
        meta = {
            "kcl_sym": {n: pretty(sp.Eq(kcl[n],0)) for n in unknown_nodes},
            "kcl_num": {}
        }
        return num, meta, eqs

    # Transient step (Backward Euler, smooth models)
    def solve_transient_step(self, dt):
        self._ensure_all_device_nodes()
        node_map = self.node_syms()
        unknown_nodes = [n for n,inf in self.nodes.items() if not inf.get("fixed",False)]
        if not unknown_nodes: return None
        unknown_syms = [node_map[n] for n in unknown_nodes]

        prev = {}
        for n,inf in self.nodes.items():
            prev[n] = float(inf["value"]) if inf.get("fixed",False) else float(self.prev_node_voltages.get(n,0.0))

        kcl = {n: sp.Integer(0) for n in unknown_nodes}
        def add_dict(d):
            for node,cur in d.items():
                if node in kcl: kcl[node] += cur

        # R / I / (ideal blocks if any)
        for dev in self.devices:
            if isinstance(dev, Resistor) or isinstance(dev, CurrentSourceDC) or isinstance(dev, CrossConductance) or isinstance(dev, PolyCurrentToGND) or isinstance(dev, ComparatorToCurrent):
                add_dict(dev.kcl_contribs(self.nodes, node_map))

        # Capacitors (Backward Euler)
        for dev in self.devices:
            if isinstance(dev, Capacitor):
                Va, Vb = node_map[dev.a], node_map[dev.b]
                Iab = sp.Float(dev.C) * ((Va - Vb) - (sp.Float(prev[dev.a]) - sp.Float(prev[dev.b]))) / sp.Float(dt)
                if dev.a in kcl: kcl[dev.a] += Iab
                if dev.b in kcl: kcl[dev.b] -= Iab

        # MOS
        for dev in self.devices:
            if isinstance(dev, NMOS):
                Id = dev.Id_expr(node_map)
                if dev.d in kcl: kcl[dev.d] += Id
                if dev.s in kcl: kcl[dev.s] -= Id
            elif isinstance(dev, PMOS):
                Id_d = dev.Id_d_leaving(node_map)
                if dev.d in kcl: kcl[dev.d] += Id_d
                if dev.s in kcl: kcl[dev.s] -= Id_d

        eqs = [sp.Eq(kcl[n], 0) for n in unknown_nodes]
        guess = [prev[n] for n in unknown_nodes]

        def _try_nsolve(g):
            root = sp.nsolve(eqs, unknown_syms, g, tol=1e-12, maxsteps=80, prec=60)
            nums = {}
            for i,n in enumerate(unknown_nodes):
                val = float(root[i])
                if not isfinite(val): val = prev[n]
                nums[n] = val
            return nums

        try:
            num = _try_nsolve(guess)
        except Exception:
            try:
                guess2 = [g + (0.02*(i%3-1)) for i,g in enumerate(guess)]
                num = _try_nsolve(guess2)
            except Exception:
                return None

        self.t += dt
        self._snapshot_voltages(num)
        meta = {"kcl_sym": {n: pretty(sp.Eq(kcl[n],0)) for n in unknown_nodes}}
        return num, meta, eqs

    # SPICE export (ideal blocks noted but skipped)
    def export_netlist(self, path):
        lines=[]
        lines.append("* Auto-generated netlist (mosfet10)")
        vdd=self.nodes["VDD"]["value"]
        lines.append(f"VDD VDD 0 {vdd}")
        model_cards=set()
        for d in self.devices:
            if isinstance(d,NMOS):
                m=f"MOD_{d.name}"; model_cards.add(f".model {m} NMOS (VTO={d.Vt} KP={d.k} LAMBDA={d.lam})")
                lines.append(f"M{d.name} {d.d} {d.g} {d.s} 0 {m} W=1u L=1u")
            elif isinstance(d,PMOS):
                m=f"MOD_{d.name}"; model_cards.add(f".model {m} PMOS (VTO={-d.Vt} KP={d.k} LAMBDA={d.lam})")
                lines.append(f"M{d.name} {d.d} {d.g} {d.s} 0 {m} W=1u L=1u")
        for d in self.devices:
            if isinstance(d,Resistor):      lines.append(f"R{d.name} {d.a} {d.b} {d.R}")
            elif isinstance(d,CurrentSourceDC): lines.append(f"I{d.name} {d.a} {d.b} {d.I}")
            elif isinstance(d,Capacitor):   lines.append(f"C{d.name} {d.a} {d.b} {d.C}")
            elif isinstance(d,PolyCurrentToGND): lines.append(f"* POLY {d.name} on {d.n} (not exported)")
            elif isinstance(d,CrossConductance): lines.append(f"* XG   {d.name}: I({d.a})=g·V({d.ctrl}) (not exported)")
            elif isinstance(d,ComparatorToCurrent): lines.append(f"* CMP  {d.name} tanh comparator (not exported)")
        lines.append(".op"); lines += list(model_cards); lines.append(".end")
        with open(path,"w") as f: f.write("\n".join(lines))

# ============================== Panels ==============================
class DraggablePanel(tk.Toplevel):
    def __init__(self, master, title, fg="#eee"):
        super().__init__(master); self.title(title); self.configure(bg="#151515")
        self.text = tk.Text(self, width=70, height=17, bg="#0d0d0d", fg=fg, insertbackground=fg,
                            font=("Consolas", 11), bd=0, highlightthickness=0)
        self.text.pack(fill="both", expand=True)
        self._ox=self._oy=0
        self.bind("<ButtonPress-1>", self._start); self.bind("<B1-Motion>", self._move)
    def _start(self,e): self._ox,self._oy=e.x,e.y
    def _move(self,e): self.geometry(f"+{self.winfo_x()+e.x-self._ox}+{self.winfo_y()+e.y-self._oy}")
    def set(self,s): self.text.config(state="normal"); self.text.delete("1.0","end"); self.text.insert("end",s); self.text.config(state="disabled")

# ============================== Editor ==============================
class Editor(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MOS Equation Lab — mosfet10")
        self.geometry("1280x800"); self.configure(bg="#101010")
        self.circ = Circuit()

        # transient ui vars
        self.running=False
        self.var_dt   = tk.DoubleVar(value=1e-4)
        self.var_tmax = tk.DoubleVar(value=5e-3)
        self.var_csv  = tk.StringVar(value="")

        self._build_left()
        self._build_canvas()
        self._build_panels()

        self.selected_node=None; self.selected_dev=None; self.dragging_node=None
        self.bind("<Delete>", lambda e: self.delete_selected())

        # Load Wijekoon–Dudek neuron on startup
        self.add_template_wijekoon2008()

    # ---------- UI build ----------
    def _build_left(self):
        sidebar = ScrollableFrame(self, width=360, height=760, bg="#171717")
        sidebar.pack(side="left", fill="y", padx=6, pady=6)
        p = sidebar.inner

        tk.Label(p, text="VDD (V)", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        self.var_vdd = tk.DoubleVar(value=self.circ.nodes["VDD"]["value"])
        tk.Spinbox(p, textvariable=self.var_vdd, from_=0.0, to=20.0, increment=0.01, width=12, font=("Consolas",11),
                   command=self.on_vdd).pack(anchor="w")
        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)

        tk.Label(p, text="Add:", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        tk.Button(p, text="+ Node", command=self.add_node).pack(fill="x", pady=2)
        tk.Button(p, text="+ Resistor", command=lambda:self.add_dev("R")).pack(fill="x", pady=2)
        tk.Button(p, text="+ NMOS", command=lambda:self.add_dev("NMOS")).pack(fill="x", pady=2)
        tk.Button(p, text="+ PMOS", command=lambda:self.add_dev("PMOS")).pack(fill="x", pady=2)
        tk.Button(p, text="+ CurrentSrc (DC)", command=lambda:self.add_dev("I")).pack(fill="x", pady=2)
        tk.Button(p, text="+ Capacitor", command=lambda:self.add_dev("C")).pack(fill="x", pady=2)
        tk.Button(p, text="+ Comparator (±I @ OUT)", command=lambda:self.add_dev("CMP")).pack(fill="x", pady=2)

        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        tk.Button(p, text="Solve DC", command=self.refresh).pack(fill="x", pady=3)
        tk.Button(p, text="Export netlist…", command=self.export_netlist).pack(fill="x", pady=3)

        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        tk.Button(p, text="Delete Selected", command=self.delete_selected).pack(fill="x", pady=3)
        tk.Button(p, text="Cleanup Unused Nodes", command=self.cleanup_unused_nodes).pack(fill="x", pady=3)

        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        tk.Label(p, text="Transient", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        r1=tk.Frame(p,bg="#171717"); r1.pack(fill="x", pady=2)
        tk.Label(r1, text="dt (s)", fg="#ddd", bg="#171717").pack(side="left")
        tk.Entry(r1, textvariable=self.var_dt, width=10).pack(side="left", padx=4)
        r2=tk.Frame(p,bg="#171717"); r2.pack(fill="x", pady=2)
        tk.Label(r2, text="t_max (s)", fg="#ddd", bg="#171717").pack(side="left")
        tk.Entry(r2, textvariable=self.var_tmax, width=10).pack(side="left", padx=4)
        r3=tk.Frame(p,bg="#171717"); r3.pack(fill="x", pady=2)
        tk.Label(r3, text="CSV path", fg="#ddd", bg="#171717").pack(side="left")
        tk.Entry(r3, textvariable=self.var_csv, width=18).pack(side="left", padx=4)
        tk.Button(p, text="Step dt", command=self.transient_step).pack(fill="x", pady=2)
        tk.Button(p, text="Run", command=self.transient_run).pack(fill="x", pady=2)
        tk.Button(p, text="Stop", command=self.transient_stop).pack(fill="x", pady=2)

        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        tk.Button(p, text="Template: Wijekoon–Dudek 2008", command=self.add_template_wijekoon2008).pack(fill="x", pady=4)
        tk.Button(p, text="Template: Izhikevich (ideal)", command=self.add_template_izhikevich).pack(fill="x", pady=4)
        tk.Button(p, text="Template: Izhikevich (NMOS-only v)", command=self.add_template_izhikevich_nmos).pack(fill="x", pady=4)

        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)
        tk.Label(p, text="Inspector", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        self.inspector = tk.Frame(p, bg="#171717"); self.inspector.pack(fill="x", pady=6)
        self.ins_text = tk.Text(self.inspector, width=30, height=18, bg="#0d0d0d", fg="#ccc",
                                insertbackground="#ccc", font=("Consolas", 10), bd=0, highlightthickness=0)
        self.ins_text.pack(fill="both", expand=True)

    def _build_canvas(self):
        self.canvas = tk.Canvas(self, width=900, height=780, bg="#0d0d0d", highlightthickness=0)
        self.canvas.pack(side="right", fill="both", expand=True, padx=6, pady=6)
        self.canvas.bind("<Button-1>", self.on_canvas_click)
        self.canvas.bind("<B1-Motion>", self.on_canvas_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_canvas_release)

    def _build_panels(self):
        self.panel_kcl = DraggablePanel(self, "Node equations (KCL)", fg="#cfe"); self.panel_kcl.geometry("+160+60")
        self.panel_dev = DraggablePanel(self, "Device notes", fg="#cfc"); self.panel_dev.geometry("+760+60")

    # ---------- Common ops ----------
    def clear_circuit(self):
        self.circ.nodes = {
            "VDD": {"x":140,"y":80,"fixed":True,"value":3.3},
            "GND": {"x":140,"y":520,"fixed":True,"value":0.0},
        }
        self.circ.devices = []
        self.circ.t = 0.0
        self.circ.prev_node_voltages = {}
        self.selected_node = None
        self.selected_dev = None
        self.refresh()

    def on_vdd(self):
        self.circ.nodes["VDD"]["value"] = float(self.var_vdd.get()); self.refresh()

    def add_node(self):
        name=f"N{len(self.circ.nodes)+1}"
        self.circ.ensure_node(name, x=220, y=200)
        self.refresh()

    def add_dev(self, kind):
        nid = uuid.uuid4().hex[:4].upper()
        if kind=="R":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2")
            self.circ.devices.append(Resistor(f"R{nid}","N1","N2",10000))
        elif kind=="NMOS":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2")
            self.circ.devices.append(NMOS(f"M{nid}","N2","N1","GND",2e-4,0.6,0.0))
        elif kind=="PMOS":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2")
            self.circ.devices.append(PMOS(f"MP{nid}","N2","N1","VDD",2e-4,0.6,0.0))
        elif kind=="I":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2")
            self.circ.devices.append(CurrentSourceDC(f"I{nid}","N1","N2",1e-3))
        elif kind=="C":
            self.circ.ensure_node("N1"); self.circ.devices.append(Capacitor(f"C{nid}","N1","GND",1e-9))
        elif kind=="CMP":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2"); self.circ.ensure_node("N3")
            self.circ.devices.append(ComparatorToCurrent(f"CMP{nid}","N1","N2","N3",1e-3,0.0,100.0))
        self.refresh()

    # ---------- Canvas interactions ----------
    def nodes_pos(self):
        pos={n:(inf["x"],inf["y"]) for n,inf in self.circ.nodes.items()}
        for d in self.circ.devices:
            if hasattr(d,"pins"):
                for _,node in d.pins:
                    if node not in pos and node not in ("GND","VDD"):
                        self.circ.ensure_node(node, x=450, y=300); pos[node]=(450,300)
        return pos
    def find_node_at(self,x,y,r=10):
        for n,inf in self.circ.nodes.items():
            nx,ny=inf["x"],inf["y"]
            if (x-nx)**2 + (y-ny)**2 <= r*r: return n
        return None
    def pick_device_near(self,x,y):
        pos=self.nodes_pos()
        best=None; mind=1e12
        for d in self.circ.devices:
            ps=[]
            if isinstance(d,(Resistor,CurrentSourceDC,Capacitor)): ps=[pos.get(d.a,(0,0)),pos.get(d.b,(0,0))]
            elif isinstance(d,(NMOS,PMOS)): ps=[pos.get(d.d,(0,0)),pos.get(d.g,(0,0)),pos.get(d.s,(0,0))]
            elif isinstance(d,PolyCurrentToGND): ps=[pos.get(d.n,(0,0))]
            elif isinstance(d,CrossConductance): ps=[pos.get(d.a,(0,0)),pos.get(d.ctrl,(0,0))]
            elif isinstance(d,ComparatorToCurrent): ps=[pos.get(d.plus,(0,0)),pos.get(d.minus,(0,0)),pos.get(d.out,(0,0))]
            if not ps: continue
            cx=sum(p[0] for p in ps)/len(ps); cy=sum(p[1] for p in ps)/len(ps)
            d2=(cx-x)**2+(cy-y)**2
            if d2<mind: mind, best=d2, d
        return best
    def on_canvas_click(self,e):
        n=self.find_node_at(e.x,e.y)
        if n:
            self.selected_node=n; self.selected_dev=None; self.dragging_node=n; self.update_inspector(); return
        self.selected_node=None; self.selected_dev=self.pick_device_near(e.x,e.y); self.dragging_node=None
        self.update_inspector()
    def on_canvas_drag(self,e):
        if self.dragging_node:
            inf=self.circ.nodes[self.dragging_node]; inf["x"],inf["y"]=e.x,e.y; self.refresh(redraw_only=True)
    def on_canvas_release(self,e): self.dragging_node=None

    # ---------- Render & panels ----------
    def refresh(self, redraw_only=False):
        C=self.canvas; C.delete("all")
        pos=self.nodes_pos()
        C.create_text(self.circ.nodes["VDD"]["x"]-30, self.circ.nodes["VDD"]["y"], text="VDD", fill="#58a6ff", font=("Consolas",12))
        for d in self.circ.devices: d.draw(C,pos)
        for name,inf in self.circ.nodes.items():
            x,y=inf["x"],inf["y"]; col="#58a6ff" if name=="VDD" else ("#888" if name=="GND" else "#fff")
            C.create_oval(x-5,y-5,x+5,y+5, fill=col, outline="")
            C.create_text(x+8,y-12, text=name, fill=col, font=("Consolas",10), anchor="w")
        if redraw_only: return

        res = self.circ.solve_dc()
        if not res:
            self.panel_kcl.set("No DC solution."); self.panel_dev.set("No DC solution."); return
        sol, meta, _ = res
        s1=["Node voltages:"]+[f"  {n} = {v:.9f} V" for n,v in sol.items()]
        s1.append("\nKCL (symbolic):")
        for n,txt in (meta.get("kcl_sym") or {}).items(): s1.append(f"[{n}]\n{txt}")
        self.panel_kcl.set("\n".join(s1) if len(s1)>1 else "No active equations.")
        self.update_inspector()

    def update_inspector(self):
        t=[]
        if self.selected_node:
            inf=self.circ.nodes[self.selected_node]; t.append(f"[Node] {self.selected_node}")
            t.append(f"  fixed: {inf.get('fixed',False)}")
            if inf.get("fixed",False): t.append(f"  value: {inf['value']} V")
        elif self.selected_dev:
            d=self.selected_dev; t.append(f"[{d.kind}] {d.name}")
            if isinstance(d,Resistor): t+= [f"  a:{d.a}", f"  b:{d.b}", f"  R:{d.R:g}Ω"]
            elif isinstance(d,CurrentSourceDC): t+= [f"  a:{d.a}", f"  b:{d.b}", f"  I:{d.I:g}A"]
            elif isinstance(d,Capacitor): t+= [f"  a:{d.a}", f"  b:{d.b}", f"  C:{d.C:g}F"]
            elif isinstance(d,NMOS): t+= [f"  d:{d.d} g:{d.g} s:{d.s}", f"  k:{d.k:g}  Vt:{d.Vt:g}  λ:{d.lam:g}"]
            elif isinstance(d,PMOS): t+= [f"  d:{d.d} g:{d.g} s:{d.s}", f"  k:{d.k:g}  |Vt|:{d.Vt:g}  λ:{d.lam:g}"]
            elif isinstance(d,ComparatorToCurrent): t+= [f"  +:{d.plus} -:{d.minus} out:{d.out}", f"  Iout:{d.Iout:g}  off:{d.offset:g}  gain:{d.gain:g}"]
            elif isinstance(d,PolyCurrentToGND): t+= [f"  node:{d.n} g1:{d.g1:g} g2:{d.g2:g} I0:{d.I0:g}"]
            elif isinstance(d,CrossConductance): t+= [f"  a:{d.a} ctrl:{d.ctrl} g:{d.g:g}"]
        else:
            t.append("Click a node to drag.\nClick a device to edit.")
        self.ins_text.config(state="normal"); self.ins_text.delete("1.0","end"); self.ins_text.insert("end","\n".join(t)); self.ins_text.config(state="disabled")
        # simple editors could be added here (omitted for brevity; you already have versions)

    # ---------- Delete & Cleanup ----------
    def delete_selected(self):
        if self.selected_dev:
            try: self.circ.devices.remove(self.selected_dev)
            except ValueError: pass
            self.selected_dev=None; self.refresh(); return
        if self.selected_node:
            n=self.selected_node; inf=self.circ.nodes.get(n,{})
            if inf.get("fixed",False): messagebox.showwarning("Delete node", f"'{n}' is fixed."); return
            used=[]
            for d in self.circ.devices:
                if hasattr(d,"pins"):
                    if n in [node for _,node in d.pins]: used.append(d.name)
            if used: messagebox.showwarning("Delete node", f"Node '{n}' used by: {', '.join(used)}"); return
            self.circ.nodes.pop(n,None); self.selected_node=None; self.refresh()
    def cleanup_unused_nodes(self):
        ref=set()
        for d in self.circ.devices:
            if hasattr(d,"pins"):
                for _,node in d.pins: ref.add(node)
        to_del=[n for n,inf in self.circ.nodes.items() if not inf.get("fixed",False) and n not in ref]
        for n in to_del: self.circ.nodes.pop(n,None)
        messagebox.showinfo("Cleanup", f"Removed {len(to_del)} unused node(s).")
        self.refresh()

    # ---------- Transient ----------
    def transient_step(self):
        dt=float(self.var_dt.get())
        res=self.circ.solve_transient_step(dt)
        if res is None:
            messagebox.showwarning("Transient","No solution this step."); return
        sol, meta, _ = res

        # Spike/reset (support both V/U and v/u)
        vkey = "V" if "V" in sol else ("v" if "v" in sol else None)
        ukey = "U" if "U" in sol else ("u" if "u" in sol else None)
        if hasattr(self,"izh") and vkey and ukey:
            Vth = self.izh.get("Vth", 1.70); Vc = self.izh.get("Vc", 0.60); D = self.izh.get("D", 0.05)
            if sol[vkey] >= Vth:
                self.circ.prev_node_voltages[vkey] = float(Vc)
                self.circ.prev_node_voltages[ukey] = float(sol[ukey] + D)
                sol[vkey] = float(Vc); sol[ukey] = float(sol[ukey] + D)

        path=self.var_csv.get().strip()
        if path:
            headers=["t"]+list(sol.keys()); write_header=not os.path.exists(path)
            with open(path,"a") as f:
                if write_header: f.write(",".join(headers)+"\n")
                f.write(",".join([f"{self.circ.t}"]+[f"{sol[n]}" for n in sol.keys()])+"\n")
        self.refresh()

    def transient_run(self):
        if self.running: return
        self.running=True
        tmax=float(self.var_tmax.get())
        def loop():
            if not self.running or self.circ.t>=tmax:
                self.running=False; return
            self.transient_step()
            self.after(1, loop)
        loop()
    def transient_stop(self): self.running=False

    # ---------- Export ----------
    def export_netlist(self):
        path=filedialog.asksaveasfilename(defaultextension=".cir", filetypes=[("Spice netlist",".cir"),("All files","*.*")])
        if not path: return
        try: self.circ.export_netlist(path)
        except Exception as e:
            messagebox.showerror("Export failed", str(e)); return
        messagebox.showinfo("Export", f"Wrote:\n{path}")

    # ---------- Templates ----------
    def add_template_wijekoon2008(self):
        """Wijekoon & Dudek 2008 compact silicon neuron; spike/reset handled in software."""
        self.clear_circuit()
        c=self.circ
        # nodes
        c.ensure_node("V",  x=420, y=220)
        c.ensure_node("U",  x=420, y=360)
        c.ensure_node("B",  x=260, y=200)
        c.ensure_node("VB", x=260, y=140)
        c.ensure_node("IN", x=200, y=240)
        # caps
        c.devices += [Capacitor("Cv","V","GND", 0.1e-12), Capacitor("Cu","U","GND", 1.0e-12)]
        # small bias divider
        c.devices += [Resistor("Rvb1","VDD","VB",200e3), Resistor("Rvb2","VB","GND",200e3)]
        # M1–M3 positive feedback/mirror (scaled for this model)
        k1,Vt1,lam1 = 2e-4, 0.5, 0.01
        c.devices += [NMOS("M1","B","V","GND", k1,Vt1,lam1),
                      NMOS("M2","B","B","GND", k1,Vt1,lam1),
                      NMOS("M3","V","B","GND", k1,Vt1,lam1)]
        # M4 leak controlled by U
        c.devices.append(NMOS("M4","V","U","GND", 2e-4, 0.5, 0.0))
        # M6 sink at U
        c.devices.append(NMOS("M6","U","U","GND", 2e-4, 0.5, 0.0))
        # M7 scaled mirror from V to U
        c.devices += [NMOS("M7ref","B","V","GND", 2e-5, 0.5, 0.01),
                      NMOS("M7","U","B","GND",    2e-5, 0.5, 0.01)]
        # input
        c.devices.append(Resistor("Rin","IN","V", 100e3))
        # tiny linear conductance to aid stability
        c.devices.append(NMOS("Ml","V","VB","GND", 1e-5, 0.3, 0.0))
        # reset params
        self.izh = {"Vth":1.70, "Vc":0.60, "D":0.05}
        c.prev_node_voltages.update({"V":0.3,"U":0.2})
        self.selected_node="V"
        self.panel_dev.set("Wijekoon–Dudek 2008 neuron loaded. Transient spike/reset active (V>=Vth → V:=Vc, U+=D).")
        self.refresh()

    def add_template_izhikevich(self):
        # Ideal math template (kept as option)
        self.clear_circuit()
        c=self.circ
        c.ensure_node("v", x=360,y=220); c.ensure_node("u", x=360,y=360); c.ensure_node("Iin", x=200,y=220)
        self.izh={"Vth":30.0,"Vc":-65.0,"D":8.0}  # map to v/u keys below
        c.devices += [Capacitor("Cv","v","GND",1.0), Capacitor("Cu","u","GND",1.0)]
        c.devices.append(PolyCurrentToGND("Pv","v",g1=5.0,g2=0.04,I0=140.0))
        c.devices.append(CrossConductance("Gvu","v","u",g=1.0))
        c.devices.append(CrossConductance("Guv","u","v",g=0.02*0.2))
        c.devices.append(CrossConductance("Guu","u","u",g=-0.02))
        c.prev_node_voltages.update({"v":-65.0,"u":-13.0})
        self.selected_node="v"; self.panel_dev.set("Ideal Izhikevich template ready."); self.refresh()

    def add_template_izhikevich_nmos(self):
        # NMOS-only dv/dt approximation (safe start settings)
        self.clear_circuit()
        c=self.circ
        c.ensure_node("v", x=360,y=220); c.ensure_node("u", x=360,y=360)
        c.ensure_node("VB", x=200,y=160); c.ensure_node("Xref", x=560,y=140)
        c.ensure_node("Xu", x=560,y=300); c.ensure_node("Iin", x=200,y=220)
        c.devices += [Capacitor("Cv","v","GND",1.0), Capacitor("Cu","u","GND",1.0)]
        c.devices += [Resistor("Rb1","VDD","VB",200e3), Resistor("Rb2","VB","GND",200e3)]
        c.devices.append(NMOS("MQ","v","v","GND",2e-4,-0.6,0.01))
        c.devices.append(NMOS("ML","v","VB","GND",2e-4,-0.2,0.0))
        c.devices += [Resistor("Rref","VDD","Xref",100e3), NMOS("Mref","Xref","Xref","GND",2e-4,-0.6,0.01),
                      NMOS("Mbias","v","Xref","GND",2e-4,-0.6,0.01)]
        c.devices += [Resistor("Ru","u","GND",100e3), NMOS("Mu_ref","Xu","Xu","GND",2e-4,-0.5,0.01),
                      NMOS("Mu_sink","v","Xu","GND",2e-4,-0.5,0.01)]
        c.devices.append(Resistor("Rin","Iin","v",100e3))
        c.prev_node_voltages.update({"v":-65.0,"u":-13.0})
        self.selected_node="v"; self.panel_dev.set("NMOS-only dv/dt: MQ≈v², ML≈g·v, Mbias const, Mu_sink ~ u."); self.refresh()

# ============================== Run ==============================
if __name__ == "__main__":
    Editor().mainloop()
