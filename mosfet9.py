#!/usr/bin/env python3
# mosfet8.py
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import sympy as sp
import itertools, uuid, os
from math import isfinite

# ---------- small utils ----------
def sym(name): return sp.Symbol(name, real=True)
def pretty(e): return sp.pretty(e)

# ---------- scrollable sidebar ----------
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

# ---------- device base ----------
class Device:
    kind = "Device"
    def __init__(self, name):
        self.name = name
        self.x, self.y = 300, 300
    @property
    def pins(self): return []
    def kcl_contribs(self, node_syms, node_to_sym): return {}
    def draw(self, C, nodes_pos): pass

# ---------- primitive devices ----------
class Resistor(Device):
    kind = "R"
    def __init__(self, name, a="N1", b="N2", R=1000.0):
        super().__init__(name); self.a, self.b, self.R = a, b, float(R)
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Va, Vb = node_to_sym[self.a], node_to_sym[self.b]
        Iab = (Va - Vb)/self.R
        return {self.a: Iab, self.b: -Iab}
    def draw(self, C, nodes_pos):
        ax, ay = nodes_pos.get(self.a, (self.x-40, self.y))
        bx, by = nodes_pos.get(self.b, (self.x+40, self.y))
        C.create_line(ax, ay, bx, by, fill="#aaa", width=2)
        C.create_text((ax+bx)/2, (ay+by)/2-12, text=f"{self.name}={self.R:g}Ω", fill="#aaa", font=("Consolas",10))

class CurrentSourceDC(Device):
    kind = "I"
    def __init__(self, name, a="N1", b="N2", I=1e-3):
        super().__init__(name); self.a, self.b, self.I = a, b, float(I)
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def kcl_contribs(self, node_syms, node_to_sym):
        return {self.a: sp.Float(self.I), self.b: -sp.Float(self.I)}
    def draw(self, C, nodes_pos):
        ax, ay = nodes_pos.get(self.a, (self.x-40, self.y))
        bx, by = nodes_pos.get(self.b, (self.x+40, self.y))
        C.create_line(ax, ay, bx, by, fill="#cc8", width=2)
        C.create_text((ax+bx)/2, (ay+by)/2-12, text=f"{self.name}={self.I:g}A", fill="#cc8", font=("Consolas",10))

class Capacitor(Device):
    kind = "C"
    def __init__(self, name, a="N1", b="N2", Cval=1e-9):
        super().__init__(name); self.a, self.b, self.C = a, b, float(Cval)
    @property
    def pins(self): return [("a", self.a), ("b", self.b)]
    def kcl_contribs(self, node_syms, node_to_sym):
        return {}  # DC open; transient handled in solver
    def draw(self, C, nodes_pos):
        ax, ay = nodes_pos.get(self.a, (self.x-40, self.y))
        bx, by = nodes_pos.get(self.b, (self.x+40, self.y))
        mx, my = (ax+bx)/2, (ay+by)/2
        C.create_line(ax, ay, mx-6, my, fill="#8ac", width=2)
        C.create_line(mx-6, my-12, mx-6, my+12, fill="#8ac", width=2)
        C.create_line(mx+6, my-12, mx+6, my+12, fill="#8ac", width=2)
        C.create_line(mx+6, my, bx, by, fill="#8ac", width=2)
        C.create_text(mx, my-18, text=f"{self.name}={self.C:g}F", fill="#8ac", font=("Consolas",10))

# ---------- MOSFETs ----------
def nmos_Id(Vgs, Vds, k, Vt, lam, region):
    Id_lin = k*((Vgs - Vt)*Vds - sp.Rational(1,2)*Vds**2)
    Id_sat = sp.Rational(1,2)*k*(Vgs - Vt)**2
    if lam and lam != 0: Id_sat = Id_sat*(1 + lam*Vds)
    if region == "cutoff": return sp.Integer(0)
    if region == "linear": return Id_lin
    if region == "sat":    return Id_sat
    raise ValueError("region")

class NMOS(Device):
    kind = "NMOS"
    def __init__(self, name, d="N1", g="N2", s="GND", k=2e-4, Vt=0.6, lam=0.0):
        super().__init__(name)
        self.d, self.g, self.s = d, g, s
        self.k, self.Vt, self.lam = float(k), float(Vt), float(lam)

    @property
    def pins(self): 
        return [("d", self.d), ("g", self.g), ("s", self.s)]

    def Id_expr(self, node_to_sym, region):
        Vd = node_to_sym[self.d]
        Vg = node_to_sym[self.g]
        Vs = node_to_sym[self.s]
        Vgs = Vg - Vs
        Vds = Vd - Vs
        Id  = nmos_Id(Vgs, Vds, sp.Float(self.k), sp.Float(self.Vt), sp.Float(self.lam), region)
        return Id, Vgs, Vds

    def draw(self, C, nodes_pos):
        dx,dy = nodes_pos.get(self.d, (self.x, self.y-30))
        gx,gy = nodes_pos.get(self.g, (self.x-40, self.y))
        sx,sy = nodes_pos.get(self.s, (self.x, self.y+30))
        cx,cy = (dx+gx+sx)/3, (dy+gy+sy)/3; w,h=90,60
        C.create_rectangle(cx-w/2, cy-h/2, cx+w/2, cy+h/2, outline="#7fd", width=2)
        C.create_text(cx, cy-h/2-12, text=f"{self.name} (NMOS)", fill="#7fd", font=("Consolas",10))
        C.create_line(dx,dy, cx,cy-h/2, fill="#aaa", width=2); C.create_oval(dx-3,dy-3,dx+3,dy+3,fill="#fff",outline="")
        C.create_text(dx+6,dy-10,text="D",fill="#aaa",font=("Consolas",9),anchor="w")
        C.create_line(gx,gy, cx- w/2,cy, fill="#aaa", width=2); C.create_oval(gx-3,gy-3,gx+3,gy+3,fill="#fff",outline="")
        C.create_text(gx+6,gy-10,text="G",fill="#aaa",font=("Consolas",9),anchor="w")
        C.create_line(sx,sy, cx,cy+h/2, fill="#aaa", width=2); C.create_oval(sx-3,sy-3,sx+3,sy+3,fill="#fff",outline="")
        C.create_text(sx+6,sy-10,text="S",fill="#aaa",font=("Consolas",9),anchor="w")

def pmos_Id_mag(Vsg, Vsd, k, Vt, lam, region):
    Id_lin = k*((Vsg - Vt)*Vsd - sp.Rational(1,2)*Vsd**2)
    Id_sat = sp.Rational(1,2)*k*(Vsg - Vt)**2
    if lam and lam != 0: Id_sat = Id_sat*(1 + lam*Vsd)
    if region == "cutoff": return sp.Integer(0)
    if region == "linear": return Id_lin
    if region == "sat":    return Id_sat
    raise ValueError("region")

class PMOS(Device):
    kind = "PMOS"
    def __init__(self, name, d="N1", g="N2", s="VDD", k=2e-4, Vt=0.6, lam=0.0):
        super().__init__(name); self.d,self.g,self.s=d,g,s; self.k=float(k); self.Vt=float(Vt); self.lam=float(lam)
    @property
    def pins(self): return [("d", self.d), ("g", self.g), ("s", self.s)]
    def Id_d_leaving(self, node_to_sym, region):
        Vd,Vg,Vs = node_to_sym[self.d], node_to_sym[self.g], node_to_sym[self.s]
        Vsg, Vsd = Vs-Vg, Vs-Vd
        Idmag = pmos_Id_mag(Vsg, Vsd, sp.Float(self.k), sp.Float(self.Vt), sp.Float(self.lam), region)
        return -Idmag, Vsg, Vsd  # leaving drain
    def draw(self, C, nodes_pos):
        dx,dy = nodes_pos.get(self.d, (self.x, self.y-30))
        gx,gy = nodes_pos.get(self.g, (self.x-40, self.y))
        sx,sy = nodes_pos.get(self.s, (self.x, self.y+30))
        cx,cy = (dx+gx+sx)/3, (dy+gy+sy)/3; w,h=90,60
        C.create_rectangle(cx-w/2, cy-h/2, cx+w/2, cy+h/2, outline="#f8a", width=2)
        C.create_text(cx, cy-h/2-12, text=f"{self.name} (PMOS)", fill="#f8a", font=("Consolas",10))
        C.create_line(dx,dy, cx,cy-h/2, fill="#aaa", width=2); C.create_oval(dx-3,dy-3,dx+3,dy+3,fill="#fff",outline="")
        C.create_text(dx+6,dy-10,text="D",fill="#aaa",font=("Consolas",9),anchor="w")
        C.create_line(gx,gy, cx- w/2,cy, fill="#aaa", width=2); C.create_oval(gx-3,gy-3,gx+3,gy+3,fill="#fff",outline="")
        C.create_text(gx+6,gy-10,text="G",fill="#aaa",font=("Consolas",9),anchor="w")
        C.create_line(sx,sy, cx,cy+h/2, fill="#aaa", width=2); C.create_oval(sx-3,sy-3,sx+3,sy+3,fill="#fff",outline="")
        C.create_text(sx+6,sy-10,text="S",fill="#aaa",font=("Consolas",9),anchor="w")

# ---------- idealized helpers ----------
class PolyCurrentToGND(Device):
    kind = "POLY"
    def __init__(self, name, n="N1", g1=0.0, g2=0.0, I0=0.0):
        super().__init__(name); self.n,self.g1,self.g2,self.I0 = n,float(g1),float(g2),float(I0)
    @property
    def pins(self): return [("n", self.n)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Vn = node_to_sym[self.n]
        return { self.n: sp.Float(self.g1)*Vn + sp.Float(self.g2)*Vn**2 + sp.Float(self.I0) }
    def draw(self, C, nodes_pos):
        x,y = nodes_pos.get(self.n, (self.x, self.y))
        C.create_oval(x-12,y-12,x+12,y+12, outline="#9cf", width=2)
        C.create_text(x, y-22, text=f"{self.name}\nI=g1·V+g2·V²+I0", fill="#9cf", font=("Consolas",9), anchor="s")

class CrossConductance(Device):
    kind = "XG"
    def __init__(self, name, a="N1", ctrl="N2", g=0.0):
        super().__init__(name); self.a,self.ctrl,self.g=a,ctrl,float(g)
    @property
    def pins(self): return [("a", self.a), ("ctrl", self.ctrl)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Va = node_to_sym[self.a]; Vc = node_to_sym[self.ctrl]
        return { self.a: sp.Float(self.g)*Vc }
    def draw(self, C, nodes_pos):
        ax,ay = nodes_pos.get(self.a,(self.x,self.y))
        cx,cy = nodes_pos.get(self.ctrl,(self.x+60,self.y))
        C.create_line(ax,ay,cx,cy, fill="#fc9", width=2, dash=(4,2))
        C.create_text((ax+cx)/2, (ay+cy)/2 - 12, text=f"{self.name}: I({self.a})=g·V({self.ctrl})", fill="#fc9", font=("Consolas",9))

class ComparatorToCurrent(Device):
    kind = "CMP"
    def __init__(self, name, plus="N1", minus="N2", out="N3", Iout=1e-3, offset=0.0):
        super().__init__(name); self.plus,self.minus,self.out=plus,minus,out; self.Iout=float(Iout); self.offset=float(offset)
    @property
    def pins(self): return [("plus", self.plus), ("minus", self.minus), ("out", self.out)]
    def kcl_contribs(self, node_syms, node_to_sym):
        Vp = node_to_sym[self.plus]; Vm = node_to_sym[self.minus]
        cond = sp.Ge(Vp - Vm - sp.Float(self.offset), 0)
        I = sp.Piecewise((sp.Float(self.Iout), cond), (-sp.Float(self.Iout), True))
        return { self.out: I }
    def draw(self, C, nodes_pos):
        px,py = nodes_pos.get(self.plus,(self.x-30,self.y-10))
        mx,my = nodes_pos.get(self.minus,(self.x-30,self.y+10))
        ox,oy = nodes_pos.get(self.out,(self.x+50,self.y))
        C.create_polygon(self.x-24,self.y-20,self.x-24,self.y+20,self.x+12,self.y, outline="#ffa", fill="", width=2)
        C.create_text(self.x-10,self.y-28,text=f"{self.name}",fill="#ffa",font=("Consolas",10))
        C.create_line(px,py, self.x-24,self.y-10, fill="#ffa", width=2)
        C.create_line(mx,my, self.x-24,self.y+10, fill="#ffa", width=2)
        C.create_line(self.x+12,self.y, ox,oy, fill="#ffa", width=2)
        C.create_text((ox+self.x+12)/2, oy-12, text=f"±{self.Iout:g}A", fill="#ffa", font=("Consolas",9))

# ---------- circuit core ----------
class Circuit:
    def __init__(self):
        self.nodes = {
            "VDD": {"x":140,"y":80,"fixed":True,"value":3.3},
            "GND": {"x":140,"y":520,"fixed":True,"value":0.0},
            "X":   {"x":360,"y":260,"fixed":False},
            "Y":   {"x":560,"y":260,"fixed":False},
        }
        self.devices = [
            Resistor("Rb","VDD","X",10000.0),
            Resistor("Rload","VDD","Y",20000.0),
            NMOS("M1","X","X","GND",2e-4,0.6,0.0),
            NMOS("M2","Y","X","GND",2e-4,0.6,0.0),
        ]
        self.t = 0.0
        self.prev_node_voltages = {}

    def ensure_node(self, name, x=220, y=200):
        if name not in self.nodes:
            self.nodes[name] = {"x":x,"y":y,"fixed":False}

    def _ensure_all_device_nodes(self):
        for d in self.devices:
            if hasattr(d,"pins"):
                for _,node in d.pins:
                    if node not in ("GND","VDD") and node not in self.nodes:
                        self.ensure_node(node)

    def node_syms(self):
        m={}
        for n,inf in self.nodes.items():
            m[n] = sp.Float(inf["value"]) if inf.get("fixed",False) else sym(f"V_{n}")
        return m

    def _snapshot_voltages(self, solved):
        self.prev_node_voltages = {n: float(solved.get(n,0.0))
                                   for n,inf in self.nodes.items() if not inf.get("fixed",False)}

    def solve_dc(self):
        self._ensure_all_device_nodes()
        node_map = self.node_syms()
        unknown_nodes = [n for n,inf in self.nodes.items() if not inf.get("fixed",False)]
        unknown_syms  = [node_map[n] for n in unknown_nodes]

        def empty_kcl(): return {n: sp.Integer(0) for n in unknown_nodes}
        base_kcl = empty_kcl()
        any_contrib = {n: False for n in unknown_nodes}

        for d in self.devices:
            if isinstance(d,(Resistor,CurrentSourceDC)):
                contrib = d.kcl_contribs(self.nodes, node_map)
                for node,cur in contrib.items():
                    if node in base_kcl:
                        base_kcl[node]+=cur; any_contrib[node]=True

        mos_list = [d for d in self.devices if isinstance(d,(NMOS,PMOS))]
        region_opts=["sat","linear","cutoff"]

        def _empty_result(eqs):
            return {}, {"regions":{}, "kcl_sym":{n:pretty(sp.Eq(base_kcl[n],0)) for n in unknown_nodes}, "kcl_num":{}}, eqs

        if not mos_list and not any(any_contrib.values()):
            eqs=[sp.Eq(base_kcl[n],0) for n in unknown_nodes]
            return _empty_result(eqs)

        if not mos_list:
            eqs=[sp.Eq(base_kcl[n],0) for n in unknown_nodes]
            try: sols = sp.solve(eqs, unknown_syms, dict=True) or []
            except Exception: sols=[]
            if not sols: return _empty_result(eqs)
            s=sols[0]; num={}
            for n in unknown_nodes:
                symn=node_map[n]
                num[n]= float(s[symn]) if symn in s else float(self.prev_node_voltages.get(n,0.0))
            self._snapshot_voltages(num)
            meta={"regions":{}, "kcl_sym":{n:pretty(sp.Eq(base_kcl[n],0)) for n in unknown_nodes}, "kcl_num":{}}
            return num, meta, eqs

        for choice in itertools.product(region_opts, repeat=len(mos_list)):
            kcl={k:v for k,v in base_kcl.items()}
            preds=[]
            for dev,rgn in zip(mos_list, choice):
                if isinstance(dev,NMOS):
                    Id,Vgs,Vds = dev.Id_expr(node_map, rgn)
                    if dev.d in kcl: kcl[dev.d]+=Id; any_contrib[dev.d]=True
                    if dev.s in kcl: kcl[dev.s]-=Id; any_contrib[dev.s]=True
                    vt=dev.Vt
                    preds += [sp.Le(Vgs,vt)] if rgn=="cutoff" else (
                             [sp.Ge(Vgs,vt), sp.Le(Vds, Vgs-vt)] if rgn=="linear" else
                             [sp.Ge(Vgs,vt), sp.Ge(Vds, Vgs-vt)])
                else:
                    Id_d,Vsg,Vsd = dev.Id_d_leaving(node_map, rgn)
                    if dev.d in kcl: kcl[dev.d]+=Id_d; any_contrib[dev.d]=True
                    if dev.s in kcl: kcl[dev.s]-=Id_d; any_contrib[dev.s]=True
                    vt=dev.Vt
                    preds += [sp.Le(Vsg,vt)] if rgn=="cutoff" else (
                             [sp.Ge(Vsg,vt), sp.Le(Vsd, Vsg-vt)] if rgn=="linear" else
                             [sp.Ge(Vsg,vt), sp.Ge(Vsd, Vsg-vt)])

            eqs=[sp.Eq(kcl[n],0) for n in unknown_nodes]
            try: sols= sp.solve(eqs, unknown_syms, dict=True) or []
            except Exception: sols=[]
            for s in sols:
                num={}
                ok=True
                for n in unknown_nodes:
                    symn=node_map[n]
                    if symn in s:
                        try: num[n]=float(s[symn])
                        except Exception: ok=False; break
                    else:
                        num[n]=float(self.prev_node_voltages.get(n,0.0))
                if not ok or not all(isfinite(v) for v in num.values()): continue
                subs={node_map[k]:v for k,v in num.items()}
                for name,inf in self.nodes.items():
                    if inf.get("fixed",False): subs[node_map.get(name,sym("dummy"))]=float(inf["value"])
                good=True
                for p in preds:
                    try:
                        if not bool(p.subs(subs)): good=False; break
                    except Exception:
                        good=False; break
                if not good: continue
                self._snapshot_voltages(num)
                kcl_sym={n:pretty(sp.Eq(sp.simplify(kcl[n]),0)) for n in unknown_nodes}
                kcl_num={n:pretty(sp.Eq(sp.nsimplify(kcl[n].subs(subs),rational=False),0)) for n in unknown_nodes}
                meta={"regions":dict(zip([m.name for m in mos_list], choice)), "kcl_sym":kcl_sym, "kcl_num":kcl_num}
                return num, meta, eqs

        return _empty_result([sp.Eq(base_kcl[n],0) for n in unknown_nodes])

    def solve_transient_step(self, dt):
        self._ensure_all_device_nodes()
        node_map = self.node_syms()
        unknown_nodes=[n for n,inf in self.nodes.items() if not inf.get("fixed",False)]
        unknown_syms=[node_map[n] for n in unknown_nodes]
        prev={}
        for n,inf in self.nodes.items():
            prev[n]= float(inf["value"]) if inf.get("fixed",False) else float(self.prev_node_voltages.get(n,0.0))

        def empty_kcl(): return {n: sp.Integer(0) for n in unknown_nodes}
        base_kcl=empty_kcl()
        for d in self.devices:
            if isinstance(d,(Resistor,CurrentSourceDC)):
                contrib=d.kcl_contribs(self.nodes, node_map)
                for node,cur in contrib.items():
                    if node in base_kcl: base_kcl[node]+=cur
        for d in self.devices:
            if isinstance(d,Capacitor):
                a,b,Cv=d.a,d.b,d.C
                Va,Vb=node_map[a],node_map[b]
                Iab = sp.Float(Cv) * ((Va - Vb) - (sp.Float(prev[a]) - sp.Float(prev[b])))/sp.Float(dt)
                if a in base_kcl: base_kcl[a]+=Iab
                if b in base_kcl: base_kcl[b]-=Iab
        for d in self.devices:
            if isinstance(d,PolyCurrentToGND) or isinstance(d,CrossConductance) or isinstance(d,ComparatorToCurrent):
                contrib=d.kcl_contribs(self.nodes, node_map)
                for node,cur in contrib.items():
                    if node in base_kcl: base_kcl[node]+=cur

        mos_list=[d for d in self.devices if isinstance(d,(NMOS,PMOS))]
        region_opts=["sat","linear","cutoff"] if mos_list else []

        if not mos_list:
            eqs=[sp.Eq(base_kcl[n],0) for n in unknown_nodes]
            try: sols= sp.solve(eqs, unknown_syms, dict=True) or []
            except Exception: sols=[]
            if not sols: return None
            s=sols[0]
            num={}
            for n in unknown_nodes:
                symn=node_map[n]
                num[n]= float(s[symn]) if symn in s else float(prev[n])
            self.t += dt
            self._snapshot_voltages(num)
            meta={"kcl_sym":{n:pretty(sp.Eq(sp.simplify(base_kcl[n]),0)) for n in unknown_nodes}}
            return num, meta, eqs

        for choice in itertools.product(region_opts, repeat=len(mos_list)):
            kcl={k:v for k,v in base_kcl.items()}
            preds=[]
            for dev,rgn in zip(mos_list,choice):
                if isinstance(dev,NMOS):
                    Id,Vgs,Vds=dev.Id_expr(node_map,rgn)
                    if dev.d in kcl: kcl[dev.d]+=Id
                    if dev.s in kcl: kcl[dev.s]-=Id
                    vt=dev.Vt
                    preds += [sp.Le(Vgs,vt)] if rgn=="cutoff" else (
                             [sp.Ge(Vgs,vt), sp.Le(Vds, Vgs-vt)] if rgn=="linear" else
                             [sp.Ge(Vgs,vt), sp.Ge(Vds, Vgs-vt)])
                else:
                    Id_d,Vsg,Vsd=dev.Id_d_leaving(node_map,rgn)
                    if dev.d in kcl: kcl[dev.d]+=Id_d
                    if dev.s in kcl: kcl[dev.s]-=Id_d
                    vt=dev.Vt
                    preds += [sp.Le(Vsg,vt)] if rgn=="cutoff" else (
                             [sp.Ge(Vsg,vt), sp.Le(Vsd, Vsg-vt)] if rgn=="linear" else
                             [sp.Ge(Vsg,vt), sp.Ge(Vsd, Vsg-vt)])
            eqs=[sp.Eq(kcl[n],0) for n in unknown_nodes]
            try: sols= sp.solve(eqs, unknown_syms, dict=True) or []
            except Exception: sols=[]
            for s in sols:
                num={}
                for n in unknown_nodes:
                    symn=node_map[n]
                    num[n]= float(s[symn]) if symn in s else float(prev[n])
                if not all(isfinite(v) for v in num.values()): continue
                subs={node_map[k]:v for k,v in num.items()}
                for name,inf in self.nodes.items():
                    if inf.get("fixed",False): subs[node_map.get(name,sym("dummy"))]=float(inf["value"])
                good=True
                for p in preds:
                    try:
                        if not bool(p.subs(subs)): good=False; break
                    except Exception:
                        good=False; break
                if not good: continue
                self.t += dt
                self._snapshot_voltages(num)
                meta={"regions": dict(zip([m.name for m in mos_list], choice)),
                      "kcl_sym": {n: pretty(sp.Eq(sp.simplify(kcl[n]),0)) for n in unknown_nodes}}
                return num, meta, eqs
        return None

    def export_netlist(self, path):
        lines=[]
        lines.append("* Auto-generated netlist (mosfet8)")
        vdd=self.nodes["VDD"]["value"]
        lines.append(f"VDD VDD 0 {vdd}")
        model_cards=[]
        for d in self.devices:
            if isinstance(d,NMOS):
                mname=f"MOD_{d.name}"
                model_cards.append(f".model {mname} NMOS (VTO={d.Vt} KP={d.k} LAMBDA={d.lam})")
                lines.append(f"M{d.name} {d.d} {d.g} {d.s} 0 {mname} W=1u L=1u")
            elif isinstance(d,PMOS):
                mname=f"MOD_{d.name}"
                model_cards.append(f".model {mname} PMOS (VTO={-d.Vt} KP={d.k} LAMBDA={d.lam})")
                lines.append(f"M{d.name} {d.d} {d.g} {d.s} 0 {mname} W=1u L=1u")
        for d in self.devices:
            if isinstance(d,Resistor):
                lines.append(f"R{d.name} {d.a} {d.b} {d.R}")
            elif isinstance(d,CurrentSourceDC):
                lines.append(f"I{d.name} {d.a} {d.b} {d.I}")
            elif isinstance(d,Capacitor):
                lines.append(f"C{d.name} {d.a} {d.b} {d.C}")
            elif isinstance(d,PolyCurrentToGND):
                lines.append(f"* POLY {d.name} on {d.n}: I = {d.g1}*V + {d.g2}*V^2 + {d.I0}  (not exported)")
            elif isinstance(d,CrossConductance):
                lines.append(f"* XG {d.name} : I({d.a}) = {d.g} * V({d.ctrl})  (not exported)")
            elif isinstance(d,ComparatorToCurrent):
                lines.append(f"* CMP {d.name} : OUT ±{d.Iout} A vs (plus-minus-{d.offset})  (not exported)")
        lines.append(".op"); lines += model_cards; lines.append(".end")
        with open(path,"w") as f: f.write("\n".join(lines))

# ---------- Panels ----------
class DraggablePanel(tk.Toplevel):
    def __init__(self, master, title, fg="#eee"):
        super().__init__(master); self.title(title); self.configure(bg="#151515")
        self.text = tk.Text(self, width=70, height=16, bg="#0d0d0d", fg=fg,
                            insertbackground=fg, font=("Consolas", 11), bd=0, highlightthickness=0)
        self.text.pack(fill="both", expand=True)
        self._ox=self._oy=0
        self.bind("<ButtonPress-1>", self._start); self.bind("<B1-Motion>", self._move)
    def _start(self,e): self._ox,self._oy=e.x,e.y
    def _move(self,e): self.geometry(f"+{self.winfo_x()+e.x-self._ox}+{self.winfo_y()+e.y-self._oy}")
    def set(self,s): self.text.config(state="normal"); self.text.delete("1.0","end"); self.text.insert("end",s); self.text.config(state="disabled")

# ---------- Editor ----------
class Editor(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MOS Equation Lab — mosfet8")
        self.geometry("1280x800"); self.configure(bg="#101010")
        self.circ = Circuit()

        # transient vars BEFORE building UI
        self.running=False
        self.var_dt   = tk.DoubleVar(value=1e-3)
        self.var_tmax = tk.DoubleVar(value=0.5)
        self.var_csv  = tk.StringVar(value="")

        self._build_left()
        self._build_canvas()
        self._build_panels()

        self.selected_node=None; self.selected_dev=None; self.dragging_node=None
        self.bind("<Delete>", lambda e: self.delete_selected())
        self.refresh()

    # ----- UI build
    def _build_left(self):
        sidebar = ScrollableFrame(self, width=360, height=760, bg="#171717")
        sidebar.pack(side="left", fill="y", padx=6, pady=6)
        p = sidebar.inner

        tk.Label(p, text="VDD (V)", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        self.var_vdd = tk.DoubleVar(value=self.circ.nodes["VDD"]["value"])
        tk.Spinbox(p, textvariable=self.var_vdd, from_=0.0, to=20.0, increment=0.01,
                   width=12, font=("Consolas",11), command=self.on_vdd).pack(anchor="w")
        ttk.Separator(p, orient="horizontal").pack(fill="x", pady=6)

        tk.Label(p, text="Add:", fg="#ddd", bg="#171717", font=("Consolas", 11)).pack(anchor="w")
        tk.Button(p, text="+ Node", command=self.add_node).pack(fill="x", pady=2)
        tk.Button(p, text="+ Resistor", command=lambda:self.add_dev("R")).pack(fill="x", pady=2)
        tk.Button(p, text="+ NMOS", command=lambda:self.add_dev("NMOS")).pack(fill="x", pady=2)
        tk.Button(p, text="+ Depletion NMOS", command=lambda:self.add_dev("NMOS_DEP")).pack(fill="x", pady=2)
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
        tk.Button(p, text="Template: Izhikevich (ideal)", command=self.add_template_izhikevich).pack(fill="x", pady=4)
        tk.Button(p, text="Template: Izhikevich (NMOS-only v)", command=self.add_template_izhikevich_nmos).pack(fill="x", pady=4)
        tk.Button(p, text="Template: u via PMOS mirror", command=self.add_template_u_with_PMOS_mirror).pack(fill="x", pady=4)

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

    # ----- palette
    def on_vdd(self):
        self.circ.nodes["VDD"]["value"]=float(self.var_vdd.get()); self.refresh()

    def add_node(self):
        name=f"N{len(self.circ.nodes)+1}"; self.circ.ensure_node(name, x=220, y=200); self.refresh()

    def add_dev(self, kind):
        nid = uuid.uuid4().hex[:4].upper()
        if kind=="R":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2")
            self.circ.devices.append(Resistor(f"R{nid}","N1","N2",10000))
        elif kind=="NMOS":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2")
            self.circ.devices.append(NMOS(f"M{nid}","N2","N1","GND",2e-4,0.6,0.0))
        elif kind=="NMOS_DEP":
            self.circ.ensure_node("N1"); self.circ.ensure_node("N2")
            self.circ.devices.append(NMOS(f"MD{nid}","N2","N1","GND",2e-4,-0.6,0.0))
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
            self.circ.devices.append(ComparatorToCurrent(f"CMP{nid}","N1","N2","N3",1e-3,0.0))
        self.refresh()

    # ----- canvas interactions
    def nodes_pos(self):
        pos={n:(inf["x"],inf["y"]) for n,inf in self.circ.nodes.items()}
        for d in self.circ.devices:
            if hasattr(d,"pins"):
                for _,node in d.pins:
                    if node not in pos and node not in ("GND","VDD"):
                        self.circ.ensure_node(node, x=450, y=300)
                        pos[node]=(450,300)
        return pos

    def find_node_at(self,x,y,r=10):
        for n,inf in self.circ.nodes.items():
            nx,ny=inf["x"],inf["y"]
            if (x-nx)**2 + (y-ny)**2 <= r*r: return n
        return None

    def pick_device_near(self, x, y):
        nodes_pos=self.nodes_pos()
        best=None; mind=1e12
        for d in self.circ.devices:
            ps=[]
            if isinstance(d,(Resistor,CurrentSourceDC,Capacitor)):
                ps=[nodes_pos.get(d.a,(0,0)), nodes_pos.get(d.b,(0,0))]
            elif isinstance(d,(NMOS,PMOS)):
                ps=[nodes_pos.get(d.d,(0,0)), nodes_pos.get(d.g,(0,0)), nodes_pos.get(d.s,(0,0))]
            elif isinstance(d,(PolyCurrentToGND,)):
                ps=[nodes_pos.get(d.n,(0,0))]
            elif isinstance(d,(CrossConductance,)):
                ps=[nodes_pos.get(d.a,(0,0)), nodes_pos.get(d.ctrl,(0,0))]
            elif isinstance(d,(ComparatorToCurrent,)):
                ps=[nodes_pos.get(d.plus,(0,0)), nodes_pos.get(d.minus,(0,0)), nodes_pos.get(d.out,(0,0))]
            if not ps: continue
            cx=sum(p[0] for p in ps)/len(ps); cy=sum(p[1] for p in ps)/len(ps)
            dist=(cx-x)**2 + (cy-y)**2
            if dist<mind: mind, best=dist, d
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

    # ----- render & panels
    def refresh(self, redraw_only=False):
        C=self.canvas; C.delete("all")
        pos=self.nodes_pos()
        C.create_text(self.circ.nodes["VDD"]["x"]-30, self.circ.nodes["VDD"]["y"], text="VDD", fill="#58a6ff", font=("Consolas",12))
        for d in self.circ.devices: d.draw(C, pos)
        for name,inf in self.circ.nodes.items():
            x,y=inf["x"],inf["y"]; color="#58a6ff" if name=="VDD" else ("#888" if name=="GND" else "#fff")
            C.create_oval(x-5,y-5,x+5,y+5, fill=color, outline=""); C.create_text(x+8,y-12,text=name, fill=color, font=("Consolas",10), anchor="w")
        if redraw_only: return

        res = self.circ.solve_dc()
        if not res:
            self.panel_kcl.set("No DC solution."); self.panel_dev.set("No DC solution.")
            return
        sol, meta, eqs = res
        s1=["Node voltages:"]+[f"  {n} = {v:.9f} V" for n,v in sol.items()]
        s1.append("\nKCL (symbolic):")
        for n,txt in (meta.get("kcl_sym") or {}).items(): s1.append(f"[{n}]\n{txt}")
        s1.append("\nKCL (numeric):")
        for n,txt in (meta.get("kcl_num") or {}).items(): s1.append(f"[{n}]\n{txt}")
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
            if isinstance(d,Resistor):
                t+= [f"  a:{d.a}", f"  b:{d.b}", f"  R:{d.R:g}Ω"]
            elif isinstance(d,CurrentSourceDC):
                t+= [f"  a:{d.a}", f"  b:{d.b}", f"  I:{d.I:g}A"]
            elif isinstance(d,Capacitor):
                t+= [f"  a:{d.a}", f"  b:{d.b}", f"  C:{d.C:g}F"]
            elif isinstance(d,NMOS):
                t+= [f"  d:{d.d} g:{d.g} s:{d.s}", f"  k:{d.k:g} A/V²  Vt:{d.Vt:g} V  λ:{d.lam:g} 1/V"]
            elif isinstance(d,PMOS):
                t+= [f"  d:{d.d} g:{d.g} s:{d.s}", f"  k:{d.k:g} A/V²  |Vt|:{d.Vt:g} V  λ:{d.lam:g} 1/V"]
            elif isinstance(d,PolyCurrentToGND):
                t+= [f"  node:{d.n}", f"  g1:{d.g1:g}", f"  g2:{d.g2:g}", f"  I0:{d.I0:g} A"]
            elif isinstance(d,CrossConductance):
                t+= [f"  a:{d.a}", f"  ctrl:{d.ctrl}", f"  g:{d.g:g} A/V"]
            elif isinstance(d,ComparatorToCurrent):
                t+= [f"  +:{d.plus}", f"  -:{d.minus}", f"  out:{d.out}", f"  Iout:{d.Iout:g} A", f"  offset:{d.offset:g} V"]
        else:
            t.append("Click a node to drag.\nClick a device to edit.")
        self.ins_text.config(state="normal"); self.ins_text.delete("1.0","end"); self.ins_text.insert("end","\n".join(t)); self.ins_text.config(state="disabled")

        # small property editors
        for w in self.inspector.pack_slaves():
            if isinstance(w, tk.Frame) and getattr(w, "_is_editor", False): w.destroy()
        if self.selected_dev:
            ed=tk.Frame(self.inspector,bg="#171717"); ed._is_editor=True; ed.pack(fill="x", pady=4)
            nodes=list(self.circ.nodes.keys())
            def add_dd(lbl,get,setter):
                row=tk.Frame(ed,bg="#171717"); row.pack(fill="x", pady=1)
                tk.Label(row,text=lbl,fg="#ddd",bg="#171717",font=("Consolas",10)).pack(side="left")
                var=tk.StringVar(value=get()); om=ttk.OptionMenu(row,var,get(),*nodes,command=lambda v:setter(v) or self.refresh())
                om.pack(side="left")
            def add_num(lbl,get,setter):
                row=tk.Frame(ed,bg="#171717"); row.pack(fill="x", pady=1)
                tk.Label(row,text=lbl,fg="#ddd",bg="#171717",font=("Consolas",10)).pack(side="left")
                var=tk.DoubleVar(value=get()); sb=tk.Spinbox(row,textvariable=var,from_=-1e6,to=1e6,increment=0.001,
                                                             width=12,font=("Consolas",10),
                                                             command=lambda:setter(var.get()) or self.refresh())
                sb.pack(side="left")
            d=self.selected_dev
            if isinstance(d,Resistor):
                add_dd("a", lambda:d.a, lambda v:setattr(d,"a",v))
                add_dd("b", lambda:d.b, lambda v:setattr(d,"b",v))
                add_num("R (Ω)", lambda:d.R, lambda v:setattr(d,"R",float(v)))
            elif isinstance(d,CurrentSourceDC):
                add_dd("a", lambda:d.a, lambda v:setattr(d,"a",v))
                add_dd("b", lambda:d.b, lambda v:setattr(d,"b",v))
                add_num("I (A)", lambda:d.I, lambda v:setattr(d,"I",float(v)))
            elif isinstance(d,Capacitor):
                add_dd("a", lambda:d.a, lambda v:setattr(d,"a",v))
                add_dd("b", lambda:d.b, lambda v:setattr(d,"b",v))
                add_num("C (F)", lambda:d.C, lambda v:setattr(d,"C",float(v)))
            elif isinstance(d,NMOS):
                add_dd("d", lambda:d.d, lambda v:setattr(d,"d",v))
                add_dd("g", lambda:d.g, lambda v:setattr(d,"g",v))
                add_dd("s", lambda:d.s, lambda v:setattr(d,"s",v))
                add_num("k (A/V²)", lambda:d.k, lambda v:setattr(d,"k",float(v)))
                add_num("Vt (V)",   lambda:d.Vt,lambda v:setattr(d,"Vt",float(v)))
                add_num("λ (1/V)",  lambda:d.lam,lambda v:setattr(d,"lam",float(v)))
            elif isinstance(d,PMOS):
                add_dd("d", lambda:d.d, lambda v:setattr(d,"d",v))
                add_dd("g", lambda:d.g, lambda v:setattr(d,"g",v))
                add_dd("s", lambda:d.s, lambda v:setattr(d,"s",v))
                add_num("k (A/V²)", lambda:d.k, lambda v:setattr(d,"k",float(v)))
                add_num("|Vt| (V)", lambda:d.Vt,lambda v:setattr(d,"Vt",float(v)))
                add_num("λ (1/V)",  lambda:d.lam,lambda v:setattr(d,"lam",float(v)))
            elif isinstance(d,PolyCurrentToGND):
                add_dd("node", lambda:d.n, lambda v:setattr(d,"n",v))
                add_num("g1 (A/V)", lambda:d.g1, lambda v:setattr(d,"g1",float(v)))
                add_num("g2 (A/V²)",lambda:d.g2, lambda v:setattr(d,"g2",float(v)))
                add_num("I0 (A)",   lambda:d.I0, lambda v:setattr(d,"I0",float(v)))
            elif isinstance(d,CrossConductance):
                add_dd("a", lambda:d.a, lambda v:setattr(d,"a",v))
                add_dd("ctrl", lambda:d.ctrl, lambda v:setattr(d,"ctrl",v))
                add_num("g (A/V)", lambda:d.g, lambda v:setattr(d,"g",float(v)))
            elif isinstance(d,ComparatorToCurrent):
                add_dd("plus",  lambda:d.plus,  lambda v:setattr(d,"plus",v))
                add_dd("minus", lambda:d.minus, lambda v:setattr(d,"minus",v))
                add_dd("out",   lambda:d.out,   lambda v:setattr(d,"out",v))
                add_num("Iout (A)", lambda:d.Iout, lambda v:setattr(d,"Iout",float(v)))
                add_num("offset (V)", lambda:d.offset, lambda v:setattr(d,"offset",float(v)))

    # ----- delete & cleanup
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
            if used:
                messagebox.showwarning("Delete node", f"Node '{n}' used by: {', '.join(used)}"); return
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

    # ----- transient
    def transient_step(self):
        dt=float(self.var_dt.get())
        res=self.circ.solve_transient_step(dt)
        if res is None:
            messagebox.showwarning("Transient","No solution this step."); return
        sol, meta, _ = res
        if hasattr(self,"izh"):
            v=sol.get("v"); u=sol.get("u")
            if v is not None and u is not None and v >= self.izh["v_spike"]:
                c=self.izh["c"]; d=self.izh["d"]
                self.circ.prev_node_voltages["v"]=float(c)
                self.circ.prev_node_voltages["u"]=float(u+d)
                sol["v"]=float(c); sol["u"]=float(u+d)
        path=self.var_csv.get().strip()
        if path:
            headers=["t"]+list(sol.keys())
            write_header=not os.path.exists(path)
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

    # ----- export
    def export_netlist(self):
        path=filedialog.asksaveasfilename(defaultextension=".cir", filetypes=[("Spice netlist",".cir"),("All","*.*")])
        if not path: return
        try: self.circ.export_netlist(path)
        except Exception as e:
            messagebox.showerror("Export failed", str(e)); return
        messagebox.showinfo("Export", f"Wrote:\n{path}")

    # ----- templates
    def add_template_izhikevich(self):
        c=self.circ
        c.ensure_node("v", x=360,y=220); c.ensure_node("u", x=360,y=360); c.ensure_node("Iin", x=200,y=220)
        self.izh={"Cv":1.0,"Cu":1.0,"a":0.02,"b":0.2,"c":-65.0,"d":8.0,"v_spike":30.0,"g2":0.04,"g1":5.0,"I0":140.0,"Iin_gain":1.0}
        c.devices += [Capacitor("Cv","v","GND",self.izh["Cv"]), Capacitor("Cu","u","GND",self.izh["Cu"])]
        c.devices.append(PolyCurrentToGND("Pv","v",g1=self.izh["g1"],g2=self.izh["g2"],I0=self.izh["I0"]))
        c.devices.append(CrossConductance("Gvu","v","u",g=1.0))
        c.devices.append(CrossConductance("GvI","v","Iin",g=self.izh["Iin_gain"]))
        c.devices.append(CrossConductance("Guv","u","v",g=self.izh["a"]*self.izh["b"]))
        c.devices.append(CrossConductance("Guu","u","u",g=-self.izh["a"]))
        c.prev_node_voltages.update({"v":-65.0,"u":self.izh["b"]*(-65.0)})
        self.selected_node="v"; self.panel_dev.set("Ideal Izhikevich template ready."); self.refresh()

    def add_template_izhikevich_nmos(self):
        c=self.circ
        c.ensure_node("v", x=360,y=220); c.ensure_node("u", x=360,y=360)
        c.ensure_node("VB", x=200,y=160); c.ensure_node("Xref", x=560,y=140)
        c.ensure_node("Xu", x=560,y=300); c.ensure_node("Iin", x=200,y=220)
        Cv,Cu=1.0,1.0
        kq,Vtq,lamq=2e-4,-0.6,0.01
        kl,Vtl,laml=2e-4,-0.2,0.0
        Rb1,Rb2=200000.0,200000.0
        Rref=100000.0; kbias,Vtbias,lambias=2e-4,-0.6,0.01
        Ru=100000.0; ku,Vtu,lamu=2e-4,-0.5,0.01; gain_u=1.0
        Rin=100000.0
        c.devices += [Capacitor("Cv","v","GND",Cv), Capacitor("Cu","u","GND",Cu)]
        c.devices += [Resistor("Rb1","VDD","VB",Rb1), Resistor("Rb2","VB","GND",Rb2)]
        c.devices.append(NMOS("MQ","v","v","GND",kq,Vtq,lamq))
        c.devices.append(NMOS("ML","v","VB","GND",kl,Vtl,laml))
        c.devices += [Resistor("Rref","VDD","Xref",Rref), NMOS("Mref","Xref","Xref","GND",kbias,Vtbias,lambias),
                      NMOS("Mbias","v","Xref","GND",kbias,Vtbias,lambias)]
        c.devices += [Resistor("Ru","u","GND",Ru), NMOS("Mu_ref","Xu","Xu","GND",ku,Vtu,lamu),
                      NMOS("Mu_sink","v","Xu","GND",ku*gain_u,Vtu,lamu)]
        c.devices.append(Resistor("Rin","Iin","v",Rin))
        c.prev_node_voltages.update({"v":-65.0,"u":-13.0})
        self.selected_node="v"
        self.panel_dev.set("NMOS-only dv/dt: MQ≈v², ML≈g·v, Mbias const, Mu_sink ~ u.")
        self.refresh()

    def add_template_u_with_PMOS_mirror(self):
        c=self.circ
        c.ensure_node("u", x=360,y=360); c.ensure_node("v", x=360,y=220); c.ensure_node("Xp", x=620,y=360)
        Cu=1.0; Rp=100000.0; kp,Vtp,lamp=2e-4,0.6,0.01
        kn_ctrl,Vtn_ctrl,lamn_ctrl=2e-4,-0.2,0.0
        kn_u,Vtn_u,lamn_u=2e-4,-0.2,0.0
        c.devices.append(Capacitor("Cu","u","GND",Cu))
        c.devices += [Resistor("Rp","VDD","Xp",Rp), PMOS("Mpref","Xp","Xp","VDD",kp,Vtp,lamp), PMOS("Mpout","u","Xp","VDD",kp,Vtp,lamp)]
        c.devices.append(NMOS("Mctrl","Xp","v","GND",kn_ctrl,Vtn_ctrl,lamn_ctrl))
        c.devices.append(NMOS("Mu_sink","u","u","GND",kn_u,Vtn_u,lamn_u))
        self.selected_node="u"; self.panel_dev.set("u via PMOS mirror: Mpout sources into u; Mu_sink sinks ~ u."); self.refresh()

# ---------- run ----------
if __name__ == "__main__":
    Editor().mainloop()
