"""Figure 1 (machine overview): draws panels A-C and stacks them above panels D-G.

A: three ways to keep bacteria growing for a long time (simulated 24 h of cell density for daily
   1:100 dilutions, a turbidostat that loses control to a biofilm, and the two-chamber
   turbidostat), with small schematics.
B: timeline of the experiment (0-60,000 generations): 4 h transfer period until ~19,000
   generations, then 2 h; dots = time points sampled in panels D-G.
C: zoom on ~19,000 generations in lineage D: OD traces with a 4 h transfer period (biofilm) and
   after switching to 2 h (../datasets/adhesion/biofilm_OD_traces.json).

Panels D-G (output/adhesion_DG.svg) are made by make_collage.ipynb, which must run first.
Writes output/adhesion.svg, output/adhesion.png and output/PDF/adhesion.pdf;
the PNG and PDF are rendered with rsvg-convert (librsvg), with Noto Sans installed.

Usage (from the figures/ folder): python make_adhesion.py
"""
import json, math, os, re, subprocess
import xml.etree.ElementTree as ET
import numpy as np

DATA = "../datasets/adhesion/biofilm_OD_traces.json"
DG = "output/adhesion_DG.svg"
OUT = "output/adhesion"
FIG_W = 504                    # pt (7 in)
W_A = 282                      # panel A width; panels B and C fill the rest
CAP_EM = 0.714                 # Noto Sans cap height (em)
TOP = 12.5 - CAP_EM * 6        # top of the content of the top row ("24 hours", "Generations")
MM = 72 / 25.4                 # 1 mm in pt

# shared style
GLASS, CAP, DARK = "#e6e6e6", "#cccccc", "#666666"
TEAL = "#3996ac"                     # medium inflow / waste needles
ORANGE, BLUE = "#f24a14", "#3b9ab2"  # chambers L and R
TSB = "#ffc8a3"                      # medium
CURVE = "#333333"
NOTE = "#666666"
LABEL = "#404040"                    # annotation text ("lag", "cleaning", ...)
BASE = "#b2b2b2"
SW = 0.55                            # curve stroke width
AX_W = 0.64                          # axis line width (as in panels D-G)


class Canvas:
    """Collects svg elements."""

    def __init__(self):
        self.out = []
        self._ids = 0

    def add(self, s):
        self.out.append(s)

    def uid(self, prefix):
        self._ids += 1
        return f"{prefix}{self._ids}"

    def text(self, x, y, s, size=6, fill=LABEL, anchor="start", bold=False, rotate=None, lines=None, leading=1.15):
        """Text at baseline y. lines: several lines (a tuple is set in regular weight)."""
        a = ' font-weight="bold"' if bold else ""
        a += f' transform="rotate({rotate} {x:.2f} {y:.2f})"' if rotate is not None else ""
        if lines:
            s = "".join(f'<tspan x="{x:.2f}" dy="{0 if i == 0 else size * leading:.2f}"'
                        + (' font-weight="normal"' if isinstance(l, tuple) else "")
                        + f'>{l[0] if isinstance(l, tuple) else l}</tspan>' for i, l in enumerate(lines))
        self.add(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}"{a}>{s}</text>')

    def line(self, x0, y0, x1, y1, color=DARK, sw=AX_W, dash=None):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<path d="M {x0:.2f} {y0:.2f} L {x1:.2f} {y1:.2f}" stroke="{color}" stroke-width="{sw}" fill="none"{da}/>')

    def arrowhead(self, x, y, ang, color=TEAL, size=3.0):
        """Filled triangle with its tip at (x, y), pointing along ang (radians)."""
        a1, a2 = ang + math.radians(150), ang - math.radians(150)
        p = [(x, y), (x + size * math.cos(a1), y + size * math.sin(a1)), (x + size * math.cos(a2), y + size * math.sin(a2))]
        self.add('<path d="M ' + " L ".join(f"{a:.2f} {b:.2f}" for a, b in p) + f' Z" fill="{color}"/>')

    def path(self, pts, color, sw=SW, fill_to=None, alpha=0.1):
        """Polyline; with fill_to (a baseline y), the area under it is shaded in the same colour."""
        d = "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in pts)
        if fill_to is not None:
            self.add(f'<path d="{d} L {pts[-1][0]:.2f} {fill_to:.2f} L {pts[0][0]:.2f} {fill_to:.2f} Z" fill="{color}" '
                     f'fill-opacity="{alpha}" stroke="none"/>')
        self.add(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{sw}" stroke-linejoin="round" '
                 f'stroke-linecap="round"/>')


# ================================================================================================
# Panel A
# ================================================================================================
DT = 1 / 3600                  # h (1 s)
T = np.arange(0, 24 + DT / 2, DT)
MU = math.log(2) / (20 / 60)   # TSB, 20 min doubling
OD_MAX = 0.375                 # plateau OD (daily dilution and biofilm)


def mu(od, K=OD_MAX, od_slow=0.1):
    """Per-capita growth rate in TSB: exponential (20 min doubling) up to OD 0.1, then a
    Gompertz-type slowdown, mu ~ ln(K/OD), reaching zero at the carrying capacity K."""
    if od <= od_slow:
        return MU
    return MU * max(math.log(K / od), 0.0) / math.log(K / od_slow)


def daily_dilution(t_dil=5.0, lag=1.5, td=55 / 60, dil=100, n_slow=0.4):
    """LTEE: 1:100 daily dilution in DM25 (density normalised to 1 at saturation).
    After a ~1.5 h lag, cells grow at the ancestor's rate (~55 min doubling); as glucose runs
    out, growth slows Gompertz-style (per-capita rate ~ ln(K/N)), giving a smooth plateau."""
    m = math.log(2) / td
    def rate(n):
        return m if n <= n_slow else m * max(math.log(1 / n), 0.0) / math.log(1 / n_slow)
    def grow(n, hours):
        for _ in range(int(hours / DT)):
            n += rate(n) * n * DT
        return n
    n = grow(1 / dil, 24 - lag)               # density at the start of the window (previous day)
    N = np.empty_like(T)
    for i, t in enumerate(T):
        if abs(t - t_dil) < DT / 2:
            n /= dil
        if t < t_dil or t > t_dil + lag:
            n += rate(n) * n * DT
        N[i] = n
    return N


def turbidostat(od0=0.003, thr=0.1, min_gap=15 / 60, t_bio=14.0, attach=0.005, K=OD_MAX, g_bio=0.8, shed=0.6):
    """Dilutes the culture 1/2 when the OD reaches 0.1, at most once every 15 min.
    From t_bio, planktonic cells increasingly attach to the walls (attachment rate ramps up
    linearly, `attach` per hour per hour). Attached cells keep growing and count in the OD but
    are not removed by dilution; while it grows, the biofilm sheds cells into the medium. The
    effective dilution ratio (OD after / OD before) climbs gradually towards 1: the machine's
    dilutions stay visible but shrink, and only fade out near the end as the OD nears K.
    Growth: 20 min doubling (biofilm: g_bio x slower), Gompertz-type slowdown above OD 0.1."""
    P, B = od0, 0.0
    OD = np.zeros_like(T)
    last = -1.0
    for i, t in enumerate(T):
        r = mu(P + B, K=K)
        a = attach * max(t - t_bio, 0.0)
        sh = shed * r / MU * B
        P += ((r - a) * P + sh) * DT
        B += (g_bio * r * B + a * P - sh) * DT
        if P + B >= thr and t - last >= min_gap - 1e-9:
            P /= 2
            last = t
        OD[i] = P + B
    return OD


def two_chamber(thr=0.1, period=4.0, transfer=10):
    """Chambers L and R alternate: each grows for 4 h (turbidostat, 1/2 dilutions at OD 0.1),
    then 1/10 of it seeds the other chamber and it is cleaned. NaN = chamber not in use."""
    L = np.full_like(T, np.nan); R = np.full_like(T, np.nan)
    od = thr / transfer
    for i, t in enumerate(T):
        k = int(t // period + 1e-9)
        if i and int(T[i - 1] // period + 1e-9) != k:
            od /= transfer
        od += mu(od) * od * DT
        if od >= thr:
            od /= 2
        (L if k % 2 == 0 else R)[i] = od
    return L, R


def build_panel_A():
    """Returns (canvas, vertical shift that brings the top of "24 hours" to TOP, bottom of the
    lowest reactor after that shift)."""
    c = Canvas()
    X0, X1 = 110, 278              # curve column = 24 h
    sx = lambda t: X0 + (X1 - X0) * t / 24
    LOGO = 91                      # schematics column centre
    TXT = 10.08                    # titles left edge (aligned with the "A" tag)
    TITLE_SIZE = 8.75
    ROWS = {"daily": (35, 68), "turbido": (87, 120), "L": (137, 162), "R": (175, 200)}
    OD_SCALE = (68 - 35) / OD_MAX  # pt per OD unit, shared by all density curves
    AREA_ALPHA = {CURVE: 0.05, ORANGE: 0.1, BLUE: 0.1}
    TRANSFER = dict(color=NOTE, sw=0.5, head=2.3, dash="1 1")   # dotted transfer arrows

    def path(pts, color, fill_to=None):
        c.path(pts, color, fill_to=fill_to, alpha=AREA_ALPHA[color])

    def arc_arrow(p0, p1, bulge, side=1, color=TEAL, sw=0.7, head=3.0, dash=None):
        """Arc-shaped arrow from p0 to p1, symmetric about the perpendicular bisector of the chord.
        bulge: how far the arc bows out (pt); side: +1/-1 picks which side of the chord."""
        dx, dy = p1[0] - p0[0], p1[1] - p0[1]
        L = math.hypot(dx, dy)
        nx, ny = -dy / L * side, dx / L * side
        k = 1.333 * bulge                       # cubic control offset giving a peak of ~bulge
        c0 = (p0[0] + nx * k + dx * 0.1, p0[1] + ny * k + dy * 0.1)
        c1 = (p1[0] + nx * k - dx * 0.1, p1[1] + ny * k - dy * 0.1)
        ang = math.atan2(p1[1] - c1[1], p1[0] - c1[0])
        q = (p1[0] - 0.8 * head * math.cos(ang), p1[1] - 0.8 * head * math.sin(ang))
        da = f' stroke-dasharray="{dash}"' if dash else ""
        c.add(f'<path d="M {p0[0]:.2f} {p0[1]:.2f} C {c0[0]:.2f} {c0[1]:.2f} {c1[0]:.2f} {c1[1]:.2f} {q[0]:.2f} {q[1]:.2f}" '
              f'fill="none" stroke="{color}" stroke-width="{sw}"{da}/>')
        c.arrowhead(p1[0], p1[1], ang, color, head)
        return 0.75 * k                         # max distance of the arc from the chord

    def baseline(y):
        """Zero line under each density curve: deliberately faint."""
        c.add(f'<path d="M {X0:.2f} {y:.2f} L {X1:.2f} {y:.2f}" stroke="#d9d9d9" stroke-width="0.25" fill="none"/>')

    def drop(x, y, r=0.9, fill=TSB):
        """Small liquid drop hanging below (x, y)."""
        c.add(f'<path d="M {x:.2f} {y:.2f} C {x+0.3:.2f} {y+0.9:.2f} {x+r:.2f} {y+1.3:.2f} {x+r:.2f} {y+2.0:.2f} '
              f'A {r:.2f} {r:.2f} 0 0 1 {x-r:.2f} {y+2.0:.2f} C {x-r:.2f} {y+1.3:.2f} {x-0.3:.2f} {y+0.9:.2f} {x:.2f} {y:.2f} Z" fill="{fill}"/>')

    # Vessels: light grey glass without outline, liquid inset by the glass thickness, rounded
    # bottoms, small rounded top corners, darker cap of the same width.
    GLASS_T = 0.55                # glass thickness (pt)
    CORNER = 0.6                  # radius of the small top corners (pt)
    REACTOR_CORNER = 2 * CORNER   # reactors: tube top corners and caps
    REACTOR_CAP_H = 2.2 * 1.75    # reactor lid thickness (pt)
    REACTOR_TUBE_H = 24.0         # tube length below the lid (1:3 aspect)
    REACTOR_H = REACTOR_CAP_H + REACTOR_TUBE_H

    def flask_path(cx, by, nw, nh, bw, bh, r, top_r=CORNER):
        top = by - bh - nh
        return (f"M {cx-nw:.2f} {top+top_r:.2f} Q {cx-nw:.2f} {top:.2f} {cx-nw+top_r:.2f} {top:.2f} "
                f"L {cx+nw-top_r:.2f} {top:.2f} Q {cx+nw:.2f} {top:.2f} {cx+nw:.2f} {top+top_r:.2f} "
                f"L {cx+nw:.2f} {top+nh:.2f} L {cx+bw-r*0.45:.2f} {by-r*0.9:.2f} "
                f"Q {cx+bw:.2f} {by:.2f} {cx+bw-r*1.2:.2f} {by:.2f} L {cx-bw+r*1.2:.2f} {by:.2f} "
                f"Q {cx-bw:.2f} {by:.2f} {cx-bw+r*0.45:.2f} {by-r*0.9:.2f} L {cx-nw:.2f} {top+nh:.2f} Z")

    def flask_outline(cx, by, nw, nh, bw, bh, r, lip=0.5, rim_h=1.1):
        """Outer shape of an Erlenmeyer: conical body, straight neck and a slight rim at the top."""
        top = by - bh - nh
        return (f"M {cx-nw-lip:.2f} {top+0.45:.2f} Q {cx-nw-lip:.2f} {top:.2f} {cx-nw-lip+0.45:.2f} {top:.2f} "
                f"L {cx+nw+lip-0.45:.2f} {top:.2f} Q {cx+nw+lip:.2f} {top:.2f} {cx+nw+lip:.2f} {top+0.45:.2f} "
                f"L {cx+nw+lip:.2f} {top+rim_h-0.35:.2f} Q {cx+nw+lip:.2f} {top+rim_h:.2f} {cx+nw:.2f} {top+rim_h+0.25:.2f} "
                f"L {cx+nw:.2f} {top+nh:.2f} L {cx+bw-r*0.45:.2f} {by-r*0.9:.2f} "
                f"Q {cx+bw:.2f} {by:.2f} {cx+bw-r*1.2:.2f} {by:.2f} L {cx-bw+r*1.2:.2f} {by:.2f} "
                f"Q {cx-bw:.2f} {by:.2f} {cx-bw+r*0.45:.2f} {by-r*0.9:.2f} L {cx-nw:.2f} {top+nh:.2f} "
                f"L {cx-nw:.2f} {top+rim_h+0.25:.2f} Q {cx-nw-lip:.2f} {top+rim_h:.2f} {cx-nw-lip:.2f} {top+rim_h-0.35:.2f} Z")

    FLASK = dict(nw=2.4, nh=5.0, bw=7.5, bh=11.5, r=1.6)
    FLASK_H = FLASK["bh"] + FLASK["nh"]

    def flask(cx, by, fill_frac, liquid):
        """Open Erlenmeyer flask, base centre (cx, by). Returns the y of its rim."""
        nw, nh, bw, bh, r = (FLASK[k] for k in ("nw", "nh", "bw", "bh", "r"))
        top = by - bh - nh
        c.add(f'<path d="{flask_outline(cx, by, nw, nh, bw, bh, r)}" fill="{GLASS}"/>')
        g = GLASS_T
        inner = flask_path(cx, by - g, nw - g, nh, bw - g * 1.15, bh - g, r * 0.8)
        cid = c.uid("flask")
        c.add(f'<clipPath id="{cid}"><path d="{inner}"/></clipPath>')
        ly = by - bh * fill_frac
        c.add(f'<rect x="{cx-bw:.2f}" y="{ly:.2f}" width="{2*bw:.2f}" height="{by-ly:.2f}" fill="{liquid}" clip-path="url(#{cid})"/>')
        return top

    PIP_SCALE = 0.75
    PIP_TILT = 15.0                  # degrees; the top of the pipette leans left, pivoting on the tip
    PIP_H = 21.5                     # unscaled height of the silhouette (tip to plunger button)

    def pipette(cx, tip_y, scale=PIP_SCALE, tilt=PIP_TILT):
        """Flat micropipette silhouette (same grey as the lids), tip at (cx, tip_y).
        Drawn with the tip at the origin, then scaled and tilted about the tip. Returns the
        (x, y) of the pipette's mid-point (on its axis)."""
        y_noz, y_body, t = -6.0, -10.2, -18.7     # tip/nozzle, nozzle/body, top of body
        sil = (f"M -0.25 0 L -0.9 {y_noz} L -1.3 {y_body} "
               f"L -2.1 {y_body} Q -2.6 {y_body} -2.6 {y_body-0.6} L -2.6 {t+2.0} "
               # finger hook: sticks out to the left and curls down, towards the tip
               f"L -3.4 {t+2.0} Q -3.9 {t+2.0} -3.95 {t+2.6} L -4.0 {t+3.6} Q -4.45 {t+4.1} -4.9 {t+3.6} "
               f"L -4.95 {t+2.0} Q -4.95 {t+0.6} -3.6 {t+0.6} L -2.6 {t+0.6} Q -2.6 {t} -2.0 {t} "
               # plunger rod + button
               f"L -0.6 {t} L -0.6 {t-1.6} L -1.5 {t-1.6} Q -1.7 {t-1.6} -1.7 {t-1.9} L -1.7 {t-2.5} "
               f"Q -1.7 {t-2.8} -1.4 {t-2.8} L 1.4 {t-2.8} Q 1.7 {t-2.8} 1.7 {t-2.5} L 1.7 {t-1.9} "
               f"Q 1.7 {t-1.6} 1.5 {t-1.6} L 0.6 {t-1.6} L 0.6 {t} L 2.0 {t} Q 2.6 {t} 2.6 {t+0.6} "
               f"L 2.6 {y_body-0.6} Q 2.6 {y_body} 2.1 {y_body} L 1.3 {y_body} L 0.9 {y_noz} L 0.25 0 Z")
        c.add(f'<path d="{sil}" fill="{CAP}" transform="translate({cx:.2f} {tip_y:.2f}) rotate({-tilt}) scale({scale})"/>')
        a = math.radians(tilt)
        h = PIP_H / 2 * scale
        return cx - h * math.sin(a), tip_y - h * math.cos(a)

    def quarter_arrow(p0, p1, color=TEAL, sw=0.7, head=3.0):
        """Quarter-circle arrow: leaves p0 going straight up and ends at p1 going right
        (requires p1 = p0 + (R, -R)). The arrowhead starts where the arc ends, pointing right."""
        R = p1[0] - p0[0]
        c.add(f'<path d="M {p0[0]:.2f} {p0[1]:.2f} A {R:.2f} {R:.2f} 0 0 1 {p1[0]:.2f} {p1[1]:.2f}" fill="none" '
              f'stroke="{color}" stroke-width="{sw}"/>')
        base = head * math.cos(math.radians(30))          # tip-to-base length of the arrowhead
        c.arrowhead(p1[0] + base, p1[1], 0.0, color, head)

    def title(lines, cy):
        """Row title, vertically centred on cy (cap height 0.72 em, leading 1.15)."""
        lead = TITLE_SIZE * 1.15
        b1 = cy - ((len(lines) - 1) * lead - 0.72 * TITLE_SIZE) / 2
        c.text(TXT, b1, None, size=TITLE_SIZE, fill="#333333", bold=True, lines=lines)

    def tube_path(cx, top, by, w, top_r=CORNER):
        r = w / 2
        return (f"M {cx-r:.2f} {top+top_r:.2f} Q {cx-r:.2f} {top:.2f} {cx-r+top_r:.2f} {top:.2f} "
                f"L {cx+r-top_r:.2f} {top:.2f} Q {cx+r:.2f} {top:.2f} {cx+r:.2f} {top+top_r:.2f} "
                f"L {cx+r:.2f} {by-r:.2f} A {r:.2f} {r:.2f} 0 0 1 {cx-r:.2f} {by-r:.2f} Z")

    def reactor(cx, ty, h=REACTOR_TUBE_H, letter=None, liquid=TSB, lw=0.65):
        """Test-tube reactor (1:3), cap top at ty. Medium comes in through the cap (drop at the tip
        of the needle); the waste needle sits at the liquid surface. Returns geometry."""
        w, cap_h, ext = 8.0, REACTOR_CAP_H, 7.0
        top = ty + cap_h * 0.45                 # the cap sits on the rim of the tube
        by = ty + cap_h + h
        c.add(f'<path d="{tube_path(cx, top, by, w, top_r=REACTOR_CORNER)}" fill="{GLASS}"/>')
        g = GLASS_T
        inner = tube_path(cx, top + g, by - g, w - 2 * g, top_r=REACTOR_CORNER - g)
        cid = c.uid("tube")
        c.add(f'<clipPath id="{cid}"><path d="{inner}"/></clipPath>')
        ly = ty + cap_h + h * 0.35
        c.add(f'<rect x="{cx-w/2:.2f}" y="{ly:.2f}" width="{w:.2f}" height="{by-ly:.2f}" fill="{liquid}" clip-path="url(#{cid})"/>')
        c.add(f'<rect x="{cx-w/2:.2f}" y="{ty:.2f}" width="{w:.2f}" height="{cap_h:.2f}" rx="{REACTOR_CORNER}" fill="{CAP}"/>')
        xi, xo = cx - 2.0, cx + 2.0
        # medium in: needle stops at the ceiling, with a drop of medium at its tip
        c.add(f'<path d="M {xi:.2f} {ty-ext:.2f} L {xi:.2f} {ty+cap_h+0.2:.2f}" stroke="{TEAL}" stroke-width="{lw}" fill="none"/>')
        drop(xi, ty + cap_h + 0.2, r=0.85)
        c.arrowhead(xi, ty - 0.9, math.pi / 2)
        # waste out: needle down to the liquid surface
        c.add(f'<path d="M {xo:.2f} {ty-ext:.2f} L {xo:.2f} {ly:.2f}" stroke="{TEAL}" stroke-width="{lw}" fill="none"/>')
        c.arrowhead(xo, ty - ext - 0.2, -math.pi / 2)
        if letter:
            # centre the glyph's ink (Noto Sans: L and R sit ~0.04 em right of their advance centre,
            # cap height 0.714 em) on the middle of the tube's rounded bottom
            size = 4.5
            ink_dx = {"L": 0.036, "R": 0.041}[letter] * size
            c.text(cx - ink_dx, by - w / 2 + 0.357 * size, letter, size=size, fill="#000000", anchor="middle")
        return dict(w=w, by=by, top_y=ty - ext - 0.2)

    # ---------------------------------------------------------------- header: 24 h bracket
    yb = 26
    c.add(f'<path d="M {X0:.2f} {yb+3:.2f} L {X0:.2f} {yb:.2f} L {X1:.2f} {yb:.2f} L {X1:.2f} {yb+3:.2f}" '
          f'stroke="{NOTE}" stroke-width="{SW}" fill="none"/>')
    c.text((X0 + X1) / 2, yb - 2.5, "24 hours", size=7, fill="#000000", anchor="middle")

    # ---------------------------------------------------------------- 1. daily dilutions
    y0, y1 = ROWS["daily"]
    t_dil, lag = 5.0, 1.5
    c.add('<g id="daily-dilutions">')
    title(["Daily", "dilutions", ("(e.g., the",), ("Lenski LTEE)",)], (y0 + y1) / 2)
    N = daily_dilution(t_dil, lag)
    baseline(y1)
    yy = lambda v: y1 - OD_SCALE * OD_MAX * v     # N is normalised to the plateau (OD_MAX)
    pre = T < t_dil
    path([(sx(a), yy(b)) for a, b in zip(T[pre][::10], N[pre][::10])] + [(sx(t_dil), yy(N[pre][-1]))], CURVE, fill_to=y1)
    post = ~pre
    path([(sx(t_dil), yy(N[pre][-1]))] + [(sx(a), yy(b)) for a, b in zip(T[post][::10], N[post][::10])]
         + [(sx(T[-1]), yy(N[-1]))], CURVE, fill_to=y1)
    # schematic: saturated culture -> micropipette -> fresh medium. A drop hangs just above the
    # mouth of flask 2 and the pipette tip touches it. The pipetting arrow is a quarter circle that
    # starts just above the centre of flask 1 and ends, pointing right, at the pipette's mid-height.
    fx1, fx2 = LOGO - 9, LOGO + 9
    DROP_H = 2.0 + 0.7                                     # drop() height for r = 0.7
    tip_above_rim = 0.8 + DROP_H                           # 0.8 pt between the drop and the rim
    sketch_h = FLASK_H + tip_above_rim + PIP_H * PIP_SCALE
    fby = (y0 + y1) / 2 + sketch_h / 2                     # flask base, centring the sketch on the row
    rim1 = flask(fx1, fby, 0.62, TSB)
    rim2 = flask(fx2, fby, 0.62, TSB)
    tip_y = rim2 - tip_above_rim
    _, y_mid = pipette(fx2, tip_y)
    drop(fx2, tip_y, r=0.7)
    y_start = rim1 - 1.0
    quarter_arrow((fx1, y_start), (fx1 + y_start - y_mid, y_mid))
    # annotations; offsets measured to the ink (Noto Sans 6 pt: "l" left bearing 0.51 pt,
    # "n" right bearing 0.49 pt, x-height 3.22 pt)
    c.text(sx(t_dil) + 1.5 * MM - 0.51, y1 - 3.5, "lag")                 # 1.5 mm right of the dilution
    t_01 = T[np.argmax(post & (N * OD_MAX >= 0.1))]                      # where the culture reaches OD 0.1
    x_01, y_01 = sx(t_01), yy(0.1 / OD_MAX)
    W_EXPONENTIAL = 33.19                                                # width of "exponential" at 6 pt
    c.text(x_01 - W_EXPONENTIAL / 2 + 4.0, y_01 - 1.4 - 1.0 - 6 * 1.15, None, anchor="middle",
           lines=["exponential", "growth"])                              # just above the curve
    c.text(W_A - 3 * MM + 0.49, y0 + 2 * MM + 3.22, "starvation", anchor="end")   # 3 mm from the right edge, 2 mm below the plateau
    c.add("</g>")

    # ---------------------------------------------------------------- 2. turbidostat
    y0, y1 = ROWS["turbido"]
    c.add('<g id="turbidostat">')
    title(["Turbidostat"], (y0 + y1) / 2)
    OD = turbidostat()
    baseline(y1)
    yyT = lambda v: y1 - OD_SCALE * v
    path([(sx(a), yyT(b)) for a, b in zip(T[::20], OD[::20])] + [(sx(T[-1]), yyT(OD[-1]))], CURVE, fill_to=y1)
    gT = reactor(LOGO, y1 - REACTOR_H - 0.2)
    drops = np.where(np.diff(OD) < -0.002)[0]
    c.text(sx(T[drops[0]]), yyT(0.1) - 1 * MM, "exponential growth")   # starts at the first dilution
    cap_ht = 0.72 * 6                                       # cap height of 6 pt Noto Sans
    lab_top = gT["top_y"]
    lab_bot = lab_top + cap_ht + 2 * 6 * 1.05 + 1.6        # bottom of "dilutions)" incl. descender
    t_hit = T[np.argmax(yyT(OD) < lab_bot + 1.5)]           # where the rising curve reaches the label
    W_DILUTIONS = 26.22                                     # width of "dilutions)" at 6 pt (widest line)
    c.text(sx(t_hit) - 2 - W_DILUTIONS / 2, lab_top + cap_ht, None, anchor="middle",
           lines=["biofilm", "(escapes", "dilutions)"], leading=1.05)
    c.add("</g>")

    # ---------------------------------------------------------------- 3. two-chamber turbidostat
    c.add('<g id="two-chamber-turbidostat">')
    # each curve is placed so that the middle of its range (OD 0 to 0.1) is level with the middle
    # of its reactor (cap + tube)
    tyL, tyR = ROWS["L"][1] - REACTOR_H - 0.2, ROWS["R"][1] - REACTOR_H - 0.2   # reactor cap tops
    cL, cR = tyL + REACTOR_H / 2, tyR + REACTOR_H / 2
    l1 = cL + OD_SCALE * 0.05                              # curve baselines
    r1 = cR + OD_SCALE * 0.05
    title(["Two-chamber", "turbidostat", ("(this work)",)], (cL + cR) / 2 - 3.0)
    Lc, Rc = two_chamber()
    yyL = lambda v: l1 - OD_SCALE * v
    yyR = lambda v: r1 - OD_SCALE * v
    ends = {}
    for b, C, col, yf in [(l1, Lc, ORANGE, yyL), (r1, Rc, BLUE, yyR)]:
        baseline(b)
        for k in range(6):
            m = (T >= k * 4) & (T < (k + 1) * 4) & ~np.isnan(C)
            if not m.any():
                continue
            t_end = min((k + 1) * 4, 24)
            pts = [(sx(t), yf(v)) for t, v in zip(T[m][::10], C[m][::10])] + [(sx(t_end), yf(C[m][-1]))]
            path(pts, col, fill_to=b)
            ends[k] = pts[-1]
    for k in range(6):
        if (k + 1) * 4 < 24:  # transfer: dotted arrow from the end of the donor curve to the recipient
            x, y_end = ends[k]
            if k % 2 == 0:
                ya, yb2, ang = y_end + 1.5, yyR(0.002) - 1.0, math.pi / 2
            else:
                ya, yb2, ang = y_end - 1.5, yyL(0.002) + 1.0, -math.pi / 2
            q = yb2 - 0.8 * TRANSFER["head"] * math.sin(ang)
            c.add(f'<path d="M {x:.2f} {ya:.2f} L {x:.2f} {q:.2f}" stroke="{NOTE}" stroke-width="{TRANSFER["sw"]}" '
                  f'stroke-dasharray="{TRANSFER["dash"]}" fill="none"/>')
            c.arrowhead(x, yb2, ang, color=NOTE, size=TRANSFER["head"])
        idle = r1 if k % 2 == 0 else l1
        c.text(sx(k * 4 + 2), idle - OD_SCALE * 0.05 + 0.36 * 6, "cleaning", anchor="middle")   # centred at OD 0.05
    # two reactors aligned with the two curves, dotted transfer arrows on each side
    gL = reactor(LOGO, tyL, letter="L")
    gR = reactor(LOGO, tyR, letter="R")
    xl, xr = LOGO - gL["w"] / 2 - 1.8, LOGO + gL["w"] / 2 + 1.8
    yL, yR = gL["by"] - 7, gR["by"] - 7
    bow = arc_arrow((xl, yL), (xl, yR), bulge=7.5, **TRANSFER)        # L -> R (left)
    arc_arrow((xr, yR), (xr, yL), bulge=7.5, **TRANSFER)              # R -> L (right)
    c.text(xl - bow - 2.4, (yL + yR) / 2, "transfer", size=5, anchor="middle", rotate=-90)
    c.add("</g>")

    # shift everything (but the "A" tag) so that the top of "24 hours" is at TOP
    shift = (yb - 2.5 - CAP_EM * 7) - TOP
    return c, shift, gR["by"] - shift


# ================================================================================================
# Panels B and C
# ================================================================================================
def build_panels_BC(a_bottom):
    """Panels B and C in their own coordinates (x from 0 to FIG_W - W_A); the bottom of C's
    "Time (h)" title is level with a_bottom, the bottom of panel A."""
    c = Canvas()
    COLOR = {"L": ORANGE, "R": BLUE}
    E_DARK, E_LIGHT = "#4ea4b9", "#89c2d0"     # dark and light OFF (blue) colours of panel E
    TAG = dict(size=12, fill="#000000", bold=True)

    # ================================================================ B: timeline
    c.add('<g id="panel-B">')
    c.text(4.32, 11.52, "B", **TAG)
    GX0, GX1 = 16.7, 202.7
    sg = lambda g: GX0 + (GX1 - GX0) * g / 60.0          # generations (x1000) -> x
    x19 = sg(19.0)                                        # change of transfer period
    # top to bottom: "Generations" (its top level with panel A's "24 hours"), generation numbers,
    # axis with the sampled time points (panels D-G), then the timeline bar
    cap_y = TOP + CAP_EM * 7                              # "Generations" baseline (7 pt)
    num_y = cap_y + 8.5                                   # generation numbers baseline
    ax_y = num_y + 3.0                                    # axis line
    bar_t, bar_b = ax_y + 4.0, ax_y + 16.0
    APEX_Y = TOP + 37.0                                   # apex of the zoom lines
    r, tip = 1.2, 7.0
    W_NUM = 21.90                                         # width of "60,000" at 7 pt
    c.text(sg(60) + W_NUM / 2, cap_y, "Generations", size=7, fill="#000000", anchor="end")   # flush right
    c.line(sg(0), ax_y, sg(60), ax_y)
    for g in (0, 15, 30, 45, 60):
        c.add(f'<circle cx="{sg(g):.2f}" cy="{ax_y:.2f}" r="1.3" fill="#333333"/>')
        c.text(sg(g), num_y, f"{g * 1000:,}", size=7, fill="#000000", anchor="middle")
    c.add(f'<rect x="{sg(0):.2f}" y="{bar_t:.2f}" width="{x19 - sg(0) - 0.6:.2f}" height="{bar_b - bar_t:.2f}" rx="{r}" fill="{E_DARK}"/>')
    # second stretch ends in a triangle: the experiment continues beyond 60,000 generations
    xa, xb, ym = x19 + 0.6, sg(60), (bar_t + bar_b) / 2
    c.add(f'<path d="M {xa + r:.2f} {bar_t:.2f} L {xb:.2f} {bar_t:.2f} L {xb + tip:.2f} {ym:.2f} L {xb:.2f} {bar_b:.2f} '
          f'L {xa + r:.2f} {bar_b:.2f} Q {xa:.2f} {bar_b:.2f} {xa:.2f} {bar_b - r:.2f} L {xa:.2f} {bar_t + r:.2f} '
          f'Q {xa:.2f} {bar_t:.2f} {xa + r:.2f} {bar_t:.2f} Z" fill="{E_LIGHT}"/>')
    LAB = 7.46                                            # bar labels font size
    c.text((sg(0) + x19) / 2, ym + 0.36 * LAB, "4 h transfers", size=LAB, fill="#ffffff", anchor="middle")
    c.text((x19 + xb) / 2, ym + 0.36 * LAB, "2 h transfers", size=LAB, fill="#ffffff", anchor="middle")
    c.line(x19, ax_y, x19, APEX_Y, color="#7f7f7f", dash="3 3")   # same style as in D-G
    c.add("</g>")

    # ================================================================ C: zoom on 19,000 generations
    data = json.load(open(DATA))
    c.add('<g id="panel-C">')

    # vertical layout, from the bottom up
    T_TITLE = a_bottom - 1.1                  # "Time (h)" baseline: its parentheses end level with panel A
    TICK_LAB = T_TITLE - 9.5                  # x tick labels (7 pt)
    PY1 = TICK_LAB - CAP_EM * 7 - 2.0 - 2.5   # x axis = OD 0 (ticks 2.5 pt, then 2 pt gap)
    PY0 = PY1 - 63.37                         # top of the plots
    SUB_Y = PY0 - 5.0                         # baseline of the subpanel titles
    OD_TOP = max(v for _, pts in data["4h"] for _, v in pts)   # tallest peak touches the top of the plot
    T_MAX = 10.3
    PANELS = [("4h", 34.0, 118.0, "4 h transfer period"), ("2h", 128.0, 212.0, "2 h transfer period")]
    yy = lambda od: PY1 - (PY1 - PY0) * od / OD_TOP

    # zoom lines from the 19,000-gen mark to the two plots, and the title of C inside them
    FAN_END = SUB_Y - CAP_EM * 7.5 - 2.0
    fan_l, fan_r = (PANELS[0][1], FAN_END), (PANELS[1][2], FAN_END)
    c.line(x19, APEX_Y, *fan_l, color=BASE, sw=0.5)
    c.line(x19, APEX_Y, *fan_r, color=BASE, sw=0.5)
    def fan_edges(y):
        f = (y - APEX_Y) / (FAN_END - APEX_Y)
        return x19 + (fan_l[0] - x19) * f, x19 + (fan_r[0] - x19) * f
    W1, W2 = 67.83, 113.00                    # widths: line 1 regular (7.5 pt), line 2 bold (7 pt)
    # vertically: centred between the apex and the end of the fan, lowered only if it does not fit;
    # horizontally: a common centre for both lines, inside the fan with >= 1.5 pt margins
    block = CAP_EM * 7.5 + 8.6 + 1.5          # cap of line 1 .. descender of line 2
    B2 = (APEX_Y + FAN_END) / 2 + block / 2 - 1.5
    while True:
        B1 = B2 - 8.6
        lo1, hi1 = fan_edges(B1 - CAP_EM * 7.5); lo2, hi2 = fan_edges(B2 - CAP_EM * 7)
        lo = max(lo1 + W1 / 2 + 1.5, lo2 + W2 / 2 + 1.5); hi = min(hi1 - W1 / 2 - 1.5, hi2 - W2 / 2 - 1.5)
        if lo <= hi:
            break
        B2 += 0.25
        assert B2 < FAN_END - 1.0, "title does not fit in the zoom fan"
    XC = (lo + hi) / 2
    c.text(XC, B1, "19,000 generations", size=7.5, fill="#333333", anchor="middle")
    c.text(XC, B2, "D evolves &lt;4 h biofilm formation", size=7, fill="#333333", anchor="middle", bold=True)
    c.text(4.32, B1, "C", **TAG)

    # subpanel titles, with an arrow showing the change
    cx = [(p[1] + p[2]) / 2 for p in PANELS]
    for (key, x0, x1, label), x in zip(PANELS, cx):
        c.text(x, SUB_Y, label, size=7.5, fill="#333333", anchor="middle", bold=True)
    WL = 69.6                                 # width of the subpanel titles
    ya = SUB_Y - 0.36 * 7.5
    x_from, x_to = cx[0] + WL / 2 + 2.5, cx[1] - WL / 2 - 2.5
    c.line(x_from, ya, x_to - 2.0, ya)
    c.arrowhead(x_to, ya, 0.0, DARK, 2.5)

    # the two plots (axes start exactly at 0 h and OD 0)
    peaks4 = []
    for key, x0, x1, label in PANELS:
        sx = lambda t, x0=x0, x1=x1: x0 + (x1 - x0) * t / T_MAX
        for ch, pts in data[key]:              # one trace per growth period, coloured by chamber
            pts = [(max(t, 0.0), max(v, 0.0)) for t, v in pts if t <= T_MAX + 1e-6]
            c.path([(sx(t), yy(v)) for t, v in pts], COLOR[ch], fill_to=yy(0))
            tm, vm = max(pts, key=lambda p: p[1])
            if key == "2h":                    # chamber letter above the peak
                c.text(sx(tm), yy(vm) - 2.5, ch, size=5.5, fill=COLOR[ch], anchor="middle")
            else:                              # placed after the arrows (below)
                peaks4.append((ch, tm, vm))
        c.line(x0, yy(0.1), x1, yy(0.1), color="#999999", sw=0.5, dash="1.5 1.5")   # dilution threshold
        c.line(x0, PY1, x1, PY1)
        c.line(x0, PY0, x0, PY1)
        for t in range(0, 11, 2):
            c.line(sx(t), PY1, sx(t), PY1 + 2.5)
            c.text(sx(t), TICK_LAB, str(t), size=7, fill="#000000", anchor="middle")
        for od in (0.0, 0.1, 0.2):
            c.line(x0 - 2.5, yy(od), x0, yy(od))
            if key == "4h":
                c.text(x0 - 4.0, yy(od) + 0.36 * 7, f"{od:.1f}", size=7, fill="#000000", anchor="end")

    # "biofilms", with small arrows to the three cycles
    sx4 = lambda t: PANELS[0][1] + (PANELS[0][2] - PANELS[0][1]) * t / T_MAX
    W_BIO = 22.51                                           # width of "biofilms" at 6 pt
    bx, by = sx4(4.0), yy(0.215)                            # left end / baseline of the word
    c.text(bx, by, "biofilms")
    # each arrow lies on the line from the centre of the word to the point where the cycle's OD
    # first exceeds 0.11 for good (after its last dilution peak below 0.11); it starts 1 mm from
    # the word and stops 1 mm before it first meets a density curve
    wx0, wx1, wy0, wy1 = bx, bx + W_BIO, by - 0.76 * 6, by  # word box (ascenders, no descenders)
    cxw, cyw = (wx0 + wx1) / 2, (wy0 + wy1) / 2
    flat = [(sx4(max(t, 0)), yy(max(v, 0))) for _, pts in data["4h"] for t, v in pts if t <= T_MAX + 1e-6]
    segs = list(zip(flat[:-1], flat[1:]))
    def seg_dist(p, a, b):
        ax_, ay_ = a; bx_, by_ = b; dx, dy = bx_ - ax_, by_ - ay_
        L2 = dx * dx + dy * dy
        u = 0 if L2 == 0 else max(0, min(1, ((p[0] - ax_) * dx + (p[1] - ay_) * dy) / L2))
        return math.hypot(p[0] - ax_ - u * dx, p[1] - ay_ - u * dy)
    arrows = []
    for _, pts in data["4h"]:
        peaks_below = [i for i in range(1, len(pts) - 1)
                       if pts[i][1] >= pts[i - 1][1] and pts[i][1] >= pts[i + 1][1] and pts[i][1] < 0.11]
        i0 = peaks_below[-1] if peaks_below else 0
        k = next(i for i in range(i0, len(pts)) if pts[i][1] > 0.11)
        (ta, va), (tb, vb) = pts[k - 1], pts[k]
        tc = ta + (tb - ta) * (0.11 - va) / (vb - va)
        px, py = sx4(tc), yy(0.11)
        dx, dy = px - cxw, py - cyw
        D = math.hypot(dx, dy); ux, uy = dx / D, dy / D
        s0 = min(((wx1 if ux > 0 else wx0) - cxw) / ux if ux else 1e9,
                 ((wy1 if uy > 0 else wy0) - cyw) / uy if uy else 1e9) + MM   # out of the word box
        s1 = D                                              # first contact with the curve
        sv = s0
        while sv < D:
            q = (cxw + sv * ux, cyw + sv * uy)
            if min(seg_dist(q, a, b) for a, b in segs) < 0.3:
                s1 = sv; break
            sv += 0.1
        s1 -= MM
        ex, ey = cxw + s1 * ux, cyw + s1 * uy
        c.line(cxw + s0 * ux, cyw + s0 * uy, ex - 1.6 * ux, ey - 1.6 * uy, color=LABEL, sw=0.4)
        arrows.append(((cxw + s0 * ux, cyw + s0 * uy), (ex, ey)))
        c.arrowhead(ex, ey, math.atan2(uy, ux), LABEL, 1.9)
    # 4 h chamber letters: left of their peak, tops level with it, 2.6 pt from the peak, slid
    # further left if they would touch the curve or an arrow (the second one goes above its peak)
    LETTER_W = {"R": 3.42, "L": 2.88}                       # widths at 5.5 pt
    def dense(a, b, step=0.2):
        n = max(1, int(math.hypot(b[0] - a[0], b[1] - a[1]) / step))
        return [(a[0] + (b[0] - a[0]) * i / n, a[1] + (b[1] - a[1]) * i / n) for i in range(n + 1)]
    obst = [p for a, b in segs + arrows for p in dense(a, b)]
    for k, (ch, tm, vm) in enumerate(peaks4):
        if k == 1:
            c.text(sx4(tm), yy(vm) - 2.5, ch, size=5.5, fill=COLOR[ch], anchor="middle")
            continue
        w, cap, top = LETTER_W[ch], CAP_EM * 5.5, yy(vm)
        xl = sx4(tm) - 2.6 - w
        while any(xl - 0.7 <= px <= xl + w + 0.7 and top - 0.7 <= py <= top + cap + 0.7 for px, py in obst):
            xl -= 0.2
        c.text(xl, top + cap, ch, size=5.5, fill=COLOR[ch])
    c.text(cx[1], yy(0.2), "biofilm disappears", anchor="middle")
    # axis titles
    c.text((PANELS[0][1] + PANELS[1][2]) / 2, T_TITLE, "Time (h)", size=7.5, fill="#000000", anchor="middle")
    c.text(10.5, (PY0 + PY1) / 2, "OD (600 nm)", size=7.5, fill="#000000", anchor="middle", rotate=-90)
    c.add("</g>")
    return c


# ================================================================================================
# Assemble: panels A-C on top, panels D-G (from make_collage.ipynb) below
# ================================================================================================
GAP = 10.0                     # from the bottom of panel A to the top of the "D" tag


def d_tag_top(svg):
    """y of the top of the "D" tag in the svglite output of make_collage.ipynb."""
    for t in svg.iter("{http://www.w3.org/2000/svg}text"):
        if (t.text or "").strip() == "D":
            m = re.search(r"translate\(([-\d.]+),\s*([-\d.]+)\)", t.get("transform", ""))
            y = float(m.group(2)) if m else float(t.get("y"))
            size = float(re.search(r"font-size:\s*([\d.]+)px", t.get("style")).group(1))
            return y - CAP_EM * size
    raise ValueError(f"no 'D' tag in {DG}")


a, a_shift, a_bottom = build_panel_A()
bc = build_panels_BC(a_bottom)

ET.register_namespace("", "http://www.w3.org/2000/svg")
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
dg = ET.parse(DG).getroot()
dg_w, dg_h = (float(dg.get(k).replace("pt", "")) for k in ("width", "height"))
assert abs(dg_w - FIG_W) < 0.5, f"panels D-G are {dg_w} pt wide, expected {FIG_W}"
dy = a_bottom + GAP - d_tag_top(dg)
fig_h = round(dy + dg_h, 1)
dg.set("x", "0"); dg.set("y", f"{dy:.3f}")
dg.set("width", f"{dg_w}"); dg.set("height", f"{dg_h}")

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}pt" height="{fig_h}pt" viewBox="0 0 {FIG_W} {fig_h}" '
       f'font-family="\'Noto Sans\', sans-serif">\n'
       f'<rect id="background" x="0" y="0" width="{FIG_W}" height="{fig_h}" fill="#ffffff"/>\n'
       f'<g id="panel-A">\n'
       f'<text x="10.08" y="11.52" font-size="12" font-weight="bold" fill="#000000">A</text>\n'
       f'<g transform="translate(0 {-a_shift:.3f})">\n' + "\n".join(a.out) + "\n</g>\n</g>\n"
       f'<g id="panels-BC" transform="translate({W_A} 0)">\n' + "\n".join(bc.out) + "\n</g>\n"
       + ET.tostring(dg, encoding="unicode") + "\n</svg>\n")
open(OUT + ".svg", "w").write(svg)

os.makedirs("output/PDF", exist_ok=True)
subprocess.run(["rsvg-convert", "-f", "png", "--dpi-x", "300", "--dpi-y", "300", "-b", "white",
                "-o", OUT + ".png", OUT + ".svg"], check=True)
subprocess.run(["rsvg-convert", "-f", "pdf", "-o", "output/PDF/adhesion.pdf", OUT + ".svg"], check=True)
print(f"wrote {OUT}.svg/.png and output/PDF/adhesion.pdf: {FIG_W} x {fig_h} pt")
