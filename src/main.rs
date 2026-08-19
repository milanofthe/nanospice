// nanospice: a classic SPICE circuit simulator in a single source file,
// hard budget 1000 lines of code, std only.
//
// Algorithm (Berkeley SPICE structure): MNA with branch currents for V, L, E.
// Newton-Raphson with pnjlim junction limiting and gmin stepping fallback for
// the operating point. Transient with trapezoidal companion models and
// iteration count timestep control. AC as small signal linearization at the
// operating point. Dense LU with partial pivoting, real and complex.
//
// Matrix convention: unknown index 0 is a dummy ground row/column so device
// stamps never need to special case ground; the LU loops start at index 1.

use std::collections::HashMap;
use std::f64::consts::{PI, SQRT_2};

const GMIN: f64 = 1e-12;
const RELTOL: f64 = 1e-3;
const VNTOL: f64 = 1e-6;
const ABSTOL: f64 = 1e-12;
const VT: f64 = 0.02585;
const ITL_OP: usize = 100;
const ITL_TRAN: usize = 10;

fn die(msg: &str) -> ! {
    eprintln!("nanospice: {}", msg);
    std::process::exit(1);
}

// stdout writer that exits quietly when the pipe is closed (e.g. piped to head)
fn outln(args: std::fmt::Arguments) {
    use std::io::Write;
    if writeln!(std::io::stdout().lock(), "{}", args).is_err() {
        std::process::exit(0);
    }
}

macro_rules! out {
    ($($a:tt)*) => { outln(format_args!($($a)*)) };
}

// ---------------------------------------------------------------- numbers ---

#[derive(Clone, Copy)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    fn new(re: f64, im: f64) -> Cx {
        Cx { re, im }
    }
    fn zero() -> Cx {
        Cx::new(0.0, 0.0)
    }
    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    fn sub(self, o: Cx) -> Cx {
        Cx::new(self.re - o.re, self.im - o.im)
    }
    fn mul(self, o: Cx) -> Cx {
        Cx::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
    fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        Cx::new((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
}

// SPICE number: longest parseable prefix plus unit suffix, e.g. 4.7k, 100n, 1meg.
fn num(tok: &str) -> Option<f64> {
    let t = tok.trim();
    let mut end = t.len();
    while end > 0 && t[..end].parse::<f64>().is_err() {
        end -= 1;
    }
    if end == 0 {
        return None;
    }
    let v: f64 = t[..end].parse().unwrap();
    let suffix = &t[end..];
    let mult = if suffix.starts_with("meg") {
        1e6
    } else if suffix.starts_with('t') {
        1e12
    } else if suffix.starts_with('g') {
        1e9
    } else if suffix.starts_with('k') {
        1e3
    } else if suffix.starts_with('m') {
        1e-3
    } else if suffix.starts_with('u') {
        1e-6
    } else if suffix.starts_with('n') {
        1e-9
    } else if suffix.starts_with('p') {
        1e-12
    } else if suffix.starts_with('f') {
        1e-15
    } else {
        1.0
    };
    Some(v * mult)
}

fn numx(tok: &str) -> f64 {
    num(tok).unwrap_or_else(|| die(&format!("bad number '{}'", tok)))
}

// ---------------------------------------------------------------- circuit ---

#[derive(Clone)]
enum Wave {
    Sin { vo: f64, va: f64, fr: f64, td: f64, theta: f64 },
    Pulse { v1: f64, v2: f64, td: f64, tr: f64, tf: f64, pw: f64, per: f64 },
}

fn wave_val(dc: f64, wave: &Option<Wave>, t: f64) -> f64 {
    match wave {
        None => dc,
        Some(Wave::Sin { vo, va, fr, td, theta }) => {
            if t < *td {
                *vo
            } else {
                let tp = t - td;
                vo + va * (-tp * theta).exp() * (2.0 * PI * fr * tp).sin()
            }
        }
        Some(Wave::Pulse { v1, v2, td, tr, tf, pw, per }) => {
            let mut tp = t - td;
            if *per > 0.0 && tp > 0.0 {
                tp %= per;
            }
            if tp <= 0.0 {
                *v1
            } else if tp < *tr {
                v1 + (v2 - v1) * tp / tr
            } else if tp < tr + pw {
                *v2
            } else if tp < tr + pw + tf {
                v2 + (v1 - v2) * (tp - tr - pw) / tf
            } else {
                *v1
            }
        }
    }
}

#[derive(Clone)]
enum Dev {
    R { p: usize, n: usize, v: f64 },
    C { p: usize, n: usize, v: f64 },
    L { p: usize, n: usize, v: f64, br: usize },
    V { p: usize, n: usize, dc: f64, ac: f64, wave: Option<Wave>, br: usize },
    I { p: usize, n: usize, dc: f64, ac: f64, wave: Option<Wave> },
    D { p: usize, n: usize, is: f64, nf: f64 },
    M { d: usize, g: usize, s: usize, kp: f64, vt0: f64, lambda: f64, pmos: bool },
    E { p: usize, n: usize, cp: usize, cn: usize, k: f64, br: usize },
    G { p: usize, n: usize, cp: usize, cn: usize, k: f64 },
}

#[derive(Clone)]
enum Analysis {
    Op,
    Dc { src: String, start: f64, stop: f64, step: f64 },
    Tran { tstep: f64, tstop: f64 },
    Ac { dec: bool, n: usize, f1: f64, f2: f64 },
}

struct Circuit {
    devs: Vec<Dev>,
    names: Vec<String>,
    nodes: Vec<String>,
    analyses: Vec<Analysis>,
    nb: usize,
}

// ----------------------------------------------------------------- parser ---

fn tok<'a>(toks: &[&'a str], i: usize) -> &'a str {
    toks.get(i).copied().unwrap_or_else(|| die("missing token on device or analysis line"))
}

fn pval(toks: &[&str], key: &str, default: f64) -> f64 {
    for t in toks {
        if let Some(rest) = t.strip_prefix(key) {
            if let Some(v) = rest.strip_prefix('=') {
                return numx(v);
            }
        }
    }
    default
}

// Source spec: bare value or dc/ac keywords plus sin(...)/pulse(...) waveforms.
fn src_spec(toks: &[&str]) -> (f64, f64, Option<Wave>) {
    let (mut dc, mut ac, mut wave) = (0.0, 0.0, None);
    let mut i = 0;
    while i < toks.len() {
        let take = |i: &mut usize| -> Vec<f64> {
            let mut v = Vec::new();
            while *i < toks.len() {
                match num(toks[*i]) {
                    Some(x) => {
                        v.push(x);
                        *i += 1;
                    }
                    None => break,
                }
            }
            v
        };
        match toks[i] {
            "dc" => {
                i += 1;
                dc = *take(&mut i).first().unwrap_or(&0.0);
            }
            "ac" => {
                i += 1;
                ac = *take(&mut i).first().unwrap_or(&1.0);
            }
            "sin" => {
                i += 1;
                let v = take(&mut i);
                let g = |k: usize| v.get(k).copied().unwrap_or(0.0);
                wave = Some(Wave::Sin { vo: g(0), va: g(1), fr: g(2), td: g(3), theta: g(4) });
            }
            "pulse" => {
                i += 1;
                let v = take(&mut i);
                let g = |k: usize| v.get(k).copied().unwrap_or(0.0);
                let pw = if v.len() > 5 { g(5) } else { f64::MAX / 4.0 };
                wave = Some(Wave::Pulse {
                    v1: g(0),
                    v2: g(1),
                    td: g(2),
                    tr: g(3).max(1e-12),
                    tf: g(4).max(1e-12),
                    pw,
                    per: g(6),
                });
            }
            t => {
                if let Some(x) = num(t) {
                    dc = x;
                }
                i += 1;
            }
        }
    }
    (dc, ac, wave)
}

fn parse(src: &str) -> Circuit {
    let mut lines: Vec<String> = Vec::new();
    for (ln, raw) in src.lines().enumerate() {
        let line = raw.split(';').next().unwrap().trim().to_ascii_lowercase();
        if ln == 0 || line.is_empty() || line.starts_with('*') {
            continue;
        }
        if let Some(rest) = line.strip_prefix('+') {
            match lines.last_mut() {
                Some(prev) => {
                    prev.push(' ');
                    prev.push_str(rest);
                }
                None => die("continuation line without a preceding line"),
            }
        } else {
            lines.push(line);
        }
    }
    let mut map: HashMap<String, usize> = HashMap::new();
    map.insert("0".into(), 0);
    map.insert("gnd".into(), 0);
    let mut nodes = vec!["0".to_string()];
    let mut devs: Vec<Dev> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut analyses: Vec<Analysis> = Vec::new();
    for line in &lines {
        let clean: String = line.chars().map(|c| if "(),".contains(c) { ' ' } else { c }).collect();
        let toks: Vec<&str> = clean.split_whitespace().collect();
        if toks.is_empty() {
            continue;
        }
        if toks[0] == ".end" {
            break;
        }
        if let Some(card) = toks[0].strip_prefix('.') {
            match card {
                "op" => analyses.push(Analysis::Op),
                "dc" => analyses.push(Analysis::Dc {
                    src: tok(&toks, 1).into(),
                    start: numx(tok(&toks, 2)),
                    stop: numx(tok(&toks, 3)),
                    step: numx(tok(&toks, 4)),
                }),
                "tran" => analyses.push(Analysis::Tran {
                    tstep: numx(tok(&toks, 1)),
                    tstop: numx(tok(&toks, 2)),
                }),
                "ac" => {
                    let kind = tok(&toks, 1);
                    if kind != "dec" && kind != "lin" {
                        die("only .ac dec and .ac lin are supported");
                    }
                    analyses.push(Analysis::Ac {
                        dec: kind == "dec",
                        n: numx(tok(&toks, 2)) as usize,
                        f1: numx(tok(&toks, 3)),
                        f2: numx(tok(&toks, 4)),
                    });
                }
                _ => eprintln!("nanospice: ignoring card .{}", card),
            }
            continue;
        }
        let mut nid = |s: &str| -> usize {
            if let Some(&i) = map.get(s) {
                return i;
            }
            let i = nodes.len();
            map.insert(s.into(), i);
            nodes.push(s.into());
            i
        };
        let dev = match toks[0].chars().next().unwrap() {
            'r' => Dev::R { p: nid(tok(&toks, 1)), n: nid(tok(&toks, 2)), v: numx(tok(&toks, 3)) },
            'c' => Dev::C { p: nid(tok(&toks, 1)), n: nid(tok(&toks, 2)), v: numx(tok(&toks, 3)) },
            'l' => Dev::L { p: nid(tok(&toks, 1)), n: nid(tok(&toks, 2)), v: numx(tok(&toks, 3)), br: 0 },
            'v' => {
                let (dc, ac, wave) = src_spec(&toks[3..]);
                Dev::V { p: nid(tok(&toks, 1)), n: nid(tok(&toks, 2)), dc, ac, wave, br: 0 }
            }
            'i' => {
                let (dc, ac, wave) = src_spec(&toks[3..]);
                Dev::I { p: nid(tok(&toks, 1)), n: nid(tok(&toks, 2)), dc, ac, wave }
            }
            'd' => Dev::D {
                p: nid(tok(&toks, 1)),
                n: nid(tok(&toks, 2)),
                is: pval(&toks[3..], "is", 1e-14),
                nf: pval(&toks[3..], "n", 1.0),
            },
            'm' => Dev::M {
                d: nid(tok(&toks, 1)),
                g: nid(tok(&toks, 2)),
                s: nid(tok(&toks, 3)),
                kp: pval(&toks[4..], "kp", 2e-5),
                vt0: pval(&toks[4..], "vt0", pval(&toks[4..], "vto", 0.0)),
                lambda: pval(&toks[4..], "lambda", 0.0),
                pmos: toks.get(5).map_or(false, |t| *t == "pmos"),
            },
            'e' => Dev::E {
                p: nid(tok(&toks, 1)),
                n: nid(tok(&toks, 2)),
                cp: nid(tok(&toks, 3)),
                cn: nid(tok(&toks, 4)),
                k: numx(tok(&toks, 5)),
                br: 0,
            },
            'g' => Dev::G {
                p: nid(tok(&toks, 1)),
                n: nid(tok(&toks, 2)),
                cp: nid(tok(&toks, 3)),
                cn: nid(tok(&toks, 4)),
                k: numx(tok(&toks, 5)),
            },
            _ => die(&format!("unknown device '{}'", toks[0])),
        };
        devs.push(dev);
        names.push(toks[0].to_string());
    }
    let nn = nodes.len();
    let mut nb = 0;
    for d in devs.iter_mut() {
        if let Dev::L { br, .. } | Dev::V { br, .. } | Dev::E { br, .. } = d {
            *br = nn + nb;
            nb += 1;
        }
    }
    Circuit { devs, names, nodes, analyses, nb }
}

// ---------------------------------------------------------- linear algebra ---

struct Sys {
    m: usize,
    a: Vec<f64>,
    b: Vec<f64>,
}

impl Sys {
    fn new(m: usize) -> Sys {
        Sys { m, a: vec![0.0; m * m], b: vec![0.0; m] }
    }
    fn add(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.m + j] += v;
    }
    fn rhs(&mut self, i: usize, v: f64) {
        self.b[i] += v;
    }
}

fn lu_solve(mut a: Vec<f64>, mut b: Vec<f64>, m: usize) -> Option<Vec<f64>> {
    for k in 1..m {
        let mut p = k;
        for i in k + 1..m {
            if a[i * m + k].abs() > a[p * m + k].abs() {
                p = i;
            }
        }
        if a[p * m + k].abs() < 1e-300 {
            return None;
        }
        if p != k {
            for j in 1..m {
                a.swap(k * m + j, p * m + j);
            }
            b.swap(k, p);
        }
        let piv = a[k * m + k];
        for i in k + 1..m {
            let f = a[i * m + k] / piv;
            if f != 0.0 {
                for j in k + 1..m {
                    a[i * m + j] -= f * a[k * m + j];
                }
                b[i] -= f * b[k];
            }
        }
    }
    for i in (1..m).rev() {
        let mut s = b[i];
        for j in i + 1..m {
            s -= a[i * m + j] * b[j];
        }
        b[i] = s / a[i * m + i];
    }
    b[0] = 0.0;
    Some(b)
}

fn clu_solve(mut a: Vec<Cx>, mut b: Vec<Cx>, m: usize) -> Option<Vec<Cx>> {
    for k in 1..m {
        let mut p = k;
        for i in k + 1..m {
            if a[i * m + k].abs() > a[p * m + k].abs() {
                p = i;
            }
        }
        if a[p * m + k].abs() < 1e-300 {
            return None;
        }
        if p != k {
            for j in 1..m {
                a.swap(k * m + j, p * m + j);
            }
            b.swap(k, p);
        }
        for i in k + 1..m {
            let f = a[i * m + k].div(a[k * m + k]);
            if f.abs() != 0.0 {
                for j in k + 1..m {
                    a[i * m + j] = a[i * m + j].sub(f.mul(a[k * m + j]));
                }
                b[i] = b[i].sub(f.mul(b[k]));
            }
        }
    }
    for i in (1..m).rev() {
        let mut s = b[i];
        for j in i + 1..m {
            s = s.sub(a[i * m + j].mul(b[j]));
        }
        b[i] = s.div(a[i * m + i]);
    }
    b[0] = Cx::zero();
    Some(b)
}

// ---------------------------------------------------------- device models ---

fn diode_eval(vd: f64, is: f64, vt: f64, gmin: f64) -> (f64, f64) {
    let e = (vd / vt).min(200.0).exp();
    (is * (e - 1.0) + gmin * vd, is * e / vt + gmin)
}

// Classic SPICE junction voltage limiting.
fn pnjlim(vnew: f64, vold: f64, vt: f64, vcrit: f64) -> f64 {
    if vnew > vcrit && (vnew - vold).abs() > 2.0 * vt {
        if vold > 0.0 {
            let arg = 1.0 + (vnew - vold) / vt;
            if arg > 0.0 {
                vold + vt * arg.ln()
            } else {
                vcrit
            }
        } else {
            vt * (vnew / vt).ln()
        }
    } else {
        vnew
    }
}

// Level 1 square law, terminal voltages already in the effective frame.
fn mos_eval(vgs: f64, vds: f64, kp: f64, vt0: f64, lambda: f64) -> (f64, f64, f64) {
    let vov = vgs - vt0;
    if vov <= 0.0 {
        return (0.0, 0.0, 0.0);
    }
    let lam = 1.0 + lambda * vds;
    if vds < vov {
        let id = kp * (vov * vds - 0.5 * vds * vds);
        (id * lam, kp * vds * lam, kp * (vov - vds) * lam + id * lambda)
    } else {
        let id = 0.5 * kp * vov * vov;
        (id * lam, kp * vov * lam, id * lambda)
    }
}

// Effective drain/source frame: handles pmos polarity and source/drain swap.
// Returns (ed, es, current into ed, gm, gds); the stamps are sign free.
fn mos_op(dv: &Dev, x: &[f64]) -> (usize, usize, f64, f64, f64) {
    let Dev::M { d, g, s, kp, vt0, lambda, pmos } = dv else { unreachable!() };
    let sg = if *pmos { -1.0 } else { 1.0 };
    let (ed, es) = if sg * (x[*d] - x[*s]) >= 0.0 { (*d, *s) } else { (*s, *d) };
    let vgs = sg * (x[*g] - x[es]);
    let vds = sg * (x[ed] - x[es]);
    let (id, gm, gds) = mos_eval(vgs, vds, *kp, *vt0, *lambda);
    (ed, es, sg * id, gm, gds)
}

// ----------------------------------------------------------------- solver ---

enum Mode<'a> {
    Dc,
    Tran { h: f64, t: f64, st: &'a [[f64; 2]] },
}

fn assemble(ckt: &Circuit, x: &[f64], lim: &mut [f64], mode: &Mode, gmin: f64) -> Sys {
    let m = ckt.nodes.len() + ckt.nb;
    let mut s = Sys::new(m);
    let t_now = if let Mode::Tran { t, .. } = mode { *t } else { 0.0 };
    for (di, dv) in ckt.devs.iter().enumerate() {
        match dv {
            Dev::R { p, n, v } => {
                let g = 1.0 / v;
                s.add(*p, *p, g);
                s.add(*p, *n, -g);
                s.add(*n, *p, -g);
                s.add(*n, *n, g);
            }
            Dev::C { p, n, v } => {
                if let Mode::Tran { h, st, .. } = mode {
                    let geq = 2.0 * v / h;
                    let ieq = -(geq * st[di][0] + st[di][1]);
                    s.add(*p, *p, geq);
                    s.add(*p, *n, -geq);
                    s.add(*n, *p, -geq);
                    s.add(*n, *n, geq);
                    s.rhs(*p, -ieq);
                    s.rhs(*n, ieq);
                }
            }
            Dev::L { p, n, v, br } => {
                s.add(*p, *br, 1.0);
                s.add(*n, *br, -1.0);
                s.add(*br, *p, 1.0);
                s.add(*br, *n, -1.0);
                if let Mode::Tran { h, st, .. } = mode {
                    let req = 2.0 * v / h;
                    s.add(*br, *br, -req);
                    s.rhs(*br, -(req * st[di][0] + st[di][1]));
                }
            }
            Dev::V { p, n, dc, wave, br, .. } => {
                s.add(*p, *br, 1.0);
                s.add(*n, *br, -1.0);
                s.add(*br, *p, 1.0);
                s.add(*br, *n, -1.0);
                s.rhs(*br, wave_val(*dc, wave, t_now));
            }
            Dev::I { p, n, dc, wave, .. } => {
                let val = wave_val(*dc, wave, t_now);
                s.rhs(*p, -val);
                s.rhs(*n, val);
            }
            Dev::D { p, n, is, nf } => {
                let vt = nf * VT;
                let vcrit = vt * (vt / (SQRT_2 * is)).ln();
                let vd = pnjlim(x[*p] - x[*n], lim[di], vt, vcrit);
                lim[di] = vd;
                let (id, gd) = diode_eval(vd, *is, vt, gmin);
                let ieq = id - gd * vd;
                s.add(*p, *p, gd);
                s.add(*p, *n, -gd);
                s.add(*n, *p, -gd);
                s.add(*n, *n, gd);
                s.rhs(*p, -ieq);
                s.rhs(*n, ieq);
            }
            Dev::M { g, .. } => {
                let (ed, es, i0, gm, gds) = mos_op(dv, x);
                let gds = gds + gmin;
                let ieq = i0 - gm * x[*g] - gds * x[ed] + (gm + gds) * x[es];
                s.add(ed, *g, gm);
                s.add(ed, ed, gds);
                s.add(ed, es, -(gm + gds));
                s.add(es, *g, -gm);
                s.add(es, ed, -gds);
                s.add(es, es, gm + gds);
                s.rhs(ed, -ieq);
                s.rhs(es, ieq);
            }
            Dev::E { p, n, cp, cn, k, br } => {
                s.add(*p, *br, 1.0);
                s.add(*n, *br, -1.0);
                s.add(*br, *p, 1.0);
                s.add(*br, *n, -1.0);
                s.add(*br, *cp, -*k);
                s.add(*br, *cn, *k);
            }
            Dev::G { p, n, cp, cn, k } => {
                s.add(*p, *cp, *k);
                s.add(*p, *cn, -*k);
                s.add(*n, *cp, -*k);
                s.add(*n, *cn, *k);
            }
        }
    }
    s
}

fn newton(ckt: &Circuit, x: &mut Vec<f64>, lim: &mut [f64], mode: &Mode, gmin: f64, maxit: usize) -> Option<usize> {
    let nn = ckt.nodes.len();
    for it in 1..=maxit {
        let s = assemble(ckt, x, lim, mode, gmin);
        let xn = lu_solve(s.a, s.b, s.m)?;
        let mut conv = true;
        for i in 1..xn.len() {
            let tol = if i < nn { VNTOL } else { ABSTOL } + RELTOL * xn[i].abs().max(x[i].abs());
            if (xn[i] - x[i]).abs() > tol {
                conv = false;
            }
        }
        *x = xn;
        if conv && it > 1 {
            return Some(it);
        }
    }
    None
}

fn op_solve(ckt: &Circuit, x: &mut Vec<f64>, lim: &mut [f64]) {
    if newton(ckt, x, lim, &Mode::Dc, GMIN, ITL_OP).is_some() {
        return;
    }
    x.iter_mut().for_each(|v| *v = 0.0);
    lim.iter_mut().for_each(|v| *v = 0.0);
    let mut g = 1e-2;
    while g > GMIN {
        if newton(ckt, x, lim, &Mode::Dc, g, ITL_OP).is_none() {
            die("operating point did not converge (gmin stepping failed, floating node?)");
        }
        g /= 10.0;
    }
    if newton(ckt, x, lim, &Mode::Dc, GMIN, ITL_OP).is_none() {
        die("operating point did not converge");
    }
}

fn init_state(ckt: &Circuit, x: &[f64]) -> Vec<[f64; 2]> {
    ckt.devs
        .iter()
        .map(|d| match d {
            Dev::C { p, n, .. } => [x[*p] - x[*n], 0.0],
            Dev::L { br, .. } => [x[*br], 0.0],
            _ => [0.0, 0.0],
        })
        .collect()
}

fn update_state(ckt: &Circuit, x: &[f64], h: f64, st: &mut [[f64; 2]]) {
    for (k, d) in ckt.devs.iter().enumerate() {
        match d {
            Dev::C { p, n, v } => {
                let vn = x[*p] - x[*n];
                let inew = 2.0 * v / h * (vn - st[k][0]) - st[k][1];
                st[k] = [vn, inew];
            }
            Dev::L { v, br, .. } => {
                let inw = x[*br];
                let vnw = 2.0 * v / h * (inw - st[k][0]) - st[k][1];
                st[k] = [inw, vnw];
            }
            _ => {}
        }
    }
}

// ----------------------------------------------------------------- output ---

fn columns(ckt: &Circuit) -> Vec<String> {
    let mut c: Vec<String> = (1..ckt.nodes.len()).map(|i| format!("v({})", ckt.nodes[i])).collect();
    for (k, d) in ckt.devs.iter().enumerate() {
        if matches!(d, Dev::L { .. } | Dev::V { .. } | Dev::E { .. }) {
            c.push(format!("i({})", ckt.names[k]));
        }
    }
    c
}

fn fmt_row(vals: &[f64]) -> String {
    vals.iter().map(|v| format!("{:.9e}", v)).collect::<Vec<_>>().join(",")
}

// --------------------------------------------------------------- analyses ---

fn run_op(ckt: &Circuit) {
    let m = ckt.nodes.len() + ckt.nb;
    let mut x = vec![0.0; m];
    let mut lim = vec![0.0; ckt.devs.len()];
    op_solve(ckt, &mut x, &mut lim);
    out!("# op");
    out!("{}", columns(ckt).join(","));
    out!("{}", fmt_row(&x[1..]));
    out!("");
}

fn run_dc(ckt: &mut Circuit, src: &str, start: f64, stop: f64, step: f64) {
    let di = ckt
        .names
        .iter()
        .position(|n| n == src)
        .unwrap_or_else(|| die(&format!("unknown sweep source '{}'", src)));
    if step == 0.0 {
        die(".dc step must be nonzero");
    }
    let m = ckt.nodes.len() + ckt.nb;
    let mut x = vec![0.0; m];
    let mut lim = vec![0.0; ckt.devs.len()];
    out!("# dc");
    out!("{},{}", src, columns(ckt).join(","));
    let npts = ((stop - start) / step).round() as i64;
    for k in 0..=npts.max(0) {
        let val = start + step * k as f64;
        match &mut ckt.devs[di] {
            Dev::V { dc, .. } | Dev::I { dc, .. } => *dc = val,
            _ => die(".dc sweep source must be a v or i source"),
        }
        op_solve(ckt, &mut x, &mut lim);
        out!("{:.9e},{}", val, fmt_row(&x[1..]));
    }
    out!("");
}

fn run_tran(ckt: &Circuit, tstep: f64, tstop: f64) {
    if tstep <= 0.0 || tstop <= 0.0 {
        die("bad .tran parameters");
    }
    let m = ckt.nodes.len() + ckt.nb;
    let mut x = vec![0.0; m];
    let mut lim = vec![0.0; ckt.devs.len()];
    op_solve(ckt, &mut x, &mut lim);
    let mut st = init_state(ckt, &x);
    out!("# tran");
    out!("time,{}", columns(ckt).join(","));
    out!("{:.9e},{}", 0.0, fmt_row(&x[1..]));
    let hmin = tstop * 1e-10;
    let (mut t, mut h) = (0.0, tstep);
    while t < tstop - hmin {
        let hs = h.min(tstop - t);
        let mut xn = x.clone();
        match newton(ckt, &mut xn, &mut lim, &Mode::Tran { h: hs, t: t + hs, st: &st }, GMIN, ITL_TRAN) {
            Some(iters) => {
                t += hs;
                update_state(ckt, &xn, hs, &mut st);
                x = xn;
                out!("{:.9e},{}", t, fmt_row(&x[1..]));
                if iters <= ITL_TRAN / 2 {
                    h = (2.0 * h).min(tstep);
                }
            }
            None => {
                h /= 8.0;
                if h < hmin {
                    die(&format!("timestep too small at t = {:.3e}", t));
                }
            }
        }
    }
    out!("");
}

fn run_ac(ckt: &Circuit, dec: bool, n: usize, f1: f64, f2: f64) {
    if f1 <= 0.0 || f2 < f1 || n == 0 {
        die("bad .ac parameters");
    }
    let m = ckt.nodes.len() + ckt.nb;
    let mut x = vec![0.0; m];
    let mut lim = vec![0.0; ckt.devs.len()];
    op_solve(ckt, &mut x, &mut lim);
    // The real part of the AC matrix is exactly the DC Jacobian at the OP.
    let sdc = assemble(ckt, &x, &mut lim, &Mode::Dc, GMIN);
    let mut freqs = Vec::new();
    if dec {
        let mut k = 0;
        loop {
            let f = f1 * 10f64.powf(k as f64 / n as f64);
            if f > f2 * (1.0 + 1e-9) {
                break;
            }
            freqs.push(f);
            k += 1;
        }
    } else if n == 1 {
        freqs.push(f1);
    } else {
        for k in 0..n {
            freqs.push(f1 + (f2 - f1) * k as f64 / (n as f64 - 1.0));
        }
    }
    out!("# ac");
    let hdr: Vec<String> = columns(ckt)
        .iter()
        .flat_map(|c| [format!("mag({})", c), format!("ph({})", c)])
        .collect();
    out!("freq,{}", hdr.join(","));
    for f in freqs {
        let w = 2.0 * PI * f;
        let mut a: Vec<Cx> = sdc.a.iter().map(|&v| Cx::new(v, 0.0)).collect();
        let mut b = vec![Cx::zero(); m];
        for dv in &ckt.devs {
            match dv {
                Dev::C { p, n, v } => {
                    let c = w * v;
                    a[*p * m + *p].im += c;
                    a[*p * m + *n].im -= c;
                    a[*n * m + *p].im -= c;
                    a[*n * m + *n].im += c;
                }
                Dev::L { v, br, .. } => a[*br * m + *br].im -= w * v,
                Dev::V { ac, br, .. } => b[*br] = Cx::new(*ac, 0.0),
                Dev::I { p, n, ac, .. } => {
                    b[*p].re -= *ac;
                    b[*n].re += *ac;
                }
                _ => {}
            }
        }
        let xa = clu_solve(a, b, m).unwrap_or_else(|| die("singular matrix in .ac"));
        let mut row = vec![f];
        for v in &xa[1..] {
            row.push(v.abs());
            row.push(v.im.atan2(v.re) * 180.0 / PI);
        }
        out!("{}", fmt_row(&row));
    }
    out!("");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        die("usage: nanospice <netlist.cir>");
    }
    let src = std::fs::read_to_string(&args[1]).unwrap_or_else(|e| die(&format!("cannot read '{}': {}", args[1], e)));
    let mut ckt = parse(&src);
    if ckt.analyses.is_empty() {
        ckt.analyses.push(Analysis::Op);
    }
    for an in ckt.analyses.clone() {
        match an {
            Analysis::Op => run_op(&ckt),
            Analysis::Dc { src: s, start, stop, step } => run_dc(&mut ckt, &s, start, stop, step),
            Analysis::Tran { tstep, tstop } => run_tran(&ckt, tstep, tstop),
            Analysis::Ac { dec, n, f1, f2 } => run_ac(&ckt, dec, n, f1, f2),
        }
    }
}
