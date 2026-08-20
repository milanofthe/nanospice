// nanospice: a classic SPICE circuit simulator in a single source file,
// hard budget 1000 lines of code, std only.
//
// Algorithm (Berkeley SPICE structure): MNA with branch currents for V, L, E.
// Newton-Raphson with pnjlim junction limiting and gmin stepping fallback for
// the operating point. Transient with trapezoidal companion models, a
// quadratic polynomial predictor and LTE based timestep control, iteration
// count as nonconvergence fallback. AC as small signal linearization at the
// operating point. Sparse LU with partial pivoting, generic over real and
// complex via the Num trait.
//
// Matrix convention: unknown 0 is ground; stamps may address it freely, the
// sparse assembly simply drops row/column 0, so device stamps never special
// case ground and the solver loops start at index 1.

use std::collections::HashMap;
use std::f64::consts::{PI, SQRT_2};

const GMIN: f64 = 1e-12;
const RELTOL: f64 = 1e-3;
const VNTOL: f64 = 1e-6;
const ABSTOL: f64 = 1e-12;
const TRTOL: f64 = 7.0;
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

trait Num:
    Copy
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<Output = Self>
    + std::ops::Div<Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn mag(self) -> f64;
}

impl Num for f64 {
    fn zero() -> f64 { 0.0 }
    fn one() -> f64 { 1.0 }
    fn mag(self) -> f64 { self.abs() }
}

#[derive(Clone, Copy)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    fn new(re: f64, im: f64) -> Cx { Cx { re, im } }
}

impl std::ops::Add for Cx {
    type Output = Cx;
    fn add(self, o: Cx) -> Cx { Cx::new(self.re + o.re, self.im + o.im) }
}

impl std::ops::Sub for Cx {
    type Output = Cx;
    fn sub(self, o: Cx) -> Cx { Cx::new(self.re - o.re, self.im - o.im) }
}

impl std::ops::Mul for Cx {
    type Output = Cx;
    fn mul(self, o: Cx) -> Cx { Cx::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re) }
}

impl std::ops::Div for Cx {
    type Output = Cx;
    fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        Cx::new((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
}

impl Num for Cx {
    fn zero() -> Cx { Cx::new(0.0, 0.0) }
    fn one() -> Cx { Cx::new(1.0, 0.0) }
    fn mag(self) -> f64 { self.re.hypot(self.im) }
}

// SPICE number: longest parseable prefix plus unit suffix, e.g. 4.7k, 100n, 1meg.
const SUFFIXES: [(&str, f64); 9] = [
    ("meg", 1e6),
    ("t", 1e12),
    ("g", 1e9),
    ("k", 1e3),
    ("m", 1e-3),
    ("u", 1e-6),
    ("n", 1e-9),
    ("p", 1e-12),
    ("f", 1e-15),
];

fn num(tok: &str) -> Option<f64> {
    let t = tok.trim();
    if !t.is_ascii() {
        return None;
    }
    let mut end = t.len();
    while end > 0 && t[..end].parse::<f64>().is_err() {
        end -= 1;
    }
    if end == 0 {
        return None;
    }
    let v: f64 = t[..end].parse().unwrap();
    let mult = SUFFIXES.iter().find(|(s, _)| t[end..].starts_with(s)).map_or(1.0, |(_, m)| *m);
    Some(v * mult)
}

fn numx(tok: &str) -> f64 {
    num(tok).unwrap_or_else(|| die(&format!("bad number '{}'", tok)))
}

// ---------------------------------------------------------------- circuit ---

#[derive(Clone)]
enum Wave {
    Sin { vo: f64, va: f64, fr: f64, td: f64, theta: f64 },
    Pwl { pts: Vec<(f64, f64)>, per: f64 },
}

fn wave_val(dc: f64, wave: &Option<Wave>, t: f64) -> f64 {
    match wave {
        None => dc,
        Some(Wave::Sin { vo, va, fr, td, theta }) => {
            let tp = (t - td).max(0.0);
            vo + va * (-tp * theta).exp() * (2.0 * PI * fr * tp).sin()
        }
        Some(Wave::Pwl { pts, per }) => {
            let t0 = pts[0].0;
            let mut tt = t;
            if *per > 0.0 && tt > t0 {
                tt = t0 + (tt - t0) % per;
            }
            let mut v = pts[0].1;
            for w in pts.windows(2) {
                if tt >= w[1].0 {
                    v = w[1].1;
                } else {
                    if tt > w[0].0 {
                        v = w[0].1 + (w[1].1 - w[0].1) * (tt - w[0].0) / (w[1].0 - w[0].0);
                    }
                    break;
                }
            }
            v
        }
    }
}

#[derive(Clone)]
enum Dev {
    R { p: usize, n: usize, v: f64 },
    C { p: usize, n: usize, v: f64, ic: f64, m: f64 },
    L { p: usize, n: usize, v: f64, br: usize, ic: f64 },
    V { p: usize, n: usize, dc: f64, ac: f64, wave: Option<Wave>, br: usize },
    I { p: usize, n: usize, dc: f64, ac: f64, wave: Option<Wave> },
    D { p: usize, n: usize, is: f64, nf: f64 },
    M { d: usize, g: usize, s: usize, kp: f64, vt0: f64, lambda: f64, pmos: bool },
    E { p: usize, n: usize, cp: usize, cn: usize, k: f64, br: usize },
    G { p: usize, n: usize, cp: usize, cn: usize, k: f64 },
    Q { c: usize, b: usize, e: usize, is: f64, bf: f64, br: f64, pnp: bool },
}

#[derive(Clone)]
enum Analysis {
    Op,
    Dc { src: String, start: f64, stop: f64, step: f64 },
    Tran { tstep: f64, tstop: f64, uic: bool },
    Ac { dec: bool, n: usize, f1: f64, f2: f64 },
}

struct Circuit {
    devs: Vec<Dev>,
    names: Vec<String>,
    nodes: Vec<String>,
    analyses: Vec<Analysis>,
    prints: Vec<String>,
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
    let mut i = 0;
    while i < toks.len() {
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
                let (v1, v2, td) = (g(0), g(1), g(2));
                let (tr, tf) = (g(3).max(1e-12), g(4).max(1e-12));
                let pw = if v.len() > 5 { g(5) } else { f64::MAX / 4.0 };
                let pts = vec![(td, v1), (td + tr, v2), (td + tr + pw, v2), (td + tr + pw + tf, v1)];
                wave = Some(Wave::Pwl { pts, per: g(6) });
            }
            "pwl" => {
                i += 1;
                let v = take(&mut i);
                let pts: Vec<(f64, f64)> = v.chunks(2).filter(|c| c.len() == 2).map(|c| (c[0], c[1])).collect();
                if pts.is_empty() {
                    die("pwl needs time/value pairs");
                }
                wave = Some(Wave::Pwl { pts, per: 0.0 });
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
        let low = raw.split(';').next().unwrap().to_ascii_lowercase();
        let line = low
            .chars()
            .map(|c| if "(),".contains(c) { ' ' } else { c })
            .collect::<String>()
            .trim()
            .to_string();
        if ln == 0 || line.is_empty() || line.starts_with('*') {
            continue;
        }
        if let Some(rest) = line.strip_prefix('+') {
            let prev = lines.last_mut().unwrap_or_else(|| die("continuation line without a preceding line"));
            prev.push_str(" ");
            prev.push_str(rest);
        } else {
            lines.push(line);
        }
    }
    let mut models: HashMap<String, Vec<String>> = HashMap::new();
    for line in &lines {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.first() == Some(&".model") && toks.len() > 2 {
            models.insert(toks[1].into(), toks[2..].iter().map(|s| s.to_string()).collect());
        }
    }
    let mut map: HashMap<String, usize> = HashMap::new();
    map.insert("0".into(), 0);
    map.insert("gnd".into(), 0);
    let mut nodes = vec!["0".to_string()];
    let mut devs: Vec<Dev> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut analyses: Vec<Analysis> = Vec::new();
    let mut prints: Vec<String> = Vec::new();
    for line in &lines {
        let toks: Vec<&str> = line.split_whitespace().collect();
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
                    uic: toks.contains(&"uic"),
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
                "model" => {}
                "print" => {
                    let mut i = 1;
                    while i + 1 < toks.len() {
                        prints.push(format!("{}({})", toks[i], toks[i + 1]));
                        i += 2;
                    }
                }
                _ => eprintln!("nanospice: ignoring card .{}", card),
            }
            continue;
        }
        let expand = |from: usize| -> Vec<&str> {
            let (mut inline, mut modt) = (Vec::new(), Vec::new());
            for t in toks.iter().skip(from) {
                match models.get(*t) {
                    Some(mt) => modt.extend(mt.iter().map(|s| s.as_str())),
                    None => inline.push(*t),
                }
            }
            inline.extend(modt);
            inline
        };
        let mut nid = |s: &str| -> usize {
            if let Some(&i) = map.get(s) {
                return i;
            }
            let i = nodes.len();
            map.insert(s.into(), i);
            nodes.push(s.into());
            i
        };
        let mut caps: Vec<(usize, usize, f64, f64)> = Vec::new();
        let dev = match toks[0].chars().next().unwrap() {
            'r' => Dev::R { p: nid(tok(&toks, 1)), n: nid(tok(&toks, 2)), v: numx(tok(&toks, 3)) },
            'c' => Dev::C {
                p: nid(tok(&toks, 1)),
                n: nid(tok(&toks, 2)),
                v: numx(tok(&toks, 3)),
                ic: pval(&toks[4..], "ic", f64::NAN),
                m: pval(&toks[4..], "m", 0.0),
            },
            'l' => Dev::L {
                p: nid(tok(&toks, 1)),
                n: nid(tok(&toks, 2)),
                v: numx(tok(&toks, 3)),
                br: 0,
                ic: pval(&toks[4..], "ic", f64::NAN),
            },
            'v' => {
                let (dc, ac, wave) = src_spec(&toks[3..]);
                Dev::V { p: nid(tok(&toks, 1)), n: nid(tok(&toks, 2)), dc, ac, wave, br: 0 }
            }
            'i' => {
                let (dc, ac, wave) = src_spec(&toks[3..]);
                Dev::I { p: nid(tok(&toks, 1)), n: nid(tok(&toks, 2)), dc, ac, wave }
            }
            'd' => {
                let ext = expand(3);
                let (p, n) = (nid(tok(&toks, 1)), nid(tok(&toks, 2)));
                caps.push((p, n, pval(&ext, "cjo", 0.0), pval(&ext, "mj", 0.5)));
                Dev::D { p, n, is: pval(&ext, "is", 1e-14), nf: pval(&ext, "n", 1.0) }
            }
            // the jfet is the square law device again: kp = 2 beta and a
            // depletion threshold, so it shares the mosfet variant
            c @ ('m' | 'j') => {
                let ext = expand(if c == 'm' { 5 } else { 4 });
                let (nd, ng, ns) = (nid(tok(&toks, 1)), nid(tok(&toks, 2)), nid(tok(&toks, 3)));
                caps.push((ng, ns, pval(&ext, "cgs", 0.0), 0.0));
                caps.push((ng, nd, pval(&ext, "cgd", 0.0), 0.0));
                Dev::M {
                    d: nd,
                    g: ng,
                    s: ns,
                    kp: if c == 'm' { pval(&ext, "kp", 2e-5) } else { 2.0 * pval(&ext, "beta", 1e-4) },
                    vt0: pval(&ext, "vt0", pval(&ext, "vto", if c == 'm' { 0.0 } else { -2.0 })),
                    lambda: pval(&ext, "lambda", 0.0),
                    pmos: ext.iter().any(|t| *t == "pmos" || *t == "pjf"),
                }
            }
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
            'q' => {
                let ext = expand(4);
                let (nc, nb, ne) = (nid(tok(&toks, 1)), nid(tok(&toks, 2)), nid(tok(&toks, 3)));
                let pnp = ext.iter().any(|t| *t == "pnp");
                let (je, jc) = if pnp { ((ne, nb), (nc, nb)) } else { ((nb, ne), (nb, nc)) };
                caps.push((je.0, je.1, pval(&ext, "cje", 0.0), pval(&ext, "mj", 0.5)));
                caps.push((jc.0, jc.1, pval(&ext, "cjc", 0.0), pval(&ext, "mj", 0.5)));
                Dev::Q {
                    c: nc,
                    b: nb,
                    e: ne,
                    is: pval(&ext, "is", 1e-16),
                    bf: pval(&ext, "bf", 100.0),
                    br: pval(&ext, "br", 1.0),
                    pnp,
                }
            }
            _ => die(&format!("unknown device '{}'", toks[0])),
        };
        devs.push(dev);
        names.push(toks[0].to_string());
        // junction and gate capacitances desugar into plain capacitors
        for (a, b, cv, mg) in caps {
            if cv > 0.0 {
                devs.push(Dev::C { p: a, n: b, v: cv, ic: f64::NAN, m: mg });
                names.push(format!("c.{}", toks[0]));
            }
        }
    }
    let nn = nodes.len();
    let mut nb = 0;
    for d in devs.iter_mut() {
        if let Dev::L { br, .. } | Dev::V { br, .. } | Dev::E { br, .. } = d {
            *br = nn + nb;
            nb += 1;
        }
    }
    Circuit { devs, names, nodes, analyses, prints, nb }
}

// -------------------------------------------------------------- sparse LU ---

// Rows are kept sorted by column. Eliminated columns are removed from active
// rows, so the column k entry of any row at index >= k is always its first
// entry; pivot search and elimination need no extra index structures.
struct Sys<T> {
    rows: Vec<Vec<(usize, T)>>,
    b: Vec<T>,
}

impl<T: Num> Sys<T> {
    fn new(m: usize) -> Sys<T> {
        Sys { rows: vec![Vec::new(); m], b: vec![T::zero(); m] }
    }
    fn add(&mut self, i: usize, j: usize, v: T) {
        if i == 0 || j == 0 {
            return;
        }
        match self.rows[i].binary_search_by_key(&j, |e| e.0) {
            Ok(k) => self.rows[i][k].1 = self.rows[i][k].1 + v,
            Err(k) => self.rows[i].insert(k, (j, v)),
        }
    }
    fn rhs(&mut self, i: usize, v: T) {
        if i > 0 {
            self.b[i] = self.b[i] + v;
        }
    }
    // Conductance quad: rows p/n against columns cp/cn (cp = p, cn = n for a
    // two terminal conductance; the general form is the VCCS stamp).
    fn quad(&mut self, p: usize, n: usize, cp: usize, cn: usize, g: T) {
        self.add(p, cp, g);
        self.add(p, cn, T::zero() - g);
        self.add(n, cp, T::zero() - g);
        self.add(n, cn, g);
    }
    // Voltage-defined branch: current column and voltage row for V, L, E.
    fn branch(&mut self, p: usize, n: usize, br: usize) {
        self.add(p, br, T::one());
        self.add(n, br, T::zero() - T::one());
        self.add(br, p, T::one());
        self.add(br, n, T::zero() - T::one());
    }
    // Verilog-A style contribution: current i flows p -> n, linearized at
    // the (possibly limited) evaluation point and shifted to the raw
    // iterate; ctl lists control pairs with d i / d (x[a] - x[b]).
    fn contrib(&mut self, p: usize, n: usize, i: T, ctl: &[(usize, usize, T)], x: &[T]) {
        let mut ieq = i;
        for &(a, b, g) in ctl {
            self.quad(p, n, a, b, g);
            ieq = ieq - g * (x[a] - x[b]);
        }
        self.rhs(p, T::zero() - ieq);
        self.rhs(n, ieq);
    }
    // The dual voltage contribution for V, L, E: branch row
    // v(p) - v(n) - sum g (x[a] - x[b]) = v.
    fn vcontrib(&mut self, p: usize, n: usize, br: usize, v: T, ctl: &[(usize, usize, T)]) {
        self.branch(p, n, br);
        for &(a, b, g) in ctl {
            self.add(br, a, T::zero() - g);
            self.add(br, b, g);
        }
        self.rhs(br, v);
    }
}

// row a minus f times pivot row p, both sorted; fill-in falls out of the merge
fn merge<T: Num>(a: &[(usize, T)], p: &[(usize, T)], f: T) -> Vec<(usize, T)> {
    let mut o = Vec::with_capacity(a.len() + p.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < p.len() {
        if a[i].0 == p[j].0 {
            o.push((a[i].0, a[i].1 - f * p[j].1));
            i += 1;
            j += 1;
        } else if a[i].0 < p[j].0 {
            o.push(a[i]);
            i += 1;
        } else {
            o.push((p[j].0, T::zero() - f * p[j].1));
            j += 1;
        }
    }
    o.extend_from_slice(&a[i..]);
    for &(c, v) in &p[j..] {
        o.push((c, T::zero() - f * v));
    }
    o
}

fn solve<T: Num>(mut s: Sys<T>) -> Result<Vec<T>, usize> {
    let m = s.b.len();
    for k in 1..m {
        let (mut p, mut best) = (0, 0.0);
        for i in k..m {
            if let Some(&(c, v)) = s.rows[i].first() {
                if c == k && v.mag() > best {
                    best = v.mag();
                    p = i;
                }
            }
        }
        if best < 1e-300 {
            return Err(k);
        }
        s.rows.swap(k, p);
        s.b.swap(k, p);
        let prow = std::mem::take(&mut s.rows[k]);
        let piv = prow[0].1;
        for i in k + 1..m {
            if s.rows[i].first().map_or(false, |e| e.0 == k) {
                let f = s.rows[i][0].1 / piv;
                s.b[i] = s.b[i] - f * s.b[k];
                s.rows[i] = merge(&s.rows[i][1..], &prow[1..], f);
            }
        }
        s.rows[k] = prow;
    }
    let mut x = s.b;
    for i in (1..m).rev() {
        let mut sum = x[i];
        for &(j, v) in &s.rows[i][1..] {
            sum = sum - v * x[j];
        }
        x[i] = sum / s.rows[i][0].1;
    }
    x[0] = T::zero();
    Ok(x)
}

// ---------------------------------------------------------- device models ---

// Charge and capacitance of the graded junction: q(v) integrates
// c(v) = c0 (1 - v/vj)^-m below fc*vj, linearized above; m = 0 is the
// plain linear capacitor. vj = 1 V and fc = 0.5 are fixed.
fn ceval(v: f64, c0: f64, m: f64) -> (f64, f64) {
    let (vj, fc) = (1.0, 0.5);
    if m == 0.0 {
        return (c0 * v, c0);
    }
    let qj = |v: f64| c0 * vj / (1.0 - m) * (1.0 - (1.0 - v / vj).powf(1.0 - m));
    if v < fc * vj {
        (qj(v), c0 * (1.0 - v / vj).powf(-m))
    } else {
        let (cf, dv) = (c0 * (1.0 - fc).powf(-m), v - fc * vj);
        let slope = cf * m / (vj * (1.0 - fc));
        (qj(fc * vj) + cf * dv + 0.5 * slope * dv * dv, cf + slope * dv)
    }
}

fn diode_eval(vd: f64, is: f64, vt: f64, gmin: f64) -> (f64, f64) {
    let e = (vd / vt).min(200.0).exp();
    (is * (e - 1.0) + gmin * vd, is * e / vt + gmin)
}

// Limited junction voltage: pnjlim plus the bookkeeping shared by D and Q.
fn junction(raw: f64, lim: &mut f64, vt: f64, is: f64, limited: &mut bool) -> f64 {
    let vcrit = vt * (vt / (SQRT_2 * is)).ln();
    let v = pnjlim(raw, *lim, vt, vcrit);
    *limited |= (v - raw).abs() > VNTOL;
    *lim = v;
    v
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
    Dc { fac: f64 },
    Tran { h: f64, t: f64, st: &'a [[f64; 2]], be: bool },
}

fn assemble(ckt: &Circuit, x: &[f64], lim: &mut [[f64; 2]], mode: &Mode, gmin: f64, limited: &mut bool) -> Sys<f64> {
    let mut s = Sys::new(ckt.nodes.len() + ckt.nb);
    let t_now = if let Mode::Tran { t, .. } = mode { *t } else { 0.0 };
    let fac = if let Mode::Dc { fac } = mode { *fac } else { 1.0 };
    for (di, dv) in ckt.devs.iter().enumerate() {
        match dv {
            Dev::R { p, n, v } => s.contrib(*p, *n, (x[*p] - x[*n]) / v, &[(*p, *n, 1.0 / v)], x),
            Dev::C { p, n, v, m: mg, .. } => {
                if let Mode::Tran { h, st, be, .. } = mode {
                    let k = if *be { 1.0 } else { 2.0 };
                    let (q, c) = ceval(x[*p] - x[*n], *v, *mg);
                    let inow = k / h * (q - st[di][0]) - if *be { 0.0 } else { st[di][1] };
                    s.contrib(*p, *n, inow, &[(*p, *n, k * c / h)], x);
                }
            }
            Dev::L { p, n, v, br, .. } => match mode {
                Mode::Tran { h, st, be, .. } => {
                    let k = if *be { 1.0 } else { 2.0 };
                    let vh = k / h * st[di][0] + if *be { 0.0 } else { st[di][1] };
                    s.vcontrib(*p, *n, *br, -vh, &[(*br, 0, k * v / h)]);
                }
                _ => s.vcontrib(*p, *n, *br, 0.0, &[]),
            },
            Dev::V { p, n, dc, wave, br, .. } => {
                s.vcontrib(*p, *n, *br, fac * wave_val(*dc, wave, t_now), &[]);
            }
            Dev::I { p, n, dc, wave, .. } => {
                s.contrib(*p, *n, fac * wave_val(*dc, wave, t_now), &[], x);
            }
            Dev::D { p, n, is, nf } => {
                let vt = nf * VT;
                let vr = x[*p] - x[*n];
                let vd = junction(vr, &mut lim[di][0], vt, *is, limited);
                let (id, gd) = diode_eval(vd, *is, vt, gmin);
                s.contrib(*p, *n, id + gd * (vr - vd), &[(*p, *n, gd)], x);
            }
            Dev::M { g, .. } => {
                let (ed, es, i0, gm, gds) = mos_op(dv, x);
                s.contrib(ed, es, i0, &[(*g, es, gm), (ed, es, gds + gmin)], x);
            }
            Dev::E { p, n, cp, cn, k, br } => s.vcontrib(*p, *n, *br, 0.0, &[(*cp, *cn, *k)]),
            Dev::G { p, n, cp, cn, k } => s.contrib(*p, *n, *k * (x[*cp] - x[*cn]), &[(*cp, *cn, *k)], x),
            // Ebers-Moll transport form as three contributions: the two
            // junctions and the transport source, KCL composes the stamp.
            Dev::Q { c, b, e, is, bf, br, pnp } => {
                let sg = if *pnp { -1.0 } else { 1.0 };
                let (vbr, vcr) = (sg * (x[*b] - x[*e]), sg * (x[*b] - x[*c]));
                let vbe = junction(vbr, &mut lim[di][0], VT, *is, limited);
                let vbc = junction(vcr, &mut lim[di][1], VT, *is, limited);
                let ef = (vbe / VT).min(200.0).exp();
                let er = (vbc / VT).min(200.0).exp();
                let (gmf, gmr) = (is * ef / VT, is * er / VT);
                let (gpi, gmu) = (gmf / bf + gmin, gmr / br + gmin);
                let (dbe, dbc) = (vbr - vbe, vcr - vbc);
                let ibe = is * (ef - 1.0) / bf + gmin * vbe;
                let ibc = is * (er - 1.0) / br + gmin * vbc;
                let ict = is * (ef - er);
                s.contrib(*b, *e, sg * (ibe + gpi * dbe), &[(*b, *e, gpi)], x);
                s.contrib(*b, *c, sg * (ibc + gmu * dbc), &[(*b, *c, gmu)], x);
                s.contrib(*c, *e, sg * (ict + gmf * dbe - gmr * dbc), &[(*b, *e, gmf), (*b, *c, -gmr)], x);
            }
        }
    }
    s
}

fn branches(ckt: &Circuit) -> impl Iterator<Item = &String> {
    let is_br = |d: &&Dev| matches!(d, Dev::L { .. } | Dev::V { .. } | Dev::E { .. });
    ckt.devs.iter().zip(&ckt.names).filter(move |(d, _)| is_br(d)).map(|(_, n)| n)
}

fn unk_name(ckt: &Circuit, i: usize) -> String {
    if i < ckt.nodes.len() {
        return format!("node {}", ckt.nodes[i]);
    }
    match branches(ckt).nth(i - ckt.nodes.len()) {
        Some(n) => format!("branch {}", n),
        None => format!("unknown {}", i),
    }
}

fn tolv(i: usize, nn: usize, a: f64, b: f64) -> f64 {
    (if i < nn { VNTOL } else { ABSTOL }) + RELTOL * a.abs().max(b.abs())
}

fn newton(ckt: &Circuit, x: &mut Vec<f64>, lim: &mut [[f64; 2]], mode: &Mode, gmin: f64, maxit: usize) -> Option<usize> {
    let nn = ckt.nodes.len();
    for it in 1..=maxit {
        // convergence is denied while junction limiting is still active,
        // otherwise a clamped junction looks like a converged solution
        let mut limited = false;
        let s = assemble(ckt, x, lim, mode, gmin, &mut limited);
        let xn =
            solve(s).unwrap_or_else(|k| die(&format!("singular matrix at {}, no conduction path?", unk_name(ckt, k))));
        let conv = (1..xn.len()).all(|i| (xn[i] - x[i]).abs() <= tolv(i, nn, xn[i], x[i]));
        *x = xn;
        if conv && !limited && it > 1 {
            return Some(it);
        }
    }
    None
}

fn op_solve(ckt: &Circuit, x: &mut Vec<f64>, lim: &mut [[f64; 2]]) {
    if newton(ckt, x, lim, &Mode::Dc { fac: 1.0 }, GMIN, ITL_OP).is_some() {
        return;
    }
    // continuation fallbacks: each path is a list of (gmin, source scale)
    // stages walked with warm starts; gmin stepping, then source stepping
    let gmin_path: Vec<(f64, f64)> = (2..=12).map(|k| (10f64.powi(-k), 1.0)).collect();
    let src_path: Vec<(f64, f64)> = (1..=10).map(|k| (GMIN, k as f64 / 10.0)).collect();
    for path in [gmin_path, src_path] {
        x.iter_mut().for_each(|v| *v = 0.0);
        lim.iter_mut().for_each(|v| *v = [0.0; 2]);
        if path.iter().all(|&(g, fac)| newton(ckt, x, lim, &Mode::Dc { fac }, g, ITL_OP).is_some()) {
            return;
        }
    }
    die("operating point did not converge (gmin and source stepping failed)");
}

fn op_setup(ckt: &Circuit) -> (usize, Vec<f64>, Vec<[f64; 2]>) {
    let m = ckt.nodes.len() + ckt.nb;
    let mut x = vec![0.0; m];
    let mut lim = vec![[0.0; 2]; ckt.devs.len()];
    op_solve(ckt, &mut x, &mut lim);
    (m, x, lim)
}

// Reactive state per device: [charge, current] for C, [flux, voltage] for
// L; both duals advance by the same divided difference.
fn init_state(ckt: &Circuit, x: &[f64], uic: bool) -> Vec<[f64; 2]> {
    let pick = |ic: f64, live: f64| if uic && ic.is_finite() { ic } else { live };
    ckt.devs
        .iter()
        .map(|d| match d {
            Dev::C { p, n, v, ic, m: mg } => [ceval(pick(*ic, x[*p] - x[*n]), *v, *mg).0, 0.0],
            Dev::L { v, br, ic, .. } => [v * pick(*ic, x[*br]), 0.0],
            _ => [0.0, 0.0],
        })
        .collect()
}

fn update_state(ckt: &Circuit, x: &[f64], h: f64, st: &mut [[f64; 2]], be: bool) {
    let (f, m) = if be { (1.0, 0.0) } else { (2.0, 1.0) };
    for (k, d) in ckt.devs.iter().enumerate() {
        let q = match d {
            Dev::C { p, n, v, m: mg, .. } => ceval(x[*p] - x[*n], *v, *mg).0,
            Dev::L { v, br, .. } => v * x[*br],
            _ => continue,
        };
        st[k] = [q, f / h * (q - st[k][0]) - m * st[k][1]];
    }
}

// ----------------------------------------------------------------- output ---

fn columns(ckt: &Circuit) -> Vec<String> {
    let v = (1..ckt.nodes.len()).map(|i| format!("v({})", ckt.nodes[i]));
    v.chain(branches(ckt).map(|n| format!("i({})", n))).collect()
}

// Printed columns and their unknown indices; empty .print means everything.
fn selection(ckt: &Circuit) -> (Vec<String>, Vec<usize>) {
    let all = columns(ckt);
    for p in &ckt.prints {
        if !all.contains(p) {
            eprintln!("nanospice: unknown .print item {}", p);
        }
    }
    let (mut names, mut idx) = (Vec::new(), Vec::new());
    for (i, c) in all.into_iter().enumerate() {
        if ckt.prints.is_empty() || ckt.prints.contains(&c) {
            names.push(c);
            idx.push(i + 1);
        }
    }
    (names, idx)
}

fn pick(x: &[f64], idx: &[usize]) -> Vec<f64> {
    idx.iter().map(|&i| x[i]).collect()
}

fn fmt_row(vals: &[f64]) -> String {
    vals.iter().map(|v| format!("{:.9e}", v)).collect::<Vec<_>>().join(",")
}

// --------------------------------------------------------------- analyses ---

fn run_op(ckt: &Circuit) {
    let (_, x, _) = op_setup(ckt);
    let (names, idx) = selection(ckt);
    out!("# op");
    out!("{}", names.join(","));
    out!("{}", fmt_row(&pick(&x, &idx)));
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
    let mut lim = vec![[0.0; 2]; ckt.devs.len()];
    let (names, idx) = selection(ckt);
    out!("# dc");
    out!("{},{}", src, names.join(","));
    let npts = ((stop - start) / step).round() as i64;
    for k in 0..=npts.max(0) {
        let val = start + step * k as f64;
        match &mut ckt.devs[di] {
            // the sweep owns the source: any transient wave is dropped
            Dev::V { dc, wave, .. } | Dev::I { dc, wave, .. } => {
                *dc = val;
                *wave = None;
            }
            _ => die(".dc sweep source must be a v or i source"),
        }
        op_solve(ckt, &mut x, &mut lim);
        out!("{:.9e},{}", val, fmt_row(&pick(&x, &idx)));
    }
    out!("");
}

// Waveform corner breakpoints: transient steps must land on them.
fn breakpoints(ckt: &Circuit, tstop: f64) -> Vec<f64> {
    let mut bp = Vec::new();
    for d in &ckt.devs {
        if let Dev::V { wave: Some(Wave::Pwl { pts, per }), .. }
        | Dev::I { wave: Some(Wave::Pwl { pts, per }), .. } = d
        {
            let mut off = 0.0;
            while pts[0].0 + off < tstop {
                for &(tc, _) in pts {
                    let c = tc + off;
                    if c > 0.0 && c < tstop {
                        bp.push(c);
                    }
                }
                if *per <= 0.0 {
                    break;
                }
                off += per;
            }
        }
    }
    bp.sort_by(f64::total_cmp);
    bp.dedup();
    bp
}

fn run_tran(ckt: &Circuit, tstep: f64, tstop: f64, uic: bool) {
    if tstep <= 0.0 || tstop <= 0.0 {
        die("bad .tran parameters");
    }
    let m = ckt.nodes.len() + ckt.nb;
    let (mut x, mut lim) = (vec![0.0; m], vec![[0.0; 2]; ckt.devs.len()]);
    if !uic {
        op_solve(ckt, &mut x, &mut lim);
    }
    let nn = ckt.nodes.len();
    let mut st = init_state(ckt, &x, uic);
    let (names, idx) = selection(ckt);
    out!("# tran");
    out!("time,{}", names.join(","));
    out!("{:.9e},{}", 0.0, fmt_row(&pick(&x, &idx)));
    let hmin = tstop * 1e-12;
    // start small like classic SPICE; the LTE controller grows the step fast
    let (mut t, mut h) = (0.0, tstep / 100.0);
    // history for the quadratic predictor, most recent first: d1 = step from
    // old[0] to the current x, d2 = step from old[1] to old[0]
    let mut old: Vec<(Vec<f64>, f64)> = Vec::new();
    let bp = breakpoints(ckt, tstop);
    let mut bpi = 0;
    while t < tstop - hmin {
        while bpi < bp.len() && bp[bpi] <= t + hmin {
            bpi += 1;
        }
        let mut hs = h.min(tstop - t);
        let on_bp = bpi < bp.len() && t + hs >= bp[bpi] - hmin;
        if on_bp {
            hs = bp[bpi] - t;
        }
        let pred: Option<Vec<f64>> = (old.len() == 2).then(|| {
            let (d1, d2) = (old[0].1, old[1].1);
            let l0 = (hs + d1) * (hs + d1 + d2) / (d1 * (d1 + d2));
            let l1 = -hs * (hs + d1 + d2) / (d1 * d2);
            let l2 = hs * (hs + d1) / ((d1 + d2) * d2);
            (0..m).map(|i| l0 * x[i] + l1 * old[0].0[i] + l2 * old[1].0[i]).collect()
        });
        // the predictor doubles as the Newton starting point; the two
        // predictor-less bootstrap steps run damped backward euler
        let be = old.len() < 2;
        let mut xn = pred.clone().unwrap_or_else(|| x.clone());
        let mode = Mode::Tran { h: hs, t: t + hs, st: &st, be };
        let Some(iters) = newton(ckt, &mut xn, &mut lim, &mode, GMIN, ITL_TRAN) else {
            h = hs / 8.0;
            if h < hmin {
                die(&format!("timestep too small at t = {:.3e}", t));
            }
            continue;
        };
        let mut hnew = if iters <= ITL_TRAN / 2 { (2.0 * hs).min(tstep) } else { hs };
        if let Some(xp) = &pred {
            // Milne style LTE estimate: the corrector minus predictor gap is
            // predictor error minus trapezoid error, both proportional to the
            // third derivative, so scale the gap back to the trapezoid share.
            let (d1, d2) = (old[0].1, old[1].1);
            let etrap = hs * hs * hs / 12.0;
            let fac = etrap / (hs * (hs + d1) * (hs + d1 + d2) / 6.0 + etrap);
            let mut r: f64 = 0.0;
            for i in 1..m {
                r = r.max(fac * (xn[i] - xp[i]).abs() / (TRTOL * tolv(i, nn, xn[i], x[i])));
            }
            let scale = 0.9 / r.max(1e-8).cbrt();
            if r > 1.0 && hs > 20.0 * hmin {
                h = hs * scale.clamp(0.125, 0.9);
                continue;
            }
            hnew = (hs * scale.clamp(0.3, 2.0)).min(tstep);
        }
        t += hs;
        update_state(ckt, &xn, hs, &mut st, be);
        old.insert(0, (x, hs));
        old.truncate(2);
        x = xn;
        out!("{:.9e},{}", t, fmt_row(&pick(&x, &idx)));
        if on_bp {
            // waveform derivative is discontinuous here: restart small and
            // drop the predictor history
            old.clear();
            h = tstep / 100.0;
        } else {
            h = hnew;
        }
    }
    out!("");
}

fn run_ac(ckt: &Circuit, dec: bool, n: usize, f1: f64, f2: f64) {
    if f1 <= 0.0 || f2 < f1 || n == 0 {
        die("bad .ac parameters");
    }
    let (m, x, mut lim) = op_setup(ckt);
    // The real part of the AC matrix is exactly the DC Jacobian at the OP.
    let sdc = assemble(ckt, &x, &mut lim, &Mode::Dc { fac: 1.0 }, GMIN, &mut false);
    let freqs: Vec<f64> = if dec {
        (0..).map(|k| f1 * 10f64.powf(k as f64 / n as f64)).take_while(|f| *f <= f2 * (1.0 + 1e-9)).collect()
    } else if n == 1 {
        vec![f1]
    } else {
        (0..n).map(|k| f1 + (f2 - f1) * k as f64 / (n as f64 - 1.0)).collect()
    };
    let (names, idx) = selection(ckt);
    out!("# ac");
    let hdr: Vec<String> = names
        .iter()
        .flat_map(|c| [format!("mag({})", c), format!("ph({})", c)])
        .collect();
    out!("freq,{}", hdr.join(","));
    for f in freqs {
        let w = 2.0 * PI * f;
        let mut sc = Sys::<Cx> {
            rows: sdc.rows.iter().map(|r| r.iter().map(|&(j, v)| (j, Cx::new(v, 0.0))).collect()).collect(),
            b: vec![Cx::new(0.0, 0.0); m],
        };
        for dv in &ckt.devs {
            match dv {
                Dev::C { p, n, v, m: mg, .. } => {
                    sc.quad(*p, *n, *p, *n, Cx::new(0.0, w * ceval(x[*p] - x[*n], *v, *mg).1))
                }
                Dev::L { v, br, .. } => sc.add(*br, *br, Cx::new(0.0, -w * v)),
                Dev::V { ac, br, .. } => sc.rhs(*br, Cx::new(*ac, 0.0)),
                Dev::I { p, n, ac, .. } => {
                    sc.rhs(*p, Cx::new(-*ac, 0.0));
                    sc.rhs(*n, Cx::new(*ac, 0.0));
                }
                _ => {}
            }
        }
        let xa = solve(sc)
            .unwrap_or_else(|k| die(&format!("singular matrix in .ac at {}", unk_name(ckt, k))));
        let mut row = vec![f];
        for &i in &idx {
            row.push(xa[i].mag());
            row.push(xa[i].im.atan2(xa[i].re) * 180.0 / PI);
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
            Analysis::Tran { tstep, tstop, uic } => run_tran(&ckt, tstep, tstop, uic),
            Analysis::Ac { dec, n, f1, f2 } => run_ac(&ckt, dec, n, f1, f2),
        }
    }
}
