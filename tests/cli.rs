use std::process::Command;

fn run(name: &str, netlist: &str) -> Vec<Vec<String>> {
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    std::fs::write(&path, netlist).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_nanospice")).arg(&path).output().unwrap();
    assert!(out.status.success(), "stderr: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split(',').map(|s| s.to_string()).collect())
        .collect()
}

fn table(rows: &[Vec<String>], tag: &str) -> (Vec<String>, Vec<Vec<f64>>) {
    let i = rows.iter().position(|r| r[0] == format!("# {}", tag)).unwrap();
    let header = rows[i + 1].clone();
    let mut data = Vec::new();
    for r in &rows[i + 2..] {
        if r[0].starts_with('#') {
            break;
        }
        data.push(r.iter().map(|s| s.parse::<f64>().unwrap()).collect());
    }
    (header, data)
}

fn col(h: &[String], name: &str) -> usize {
    h.iter().position(|c| c == name).unwrap_or_else(|| panic!("no column {} in {:?}", name, h))
}

#[test]
fn op_divider() {
    let rows = run("div.cir", "divider\nv1 in 0 dc 10\nr1 in out 1k\nr2 out 0 1k\n.op\n.end\n");
    let (h, d) = table(&rows, "op");
    assert!((d[0][col(&h, "v(out)")] - 5.0).abs() < 1e-9);
    assert!((d[0][col(&h, "i(v1)")] + 5e-3).abs() < 1e-9);
}

#[test]
fn op_diode() {
    let rows = run("diode.cir", "diode\nv1 a 0 dc 5\nr1 a b 1k\nd1 b 0\n.op\n.end\n");
    let (h, d) = table(&rows, "op");
    let vb = d[0][col(&h, "v(b)")];
    assert!(vb > 0.55 && vb < 0.8, "vd = {}", vb);
    let ir = (5.0 - vb) / 1e3;
    let id = 1e-14 * ((vb / 0.02585).exp() - 1.0);
    assert!((ir - id).abs() / ir < 5e-2, "kcl mismatch: {} vs {}", ir, id);
}

#[test]
fn op_mos() {
    let rows = run(
        "mos.cir",
        "mos common source\nvdd vdd 0 dc 5\nvg g 0 dc 2\nrd vdd d 2k\nm1 d g 0 0 nmos kp=2e-3 vt0=1\n.op\n.end\n",
    );
    let (h, d) = table(&rows, "op");
    // id = 0.5 * kp * (vgs - vt0)^2 = 1 mA, v(d) = 5 - 2k * 1mA = 3 V
    assert!((d[0][col(&h, "v(d)")] - 3.0).abs() < 1e-4);
}

#[test]
fn op_model_card() {
    // n=2 diode via .model card: drop is twice the n=1 drop, so this fails
    // if the model reference is not resolved
    let rows = run(
        "model.cir",
        "model card\n.model dx d is=1e-15 n=2\nv1 a 0 dc 5\nr1 a b 1k\nd1 b 0 dx\n.op\n.end\n",
    );
    let (h, d) = table(&rows, "op");
    let vb = d[0][col(&h, "v(b)")];
    assert!(vb > 1.3 && vb < 1.7, "vd = {}", vb);
    let ir = (5.0 - vb) / 1e3;
    let id = 1e-15 * ((vb / (2.0 * 0.02585)).exp() - 1.0);
    assert!((ir - id).abs() / ir < 5e-2, "kcl mismatch: {} vs {}", ir, id);
}

#[test]
fn op_bjt() {
    // forward active npn: base current set by rb, check ic = bf * ib
    let npn = run(
        "bjt_npn.cir",
        "npn ce\nvcc vcc 0 dc 5\nrb vcc b 430k\nrc vcc c 1k\nq1 c b 0 npn is=1e-16 bf=100\n.op\n.end\n",
    );
    let (h, d) = table(&npn, "op");
    let (vb, vc) = (d[0][col(&h, "v(b)")], d[0][col(&h, "v(c)")]);
    assert!(vb > 0.7 && vb < 0.85, "vb = {}", vb);
    let ib = (5.0 - vb) / 430e3;
    let ic = (5.0 - vc) / 1e3;
    assert!((ic - 100.0 * ib).abs() / ic < 1e-3, "beta mismatch: ic={} ib={}", ic, ib);
    // pnp mirror image of the same circuit must be exactly symmetric
    let pnp = run(
        "bjt_pnp.cir",
        "pnp ce\nvee vee 0 dc -5\nrb vee b 430k\nrc vee c 1k\nq1 c b 0 pnp is=1e-16 bf=100\n.op\n.end\n",
    );
    let (h2, d2) = table(&pnp, "op");
    assert!((d2[0][col(&h2, "v(c)")] + vc).abs() < 1e-6);
    assert!((d2[0][col(&h2, "v(b)")] + vb).abs() < 1e-6);
}

#[test]
fn op_jfet() {
    // n-jfet at vgs = 0 with default beta and vto: id = beta * vto^2 = 0.4 mA
    let rows = run("jfet.cir", "jfet\nvdd vdd 0 dc 10\nrd vdd d 5k\nj1 d 0 0\n.op\n.end\n");
    let (h, d) = table(&rows, "op");
    assert!((d[0][col(&h, "v(d)")] - 8.0).abs() < 1e-4);
}

#[test]
fn dc_sweep() {
    let rows = run("sweep.cir", "sweep\nv1 in 0 dc 0\nr1 in out 1k\nr2 out 0 1k\n.dc v1 0 10 1\n.end\n");
    let (h, d) = table(&rows, "dc");
    assert_eq!(d.len(), 11);
    for r in &d {
        assert!((r[col(&h, "v(out)")] - r[0] / 2.0).abs() < 1e-9);
    }
}

#[test]
fn tran_rc() {
    let rows = run(
        "rc.cir",
        "rc step\nv1 in 0 pulse 0 1 0 1n 1n 1 2\nr1 in out 1k\nc1 out 0 1u\n.tran 10u 5m\n.end\n",
    );
    let (h, d) = table(&rows, "tran");
    assert!(d.len() > 100);
    let c = col(&h, "v(out)");
    for r in &d {
        let refv = 1.0 - (-r[0] / 1e-3).exp();
        assert!((r[c] - refv).abs() < 1e-2, "t={} v={} ref={}", r[0], r[c], refv);
    }
}

#[test]
fn tran_rl() {
    let rows = run(
        "rl.cir",
        "rl step\nv1 in 0 pulse 0 1 0 1n 1n 1 2\nr1 in a 1\nl1 a 0 1m\n.tran 10u 5m\n.end\n",
    );
    let (h, d) = table(&rows, "tran");
    let c = col(&h, "i(l1)");
    for r in &d {
        let refi = 1.0 - (-r[0] / 1e-3).exp();
        assert!((r[c] - refi).abs() < 1e-2, "t={} i={} ref={}", r[0], r[c], refi);
    }
}

#[test]
fn tran_lte_fast_tau() {
    // tau = 10u is ten times smaller than tstep, so a fixed tstep integration
    // would ring badly; the LTE controller has to resolve the edge on its own
    let rows = run(
        "lte.cir",
        "lte\nv1 in 0 pulse 0 1 0 1n 1n 1 2\nr1 in out 1k\nc1 out 0 10n\n.tran 100u 1m\n.end\n",
    );
    let (h, d) = table(&rows, "tran");
    let c = col(&h, "v(out)");
    for r in &d {
        let refv = 1.0 - (-r[0] / 1e-5).exp();
        assert!((r[c] - refv).abs() < 2e-2, "t={} v={} ref={}", r[0], r[c], refv);
    }
}

#[test]
fn op_ladder() {
    // 201 equal resistors in series exercise the sparse solver and pivoting
    let mut nl = String::from("ladder\nv1 n0 0 dc 10\n");
    for k in 0..200 {
        nl.push_str(&format!("r{} n{} n{} 1\n", k, k, k + 1));
    }
    nl.push_str("rload n200 0 1\n.op\n.end\n");
    let rows = run("ladder.cir", &nl);
    let (h, d) = table(&rows, "op");
    assert!((d[0][col(&h, "v(n100)")] - 10.0 * 101.0 / 201.0).abs() < 1e-9);
    assert!((d[0][col(&h, "i(v1)")] + 10.0 / 201.0).abs() < 1e-9);
}

#[test]
fn tran_uic_lc() {
    // undriven lc tank started by an inductor initial current: amplitude is
    // i0 * sqrt(l/c) = 31.62 mV, only reachable with uic (the op is all zero)
    let rows = run(
        "lc.cir",
        "lc tank\nc1 a 0 1u\nl1 a 0 1m ic=1m\n.tran 5u 200u uic\n.end\n",
    );
    let (h, d) = table(&rows, "tran");
    let c = col(&h, "v(a)");
    let peak = d.iter().map(|r| r[c].abs()).fold(0.0, f64::max);
    assert!((peak - 0.03162).abs() < 1e-3, "peak = {}", peak);
}

#[test]
fn ac_rc() {
    let rows = run(
        "acrc.cir",
        "ac rc corner\nv1 in 0 dc 0 ac 1\nr1 in out 1k\nc1 out 0 1u\n.ac lin 1 159.1549431 159.1549431\n.end\n",
    );
    let (h, d) = table(&rows, "ac");
    // at f = 1/(2 pi R C): |H| = 1/sqrt(2), phase = -45 deg
    assert!((d[0][col(&h, "mag(v(out))")] - 0.70711).abs() < 1e-3);
    assert!((d[0][col(&h, "ph(v(out))")] + 45.0).abs() < 0.1);
}

#[test]
fn print_selection() {
    let rows = run(
        "print.cir",
        "print\nv1 in 0 dc 10\nr1 in out 1k\nr2 out 0 1k\n.print v(out) i(v1)\n.op\n.end\n",
    );
    let (h, d) = table(&rows, "op");
    assert_eq!(h, vec!["v(out)", "i(v1)"]);
    assert!((d[0][0] - 5.0).abs() < 1e-9);
    assert!((d[0][1] + 5e-3).abs() < 1e-9);
}

fn run_fail(name: &str, netlist: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    std::fs::write(&path, netlist).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_nanospice")).arg(&path).output().unwrap();
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(!out.status.success(), "expected failure for {}", name);
    assert!(err.contains("nanospice:"), "no clean error for {}: {}", name, err);
    assert!(!err.contains("panicked"), "panicked on {}: {}", name, err);
    err
}

#[test]
fn clean_errors() {
    run_fail("f1.cir", "t\nzz1 a b 1\n.op\n.end\n");
    run_fail("f2.cir", "t\nr1 a 0 abc\n.op\n.end\n");
    run_fail("f3.cir", "t\nr1 a 0\n.op\n.end\n");
    run_fail("f4.cir", "t\nv1 a 0 dc 5\nr1 a 0 1k\n.dc vx 0 1 0.1\n.end\n");
    run_fail("f5.cir", "t\nr1 a 0 1\u{b5}\n.op\n.end\n");
    run_fail("f6.cir", "t\n.dc\n.end\n");
    // a node reachable only through a capacitor has no dc path; the error
    // must name it
    let e = run_fail("f7.cir", "t\nv1 a 0 dc 1\nr1 a 0 1k\nc1 a b 1u\n.op\n.end\n");
    assert!(e.contains("node b"), "missing node name: {}", e);
}

#[test]
fn op_random_ladder() {
    // lcg-randomized resistor ladder checked against a thevenin reduction
    let mut seed: u64 = 42;
    let mut rnd = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        1.0 + (seed >> 33) as f64 / 4.3e9 * 99.0
    };
    let n = 30;
    let (mut rs, mut rp) = (Vec::new(), Vec::new());
    let mut nl = String::from("random ladder\nv1 n0 0 dc 10\n");
    for i in 0..n {
        let (a, b) = (rnd(), rnd());
        nl.push_str(&format!("rs{} n{} n{} {}\n", i, i, i + 1, a));
        nl.push_str(&format!("rp{} n{} 0 {}\n", i, i + 1, b));
        rs.push(a);
        rp.push(b);
    }
    nl.push_str(".op\n.end\n");
    let mut req = vec![0.0; n];
    for i in (0..n).rev() {
        let right = if i + 1 < n { rs[i + 1] + req[i + 1] } else { f64::INFINITY };
        req[i] = 1.0 / (1.0 / rp[i] + 1.0 / right);
    }
    let rows = run("randlad.cir", &nl);
    let (h, d) = table(&rows, "op");
    let mut v = 10.0;
    for i in 0..n {
        v *= req[i] / (rs[i] + req[i]);
        let got = d[0][col(&h, &format!("v(n{})", i + 1))];
        assert!((got - v).abs() < 1e-6 * v + 1e-15, "node {}: {} vs {}", i + 1, got, v);
    }
}

#[test]
fn tran_lc_energy() {
    // 20 periods: the trapezoid rule must not damp the oscillation
    let rows = run("lce.cir", "lc energy\nc1 a 0 1u\nl1 a 0 1m ic=1m\n.tran 5u 4m uic\n.end\n");
    let (h, d) = table(&rows, "tran");
    let c = col(&h, "v(a)");
    let late = d.iter().filter(|r| r[0] > 3.6e-3).map(|r| r[c].abs()).fold(0.0, f64::max);
    assert!((late - 0.03162).abs() < 1.5e-3, "late peak = {}", late);
}

#[test]
fn op_latch() {
    // cross coupled nmos latch, a classic hard operating point
    let rows = run(
        "latch.cir",
        "latch\nvdd vdd 0 dc 5\nr1 vdd a 10k\nr2 vdd b 10k\nm1 a b 0 0 nmos kp=1e-3 vt0=1\nm2 b a 0 0 nmos kp=1e-3 vt0=1\n.op\n.end\n",
    );
    let (h, d) = table(&rows, "op");
    for name in ["v(a)", "v(b)"] {
        let v = d[0][col(&h, name)];
        assert!((-1e-6..=5.0 + 1e-6).contains(&v), "{} = {}", name, v);
    }
}

#[test]
fn loc_budget() {
    let src = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/src/main.rs")).unwrap();
    let n = src
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with("//")
        })
        .count();
    assert!(n <= 1000, "src/main.rs has {} LOC, budget is 1000", n);
}
