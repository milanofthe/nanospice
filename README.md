# nanospice

A classic SPICE circuit simulator in one Rust source file, capped at 1000
lines of code. No dependencies, std only. A test counts the nonblank,
noncomment lines of src/main.rs and fails above 1000. Current count: 992.

The repository is educational. The report in report/ derives every algorithm
in the simulator, explains the design decisions, and maps both to the code
section by section.

## Build and run

    cargo build --release
    target/release/nanospice circuit.cir

Output is CSV on stdout, one block per analysis. Example netlists are in
circuits/.

## Scope

Analyses: .op, .dc source sweep, .tran with adaptive timestep and uic, .ac
dec/lin. Cards: .model, .print, .end.

| Card | Device |
|---|---|
| `Rxxx p n value` | resistor |
| `Cxxx p n value [ic=v] [m=0]` | capacitor, graded junction for m > 0 |
| `Lxxx p n value [ic=i]` | inductor |
| `Vxxx p n [dc v] [ac mag] [sin(vo va f td theta)] [pulse(v1 v2 td tr tf pw per)] [pwl(t1 v1 t2 v2 ...)]` | voltage source |
| `Ixxx p n ...` | current source, same spec as V |
| `Dxxx p n [model] [is=1e-14] [n=1] [cjo=0] [mj=0.5]` | diode |
| `Mxxx d g s b nmos\|pmos [model] [kp=2e-5] [vt0=0] [lambda=0] [cgs=0] [cgd=0]` | MOSFET level 1 |
| `Jxxx d g s njf\|pjf [model] [beta=1e-4] [vto=-2] [lambda=0] [cgs=0] [cgd=0]` | JFET, square law |
| `Qxxx c b e npn\|pnp [model] [is=1e-16] [bf=100] [br=1] [cje=0] [cjc=0] [mj=0.5]` | BJT, Ebers-Moll |
| `Exxx p n cp cn gain` | VCVS |
| `Gxxx p n cp cn gm` | VCCS |

The first line is the title. `*` starts a comment line, `;` a trailing
comment, `+` continues the previous line. Node 0 or gnd is ground. Unit
suffixes: t g meg k m u n p f. Instance parameters override .model
parameters. Branch current i(vx) is measured into the positive terminal.
With uic the operating point is skipped and ic= values seed the state; the
t=0 row is then the zero vector. `.print v(out) i(v1)` selects output
columns, default is everything.

The MOSFET bulk node is parsed and ignored; vt0 is the threshold magnitude
for pmos as well. Junction capacitances desugar into internal graded
capacitors with C(v) = cjo (1 - v)^-m and the standard fc = 0.5
linearization in forward bias; gate capacitances are linear. Not
supported: subcircuits, .param, noise analysis.

## Algorithms

MNA with branch currents for V, L and E. All device stamps go through two
Verilog-A style contribution primitives, one for currents and one for
voltage-defined branches. Newton-Raphson with pnjlim
junction limiting; gmin stepping and source stepping as operating point
fallbacks.
Transient: trapezoidal companion models, quadratic predictor, LTE timestep
control, waveform breakpoints with damped backward Euler restart steps. AC: small-signal linearization at the operating
point. Linear solver: sparse LU with partial pivoting, generic over real and
complex. Derivations are in report/nanospice.pdf.

## Tests

`cargo test` runs 23 integration tests in tests/cli.rs against the built
binary: analytic references (RC and RL step response, RC corner frequency,
LC amplitude and energy conservation, MOSFET, JFET and BJT bias points), an
npn/pnp symmetry check, a randomized resistor ladder verified against a
Thevenin reduction computed in the test, error handling fuzz cases, and the
LOC budget guard. Tests and comments do not count toward the budget.

## Report

    tectonic report/nanospice.tex

The prebuilt PDF is committed at report/nanospice.pdf.

## License

MIT.
