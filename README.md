# nanospice

A classic SPICE circuit simulator in a single Rust source file with a hard
budget of 1000 lines of code. No dependencies, std only. The budget is
enforced by a test (`loc_budget` in `tests/cli.rs`), counting nonblank,
noncomment lines of `src/main.rs`.

## Build and run

    cargo build --release
    target/release/nanospice circuit.cir

Output is CSV on stdout, one block per analysis, each prefixed with a
`# <analysis>` comment line.

## Supported

Analyses: `.op`, `.dc src start stop step`, `.tran tstep tstop`,
`.ac dec|lin n fstart fstop`.

Devices:

| Card | Device |
|---|---|
| `Rxxx p n value` | resistor |
| `Cxxx p n value` | capacitor |
| `Lxxx p n value` | inductor |
| `Vxxx p n [dc v] [ac mag] [sin(vo va f td theta)] [pulse(v1 v2 td tr tf pw per)]` | voltage source |
| `Ixxx p n ...` | current source, same spec as V |
| `Dxxx p n [is=1e-14] [n=1]` | diode, Shockley |
| `Mxxx d g s b nmos\|pmos [kp=2e-5] [vt0=0] [lambda=0]` | MOSFET level 1, bulk ignored, vt0 is the threshold magnitude for pmos too |
| `Exxx p n cp cn gain` | VCVS |
| `Gxxx p n cp cn gm` | VCCS |

Netlist rules: the first line is a title, `*` starts a comment line, `;` starts
a trailing comment, `+` continues the previous line, node `0` or `gnd` is
ground, unit suffixes `t g meg k m u n p f`.

## Algorithm

Classic Berkeley SPICE structure: MNA with branch currents for V, L and E,
Newton-Raphson with pnjlim junction limiting and gmin stepping for the
operating point, trapezoidal integration with companion models and iteration
count timestep control for transient, small signal AC linearized at the
operating point. Dense LU with partial pivoting instead of sparse Markowitz is
the one deliberate departure from the classic implementation, which caps
practical circuit size at a few hundred nodes.

## Not supported

Subcircuits, `.model` cards, `.param`, BJT, noise analysis, sparse matrices.
That is the price of the budget.

## Examples and tests

Example netlists live in `circuits/`. The integration tests in `tests/cli.rs`
check the solver against analytic references (RC and RL step response, RC
corner frequency, diode KCL consistency, MOSFET bias point) and enforce the
LOC budget:

    cargo test
