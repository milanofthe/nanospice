#!/usr/bin/env python3
# Regenerates the .dat files for the report and talk plots. Run from the
# repository root after a release build:
#   python3 report/data/gen.py [path/to/nanospice]
import os
import subprocess
import sys
import tempfile
import time

BIN = sys.argv[1] if len(sys.argv) > 1 else "target/release/nanospice"
OUT = os.path.dirname(os.path.abspath(__file__))


def run(netlist):
    with tempfile.NamedTemporaryFile("w", suffix=".cir", delete=False) as f:
        f.write(netlist)
        path = f.name
    out = subprocess.run([BIN, path], capture_output=True, text=True, check=True).stdout
    os.unlink(path)
    return out


def table(out, tag):
    lines = [l for l in out.splitlines() if l.strip()]
    i = next(k for k, l in enumerate(lines) if l == "# " + tag)
    data = []
    for l in lines[i + 2:]:
        if l.startswith("#"):
            break
        data.append([float(x) for x in l.split(",")])
    return data


def write(name, header, rows):
    with open(os.path.join(OUT, name), "w") as f:
        f.write(header + "\n")
        for r in rows:
            f.write(" ".join(f"{x:.9e}" for x in r) + "\n")
    print(name, len(rows), "rows")


# rc step response, t in ms
d = table(run("rc\nv1 in 0 pulse 0 1 0 1n 1n 1 2\nr1 in out 1k\nc1 out 0 1u\n"
              ".print v(out)\n.tran 10u 5m\n.end\n"), "tran")
write("rc.dat", "t v", [[r[0] * 1e3, r[1]] for r in d])

# cmos inverter dc transfer
d = table(run("cmos inverter\nvdd vdd 0 dc 5\nvin in 0 dc 0\n"
              "m1 out in 0 0 nmos kp=2e-4 vt0=1 lambda=0.05\n"
              "m2 out in vdd vdd pmos kp=2e-4 vt0=1 lambda=0.05\n"
              ".print v(out)\n.dc vin 0 5 0.05\n.end\n"), "dc")
write("inverter.dat", "vin vout", d)

# common emitter amplifier, gain in dB over frequency
d = table(run("common emitter amplifier\nvcc vcc 0 dc 12\nvin in 0 dc 0 ac 1\n"
              "cin in b 10u\nr1 vcc b 47k\nr2 b 0 10k\nrc vcc c 4.7k\nre e 0 1k\n"
              "ce e 0 100u\nq1 c b e npn is=1e-15 bf=200\n"
              ".print v(c)\n.ac dec 10 10 10meg\n.end\n"), "ac")
import math
write("ce_bode.dat", "f db", [[r[0], 20.0 * math.log10(r[1])] for r in d])

# ring oscillator, t in ns
d = table(run("cmos ring oscillator\n"
              ".model n nmos kp=1e-3 vt0=1 cgs=5p cgd=2p\n"
              ".model p pmos kp=1e-3 vt0=1 cgs=5p cgd=2p\n"
              "vdd vdd 0 dc 5\nm1 b a 0 0 n\nm2 b a vdd vdd p\n"
              "m3 c b 0 0 n\nm4 c b vdd vdd p\nm5 a c 0 0 n\nm6 a c vdd vdd p\n"
              "i1 0 a pulse(0 1m 0 1n 1n 20n)\n.print v(a)\n.tran 1n 500n\n.end\n"),
          "tran")
write("ring.dat", "t v", [[r[0] * 1e9, r[1]] for r in d])

# lc tank from initial conditions, t in us, v in mV
d = table(run("lc tank\nc1 a 0 1u\nl1 a 0 1m ic=1m\n.print v(a)\n"
              ".tran 5u 600u uic\n.end\n"), "tran")
write("lc.dat", "t v", [[r[0] * 1e6, r[1] * 1e3] for r in d])

# operating point scaling on rc ladders, wall clock of the whole process
rows = []
for n in [100, 300, 1000, 3000, 10000, 30000]:
    nl = ["ladder", "v1 n0 0 dc 10"]
    for i in range(n):
        nl.append(f"rs{i} n{i} n{i+1} 1")
        nl.append(f"rp{i} n{i+1} 0 1")
    nl += [".print v(n1)", ".op", ".end", ""]
    with tempfile.NamedTemporaryFile("w", suffix=".cir", delete=False) as f:
        f.write("\n".join(nl))
        path = f.name
    best = min(
        (lambda t0: (subprocess.run([BIN, path], capture_output=True, check=True),
                     time.perf_counter() - t0)[1])(time.perf_counter())
        for _ in range(5 if n <= 3000 else 3)
    )
    os.unlink(path)
    rows.append([n, best * 1e3])
    print("bench", n, f"{best*1e3:.2f} ms")
write("bench.dat", "n ms", rows)
