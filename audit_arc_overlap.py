"""Detect sustained, excessive parallel bead overlap in arc-fill G-code.

This complements centerline crossing checks. Ordinary side bonding and short
terminal contacts are allowed; support and full toolhead clearance are separate
checks. Widths are commanded preview widths, not a physical extrusion model.
"""
import argparse
from collections import defaultdict
import json
import math

from audit_gcode_moves import WORDS, linearized_gcode


def audit_overlap(stream, target_z, maximum_overlap=.6, chord_error=.001):
    if not math.isfinite(target_z) or not 0 <= maximum_overlap < 1:
        raise ValueError("Expected finite Z and overlap fraction in [0, 1)")
    position = dict.fromkeys("XYZ", 0.)
    width, role = .45, ""
    pending, failures = [], []
    grid = defaultdict(list)
    max_width = 0.
    arc_paths = 0

    def finish():
        nonlocal pending, max_width, arc_paths
        if not pending:
            return
        if pending[0][4] in ("Arc bridge", "Arc overhang"):
            arc_paths += 1
            run = longest = 0.
            other_source = longest_source = None
            for a, b, bead, source, _ in pending:
                dx, dy = b[0]-a[0], b[1]-a[1]
                norm = dx*dx+dy*dy
                length = math.sqrt(norm)
                count = max(1, math.ceil(length / .025))
                reach = .5 * (bead + max_width)
                for sample in range(count):
                    t = (sample+.5)/count
                    x, y = a[0]+t*dx, a[1]+t*dy
                    duplicate = False
                    for gx in range(math.floor(x-reach), math.floor(x+reach)+1):
                        for gy in range(math.floor(y-reach), math.floor(y+reach)+1):
                            for c, d, other_bead, prior_source, _ in grid[gx, gy]:
                                vx, vy = d[0]-c[0], d[1]-c[1]
                                other_norm = vx*vx+vy*vy
                                if (dx*vx+dy*vy)**2 < .82 * norm * other_norm:
                                    continue
                                u = max(0., min(1., ((x-c[0])*vx+(y-c[1])*vy)/other_norm))
                                clearance = .5*(bead+other_bead)-maximum_overlap*min(bead, other_bead)
                                if math.hypot(x-c[0]-u*vx, y-c[1]-u*vy) < clearance:
                                    duplicate = True
                                    other_source = prior_source
                                    break
                            if duplicate:
                                break
                        if duplicate:
                            break
                    run = run+length/count if duplicate else 0.
                    if run > longest:
                        longest, longest_source = run, other_source
            if longest > .65 * pending[0][2]:
                failures.append(dict(line=pending[0][3], prior_line=longest_source,
                                     overlap_run_mm=longest, width_mm=pending[0][2]))
        # A path cannot count itself as previously deposited overlap.
        for segment in pending:
            a, b, bead, _, _ = segment
            max_width = max(max_width, bead)
            for gx in range(math.floor(min(a[0], b[0])), math.floor(max(a[0], b[0]))+1):
                for gy in range(math.floor(min(a[1], b[1])), math.floor(max(a[1], b[1]))+1):
                    grid[gx, gy].append(segment)
        pending = []

    for source, line in linearized_gcode(stream, chord_error, exclude_cleaner_extrusion=True):
        if line.startswith(";TYPE:"):
            finish()
            role = line.strip()[6:]
        if line.startswith(";WIDTH:"):
            width = float(line.strip()[7:])
        opcode = line.split(maxsplit=1)[0] if line.strip() else ""
        if opcode not in ("G0", "G1", "G92"):
            continue
        words = {k: float(v) for k, v in WORDS.findall(line.split(";")[0])}
        a = tuple(position.values())
        position.update((k, words[k]) for k in "XYZ" if k in words)
        b = tuple(position.values())
        if (opcode == "G92" or words.get("E", 0) <= 0 or math.dist(a[:2], b[:2]) < 1e-9
                or abs(a[2]-target_z) > .001 or abs(b[2]-target_z) > .001):
            finish()
        else:
            pending.append((a, b, width, source, role))
    finish()
    return dict(z=target_z, arc_paths=arc_paths, maximum_overlap=maximum_overlap, overlaps=failures)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gcode")
    parser.add_argument("--z", required=True, type=float)
    args = parser.parse_args()
    with open(args.gcode) as stream:
        result = audit_overlap(stream, args.z)
    print(json.dumps(result, indent=2))
    raise SystemExit(bool(result["overlaps"]))
