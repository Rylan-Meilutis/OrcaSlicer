"""Audit arc contact against extrusion deposited before each path.

This is a geometric contact diagnostic, not a cooling/adhesion simulation.
Both ends and the entire body are checked before adding ANY segment of the
current path. Endpoint-only checking misses outer arcs printed before inner arcs.
"""
import argparse
from collections import defaultdict
import json
import math

from audit_gcode_moves import WORDS, linearized_gcode


def audit_support(stream, target_z, chord_error=0.0001, max_support_drop=0.25, check_bodies=True,
                  min_overlap=0.):
    if not math.isfinite(target_z) or not math.isfinite(max_support_drop) or max_support_drop <= 0:
        raise ValueError("Z and maximum support drop must be finite; drop must be positive")
    if not math.isfinite(min_overlap) or min_overlap < 0:
        raise ValueError("Minimum overlap must be finite and nonnegative")
    grid = defaultdict(list)
    position = dict.fromkeys("XYZ", 0.)
    width, role = 0.45, ""
    pending = []
    failures = []
    body_failures = []
    short_paths = []
    path_count = 0

    def finish():
        nonlocal pending, path_count
        if not pending:
            return
        if pending[0][4] in ("Arc bridge", "Arc overhang") and abs(pending[0][0][2] - target_z) < .001:
            path_count += 1
            path_length = sum(math.dist(a, b) for a, b, *_ in pending)
            if path_length < 1.5 * pending[0][2]:
                # Informational: a tiny grounded path may be valid. Geometry
                # and deposited-prefix contact must still be checked separately.
                short_paths.append({"line": pending[0][3], "length_mm": path_length,
                                    "width_mm": pending[0][2]})
            for part in (("start", "end", "body") if check_bodies else ("start", "end")):
                reverse = part == "end"
                lead = (path_length
                        if part == "body" else .75 * pending[0][2])
                unsupported = 0.
                first_unsupported = None
                for a, b, bead, source, _ in (reversed(pending) if reverse else pending):
                    if reverse:
                        a, b = b, a
                    length = math.dist(a, b)
                    run = min(lead, length)
                    count = max(1, math.ceil(run / .05))
                    for index in range(count + 1):
                        t = run * index / count / length
                        p = tuple(a[k] + t * (b[k] - a[k]) for k in range(3))
                        contact = False
                        for x in range(math.floor(p[0] - 1), math.floor(p[0] + 1) + 1):
                            for y in range(math.floor(p[1] - 1), math.floor(p[1] + 1) + 1):
                                for c, d, other_width, first, last in grid[x, y]:
                                    # Same layer or immediate underlying course.
                                    if not -.001 <= p[2] - d[2] <= max_support_drop:
                                        continue
                                    delta = (d[0] - c[0], d[1] - c[1])
                                    norm = delta[0] ** 2 + delta[1] ** 2
                                    u = 0 if norm == 0 else sum((p[k] - c[k]) * delta[k] for k in range(2)) / norm
                                    # End-cap proximity is not lateral bonding:
                                    # it permits a succession of free tips to
                                    # walk outward with no bead beside them.
                                    if part != "body" and abs(p[2] - d[2]) <= .001 and ((first and u < 0) or (last and u > 1)):
                                        continue
                                    u = max(0, min(1, u))
                                    reach = .5 * (bead + other_width) - min_overlap + .01
                                    if math.hypot(p[0] - c[0] - u * delta[0], p[1] - c[1] - u * delta[1]) <= reach:
                                        contact = True
                                        break
                                if contact:
                                    break
                            if contact:
                                break
                        if not contact:
                            unsupported += run / (count + 1)
                            if first_unsupported is None:
                                first_unsupported = p
                    lead -= run
                    if lead < 1e-8:
                        break
                if unsupported > .01:
                    length_key = "unsupported_body_mm" if part == "body" else "unsupported_lead_mm"
                    (body_failures if part == "body" else failures).append({"line": pending[0][3], "end": part, length_key: unsupported,
                                     "first_unsupported_xyz": first_unsupported,
                                     "xy": (pending[-1][1] if reverse else pending[0][0])[:2]})
        closed = math.dist(pending[0][0], pending[-1][1]) < .005
        for index, (a, b, bead, _, _) in enumerate(pending):
            for x in range(math.floor(min(a[0], b[0])), math.floor(max(a[0], b[0])) + 1):
                for y in range(math.floor(min(a[1], b[1])), math.floor(max(a[1], b[1])) + 1):
                    grid[x, y].append((a, b, bead, not closed and index == 0,
                                       not closed and index == len(pending) - 1))
        pending = []

    for source, line in linearized_gcode(stream, chord_error, exclude_cleaner_extrusion=True):
        if line.startswith(";WIDTH:"):
            width = float(line.split(":")[1])
        if line.startswith(";TYPE:"):
            finish()
            role = line.strip()[6:]
        opcode = line.split(maxsplit=1)[0] if line.strip() else ""
        if opcode not in ("G0", "G1", "G92"):
            continue
        words = {key: float(value) for key, value in WORDS.findall(line.split(";")[0])}
        a = tuple(position.values())
        position.update((key, words[key]) for key in "XYZ" if key in words)
        b = tuple(position.values())
        if opcode == "G92" or words.get("E", 0) <= 0 or math.dist(a, b) < 1e-9:
            finish()
            continue
        if target_z - max_support_drop <= b[2] <= target_z + .001:
            pending.append((a, b, width, source, role))
    finish()
    return {"z": target_z, "max_support_drop": max_support_drop, "min_overlap": min_overlap,
            "arc_paths": path_count, "unsupported_ends": failures,
            "unsupported_bodies": body_failures, "short_paths": short_paths}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gcode")
    parser.add_argument("--z", type=float, required=True)
    parser.add_argument("--max-support-drop", type=float, default=0.25,
                        help="maximum supporting-layer separation; default suits 0.2 mm layers")
    parser.add_argument("--check-bodies", action="store_true", default=True,
                        help="require contact along each arc (the default)")
    parser.add_argument("--endpoints-only", action="store_false", dest="check_bodies",
                        help="legacy diagnostic; does not validate arc growth order")
    parser.add_argument("--min-overlap", type=float, default=0.,
                        help="required lateral bead overlap in mm, beyond mere contact")
    args = parser.parse_args()
    with open(args.gcode) as stream:
        result = audit_support(stream, args.z, max_support_drop=args.max_support_drop,
                               check_bodies=args.check_bodies, min_overlap=args.min_overlap)
    print(json.dumps(result, indent=2))
    raise SystemExit(bool(result["unsupported_ends"] or result["unsupported_bodies"]))
