"""Measure consecutive concentric arc families in emitted G-code.

This measures pattern coherence, not support, collisions, or model coverage.
Use alongside audit_arc_support.py and audit_arc_gcode.py.
"""
import argparse
import json
import math

from audit_gcode_moves import WORDS, linearized_gcode


def audit_families(stream, target_z, center_tolerance=.05):
    position = dict.fromkeys("XYZ", 0.)
    role, points, groups = "", [], []
    primary = None

    def finish():
        nonlocal points, primary
        if len(points) < 2:
            points = []
            return
        length = sum(math.dist(a, b) for a, b in zip(points, points[1:]))
        center, radius = None, 0.
        if len(points) >= 3:
            a, b, c = (points[i] for i in (0, len(points) // 3, 2 * len(points) // 3))
            bx, by, cx, cy = b[0]-a[0], b[1]-a[1], c[0]-a[0], c[1]-a[1]
            determinant = 2 * (bx*cy-by*cx)
            if abs(determinant) > 1e-10:
                bb, cc = bx*bx+by*by, cx*cx+cy*cy
                center = (a[0]+(cy*bb-by*cc)/determinant, a[1]+(bx*cc-cx*bb)/determinant)
                radius = math.dist(a, center)
                if any(abs(math.dist(p, center)-radius) > center_tolerance for p in points):
                    center = None
        if center is not None and primary is None:
            primary = center
        if (center is not None and groups and groups[-1]["center"] is not None
                and math.dist(center, groups[-1]["center"]) <= center_tolerance
                and radius >= groups[-1]["radii"][-1] - center_tolerance):
            group = groups[-1]
            if radius > group["radii"][-1] + center_tolerance:
                group["radii"].append(radius)
            group["paths"] += 1
            group["length_mm"] += length
        else:
            groups.append(dict(center=center, radii=[radius], paths=1, length_mm=length))
        points = []

    for _, line in linearized_gcode(stream, .001, exclude_cleaner_extrusion=True):
        if line.startswith(";TYPE:"):
            finish()
            role = line.strip()[6:]
        opcode = line.split(maxsplit=1)[0] if line.strip() else ""
        if opcode not in ("G0", "G1", "G92"):
            continue
        words = {k: float(v) for k, v in WORDS.findall(line.split(";")[0])}
        a = tuple(position.values())
        position.update((k, words[k]) for k in "XYZ" if k in words)
        b = tuple(position.values())
        if (opcode == "G92" or words.get("E", 0) <= 0 or math.dist(a, b) < 1e-9
                or role not in ("Arc bridge", "Arc overhang")
                or abs(a[2]-target_z) > .001 or abs(b[2]-target_z) > .001):
            finish()
        else:
            if not points:
                points.append(a[:2])
            points.append(b[:2])
    finish()
    recursive = [g for g in groups if g["center"] is None or primary is None
                 or math.dist(g["center"], primary) > center_tolerance]
    total = sum(g["length_mm"] for g in recursive)
    grouped = sum(g["length_mm"] for g in recursive if len(g["radii"]) >= 3)
    return dict(z=target_z, paths=sum(g["paths"] for g in groups),
                recursive_families=len(recursive),
                single_ring_families=sum(len(g["radii"]) == 1 for g in recursive),
                recursive_length_mm=total, multi_ring_length_mm=grouped,
                multi_ring_fraction=grouped/total if total else 0.)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gcode")
    parser.add_argument("--z", type=float, required=True)
    args = parser.parse_args()
    with open(args.gcode) as stream:
        print(json.dumps(audit_families(stream, args.z), indent=2))
