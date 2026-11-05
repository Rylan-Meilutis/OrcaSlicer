#!/usr/bin/env python3
"""Remove all non-planar-owned extrusion and audit the conventional remainder.

The filtered G-code retains every XYZ move so machine state remains faithful, but
removes E from spatial extrusion moves owned by the non-planar scheduler.  This
is intentionally broader than filtering preview roles: fixed-height dependency
courses may retain an Inner/Outer wall role while still belonging to the
surface-following schedule.
"""

import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


NUMBER = r"-?(?:\d+(?:\.\d*)?|\.\d+)"
AXIS_RE = re.compile(r"([XYZEF])(" + NUMBER + r")")
NONPLANAR_COMMENTS = (
    "non-planar perimeter",
    "non-planar top surface",
    "non-planar ironing",
    "surface-following feature perimeter course",
    "surface-following feature fill course",
)
EXEMPT_SUPPORT_TYPES = (
    "bridge",
    "support",
    "overhang",
)
SURFACE_REMAINDER_TYPES = (
    "inner wall",
    "outer wall",
    "top surface",
    "bottom surface",
    "brick wall",
)


def is_owned(feature: str, comment: str) -> bool:
    lower_feature = feature.lower()
    lower_comment = comment.lower()
    return (
        "non-planar" in lower_feature
        or any(marker in lower_comment for marker in NONPLANAR_COMMENTS)
    )


def segment_distance_xy(point, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    denominator = dx * dx + dy * dy
    t = 0.0 if denominator == 0.0 else max(
        0.0,
        min(1.0, ((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / denominator),
    )
    q = (a[0] + t * dx, a[1] + t * dy, a[2] + t * (b[2] - a[2]))
    return math.hypot(point[0] - q[0], point[1] - q[1]), point[2] - q[2]


def cells_for(a, b, padding=0.55):
    low_x, high_x = sorted((a[0], b[0]))
    low_y, high_y = sorted((a[1], b[1]))
    for ix in range(math.floor(low_x - padding), math.floor(high_x + padding) + 1):
        for iy in range(math.floor(low_y - padding), math.floor(high_y + padding) + 1):
            yield ix, iy


def without_extrusion(line: str) -> str:
    code, separator, comment = line.rstrip("\n").partition(";")
    code = re.sub(r"(?:^|\s)E" + NUMBER, "", code)
    code = re.sub(r"[ \t]+", " ", code).rstrip()
    suffix = (" ;" + comment) if separator else ""
    return code + suffix + " ; non-planar extrusion removed\n"


def collect_extrusion_segments(lines, owned_only=False):
    """Index extrusion before auditing so G-code order cannot hide contact."""
    position = {axis: 0.0 for axis in "XYZE"}
    relative_e = False
    feature = ""
    buckets = defaultdict(list)
    for line in lines:
        if line.startswith(";TYPE:"):
            feature = line[6:].strip()
        command, _, comment = line.partition(";")
        command = command.strip()
        if command == "M82":
            relative_e = False
            continue
        if command == "M83":
            relative_e = True
            continue
        values = {axis: float(value) for axis, value in AXIS_RE.findall(command)}
        if command.startswith("G92"):
            if "E" in values:
                position["E"] = values["E"]
            continue
        if not (command.startswith("G0") or command.startswith("G1")):
            continue
        old = position["X"], position["Y"], position["Z"]
        for axis in "XYZ":
            if axis in values:
                position[axis] = values[axis]
        new = position["X"], position["Y"], position["Z"]
        new_e = values.get("E", position["E"])
        e_delta = new_e if relative_e else new_e - position["E"]
        if not relative_e and "E" in values:
            position["E"] = new_e
        if ("E" in values and e_delta > 1e-8 and math.dist(old, new) > 0.01
                and (not owned_only or is_owned(feature, comment))):
            for cell in cells_for(old, new, 0.0):
                buckets[cell].append((old, new))
    return buckets


def main() -> int:
    if len(sys.argv) not in (2, 3, 4):
        print(
            f"usage: {sys.argv[0]} INPUT.gcode [FILTERED.gcode] [PLANAR_CONTROL.gcode]",
            file=sys.stderr,
        )
        return 2

    source = Path(sys.argv[1])
    destination = Path(sys.argv[2]) if len(sys.argv) == 3 else Path("/tmp/nonplanar-stripped.gcode")
    if len(sys.argv) == 4:
        destination = Path(sys.argv[2])
    control = Path(sys.argv[3]) if len(sys.argv) == 4 else None
    lines = source.read_text(errors="replace").splitlines(keepends=True)
    position = {axis: 0.0 for axis in "XYZE"}
    relative_e = False
    feature = ""
    buckets = defaultdict(list)
    removed_buckets = collect_extrusion_segments(lines, owned_only=True)
    control_buckets = (
        collect_extrusion_segments(control.read_text(errors="replace").splitlines(keepends=True))
        if control is not None else defaultdict(list)
    )
    counts = Counter()
    unsupported_samples = []
    run = None
    runs = []

    with destination.open("w") as output_stream:
        for line_number, line in enumerate(lines, 1):
            if line.startswith(";TYPE:"):
                if run is not None:
                    runs.append(run)
                    run = None
                feature = line[6:].strip()

            command, _, comment = line.partition(";")
            command = command.strip()
            if command == "M82":
                relative_e = False
                output_stream.write(line)
                continue
            if command == "M83":
                relative_e = True
                output_stream.write(line)
                continue
            values = {axis: float(value) for axis, value in AXIS_RE.findall(command)}
            if command.startswith("G92"):
                if "E" in values:
                    position["E"] = values["E"]
                output_stream.write(line)
                continue
            if not (command.startswith("G0") or command.startswith("G1")):
                output_stream.write(line)
                continue

            old = position["X"], position["Y"], position["Z"]
            for axis in "XYZ":
                if axis in values:
                    position[axis] = values[axis]
            new = position["X"], position["Y"], position["Z"]
            new_e = values.get("E", position["E"])
            e_delta = new_e if relative_e else new_e - position["E"]
            if not relative_e and "E" in values:
                position["E"] = new_e
            spatial = math.dist(old, new) > 0.01
            extruding = "E" in values and e_delta > 1e-8 and spatial
            owned = extruding and is_owned(feature, comment)

            if owned:
                counts[("removed", feature)] += 1
                if run is not None:
                    runs.append(run)
                    run = None
                output_stream.write(without_extrusion(line))
                continue

            output_stream.write(line)
            if not extruding:
                if spatial and run is not None:
                    runs.append(run)
                    run = None
                continue

            midpoint = tuple((old[index] + new[index]) * 0.5 for index in range(3))
            candidates = set()
            for cell in cells_for(midpoint, midpoint):
                candidates.update(buckets.get(cell, ()))
            # Benchy uses a 0.4 mm first layer. There cannot be a preceding
            # layer below it, so its model walls are grounded by definition.
            supported_below = midpoint[2] <= 0.41
            nearest = (999.0, 999.0)
            for a, b in candidates:
                distance, vertical_gap = segment_distance_xy(midpoint, a, b)
                if distance < nearest[0]:
                    nearest = distance, vertical_gap
                if distance <= 0.5 and 0.025 <= vertical_gap <= 0.32:
                    supported_below = True
                    break
            exempt = any(token in feature.lower() for token in EXEMPT_SUPPORT_TYPES)
            replacement_dependent = False
            if not supported_below and not exempt:
                removed_candidates = set()
                for cell in cells_for(midpoint, midpoint):
                    removed_candidates.update(removed_buckets.get(cell, ()))
                for a, b in removed_candidates:
                    distance, vertical_gap = segment_distance_xy(midpoint, a, b)
                    if distance <= 0.5 and abs(vertical_gap) <= 0.32:
                        replacement_dependent = True
                        break
            present_in_planar_control = False
            if not supported_below and not exempt and not replacement_dependent and control is not None:
                control_candidates = set()
                for cell in cells_for(midpoint, midpoint):
                    control_candidates.update(control_buckets.get(cell, ()))
                for a, b in control_candidates:
                    distance, vertical_gap = segment_distance_xy(midpoint, a, b)
                    if distance <= 0.05 and abs(vertical_gap) <= 0.02:
                        present_in_planar_control = True
                        break
            status = (
                "supported" if supported_below or exempt
                else "replacement-dependent" if replacement_dependent
                else "planar-control" if present_in_planar_control
                else "unsupported"
            )
            counts[(status, feature)] += 1
            if status == "unsupported" and len(unsupported_samples) < 80:
                unsupported_samples.append(
                    (line_number, feature, midpoint, nearest, line.strip()[:140])
                )
            length = math.dist(old, new)
            if run is None:
                run = {"feature": feature, "line": line_number, "length": 0.0,
                       "segments": 0, "unsupported": 0, "minimum_z": midpoint[2],
                       "maximum_z": midpoint[2]}
            run["length"] += length
            run["segments"] += 1
            run["unsupported"] += status == "unsupported"
            run["minimum_z"] = min(run["minimum_z"], midpoint[2])
            run["maximum_z"] = max(run["maximum_z"], midpoint[2])
            for cell in cells_for(old, new, 0.0):
                buckets[cell].append((old, new))

    if run is not None:
        runs.append(run)

    floating_runs = [
        item for item in runs
        if item["segments"] > 0
        and item["unsupported"] / item["segments"] >= 0.5
        and not any(token in item["feature"].lower() for token in EXEMPT_SUPPORT_TYPES)
        and any(token in item["feature"].lower() for token in SURFACE_REMAINDER_TYPES)
    ]
    print(f"filtered_gcode {destination}")
    print("removed")
    for (status, name), count in sorted(counts.items()):
        if status == "removed":
            print(f"  {name or '<untyped>'}: {count}")
    print("unsupported conventional extrusion segments")
    for (status, name), count in sorted(counts.items()):
        if status == "unsupported":
            print(f"  {name or '<untyped>'}: {count}")
    print("conventional segments intentionally dependent on removed replacement")
    for (status, name), count in sorted(counts.items()):
        if status == "replacement-dependent":
            print(f"  {name or '<untyped>'}: {count}")
    if control is not None:
        print("otherwise unsupported segments also present in planar control")
        for (status, name), count in sorted(counts.items()):
            if status == "planar-control":
                print(f"  {name or '<untyped>'}: {count}")
    print(f"floating_surface_remainder_runs {len(floating_runs)}")
    for item in sorted(floating_runs, key=lambda value: (-value["maximum_z"], value["line"]))[:60]:
        ratio = item["unsupported"] / item["segments"]
        print(
            f"  line {item['line']} {item['feature'] or '<untyped>'}: "
            f"{item['length']:.3f}mm, {ratio:.0%} unsupported, "
            f"Z{item['minimum_z']:.3f}..{item['maximum_z']:.3f}"
        )
    print("unsupported samples")
    for line_number, name, midpoint, nearest, text in unsupported_samples[:30]:
        xyz = tuple(round(value, 3) for value in midpoint)
        proximity = tuple(round(value, 3) for value in nearest)
        print(f"  line {line_number} {name}: xyz={xyz} nearest={proximity} {text}")
    return 1 if floating_runs else 0


if __name__ == "__main__":
    raise SystemExit(main())
