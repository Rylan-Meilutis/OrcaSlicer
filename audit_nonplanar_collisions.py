import math
import re
import sys
from collections import Counter, defaultdict
from audit_gcode_moves import linearized_gcode


NUMBER = r"-?(?:\d+(?:\.\d*)?|\.\d+)"
AXIS_RE = re.compile(r"([XYZEF])(" + NUMBER + r")")
WIDTH_RE = re.compile(r";WIDTH:(" + NUMBER + r")")


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def segment_intersection(first, second):
    a, b = first[0], first[1]
    c, d = second[0], second[1]
    r = (b[0] - a[0], b[1] - a[1])
    s = (d[0] - c[0], d[1] - c[1])
    denominator = cross(r, s)
    relative = (c[0] - a[0], c[1] - a[1])
    if abs(denominator) < 1e-9:
        return None
    t = cross(relative, s) / denominator
    u = cross(relative, r) / denominator
    if 1e-4 < t < 1.0 - 1e-4 and 1e-4 < u < 1.0 - 1e-4:
        return t, u
    return None


def point_segment_distance(point, a, b):
    delta = (b[0] - a[0], b[1] - a[1])
    denominator = delta[0] * delta[0] + delta[1] * delta[1]
    if denominator == 0.0:
        return math.hypot(point[0] - a[0], point[1] - a[1]), 0.0
    t = max(0.0, min(1.0, ((point[0] - a[0]) * delta[0] +
                            (point[1] - a[1]) * delta[1]) / denominator))
    nearest = (a[0] + t * delta[0], a[1] + t * delta[1])
    return math.hypot(point[0] - nearest[0], point[1] - nearest[1]), t


def segment_proximity(first, second):
    candidates = []
    for point, segment, owner_t in (
        (first[0], second, 0.0), (first[1], second, 1.0),
        (second[0], first, 0.0), (second[1], first, 1.0)):
        distance, segment_t = point_segment_distance(point, segment[0], segment[1])
        if segment is second:
            candidates.append((distance, owner_t, segment_t))
        else:
            candidates.append((distance, segment_t, owner_t))
    return min(candidates)


def position_at(segment, t):
    return tuple(segment[0][axis] + t * (segment[1][axis] - segment[0][axis])
                 for axis in range(3))


def collinear_overlap(first, second):
    """Return interior samples for a positive-length shared XY interval."""
    a, b = first
    c, d = second
    r = (b[0] - a[0], b[1] - a[1])
    length_squared = r[0] ** 2 + r[1] ** 2
    if length_squared < 1e-12:
        return None
    length = math.sqrt(length_squared)
    if any(abs(cross(r, (p[0] - a[0], p[1] - a[1]))) > 1e-6 * length
           for p in (c, d)):
        return None
    parameters = [((p[0] - a[0]) * r[0] + (p[1] - a[1]) * r[1]) /
                  length_squared for p in (c, d)]
    lo, hi = max(0.0, min(parameters)), min(1.0, max(parameters))
    if (hi - lo) * length <= 1e-6:
        return None
    t = (lo + hi) / 2
    _, u = point_segment_distance(position_at(first, t), c, d)
    return t, u


path = sys.argv[1]
maximum_z_gap = float(sys.argv[2]) if len(sys.argv) > 2 else 0.08
different_feature_courses = "--different-feature-courses-only" in sys.argv[3:]
include_nonplanar_walls = "--include-nonplanar-walls" in sys.argv[3:] or different_feature_courses
position = {"X": 0.0, "Y": 0.0, "Z": 0.0, "E": 0.0}
relative_e = False
relative_xyz = False
flushing = False
feature = ""
width = 0.45
layer = -1
segments = []
grid = defaultdict(list)
conflicts = []
counts = Counter()
vertical_relations = Counter()
extrusion_chain = 0
previous_feature_wall = None
wall_discontinuities = []

with open(path, errors="replace") as stream:
    for line_number, line in linearized_gcode(stream, exclude_cleaner_extrusion=True):
        line = line.lstrip()
        tag = line[1:].strip() if line.startswith(";") else ""
        if tag == "FLUSH_START":
            flushing = True
        elif tag == "FLUSH_END":
            flushing = False
        if line.startswith(";LAYER_CHANGE"):
            layer += 1
        if line.startswith(";TYPE:"):
            feature = line[6:].strip()
        width_match = WIDTH_RE.match(line)
        if width_match:
            width = float(width_match.group(1))
        command = line.split(";", 1)[0].strip()
        opcode = command.split()[0] if command else ""
        if opcode in ("G90", "G91"):
            relative_xyz = opcode == "G91"
            continue
        if command == "M82":
            relative_e = False
            continue
        if command == "M83":
            relative_e = True
            continue
        values = {key: float(value) for key, value in AXIS_RE.findall(command)}
        if opcode == "G92":
            for axis in "XYZE":
                if axis in values:
                    position[axis] = values[axis]
            continue
        if opcode not in ("G0", "G1"):
            continue
        old = (position["X"], position["Y"], position["Z"])
        new_e = values.get("E", position["E"])
        for axis in "XYZ":
            if axis in values:
                position[axis] = position[axis] + values[axis] if relative_xyz else values[axis]
        new = (position["X"], position["Y"], position["Z"])
        e_delta = (values.get("E", 0.0) if relative_e or relative_xyz
                   else new_e - position["E"])
        position["E"] += e_delta
        if flushing or "E" not in values or e_delta <= 1e-8 or math.dist(old, new) <= 0.01:
            if flushing or e_delta < -1e-8 or math.dist(old[:2], new[:2]) > 0.005:
                previous_feature_wall = None
            if math.dist(old, new) > 0.01 or e_delta < -1e-8:
                extrusion_chain += 1
            continue

        comment = line.split(";", 1)[1].strip() if ";" in line else ""
        nonplanar = ("non-planar" in feature.lower() or "non-planar" in comment.lower()
                     or "surface-following feature" in comment.lower())
        perimeter = "wall" in feature.lower() or "perimeter" in comment.lower()
        current = {
            "points": (old, new), "feature": feature, "comment": comment,
            "nonplanar": nonplanar, "perimeter": perimeter, "width": width,
            "line": line_number, "layer": layer,
            "chain": extrusion_chain,
        }
        course_match = re.search(r"surface-following feature perimeter course (\d+)", comment)
        current["course"] = int(course_match.group(1)) if course_match else None
        if previous_feature_wall is not None and current["course"] is not None:
            endpoint = previous_feature_wall["points"][1]
            if (current["course"] == previous_feature_wall["course"] and
                math.dist(endpoint[:2], old[:2]) <= 0.005 and
                abs(endpoint[2] - old[2]) > 0.005):
                wall_discontinuities.append((previous_feature_wall["line"], line_number,
                    current["course"], round(old[2] - endpoint[2], 4)))
        previous_feature_wall = current if current["course"] is not None else None
        if nonplanar or perimeter:
            nearby = set()
            padding = 0.1
            for x in range(math.floor(min(old[0], new[0]) - padding),
                           math.floor(max(old[0], new[0]) + padding) + 1):
                for y in range(math.floor(min(old[1], new[1]) - padding),
                               math.floor(max(old[1], new[1]) + padding) + 1):
                    nearby.update(grid.get((x, y), ()))
            for previous_index in nearby:
                previous = segments[previous_index]
                if different_feature_courses and (current["course"] is None or
                    previous["course"] is None or current["course"] == previous["course"]):
                    continue
                if not (current["nonplanar"] or previous["nonplanar"]):
                    continue
                if not (current["perimeter"] or previous["perimeter"]):
                    continue
                if current["nonplanar"] and previous["nonplanar"] and not include_nonplanar_walls:
                    continue
                intersection = segment_intersection(current["points"], previous["points"])
                kind = "crossing"
                if intersection is None:
                    overlap = collinear_overlap(current["points"], previous["points"])
                    if overlap is not None:
                        current_t, previous_t = overlap
                    else:
                        distance, current_t, previous_t = segment_proximity(
                            current["points"], previous["points"])
                        if distance > 0.06:
                            continue
                        # Joined paths and loop closures are not duplicate edges.
                        # Check positive-length overlap first so backtracking is
                        # still detected even when both endpoints are shared.
                        if (current_t in (0.0, 1.0) and previous_t in (0.0, 1.0) and
                            math.dist(position_at(current["points"], current_t),
                                      position_at(previous["points"], previous_t)) <= 1e-6):
                            continue
                    kind = "near-overlap"
                else:
                    current_t, previous_t = intersection
                current_point = position_at(current["points"], current_t)
                previous_point = position_at(previous["points"], previous_t)
                signed_z_gap = current_point[2] - previous_point[2]
                z_gap = abs(signed_z_gap)
                if z_gap > maximum_z_gap:
                    continue
                if "--below-prior-only" in sys.argv[3:] and signed_z_gap >= -0.001:
                    continue
                counts[(kind, current["feature"], previous["feature"])] += 1
                vertical_relations[(kind, "below prior extrusion" if signed_z_gap < -0.001
                                    else "above prior extrusion" if signed_z_gap > 0.001
                                    else "same height")] += 1
                if len(conflicts) < 50:
                    conflicts.append((
                        kind, current["line"], previous["line"], current["layer"],
                        round(signed_z_gap, 4), current["feature"], previous["feature"],
                        tuple(round(value, 3) for value in current_point)))

        segment_index = len(segments)
        segments.append(current)
        for x in range(math.floor(min(old[0], new[0])),
                       math.floor(max(old[0], new[0])) + 1):
            for y in range(math.floor(min(old[1], new[1])),
                           math.floor(max(old[1], new[1])) + 1):
                grid[(x, y)].append(segment_index)

print("conflict counts")
for key, value in sorted(counts.items()):
    print(key, value)
print("vertical relations")
for key, value in sorted(vertical_relations.items()):
    print(key, value)
print("samples (signed Z gap: current minus previously deposited)")
for conflict in conflicts:
    print(conflict)
print("feature-wall vertical discontinuities", len(wall_discontinuities))
for discontinuity in wall_discontinuities[:50]:
    print(discontinuity)

# Allow the same real-G-code audit to be used as a failing regression check,
# without changing the report-only behavior of existing invocations.
sys.exit(1 if ("--fail-on-conflict" in sys.argv[3:] and counts) or
             ("--fail-on-wall-discontinuity" in sys.argv[3:] and wall_discontinuities) else 0)
