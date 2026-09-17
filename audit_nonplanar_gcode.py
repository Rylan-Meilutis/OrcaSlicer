import math
import re
import sys
from collections import defaultdict, Counter
from audit_gcode_moves import linearized_gcode

number = r"-?(?:\d+(?:\.\d*)?|\.\d+)"
axis_re = re.compile(r"([XYZEF])(" + number + r")")
buckets = defaultdict(list)
counts = Counter()
feature_course_support = Counter()
samples = []
run_lengths = []
active_run = None
position = {"X": 0.0, "Y": 0.0, "Z": 0.0, "E": 0.0}
feature = ""
height = 0.2
relative_e = False
extrusion_chain = 0
sample_course = (int(sys.argv[sys.argv.index("--sample-course") + 1])
                 if "--sample-course" in sys.argv else None)
required_top = None
top_length = 0.0
if "--require-top-z-range" in sys.argv:
    start = sys.argv.index("--require-top-z-range") + 1
    required_top = tuple(map(float, sys.argv[start:start + 3]))
    if (len(required_top) != 3 or not all(map(math.isfinite, required_top)) or
            required_top[0] >= required_top[1] or required_top[2] <= 0):
        raise SystemExit("--require-top-z-range requires MIN_Z MAX_Z MIN_LENGTH_MM")

def cells_for(a, b, padding=0.45):
    lo_x, hi_x = sorted((a[0], b[0])); lo_y, hi_y = sorted((a[1], b[1]))
    for ix in range(math.floor(lo_x-padding), math.floor(hi_x+padding)+1):
        for iy in range(math.floor(lo_y-padding), math.floor(hi_y+padding)+1):
            yield (ix, iy)

def point_segment_xy(p, a, b):
    dx, dy = b[0]-a[0], b[1]-a[1]
    den = dx*dx + dy*dy
    t = 0.0 if den == 0 else max(0.0, min(1.0, ((p[0]-a[0])*dx+(p[1]-a[1])*dy)/den))
    q = (a[0]+t*dx, a[1]+t*dy, a[2]+t*(b[2]-a[2]))
    return math.hypot(p[0]-q[0], p[1]-q[1]), p[2]-q[2]

with open(sys.argv[1], errors="replace") as stream:
    for lineno, line in linearized_gcode(stream, exclude_cleaner_extrusion=True):
        if line.startswith(";HEIGHT:"):
            height = float(line.split(":", 1)[1])
        if line.startswith(";TYPE:"):
            if active_run is not None:
                run_lengths.append(active_run)
                active_run = None
            feature = line[6:].strip()
        command = line.split(";", 1)[0].strip()
        if command == "M82":
            relative_e = False
            continue
        if command == "M83":
            relative_e = True
            continue
        if command.startswith("G92"):
            values = {key: float(value) for key, value in axis_re.findall(command)}
            for axis in "XYZE":
                if axis in values:
                    position[axis] = values[axis]
            continue
        if not (line.startswith("G0") or line.startswith("G1")):
            continue
        old = (position["X"], position["Y"], position["Z"])
        values = {key: float(value) for key, value in axis_re.findall(line.split(";",1)[0])}
        new_e = values.get("E", position["E"])
        for key in "XYZ":
            if key in values: position[key] = values[key]
        new = (position["X"], position["Y"], position["Z"])
        e_delta = new_e if relative_e else new_e - position["E"]
        extruding = "E" in values and e_delta > 1e-8 and math.hypot(new[0]-old[0], new[1]-old[1]) > 0.01
        if not relative_e:
            position["E"] = new_e
        if not extruding:
            if math.dist(old, new) > 0.01 or ("E" in values and e_delta < -1e-8):
                extrusion_chain += 1
            if math.dist(old, new) > 0.01 and active_run is not None:
                run_lengths.append(active_run)
                active_run = None
            continue
        is_np = ("Non-planar" in feature or "non-planar" in line.lower() or
                 "surface-following feature" in line.lower())
        if is_np:
            segment_length = math.dist(old, new)
            if required_top and feature == "Non-planar top surface":
                low, high, _ = required_top
                dz = new[2] - old[2]
                if abs(dz) < 1e-12:
                    fraction = float(low <= old[2] <= high)
                else:
                    entry, leave = sorted(((low - old[2]) / dz, (high - old[2]) / dz))
                    fraction = max(0.0, min(1.0, leave) - max(0.0, entry))
                top_length += fraction * segment_length
            if active_run is None:
                active_run = [feature, segment_length, lineno]
            else:
                active_run[1] += segment_length
            # Walls and transitions can float at their ends too. A midpoint
            # test alone incorrectly accepted long, partially supported walls.
            sample_count = max(1, math.ceil(segment_length / 0.2))
            segment_supported = True
            for sample_index in range(sample_count):
                t = (sample_index + 0.5) / sample_count
                p = tuple(old[i] + t * (new[i] - old[i]) for i in range(3))
                candidates = set()
                for key in cells_for(p, p): candidates.update(buckets.get(key, ()))
                supported = False
                nearest = (999.0, 999.0)
                nearest_height_error = float("inf")
                for a,b,chain in candidates:
                    # Earlier segments of this same rising wall are not a
                    # deposited lower course, even if their endpoints happen
                    # to be within one bead width and height of this sample.
                    if chain == extrusion_chain:
                        continue
                    d,gap = point_segment_xy(p,a,b)
                    if d <= 0.4 and abs(gap - height) < nearest_height_error:
                        nearest=(d,gap)
                        nearest_height_error = abs(gap - height)
                    if d <= 0.4 and 0.04 <= gap <= height + 0.06:
                        supported=True; break
                segment_supported &= supported
                matching_course = (sample_course is None or
                    re.search(r"surface-following feature perimeter course " + str(sample_course) + r"\b", line))
                if not supported and matching_course and len(samples) < 8:
                    samples.append((lineno, feature, tuple(round(v,3) for v in p), tuple(round(v,3) for v in nearest), line.strip()[:130]))
            counts[(feature, "supported" if segment_supported else "unsupported")] += 1
            course = re.search(r"surface-following feature perimeter course (\d+)", line)
            if course:
                feature_course_support[(int(course.group(1)), "supported" if segment_supported else "unsupported")] += 1
        for key in cells_for(old,new,0.0): buckets[key].append((old,new,extrusion_chain))

if active_run is not None:
    run_lengths.append(active_run)

print("counts")
for key, value in sorted(counts.items()): print(key, value)
print("samples")
for sample in samples: print(sample)
print("feature course support")
for key, value in sorted(feature_course_support.items()): print(key, value)
print("runs")
run_counts = Counter()
for feature_name, length, unused_line in run_lengths:
    for threshold in (0.2, 0.4, 0.8):
        if length < threshold:
            run_counts[(feature_name, threshold)] += 1
for key, value in sorted(run_counts.items()):
    print(key, value)
for feature_name in sorted({run[0] for run in run_lengths}):
    lengths = [run[1] for run in run_lengths if run[0] == feature_name]
    print(feature_name, "count", len(lengths), "min", round(min(lengths), 3),
          "mean", round(sum(lengths) / len(lengths), 3))
print("short run samples")
for feature_name, length, line_number in sorted(
        (run for run in run_lengths if run[1] < 0.8), key=lambda run: run[1])[:40]:
    print(feature_name, round(length, 3), line_number)

missing_top = required_top is not None and top_length < required_top[2]
if required_top:
    # Presence guard, not a claim of complete mesh coverage. A clean contact
    # audit must not silently pass when the intended surface falls back.
    print("required top range", required_top[:2], "length", top_length,
          "minimum", required_top[2], "FAIL" if missing_top else "PASS")
sys.exit(1 if missing_top or ("--fail-on-unsupported" in sys.argv[2:] and
         (not counts or any(status == "unsupported" and count
                            for (_, status), count in counts.items()))) else 0)
