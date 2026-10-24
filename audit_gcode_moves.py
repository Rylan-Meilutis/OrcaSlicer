"""Canonical G0/G1 moves for the diagnostic audits, retaining source line numbers.

XY arcs are tessellated to a bounded chord error. Unsupported arc planes fail
explicitly instead of silently dropping deposited material from an audit.
"""
import math
import re

WORDS = re.compile(r"([A-Z])([-+]?(?:\d+(?:\.\d*)?|\.\d+))", re.I)


def linearized_gcode(stream, chord_error=0.005, exclude_cleaner_extrusion=False):
    if chord_error <= 0:
        raise ValueError("chord error must be positive")
    position = dict.fromkeys("XYZE", 0.0)
    relative_xyz = False
    relative_e = None
    absolute_center = False
    plane = "G17"
    in_cleaner = False
    yield 0, "G90\n"
    yield 0, "M83\n"
    for lineno, line in enumerate(stream, 1):
        command, separator, comment = line.partition(";")
        tokens = command.split()
        opcode = tokens[0].upper() if tokens else ""
        words = {key.upper(): float(value) for key, value in WORDS.findall(command)}
        if exclude_cleaner_extrusion and opcode == "G12":
            if words.get("S") == 90:
                in_cleaner = True
            elif words.get("S") == 91:
                in_cleaner = False
        if opcode in ("G90", "G91"):
            relative_xyz = opcode == "G91"
            relative_e = None
            yield lineno, "G90\n"
            continue
        if opcode in ("M82", "M83"):
            relative_e = opcode == "M83"
            yield lineno, "M83\n"
            continue
        if opcode in ("G90.1", "G91.1"):
            absolute_center = opcode == "G90.1"
            continue
        if opcode in ("G17", "G18", "G19"):
            plane = opcode
            continue
        if opcode == "G92":
            for axis in "XYZE":
                if axis in words:
                    position[axis] = words[axis]
            yield lineno, "G92 " + " ".join(
                f"{axis}{position[axis]:.9f}" for axis in "XYZ" if axis in words) + "\n"
            continue
        if opcode not in ("G0", "G1", "G2", "G3"):
            yield lineno, line
            continue
        old = position.copy()
        for axis in "XYZ":
            if axis in words:
                position[axis] = words[axis] + (old[axis] if relative_xyz else 0.)
        de = 0.
        if "E" in words:
            e_is_relative = relative_xyz if relative_e is None else relative_e
            de = words["E"] if e_is_relative else words["E"] - old["E"]
            position["E"] += de
        # INDX cleaner macros move the nozzle to the bucket outside the
        # explicit XYZ stream. Its purge is not model material, even if TYPE
        # still names the preceding wall. Keep modal/E state for the return.
        if exclude_cleaner_extrusion and in_cleaner:
            de = 0.
        samples = 1
        center = None
        if opcode in ("G2", "G3"):
            if plane != "G17":
                raise ValueError(f"line {lineno}: unsupported arc plane {plane}")
            clockwise = opcode == "G2"
            def sweep_from(cx, cy):
                start = math.atan2(old["Y"] - cy, old["X"] - cx)
                end = math.atan2(position["Y"] - cy, position["X"] - cx)
                sweep = (start - end if clockwise else end - start) % math.tau
                return start, sweep if sweep > 1e-12 else math.tau
            if "I" in words or "J" in words:
                cx = words.get("I", old["X"] if absolute_center else 0.)
                cy = words.get("J", old["Y"] if absolute_center else 0.)
                if not absolute_center:
                    cx += old["X"]
                    cy += old["Y"]
                radius = math.hypot(old["X"] - cx, old["Y"] - cy)
            elif "R" in words:
                dx, dy = position["X"] - old["X"], position["Y"] - old["Y"]
                chord = math.hypot(dx, dy)
                radius = abs(words["R"])
                if chord == 0 or chord > 2 * radius + 1e-9:
                    raise ValueError(f"line {lineno}: invalid radius arc")
                h = math.sqrt(max(0., radius * radius - 0.25 * chord * chord))
                candidates = [(0.5 * (old["X"] + position["X"]) - sign * dy * h / chord,
                               0.5 * (old["Y"] + position["Y"]) + sign * dx * h / chord)
                              for sign in (-1., 1.)]
                cx, cy = min(candidates, key=lambda c: sweep_from(*c)[1]) if words["R"] >= 0 else max(
                    candidates, key=lambda c: sweep_from(*c)[1])
            else:
                raise ValueError(f"line {lineno}: arc has no center or radius")
            if radius <= 0:
                raise ValueError(f"line {lineno}: zero-radius arc")
            end_radius = math.hypot(position["X"] - cx, position["Y"] - cy)
            if abs(end_radius - radius) > 0.05:
                raise ValueError(f"line {lineno}: inconsistent arc radii")
            start, sweep = sweep_from(cx, cy)
            step = min(math.pi / 2, 2 * math.acos(max(-1., 1. - chord_error / radius)))
            samples = max(1, math.ceil(sweep / step))
            center = cx, cy
        for idx in range(1, samples + 1):
            fraction = idx / samples
            xyz = [old[axis] + fraction * (position[axis] - old[axis]) for axis in "XYZ"]
            if center is not None and idx != samples:
                angle = start + (-1 if clockwise else 1) * sweep * fraction
                r = radius + fraction * (end_radius - radius)
                xyz[:2] = [cx + r * math.cos(angle), cy + r * math.sin(angle)]
            move = "G0" if opcode == "G0" else "G1"
            move += " " + " ".join(f"{axis}{value:.9f}" for axis, value in zip("XYZ", xyz))
            if "E" in words:
                move += f" E{de / samples:.12f}"
            if "F" in words:
                move += f" F{words['F']:.6f}"
            if separator:
                move += " ;" + comment.rstrip()
            yield lineno, move + "\n"
