"""Regression checks for the standalone G-code audit reader."""
import io
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from audit_gcode_moves import WORDS, linearized_gcode
from audit_arc_gcode import audit_crossings
from audit_arc_support import audit_support
from audit_arc_families import audit_families
from audit_arc_overlap import audit_overlap


def read_moves(gcode):
    return [(line, {key: float(value) for key, value in WORDS.findall(text)})
            for line, text in linearized_gcode(io.StringIO(gcode))
            if text.startswith(("G0 ", "G1 "))]


class AuditMovesTest(unittest.TestCase):
    def test_arc_overlap_detects_a_curved_cleanup_retrace(self):
        # KSR Z50.2: neither centerline crosses the other, but the middle of
        # the second arc follows the first within about 0.15 mm.
        source = ("G90\nM83\n;WIDTH:0.45\n;TYPE:Arc bridge\n"
                  "G0 X150.530 Y74.968 Z50.2\n"
                  "G1 X150.478 Y74.927 E.01\nG1 X150.341 Y74.850 E.01\n"
                  "G1 X150.194 Y74.794 E.01\nG1 X150.041 Y74.760 E.01\n"
                  "G1 X149.885 Y74.748 E.01\nG0 X149.529 Y74.759\n"
                  "G1 X149.681 Y74.689 E.01\nG1 X149.847 Y74.638 E.01\n"
                  "G1 X150.018 Y74.610 E.01\nG1 X150.192 Y74.605 E.01\n"
                  "G1 X150.365 Y74.623 E.01\nG1 X150.533 Y74.664 E.01\n"
                  "G1 X150.665 Y74.716 E.01\n")
        self.assertEqual(audit_crossings(io.StringIO(source))["crossings"], [])
        result = audit_overlap(io.StringIO(source), 50.2)
        self.assertEqual(len(result["overlaps"]), 1)
        self.assertGreater(result["overlaps"][0]["overlap_run_mm"], .4)

    def test_arc_overlap_detects_an_interior_retrace_without_a_crossing(self):
        source = ("G90\nM83\n;WIDTH:0.45\n;TYPE:Outer wall\n"
                  "G0 X0 Y0 Z1\nG1 X5 E1\n;TYPE:Arc bridge\n"
                  "G0 X0 Y1\nG1 X1 Y.05 E.1\nG1 X4 E.5\nG1 X5 Y1 E.1\n")
        self.assertEqual(audit_crossings(io.StringIO(source))["crossings"], [])
        result = audit_overlap(io.StringIO(source), 1.)
        self.assertEqual(len(result["overlaps"]), 1)
        self.assertGreater(result["overlaps"][0]["overlap_run_mm"], 2.9)

    def test_arc_overlap_allows_normal_bonding_and_ignores_other_layers(self):
        for z, separation in ((1., .3375), (.8, .01)):
            source = (f"G90\nM83\n;WIDTH:0.45\n;TYPE:Outer wall\nG0 X0 Y0 Z{z}\nG1 X5 E1\n"
                      f";TYPE:Arc bridge\nG0 X0 Y{separation} Z1\nG1 X5 E1\n")
            self.assertEqual(audit_overlap(io.StringIO(source), 1.)["overlaps"], [])

    def test_fitted_arc_contact_uses_a_resolved_terminal_tangent(self):
        # From KSR's narrow bridge: a coarse chord rotates the free-end
        # tangent enough to falsely reject this genuinely bonded neighbour.
        prefix = "G90\nM83\n;WIDTH:0.45\n"
        parent = (";TYPE:Outer wall\nG0 X87.47 Y129.231 Z8.2\n"
                  "G3 X90.53 Y129.231 I1.53 J.828 E.27825\n")
        child = (";TYPE:Arc bridge\nG0 X90.864 Y129.032 Z8.2\n"
                 "G2 X87.136 Y129.032 I-1.864 J1.027 E.33785\n")
        result = audit_support(io.StringIO(prefix + parent + child), 8.2,
                               check_bodies=False, min_overlap=.05)
        self.assertEqual(result["unsupported_ends"], [])
        refined = audit_support(io.StringIO(prefix + parent + child), 8.2,
                                chord_error=.000025, check_bodies=False, min_overlap=.05)
        self.assertEqual(refined["unsupported_ends"], [])
        future = audit_support(io.StringIO(prefix + child + parent), 8.2,
                               check_bodies=False, min_overlap=.05)
        self.assertEqual(len(future["unsupported_ends"]), 2)

    def test_arc_short_path_report_uses_the_extruded_width(self):
        source = ("G90\nM83\n;WIDTH:0.45\n;TYPE:Arc bridge\n"
                  "G0 X0 Y0 Z1\nG1 X.5 E.1\nG0 X1\nG1 X1.8 E.1\n")
        result = audit_support(io.StringIO(source), 1.)
        self.assertEqual(result["arc_paths"], 2)
        self.assertEqual(len(result["short_paths"]), 1)
        self.assertAlmostEqual(result["short_paths"][0]["length_mm"], .5)

    def test_arc_contact_is_not_proof_of_meaningful_overlap(self):
        prefix = "G90\nM83\n;WIDTH:0.45\n;TYPE:Outer wall\nG0 X-1 Y0 Z1\nG1 X3 E1\n"
        for separation, supported in ((.44, False), (.35, True)):
            with self.subTest(separation=separation):
                source = prefix + f";TYPE:Arc bridge\nG0 X0 Y{separation}\nG1 X2 E.1\n"
                touching = audit_support(io.StringIO(source), 1.)
                self.assertEqual(touching["unsupported_ends"], [])
                bonded = audit_support(io.StringIO(source), 1., min_overlap=.05)
                self.assertEqual(not bonded["unsupported_ends"], supported)
                self.assertEqual(not bonded["unsupported_bodies"], supported)

    def test_family_audit_counts_distinct_concentric_courses_not_repeated_paths(self):
        def arc(center, radius):
            return (f"G0 X{center+radius} Y0 Z1\n"
                    f"G3 X{center-radius} Y0 I{-radius} J0 E1\n")
        prefix = "G90\nM83\n;TYPE:Arc bridge\n" + arc(0, 1)
        cases = (
            ([(10, 1), (10, 1.4), (10, 1.8)], 1.),
            ([(10, 1), (11, 1.4), (12, 1.8)], 0.),
            ([(10, 1), (10, 1), (10, 1)], 0.),
            ([(10, 1.8), (10, 1.4), (10, 1)], 0.),
        )
        for courses, expected in cases:
            with self.subTest(courses=courses):
                result = audit_families(io.StringIO(prefix + "".join(arc(*c) for c in courses)), 1.)
                self.assertEqual(result["paths"], 4)
                self.assertAlmostEqual(result["multi_ring_fraction"], expected)

    def test_supported_arc_body_can_contact_an_existing_cap_between_anchored_ends(self):
        source = ("G90\nM83\n;WIDTH:0.45\n;TYPE:Outer wall\n"
                  "G0 X-.6 Y-1 Z1\nG1 Y1 E1\n"
                  "G0 X.6 Y-1\nG1 Y1 E1\n"
                  "G0 X0 Y-1\nG1 Y0 E1\n"
                  ";TYPE:Arc bridge\nG0 X-.6 Y.3\nG1 X.6 E1\n")
        result = audit_support(io.StringIO(source), 1.)
        self.assertEqual(result["unsupported_ends"], [])
        self.assertEqual(result["unsupported_bodies"], [])

    def test_outer_arc_cannot_use_an_inner_arc_printed_later(self):
        foundation = ("G90\nM83\n;WIDTH:0.45\n;TYPE:Outer wall\n"
                      "G0 X0 Y-5 Z1\nG1 Y5 E1\n"
                      "G0 X0 Y-1.6\nG3 X0 Y1.6 I0 J1.6 E1\n")
        inner = ";TYPE:Arc bridge\nG0 X0 Y-2\nG3 X0 Y2 I0 J2 E1\n"
        outer = ";TYPE:Arc bridge\nG0 X0 Y-2.4\nG3 X0 Y2.4 I0 J2.4 E1\n"
        bad = audit_support(io.StringIO(foundation + outer + inner), 1., check_bodies=True)
        self.assertEqual(bad["unsupported_ends"], [])
        self.assertEqual(len(bad["unsupported_bodies"]), 1)
        good = audit_support(io.StringIO(foundation + inner + outer), 1., check_bodies=True)
        self.assertEqual(good["unsupported_ends"], [])
        self.assertEqual(good["unsupported_bodies"], [])

    def test_arc_support_rejects_a_staircase_of_free_end_caps(self):
        prefix = "G90\nM83\n;WIDTH:0.45\n;TYPE:Outer wall\nG0 X0 Y-1 Z1\nG1 Y1 E1\n"
        parent = ";TYPE:Arc bridge\nG0 X0 Y0\nG1 X2 E1\n"
        child = "G0 X0 Y.3\nG1 X2.3 E1\n"
        result = audit_support(io.StringIO(prefix + parent + child), 1.)
        self.assertTrue(any(failure["xy"] == (2.3, .3)
                            for failure in result["unsupported_ends"]))

    def test_arc_support_does_not_count_future_paths_or_self_support(self):
        prefix = "G90\nM83\n;WIDTH:0.45\n;TYPE:Outer wall\nG0 X0 Y-1 Z1\nG1 Y1 E1\n"
        child = ";TYPE:Arc bridge\nG0 X0 Y0\nG1 X2 E1\n"
        parent = ";TYPE:Outer wall\nG0 X2 Y-1\nG1 Y1 E1\n"
        bad = audit_support(io.StringIO(prefix + child + parent), 1.)
        self.assertEqual(bad["arc_paths"], 1)
        self.assertEqual([failure["end"] for failure in bad["unsupported_ends"]], ["end"])
        good = audit_support(io.StringIO(prefix + parent + child), 1.)
        self.assertEqual(good["unsupported_ends"], [])

    def test_arc_support_does_not_flatten_a_raised_wall_into_an_anchor(self):
        source = ("G90\nM83\n;WIDTH:0.45\n;TYPE:Brick wall\n"
                  "G0 X0 Y-1 Z1.1\nG1 Y1 E1\n"
                  ";TYPE:Arc bridge\nG0 X0 Y0 Z1\nG1 X.4 E1\n")
        result = audit_support(io.StringIO(source), 1.)
        self.assertEqual(result["arc_paths"], 1)
        self.assertEqual(len(result["unsupported_ends"]), 2)

    def test_short_restart_cannot_borrow_contact_from_a_later_concentric_arc(self):
        # A small quarter-circle departs from a wall, but its far end has
        # nothing beside it until the neighboring larger arc is deposited.
        prefix = "G90\nM83\n;WIDTH:0.45\n;TYPE:Outer wall\nG0 X-1 Y0 Z1\nG1 X1 E1\n"
        short = ";TYPE:Arc bridge\nG0 X.7 Y0\nG3 X0 Y.7 I-.7 J0 E.1\n"
        neighbor = ";TYPE:Arc bridge\nG0 X1 Y0\nG3 X-1 Y0 I-1 J0 E.3\n"
        bad = audit_support(io.StringIO(prefix + short + neighbor), 1.)
        self.assertTrue(any(f["end"] == "end" and f["xy"] == (0., .7)
                            for f in bad["unsupported_ends"]))
        reordered = audit_support(io.StringIO(prefix + neighbor + short), 1.)
        self.assertFalse(any(f["xy"] == (0., .7) for f in reordered["unsupported_ends"]))
        # Reordering alone does not make the larger, unbacked body printable.
        self.assertTrue(reordered["unsupported_bodies"])

    def test_arc_audit_detects_short_fitted_terminal_crossing(self):
        source = ("G90\nM83\n;TYPE:Arc bridge\nG0 X133.749 Y70.899 Z50.2\n"
                  "G1 X133.716 Y70.902 E.00246\n"
                  "G2 X133.374 Y70.978 I.195 J1.681 E.0261\n"
                  "G0 X134.125 Y70.596\nG1 X134.103 Y70.602 E.00169\n"
                  "G2 X133.6 Y70.923 I.327 J1.066 E.04497\n")
        result = audit_crossings(io.StringIO(source), chord_error=.0001)
        self.assertEqual(len(result["crossings"]), 1)
        distance = min(result["crossings"][0]["endpoint_distances"])
        self.assertGreater(distance, .005)
        self.assertLess(distance, .006)

    def test_arc_crossing_audit_does_not_hide_small_endpoint_crossings(self):
        prefix = "G90\nM83\n;TYPE:Arc bridge\nG0 X0 Y0 Z1\nG1 X1 E1\n"
        for x, expected in ((0.99, 1), (1., 0), (1.01, 0)):
            with self.subTest(x=x):
                result = audit_crossings(io.StringIO(prefix + f"G0 X{x} Y-1\nG1 Y1 E1\n"))
                self.assertEqual(len(result["crossings"]), expected)

    def test_arc_crossing_audit_checks_fitted_curves_not_their_endpoint_chords(self):
        result = audit_crossings(io.StringIO(
            "G90\nM83\n;TYPE:Arc bridge\nG0 X-1 Y0 Z1\nG2 X1 Y0 I1 J0 E1\n"
            "G0 X0 Y0.5\nG1 Y1.5 E1\n"))
        self.assertEqual(len(result["crossings"]), 1)
        self.assertAlmostEqual(result["crossings"][0]["xy"][1], 1., delta=0.0005)

    def test_arc_crossing_audit_preserves_height_and_path_boundaries(self):
        result = audit_crossings(io.StringIO(
            "G90\nM83\n;TYPE:Arc bridge\nG0 X0 Y0 Z1\nG1 X1 Y1 E1\nG1 X2 Y0 E1\n"
            "G0 X0 Y1 Z2\nG1 X2 Y1 E1\n"))
        self.assertFalse(result["crossings"])

    def test_cleaner_purge_is_not_model_material_and_preserves_absolute_e(self):
        gcode = "G90\nM82\nG0 X0 Y0 Z1\nG12 S90\nG1 Y1 E5\nG12 S91\nG0 X10 Y10 Z1\nG1 X11 E6\n"
        for exclude, expected in ((False, [5., 1.]), (True, [0., 1.])):
            with self.subTest(exclude=exclude):
                values = [dict(WORDS.findall(line)) for _, line in linearized_gcode(
                    io.StringIO(gcode), exclude_cleaner_extrusion=exclude) if line.startswith("G1 ")]
                self.assertEqual([float(move["E"]) for move in values], expected)

    def test_helical_arc_preserves_geometry_extrusion_and_source_line(self):
        moves = read_moves("G90\nM83\nG0 X1 Y0 Z0\nG3 X0 Y1 Z1 I-1 J0 E2\n")
        arc = moves[1:]
        self.assertGreater(len(arc), 2)
        self.assertTrue(all(line == 4 for line, _ in arc))
        self.assertAlmostEqual(sum(p["E"] for _, p in arc), 2.)
        for idx, (_, p) in enumerate(arc, 1):
            self.assertAlmostEqual(math.hypot(p["X"], p["Y"]), 1., places=8)
            self.assertAlmostEqual(p["Z"], idx / len(arc), places=8)
        self.assertEqual(arc[-1][1]["X"], 0.)
        self.assertEqual(arc[-1][1]["Y"], 1.)
        points = [(1., 0.)] + [(p["X"], p["Y"]) for _, p in arc]
        for a, b in zip(points, points[1:]):
            self.assertLessEqual(1. - math.hypot((a[0]+b[0])/2, (a[1]+b[1])/2), 0.005)

    def test_clockwise_full_circle_and_absolute_extrusion(self):
        moves = read_moves("G90\nM82\nG92 E5\nG0 X1 Y0\nG2 X1 Y0 I-1 E7\n")
        arc = moves[1:]
        self.assertLess(arc[0][1]["Y"], 0.)
        self.assertAlmostEqual(sum(p["E"] for _, p in arc), 2.)
        self.assertEqual(arc[-1][1]["X"], 1.)

    def test_radius_minor_and_major_arcs(self):
        for direction in ("G2", "G3"):
            minor = read_moves(f"G0 X1 Y0\n{direction} X0 Y1 R1\n")
            major = read_moves(f"G0 X1 Y0\n{direction} X0 Y1 R-1\n")
            self.assertGreater(len(major), len(minor))
            self.assertEqual(major[-1][1]["X"], 0.)
            self.assertEqual(major[-1][1]["Y"], 1.)

    def test_relative_axes_reset_position_and_absolute_center(self):
        moves = read_moves("G92 X10 Y20\nG91\nM83\nG1 X1 E1\n"
                           "G90\nM82\nG90.1\nG3 X10 Y21 I10 J20 E2\n")
        self.assertEqual(moves[0][1]["X"], 11.)
        self.assertEqual(moves[0][1]["Y"], 20.)
        self.assertAlmostEqual(sum(p.get("E", 0.) for _, p in moves), 2.)
        self.assertEqual(moves[-1][1]["Y"], 21.)

    def test_invalid_arcs_fail_instead_of_disappearing(self):
        for code in ("G18\nG2 X1 I1\n", "G2 X1\n", "G2 X5 R1\n"):
            with self.assertRaises(ValueError):
                read_moves(code)


class FeatureCourseAuditTest(unittest.TestCase):
    def run_audit(self, script, gcode, *arguments):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "feature.gcode"
            path.write_text(gcode)
            return subprocess.run(
                [sys.executable, str(Path(__file__).with_name(script)), str(path), *arguments],
                capture_output=True, text=True, timeout=10)

    def test_wall_crossing_is_detected_in_either_deposition_order(self):
        for reverse in (False, True):
            for wall_z, expected in ((0.3, 1), (0.8, 0)):
                with self.subTest(reverse=reverse, wall_z=wall_z):
                    paths = [";TYPE:Non-planar top surface\nG0 X0 Y0 Z0.2\nG1 X2 Z0.4 E1\n",
                             f";TYPE:Outer wall\nG0 X1 Y-1 Z{wall_z}\nG1 Y1 E1\n"]
                    if reverse:
                        paths.reverse()
                    result = self.run_audit("audit_nonplanar_collisions.py",
                        "G90\nM83\n" + "".join(paths), "0.05",
                        "--include-nonplanar-walls", "--fail-on-conflict")
                    self.assertEqual(result.returncode, expected, result.stdout + result.stderr)

    def test_continuous_nonplanar_wall_cannot_cross_itself(self):
        for moves, expected in (("G1 X2 Y2 E1\nG1 X0 E1\nG1 X2 Y0 E1\n", 1),
                                ("G1 X2 E1\nG1 X0 E1\n", 1),
                                ("G1 X2 E1\nG1 X1 E1\n", 1),
                                ("G1 X1 E1\nG1 X2 E1\n", 0),
                                ("G1 X2 E1\nG1 Y2 E1\nG1 X0 E1\nG1 Y0 E1\n", 0)):
            with self.subTest(moves=moves):
                result = self.run_audit("audit_nonplanar_collisions.py",
                    "G90\nM83\n;TYPE:Outer wall\nG0 X0 Y0 Z0.2\n" +
                    moves.replace(" E1", " E1 ; non-planar perimeter"), "0.05",
                    "--include-nonplanar-walls", "--fail-on-conflict")
                self.assertEqual(result.returncode, expected, result.stdout + result.stderr)

    def test_shared_endpoint_is_not_a_wall_collision(self):
        result = self.run_audit("audit_nonplanar_collisions.py", """G90
M83
;TYPE:Outer wall
G0 X0 Y0 Z0.2
G1 X2 E1
;TYPE:Non-planar top surface
G1 Y2 Z0.4 E1
""", "0.05", "--include-nonplanar-walls", "--fail-on-conflict")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_required_top_range_rejects_missing_or_insufficient_surface(self):
        for role, start_z, end_z, minimum, expected in (
                ("Non-planar top surface", 1, 3, 4, 0),
                ("Non-planar top surface", 1, 3, 6, 1),
                ("Non-planar top surface", 0, 0, 1, 1),
                ("Non-planar transition", 1, 3, 1, 1)):
            with self.subTest(role=role, start_z=start_z, minimum=minimum):
                result = self.run_audit("audit_nonplanar_gcode.py", f"""G90
M83
;TYPE:{role}
G0 X0 Y0 Z{start_z}
G1 X10 Z{end_z} E1
""", "--require-top-z-range", "1.5", "2.5", str(minimum))
                self.assertEqual(result.returncode, expected, result.stdout + result.stderr)

    def test_bucket_purge_cannot_satisfy_required_roof_coverage(self):
        # Custom tool-change G-code can inherit the last model feature label.
        # A purge at the roof's Z must not make a missing roof pass validation.
        result = self.run_audit("audit_nonplanar_gcode.py", """G90
M83
;TYPE:Non-planar top surface
G0 X0 Y0 Z39
G12 S90
G1 Y20 E1
G12 S91
""", "--require-top-z-range", "38", "41", "10")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("length 0.0", result.stdout)

    def test_wall_support_checks_the_ends_not_only_the_midpoint(self):
        for extent, expected in ((10, 0), (2, 1)):
            with self.subTest(extent=extent):
                result = self.run_audit("audit_nonplanar_gcode.py", f"""G90
M83
;HEIGHT:0.2
;TYPE:Inner wall
G0 X{5 - extent / 2} Y0 Z0.2
G1 X{5 + extent / 2} E0.1
G0 X0 Z0.4
G1 X10 E0.1 ; surface-following feature perimeter course 1
""", "--fail-on-unsupported")
                self.assertEqual(result.returncode, expected, result.stdout + result.stderr)

    def test_a_rising_wall_cannot_use_itself_as_its_support(self):
        for foundation_z, expected in ((0.2, 1), (0.4, 0)):
            with self.subTest(foundation_z=foundation_z):
                result = self.run_audit("audit_nonplanar_gcode.py", f"""G90
M83
;HEIGHT:0.2
;TYPE:Inner wall
G0 X0 Y0 Z{foundation_z}
G1 X1 E0.1
G0 X0 Z0.4
G1 X0.2 E0.01
G1 X0.4 Z0.6 E0.01 ; surface-following feature perimeter course 1
""", "--fail-on-unsupported")
                self.assertEqual(result.returncode, expected, result.stdout + result.stderr)

    def test_horizontal_feature_wall_requires_support_despite_ordinary_role(self):
        for z, expected in ((0.4, 0), (0.6, 1)):
            with self.subTest(z=z):
                result = self.run_audit("audit_nonplanar_gcode.py", f"""G90
M83
;HEIGHT:0.2
;TYPE:Inner wall
G0 X0 Y0 Z0.2
G1 X2 E0.1
G0 X0 Z{z}
G1 X2 E0.1 ; surface-following feature perimeter course 1
""", "--fail-on-unsupported")
                self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
                self.assertIn("Inner wall", result.stdout)

    def test_horizontal_feature_courses_participate_in_contact_checks(self):
        for z, expected in ((0.22, 1), (0.4, 0)):
            with self.subTest(z=z):
                result = self.run_audit("audit_nonplanar_collisions.py", f"""G90
M83
;TYPE:Inner wall
G0 X0 Y0 Z0.2
G1 X2 E0.1 ; surface-following feature perimeter course 1
G0 X1 Y-1 Z{z}
G1 Y1 E0.1 ; surface-following feature perimeter course 2
""", "0.05", "--different-feature-courses-only", "--fail-on-conflict")
                self.assertEqual(result.returncode, expected, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
