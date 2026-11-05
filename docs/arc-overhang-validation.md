# Recursive arc repair

## Larger families and short-arc overlap (September 17, Downloads reproduction)

Examined `~/Downloads/ksr_fdmtest_v4_0.2mm_PLA_Prusa CORE One HF_2h59m.gcode`.
At Z50.2 it contains 542 arc paths and 206 single-ring recursive families.
The new independent `audit_arc_overlap.py` detects 11 sustained excessive
parallel bead overlaps there. These do not require a centerline crossing:
the previous crossing-only audit was insufficient. Two additional endpoint
incursions in the file are less than 0.0005 mm, below the normal 0.005 mm
crossing tolerance; they are not evidence of large transverse crossings.

The retained generator changes:

- Keep already supported open courses intact, without splitting them merely
  because their endpoints no longer touch the original foundation polygon.
- Permit supported wide recursive families to grow beyond the 15 mm seed
  radius. Every course still passes the curvature, support and region checks.
- Check short cleanup arcs for sustained interior retracing, not just retracing
  at their ends. Include fitting/rounding clearance in the overlap test.
- Revalidate native wide arcs against the actual emitted prefix after trimming.
  Rejected arcs cannot support later arcs. Choose the nearer endpoint only when
  both leads have already been validated.

Trimming interior overlap on *all* paths was rejected: it destroyed supporting
backbones and failed coverage tests. The retained interior check is limited to
cleanup paths no longer than four bead widths. Ordinary intentional side
bonding remains allowed.

Fresh same-project comparison using repository `ksr_fdmtest_v4.3mf` (not an
assumption that its settings exactly match the supplied Downloads export):

| Z50.2 measurement | Before | After |
| --- | ---: | ---: |
| Arc paths | 493 | 418 |
| Single-ring recursive families | 198 | 114 |
| Recursive length in multi-ring families | 90.75% | 93.81% |

Verified output: `build/arc-overlap-verified.BmKWYe/plate_1.gcode`.
All six arc layers (Z6.2, 8.2, 11.2, 14.2, 17.2, 50.2) pass the excessive
parallel-overlap audit and deposited-prefix endpoint support check at 0.05 mm
minimum overlap. Z50.2 also passes full-body contact. The whole output has no
proper arc crossings among 1,169,773 linearized segments at 0.0001 mm chord
error and 0.005 mm endpoint tolerance. These are geometric diagnostics, not
physical adhesion or swept-toolhead collision certification.

Short supported cleanup strokes remain (26 below 0.675 mm at Z50.2); this is
not a claim that every short path has disappeared. A 0.1 mm raster comparison
measures 1476.71 to 1475.32 mm² total commanded bead footprint at Z50.2, with
0.95 mm² of old coverage farther than 0.2 mm from new coverage. This comparison
does not establish mesh-volume completeness or prove zero holes.

The open-course and recursive-overlap regressions were run with the three old
generation policies restored: both test cases failed (43 assertions), then
the fixes were restored. All 32 Python audit tests pass, including an actual
curved overlap fixture and intentional side-bonding controls.
The final rebuilt C++ arc suite passes all 34 cases / 35,674 assertions.
App build: `build/arm64-ninja/src/OrcaSlicer.app` (not installed into
`/Applications`). The test runner logs a missing local nozzle-info cache;
the build also warns about the Homebrew assimp deployment target. Neither
prevented this local build or the tests from completing.

The enlarged 80 mm pocket retains the same coverage/coherence requirements.
An experimental assertion requiring a child radius above 15 mm was removed:
after overlap and support trimming, that particular pocket does not sustain
such a child. A larger permitted radius must not override those safety checks.
Measured improvement on the KSR export is reported above instead.

## Alternating direction and remaining fragmentation (September 17)

The final direction pass now recognizes native previous-layer foundation
polygons and actual same-height support paths. Previously it ignored those
polygons and could reverse an already serpentine arc back to the far side.
Both supported leads are checked before choosing the nearer end; raised wall
obstacles do not count as anchors. Generation order is unchanged and no
connecting extrusion is invented across gaps.

The existing travel regression was extended with native foundation inputs.
It failed before the fix (about 11.1 mm travel instead of 0.48 mm for one
consecutive pair). With the fix, all 33 arc cases / 35,669 assertions pass.
The 29 Python audit tests also pass. This is generator-level validation, not
a new KSR reslice or physical-print qualification.

**Larger-family selection remains unresolved.** Experiments comparing candidate
centers by deposited length or ring count, prioritizing broader centers or
longer frontiers, and changing family-size tiers were not retained. Broader
centers introduced coverage failures; other changes reduced coherence or left
the broad-pocket fixture unchanged. As a diagnostic, about 49% of recursive
length in that fixture belongs to consecutive families of eight or more rings.
This is not a universal printability threshold: fewer rings can still have
large radii and long supported paths. A future selection change must compare
actual path lengths/radii and model coverage as well as family coherence,
deposited-prefix support, and slicing cost. Passing support checks alone does
not establish that avoidable fragmentation has been fixed.

## Uniform narrow bridges (September 17, Z17.2 follow-up)

The narrow KSR strip now grows regularly translated arches from each supported
end instead of filling most of its length with small recursive restarts. Native
axis-aligned narrow strips select their long direction only when a real end
foundation exists. Each starter uses embedded support, a radius matching an
actual concentric seed course, and its own foundation-facing start direction.
The arch families keep that radius and regular pitch; only their central
junction and small edge residuals need recursive cleanup. Legacy callers and
the wide concentric planner retain their previous behavior.

Merely preferring six-ring children did not solve the strip. Direction-only
and partial handoff trials were also rejected. The two-ended regression caught
the reversed starter using the first family's direction, producing coverage
and coherence failures at both tested bridge angles. It now passes without
relaxing the 0.5% substantial-hole and 65% coherent-course requirements.
The coherence check recognizes regularly translated equal-radius arches as
well as expanding concentric courses; it does not count arbitrary new centers
as coherent. The fixture includes retained end walls, one/two foundations and
both perpendicular bridge angles.

Rebuilt app: `build/arm64-ninja/src/OrcaSlicer.app`.
Fresh KSR export: `build/arc-narrow-oriented.7uRAft/plate_1.gcode`.
Visual centerline inspection: `build/arc-narrow-oriented.7uRAft/paths.svg`.

- Z17.2: 162 paths instead of 735. The exported pattern has two regular arch
  runs with localized central cleanup, rather than scattered restarts along
  the whole strip. Short supported cleanup paths remain.
- All six arc layers pass the 0.05 mm endpoint-overlap diagnostic; Z50.2 also
  passes full-body contact. Intentional narrow bridge spans are not required
  to have lateral support along their complete bodies.
- No proper arc crossings in 1,171,342 linearized segments at 0.0001 mm chord
  error. This is not a complete swept-toolhead collision simulation.
- Z17.2 coarse footprint increases from 245.03 to 252.71 mm². Only 0.20 mm²
  of its previous footprint lies more than 0.2 mm from new extrusion. Z8.2,
  11.2 and 14.2 also gain coverage, with corresponding distant differences
  of 0.08, 0.17 and 0.14 mm². Z6.2 and Z50.2 footprints are unchanged.
- All 33 arc cases / 35,629 assertions and 29 Python auditor cases pass.
- Slice: 83.42 s wall / 101.61 s user while tests were also running.

Raster coverage is comparative, not exact mesh-volume certification. No
physical print has been tested; the small central cleanup is intentional and
must still be inspected in print order, not judged from length alone.

## Lower bridge anchor repair (September 17)

Native narrow bridges now use the same overlap allowance as concentric
growth. If clipping leaves a weak terminal tip, the primary arch may trim up
to one bead width from each end along its existing curve to obtain a fully
supported lead. Both leads are revalidated before committing the path. This
does not require lateral body contact across intentionally spanning bridges.
Legacy callers without native foundation geometry retain their old allowance.

An intermediate trial only tightened the allowance. Although its unit tests
passed, it lost 122.67 mm² beyond 0.2 mm of the old footprint at Z11.2; that
trial was rejected. Coverage comparison is essential alongside anchor checks.

Retained export: `build/arc-lower-trim.hyPa7V/plate_1.gcode`.
Rebuilt app: `build/arm64-ninja/src/OrcaSlicer.app`.

- All six arc layers pass the 0.05 mm overlap endpoint diagnostic. Z50.2 also
  passes full-body contact. Lower intentional bridge spans are not subjected
  to that full-body criterion. The previous export still produces 28, 61, 74
  and 73 deficient leads at Z8.2, 11.2, 14.2 and 17.2 with the same auditor.
- Zero proper crossings in 1,146,383 linearized planar segments.
- Z6.2 and Z50.2 arc footprints are unchanged. At Z8.2, 11.2, 14.2 and 17.2,
  arc footprint area changes by +0.46, -2.81, -2.39 and -0.59 mm². The old
  footprint beyond 0.2 mm of new material is 0.42, 2.50, 2.20 and 1.32 mm²
  respectively. These local coverage differences remain; the checker does
  not establish exact volumetric coverage or physical printability. Including
  all extrusion roles at Z11.2 gives 770.50 → 768.13 mm² and 2.19 mm² beyond
  that distance. No claim of zero missing model volume is made.
- All 33 arc cases pass, now with 34,967 assertions. The one-sided anchor
  matrix covers native foundations and independently buffered bonded leads
  as well as legacy callers, overlap settings and raised-wall obstacles.
- All 29 Python auditor cases pass. A new real-coordinate regression exposes
  a checker defect: 0.005 mm chord sampling can rotate a terminal tangent
  enough to reject a supported arc. Default sampling is now 0.0001 mm; the
  regression also verifies convergence at 0.000025 mm and rejects reversed
  printing order. Failures report the first unsupported sample coordinate.
- Native reslice completed in 82.58 s wall / 105.91 s user, with tests also
  running. This is not an isolated performance benchmark. No physical print
  has been tested.

## Bonded recursive restarts (September 16, late evening)

Native concentric growth now requires overlap along both leads and the body,
using the actual bead widths of retained supporting walls. Weak endpoints may
extend along their existing circle into an already supported rim, stopping as
soon as their lead is bonded and before crossing an occupied centerline. The
fill region itself is not enlarged: that trial generated unrelated rim arcs.
Recursive restarts require separate start/end runways. Families are attempted
in tiers of 12, 6, 3, 2 and 1 courses; a rejected initial radius may try the next
two radii, but an established family cannot skip an unsupported ring.

Rebuilt app: `build/arm64-ninja/src/OrcaSlicer.app`.
Fresh native KSR export: `build/arc-bond-final.Z305ov/plate_1.gcode`.

- Z50.2 passes `audit_arc_support.py --z 50.2 --min-overlap .05` with zero
  deficient ends or bodies, versus 82 ends / 73 bodies in the supplied 3h16m
  export. This is a geometric diagnostic, not a physical adhesion guarantee.
- 493 paths instead of 1,067; 90.0% of recursive length belongs to consecutive
  families of at least three radii. There remain 29 short paths; they pass the
  deposited-prefix contact/overlap checks. Short length alone is not treated
  as proof of an unsupported path.
- Zero proper arc crossings in 1,143,011 linearized planar segments, using
  0.0001 mm G2/G3 chord error.
- Z50.2 coarse footprint: 1425.62 mm² versus 1427.97 mm². Of the original
  footprint, 0.71 mm² lies more than 0.2 mm from the replacement. The SVG
  comparison shows small local changes, not large stranded pockets. This is
  not a mesh-volume completeness certificate.
- 33 arc cases / 22,039 assertions and 28 Python auditor cases pass. The new
  width regression independently checks support with 0.1, 0.45 and 0.8 mm
  retained walls. Coverage, coherence, body support and crossing tests remain.
- Slice: 102.34 s wall / 135.58 s user. This stricter implementation is slower
  than the preceding 39.64 s mere-contact implementation; optimization remains.

**Remaining limitation:** the stricter endpoint diagnostic passes Z6.2 but
flags 28, 61, 74 and 73 leads at Z8.2, 11.2, 14.2 and 17.2 respectively.
The original supplied file has exactly the same counts on those layers.
These intentionally spanning bridge layers retain their existing generation
policy; this repair must not be represented as resolving their weaker overlap
or as validating the entire KSR print. No physical print has been tested.

## Larger-family-first follow-up (September 16, evening)

The three-course preference still fragmented the KSR surface. Native
concentric growth now begins with six-course families and progressively relaxes
that minimum only after a complete stalled pass. Existing clipping, deposited
support, rollback and coverage logic is unchanged. Trying broader centers first
was also evaluated, but rejected because the family-coherence regression fell
from 84% to 59%; that center-order change is not retained.

Export: `build/arc-coherent.yPAx1Z/plate_1.gcode`.

- Z50.2: 708 paths instead of 1,067; multi-course recursive length rises from
  81.9% to 90.5%. Support audit reports zero unsupported ends or bodies.
- No proper arc crossings detected in 1,145,065 planar segments.
- Slice: 39.64 s wall / 59.45 s user on this machine.
- Z50.2 coarse footprint: 1427.91 mm² versus 1427.97 mm²; 0.47 mm² of the
  former footprint is more than 0.2 mm from its replacement. This is a local
  pattern comparison, not proof of complete model coverage.
- A new independent auditor regression checks a short quarter-circle followed
  by its larger neighbor: later contact must not qualify the short path. It
  also checks that swapping their order does not magically qualify the larger
  arc's unsupported body.

The subsequently supplied `ksr_fdmtest_v4_0.2mm_PLA_Prusa CORE One HF_3h16m.gcode`
matches the preceding 1,067-path output at Z50.2. It passes mere-contact checks,
but `audit_arc_support.py --min-overlap .05` reports **82 deficient end leads
and 73 deficient bodies**. One short arc starts at line 381353, XY146.577,
101.135, and ends at XY146.902,100.677. The larger-family export still has
47 deficient end leads and 33 deficient bodies under this stronger diagnostic.
The 0.05 mm threshold is a bonding diagnostic, not a physically calibrated
universal printability threshold.

Two stricter generation trials (full-body overlap and lead-only overlap) were
rejected: both drove the independent fixture's multi-course fraction down to
approximately 24%, below the unchanged 65% regression requirement. Neither is
retained. Simply reducing allowed support distance does not solve this problem
without disrupting recursive growth. The supplied weak-contact defect remains
open; zero mere-contact flags must not be reported as its resolution.

32 arc cases / 21,447 assertions pass for the retained larger-family change.
27 auditor cases pass, including new future-arc and touching-versus-overlapping
regressions. All six layers pass endpoint-only contact checks; Z50.2 additionally
passes full-body contact. These are explicitly weaker than the new overlap
diagnostic, which still fails. No physical print was tested.

## Concentric restart families (September 16)

Native concentric recursive growth now tries families with at least three
distinct radius courses before falling back to two-course families, then
single-course gap cleanup. Each accepted family keeps its deposited-prefix
order. Failed trial families roll back their support-index entries so later
arcs cannot borrow nonexistent support. Multi-course trials can grow beyond
the nearest model-edge distance; clipping and the existing contact/collision
checks determine which courses actually fit. Intentionally spanning narrow
bridge arches retain their existing strategy.

Rebuilt app and native KSR export: `build/arc-uniform.yCEctZ/plate_1.gcode`.

- Z50.2 recursive extrusion length in consecutive families of three or more
  radii increases from 45.8% to 81.9%, measured independently from exported
  geometry by `audit_arc_families.py`. Total arc paths decrease from 3,205 to
  1,067. Singles remain for residual gaps; this is not a guarantee that every
  clipped fragment has three neighbors.
- Zero unsupported starts/ends or bodies at Z50.2. All five lower arc layers
  pass endpoint checks (not full-body checks for intentional bridge spans).
- Zero detected proper arc crossings across 1,148,803 planar segments at
  0.0001 mm G2/G3 chord error.
- 32 arc regression cases / 21,447 assertions and 25 Python auditor cases
  pass. The family-coherence regression failed before this change: only 29%
  of recursive length belonged to multi-course families, against its 65%
  minimum. Coverage, crossing and independent body-contact tests remain.
- KSR slicing: 47.21 s wall / 67.88 s user, down from 72.49 s / 114.19 s for
  the previous body-support-validated implementation, on this machine.
- Coarse arc footprint grows from 1426.39 to 1427.97 mm² at Z50.2 and from
  1040.63 to 1043.30 mm² at Z6.2. Four intervening arc layers are unchanged.
  Local differences remain: 1.06 mm² at Z50.2 and 0.05 mm² at Z6.2 lie more
  than 0.2 mm from the replacement footprint. Visual footprint inspection
  shows scattered small differences, not the earlier large missing pockets.

The family auditor measures grouping, not support or printability. The raster
comparison is not mesh-volume certification, and no physical print was tested.

## Reopened: outer arcs precede their inner support (September 16, Z50.2)

The endpoint-only checkpoint below missed support along the arc body. In its
Z50.2 export, a 0.388 mm primary radius is followed by a 6.587 mm radius; smaller
rings are filled afterwards. The strengthened auditor identifies 90 bodies
without already-deposited contact, despite zero unsupported-end flags.

`audit_arc_support.py` now checks bodies by default. `--endpoints-only` retains
the narrower diagnostic for intentional bridge spans; it cannot certify a
concentric overhang's growth order. Narrow, translated bridge arches and
concentric overhang growth have different foundation requirements.

The new generator regression checks entire paths against independently buffered
previous paths, not the generation-time support index. A separate auditor
regression has two anchored concentric arcs: outer-first must fail and
inner-first must pass. Coverage and crossing checks remain required alongside
the contact checks. Development trials which simply lost coverage were rejected.

Verified build and KSR reslice: `build/arc-order-verified.yxZM6i/`.

- Z50.2: 3,205 arc paths, zero unsupported starts/ends and zero unsupported
  bodies against the previously deposited prefix. Body contact includes actual
  bead caps; start/end lead qualification still excludes free caps. The earlier
  overly conservative body-cap exclusion reported 122 rather than 90 failures.
- All six arc layers pass the endpoint audit. The five lower layers use the
  endpoint-only diagnostic because they include intentionally spanning narrow
  bridge arches; this is **not** a claim of full-body lateral support there.
- Zero detected arc centerline crossings in 1,163,971 planar segments using
  0.0001 mm chord error for fitted G2/G3 moves.
- 32 arc tests / 21,445 assertions and 24 Python auditor tests pass. The narrow
  bridge matrix now includes explicit native foundation masks as well as the
  legacy standalone anchor context.
- The generator checks full-body contact for native concentric fronts, retries
  stalled parents only after real deposition, and tests additional center
  offsets when the original three cannot grow a supported child. Complete arcs
  remain intact; experimental partial-body clipping and forced parent-circle
  continuation were rejected by coverage tests and removed.
- Main app rebuilt at `build/arm64-ninja/src/OrcaSlicer.app`. The stricter KSR
  slice takes 72.49 s wall / 114.19 s user, versus ~19 s for the earlier
  endpoint-only slice. Further performance work must retain the body-order and
  coverage checks rather than bypass them.

Coarse arc-footprint comparison against the previous endpoint-verified export:

| Arc Z | Before (mm²) | After (mm²) |
| --- | ---: | ---: |
| 6.2 | 1040.53 | 1040.63 |
| 8.2 | 103.37 | 103.73 |
| 11.2 | 149.83 | 150.69 |
| 14.2 | 198.13 | 198.76 |
| 17.2 | 245.18 | 245.62 |
| 50.2 | 1423.90 | 1426.39 |

Large missing pockets found during development are absent in this comparison.
Small footprint differences remain: 0.67 mm² at Z6.2 and 0.53 mm² at Z50.2
lie more than 0.2 mm from the replacement footprint on the 0.1 mm raster.
This comparison is not a mesh-volume certificate, nor physical print validation.

## Historical free-tip and recursive-frontier checkpoint (September 16)

This endpoint checkpoint is superseded by the body-order finding above. It
superseded the earlier historical results below. The generator
distinguishes the side of a deposited bead from its free end caps, including
retained wall tips. Explicit lower-layer foundations are separate from
current-layer wall-contact masks. Cached butt-ended contact corridors include
a terminal guard for fitting/rounding, without shrinking the sides of correctly
spaced neighboring arcs. Contact uses extrusion width, not the overlap-reduced
centerline pitch.

The repair also regenerates coverage rather than merely deleting unsupported
paths: original foundations and inner exposed rings can seed new families;
long narrow regions retain their validated parent-before-child order and get
multiple seed attempts; the frontier and child budgets are consistent. Rings
are generated on demand so a failed seed does not clip a whole unused family.
Recursive fragments use the tighter fitting tolerance even in narrow regions.

Final-build verification artifacts and exported G-code:
`build/arc-verified.YY4aKL/`. The preceding verification is retained in
`build/arc-final-contact.ouF4Vq/`.

- Supplied 18:53 KSR export: 86 unsupported-end flags on Z50.2.
- Fresh KSR export: 5,196 arc paths across Z6.2, 8.2, 11.2, 14.2, 17.2 and
  50.2; zero unsupported-end flags on every layer.
- Zero detected centerline crossings in 1,151,091 planar segments, auditing
  fitted G2/G3 moves at 0.0001 mm chord error.
- All 31 arc cases / 21,017 assertions and all 22 Python auditor tests pass.
  Tests now include butt-ended prefix support, zero/25% overlap, one/both
  short-end foundations, retained walls, long cantilevers, large residual
  pockets, empty explicit foundations and narrow fitted/rounded arcs.
- The new large-pocket coverage test failed during development (21.86 mm²
  missing versus a 7.5 mm² limit). It and the earlier coverage/contact tests
  were retained, not relaxed to accept deleted paths.
- Final-build repeat: 18.97 s wall / 39.54 s user (preceding run: 20.16 s).
  This is slower than the earlier
  unsafe ~10 s result, but substantially below the 80 s intermediate repair.
- Main app rebuilt at `build/arm64-ninja/src/OrcaSlicer.app`.

Coarse 0.1 mm raster footprint comparison against the supplied G-code:

| Arc Z | Supplied area (mm²) | New area (mm²) |
| --- | ---: | ---: |
| 6.2 | 1038.20 | 1040.53 |
| 8.2 | 91.47 | 103.37 |
| 11.2 | 135.17 | 149.83 |
| 14.2 | 179.53 | 198.13 |
| 17.2 | 227.19 | 245.18 |
| 50.2 | 1422.75 | 1423.90 |

These are geometric regression checks, not a guarantee of physical bridging,
cooling, adhesion, or full toolhead clearance. The footprint comparison uses
the old export, not the mesh as a perfect reference; changed bead footprints
still leave small local differences (0.15–1.70 mm² per layer more than 0.2 mm
from the replacement footprint). No physical print has been validated.

## Reopened: free-tip propagation is unresolved (September 16, 18:53 export)

The later screenshot and overwritten KSR G-code invalidate the earlier
printability conclusion. A distance-to-segment test includes a round disk at
every free endpoint. It can therefore accept an outward staircase of tips
without a previously printed bead alongside their terminal sections.

- The auditor now distinguishes free path ends from internal tessellation
  vertices and closed loops. It flags 86 endpoints on Z50.2 of the 18:53
  export; see `build/arc-endcap-validation.NqF3R9/support.json`.
- The C++ prefix regression now uses butt-ended contact corridors, retaining
  its coverage/retrace checks. Against the restored generator it fails 388
  of 7,464 assertions in that one parameterized case. It is deliberately not
  disabled or weakened. The previous green suite missed this property.
- All 22 Python auditor tests pass, including a control with a later path
  extending diagonally beyond an earlier free tip.
- Two attempted production repairs were rejected: rejecting tip-supported
  paths lost substantial coverage; trimming tips still failed contact and
  coverage checks. Those generator changes were reverted and the app rebuilt
  with its prior behavior. **The defect is not fixed.**

The growth planner must regenerate supported families/residual regions, not
only remove unsupported fragments or trim their ends. New validation must
prove lateral contact in emitted order together with retained area coverage.
The older checkpoints below remain historical measurements, not certification
that this export is printable.

## Follow-up: raised-wall anchors and short fragments (September 16)

The user's `ksr_fdmtest_v4_0.2mm_PLA_Prusa CORE One HF_3h12m.gcode`
exposed a gap in the earlier validation. At Z50.2, the deposited-prefix audit
flags trailing leads at source lines 342318 and 342360. The planner's trace
classified both as supported by an anchor polygon, while the nearest retained
wall geometry included raised brick/variable-Z paths. Flattening those walls
into XY does not establish an anchor at the arc's printing height.

The fill planner now carries separate collision and support paths. All
retained walls remain obstacles; only level walls contribute same-layer
support centerlines and anchor rims. Short-fragment support is also checked
again after retrace trimming, not only before it.

`audit_arc_support.py` checks both endpoint leads against previously emitted
extrusions, withholding the current path until both checks finish. It tests
actual G-code ordering and reconstructs G2/G3 moves. Its configurable maximum
support drop defaults to 0.25 mm for this 0.2 mm project. This is a contact
diagnostic, not a nozzle-envelope or cooling simulation. Positive controls
cover future support, self-support and raised walls. On the supplied file,
the five lower arc layers pass and Z50.2 fails.

Final verification:

- All 29 arc cases / 20,270 assertions pass, including the raised-obstacle
  variants and full-length contact for short fragments after trimming.
- All 21 Python auditor tests pass, including a fitted terminal crossing
  extracted from an intermediate export.
- Rebuilt `build/arm64-ninja/src/OrcaSlicer.app` and resliced the updated KSR
  project to `build/arc-ksr-final-support.cXVkYa/plate_1.gcode`.
- Deposited-prefix endpoint checks pass on all six arc layers (3,317 paths);
  Z50.2 contains 1,635 paths and zero flagged ends.
- Export crossing audit: zero crossings among 1,178,692 planar segments at
  0.0001 mm chord error. Clipped disk families now fit to 0.001 mm, while
  translated narrow arches retain 0.005 mm. The 0.012 mm contact guard is
  preserved; increasing it changed family selection and was rejected after
  exported-G-code checks found new failures.
- Slice time: 9.94 s wall / 20.77 s user; filament: 45.76 g.
- Coarse 0.1 mm XY raster comparison at Z50.2 measures 1,422.75 mm² of arc
  coverage versus 1,422.95 mm² in the supplied file. Changed paths lose
  6.52 mm² and gain 6.32 mm² of sampled footprint; 0.85 mm² of old footprint
  lies more than 0.2 mm from the new footprint. This comparison is not a
  gap-free-model certificate: the baseline itself contains invalid paths,
  and nozzle geometry, cooling and physical bonding remain unvalidated.

The full app rebuild initially exhausted disk space during archiving, then
succeeded on retry with lower concurrency. No project files were deleted.

## Latest checkpoint: deposited-prefix support (September 16)

The wider, one-sided 12 mm regression reproduced seven failures: short arcs
had an unsupported leading/trailing section in final print order, or retraced
an earlier bead. Narrow-only cases did not expose the failure. Generation-time
support was invalidated by later cleanup/reordering of the wider arc families.

Supported recursive wide families now retain their committed growth order,
like narrow arches: trim and validate a path before it becomes support for a
child, and require both endpoint lead sections to touch already deposited
material. Minimum useful length is checked again after child trimming.
Other explicitly supported emission paths rebuild dependencies from the
deposited prefix, deferring blocked children rather than inventing an
unsupported root. Retrace trimming includes its sampling interval allowance.

- Regression first failed against the previous implementation; all 29 arc
  cases now pass (16,538 assertions), including support order, coverage,
  curvature, retracing and fitted/rounded crossings. Coverage thresholds were
  not relaxed. An intermediate final-filter-only implementation was rejected
  because it removed unsupported fragments but left large uncovered areas.
- All 18 Python G-code auditor tests pass.
- Rebuilt the main `build/arm64-ninja/src/OrcaSlicer.app` and resliced
  `ksr_fdmtest_v4.3mf` to `build/arc-ksr-prefix.1fqAkY/plate_1.gcode`.
  Slice time: 16.39 s wall / 30.26 s user; total filament: 45.78 g.
- Export audit at 0.0001 mm chord error: zero centerline crossings across
  1,297,864 planar segments (`crossings.json` beside the export).

Endpoint support/coverage is tested at generator level; the export audit
checks crossings, not a complete physical bead-support or cooling simulation.
No physical print validation is claimed.

## Latest checkpoint: fitted-arc terminal clearance (September 16)

Contact trimming now reserves 0.012 mm normal clearance to the occupied
segment, rather than backing off 0.001 mm along an arbitrarily shallow
approach. This includes twice the shared 0.005 mm arc-fitting tolerance and
coordinate rounding allowance. Trimming happens before paths become support
for their children, as well as in the final collision sanitizer.

- The new fitted/rounded-path regression fails against the previous core
  library and passes with the repair. It exercises two resolutions and two
  coordinate origins, reconstructing rounded G2/G3 geometry.
- All 29 arc cases / 15,938 assertions pass, including coverage and bonding.
- `audit_arc_gcode.py` independently checks exported moves, tessellating arcs
  at 0.0005 mm chord error and allowing 0.005 mm endpoint contact. Its positive
  controls detect short terminal crossings and curved crossings that a single
  endpoint chord misses; all 18 Python auditor tests pass.
- With this finer audit, the previous export has two confirmed crossings
  (line pairs 111673/110852 and 111691/108828). The older coarse audit's 17
  candidates below were not all confirmed: its chord error was too large
  relative to terminal clearance.
- New export: `build/arc-ksr-clearance.3yKjAM/plate_1.gcode`, zero crossings
  among 703,967 tessellated planar segments. Filament remains 46.04 g versus
  46.01 g previously. A tighter 0.0001 mm chord-error audit also reports zero
  crossings across 1,304,330 segments. Slice time was 22.31 s wall / 41.67 s user while the
  regression suite ran concurrently, not an isolated performance comparison.
- The main arm64-ninja app and test executable were rebuilt successfully.

This is a measured repair for the supplied reproduction and regression
fixtures, not a proof of nozzle clearance or physical printability on every
model. The audit does not certify collinear overlap, bead coverage, or machine
clearance. Arc fitting tolerance alone is not a global curve-distance bound.

## Latest checkpoint: supported-frontier growth (September 16)

Narrow slots now start with a returning concentric family bounded by the
short dimension. Subsequent families grow from the deposited frontier, using
up to eight separated seed attempts and the existing parent/path budgets.
Each child is trimmed and rechecked against the committed prefix immediately
before acceptance, including siblings accepted from the same radius. This
planner preserves that order instead of subsequently trimming its parents.
The narrow-slot construction applies with recursive refinement disabled too;
that option still controls refinement of the ordinary wider root families.

The final nearest-anchor scan uses the existing spatial index, with an exact
full-search fallback when no nearby anchor is found. Coverage subtraction and
seed-window pruning remain bounded; progress and cancellation are retained.

- Expanded one-sided endpoint regression: widths 2.4, 3.2, 4.5, 5.8 and 6.8 mm.
- Arc suite: 15,810 assertions / 28 cases pass, 6.11 s wall time. This includes
  more fixtures than the previous 4.18 s run, so those suite timings are not
  a like-for-like performance comparison.
- Full KSR export: `build/arc-ksr-frontier.mIf2eG/plate_1.gcode`, 17.87 s wall /
  35.16 s user CPU versus the preceding optimized export's 17.95 / 34.42 s.
- Arc extrusion tessellation: 28,352 segments versus 28,332 previously;
  112 old segments removed and 132 added (endpoints rounded to 0.00001 mm).
- **Not a completed real-model acceptance:** a separate spatial-grid audit of
  369,272 planar extrusion segments still finds 17 proper arc-to-arc crossings
  farther than 0.005 mm from tessellated segment endpoints. Examples at Z6.2:
  source-line pairs 111610/111022, 111616/110502 and 111639/109526. The audit
  linearizes G2/G3 at 0.005 mm chord error; these need investigation against
  fitted arcs and whole-path endpoint contact, not dismissal based on the
  passing small fixtures. A follow-up check places these intersections
  0.00508–0.01166 mm from the later G2/G3 command's endpoint. This narrows the
  investigation to terminal contact/fitting/quantization, but is not yet a
  validated repair.

ContourZ's separate eight-instance exception is repaired and tested for
translated and mixed-rotation copies (9 assertions). Surface-following support
ironing now also handles multiple copies, tested with one and eight instances.
The first sidebar rebuild failed while archiving the GUI with errno 28 (disk
full). After free space rose to 32 GB, the complete app rebuilt successfully,
including the compact, collapsed-by-default Print bindings section with the
same StaticBox/Label header treatment as Printer and Project Filaments.

The sections below retain the earlier measurements and investigation history.

The supplied reproduction is `ksr_fdmtest_v4.3mf`. Its baseline is retained in
`build/arc-ksr-validation.rCVsPt/`: 150 fitted arc extrusion commands in Arc bridge
sections, including 77 with radius above 15.1 mm (maximum 32.533 mm). These are
command counts, not independent extrusion paths.

The reference implementations grow a concentric family to its boundary/radius
limit, then recursively seed more families on previously deposited arcs:

- https://github.com/stmcculloch/arc-overhang
- https://github.com/nicolai-wachenschwan/arc-overhang-prusaslicer-integration

Changes here apply the existing local 15 mm child-family policy to ordinary
root families when recursive fill is enabled. This is an engineering limit,
not a universal experimentally established printability threshold. Clipped
short paths still require curvature, anchoring, and intersection checks.

A rejected preferred seed now permits other well-separated seed candidates.
After clipping and final emission ordering, substantial residual pockets may
receive additional families, anchored to the actual surviving deposited paths.
Refinement is depth-bounded and retains the existing obstacle/support checks.
The non-recursive root behavior is unchanged.

Tests cover wide and holed surfaces, bounded radius, surviving bead coverage,
child attachment, and proper crossings. The crossing auditor has positive
controls for separate-path and self intersections. Completely straight long
paths count as infinite radius, rather than silently passing curvature checks.

These geometric tests do not certify a physical free-air print. Validate the
full KSR export and inspect its coverage and deposition order before calling
the model repaired. The earlier baseline used `-O0`; its runtime is not a
release-build performance measurement.

Latest targeted acceptance run: 1,223 assertions pass across the wide-roof,
holed-surface and exported narrow-bridge cases. An earlier broader run passed
26 of 27 cases and caught a detached narrow-arch start. Preserving that planner's
original seed policy fixes the failing case; all three affected cases were
rerun after rebuilding. Unrelated cases from the broader run were not rerun
after that final restriction.

The full-project slice in `build/arc-ksr-anchored.tD8y2V/` completed and exported
`plate_1.gcode`. It is not yet an accepted real-model result: full output
coverage and deposition order still need auditing.

## Long infill pause investigation

The unoptimized pre-performance executable takes 183.12 seconds (178.08 seconds
user CPU) for `[BoundedRadius]`, passing four assertions. Sampling the running
GUI while it displays 35% shows arc support checking, not an idle worker.

Verified changes in the rebuilt executable:

- Skip exact line tests for disjoint padded bounds inside the arc spatial grid.
- Reuse the unchanged deposited-parent support envelope across residual pockets.
- Report infill layer completion and arc growth, recursion, validation, and
  residual-pocket stages, throttled to four updates per second per object.
- Check cancellation inside the arc loops, including when UI updates are throttled.
- Build with `-DSLIC3R_UNOPTIMIZED_RELWITHDEBINFO=OFF` to retain `-O2` rather than
  the existing debug configuration's forced `-O0`. This does not enable fast-math.

New tests require reporting to leave the generated paths unchanged, cancellation
to escape from inside family generation, and infill to report beyond its initial
35% status.

Optimized build verification (September 16):

- App and both test executables rebuilt successfully, including the final
  incremental cancellation checkpoints.
- `[Progress]`: 25 assertions / 3 cases pass.
- `[ArcOverhang]`: 5,323 assertions / 28 cases pass in 4.18 seconds.
- `[BoundedRadius]`: 4 assertions pass in 1.89 seconds (1.86 seconds user CPU),
  versus 183.12 seconds / 178.08 seconds user CPU before. This measures the
  combined optimization/build changes, not an isolated algorithmic speedup.
- Full `ksr_fdmtest_v4.3mf` export completes in 17.95 seconds (34.42 seconds
  user CPU), in `build/arc-ksr-optimized.ffORf5/`. The previous unoptimized
  export ran approximately 20 minutes, based on its launch and output timestamps.
- The real slice reports arc growth, recursion, validation, residual pockets,
  and completion at 39%, rather than retaining only the initial 35% message.
- Comparing canonical tessellated Arc bridge extrusion moves with the previous
  export gives 28,332 segments in each, with zero differences in their multiset
  of endpoints (rounded to 0.00001 mm) and extrusion (0.0000001 mm). The full
  G-code is not byte-identical; ordinary wall moves also differ. This comparison
  preserves existing arc geometry, not proof of physical printability.
- G-code parser: 15 tests pass. Non-planar core: 4,272 assertions / 26 cases pass.
  The separate non-planar pipeline acceptance still fails all three wall-order
  expectations (236/239 assertions pass); the missing non-planar roof is not
  resolved by this performance repair.

## Unsupported arc-end reproduction (September 16)

Re-read the original and integrated reference descriptions and source:

- https://github.com/stmcculloch/arc-overhang
- https://github.com/nicolai-wachenschwan/arc-overhang-prusaslicer-integration
- https://raw.githubusercontent.com/nicolai-wachenschwan/arc-overhang-prusaslicer-integration/main/prusa_slicer_post_processing_script.py

The reference grows concentric families from a previous-layer start boundary,
subtracts their occupied region, and seeds recursive families on earlier arcs.
Its endpoint-extension policy is explicitly intended to improve bonding. A
recursive start away from the original body is not inherently invalid; it must
inherit an already deposited support chain.

The existing one-sided narrow-overhang regression checked only each departure
runway. Checking the trailing runway too reproduces a failure on arc zero:
start (0.225, 10), end (4.27431, 9.99927), with support only at X0–0.6.
Subsequent arcs inherit that unsupported seed. This new assertion remains
enabled and failing; the prior green arc suite did not establish both-end support.

Experiments not retained:

- Rejecting every final arc without two supported runways leaves large holes
  (4/28 cases fail, including coverage).
- Forcing a concentric backbone for every narrow region fixes the new end test
  but fails other curvature/coverage cases (3/28 cases, 95 assertions).

The production repair still needed is a supported seed/frontier planner with
bounded-curvature children, valid endpoint bonding after final clipping, and
coverage recovery. Do not disable coverage assertions or ship the rejection
experiment as a repair. Exported-path crossing validation is also still needed.

The earlier support-envelope cache now refreshes when `deposited.size()` changes:
recursive children append to `arc_prior_paths`, so the support is not invariant
after an accepted sibling. Reuse remains valid only between such additions.
