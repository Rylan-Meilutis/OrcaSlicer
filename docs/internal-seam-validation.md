# Internal seam and brick preview validation

September 16 checkpoint:

- Eligible non-scarf outer loops start along the latter half of a reserved
  inner-wall tail, connect outward, follow the outer wall, then return to the
  former half of that tail and finish at the internal starting point.
- The complete reserved tail is omitted from the ordinary inner loop, in
  both inner-first and outer-first orders. The two halves are not retraced.
- Both connectors are checked against the inner and outer wall centerlines
  and each other. Unsafe connections retain the non-extruding approach.
- Wiping follows the inner finish backward, not a stale outer-loop path.
- Scarf qualification remains in place. Internal seam reservation also
  supports constant-height brick courses at different Z, including raised
  inner-only bricks beside lower outer walls. Tails keep their original
  height, flow and brick role. Entry/exit connectors use the actual endpoint
  Z and the average cross-section of the adjoining courses, with XYZ length.
  Their rise/run is limited to 1:2; wall intersection checks still apply.
  Arbitrarily contoured/non-planar loops remain excluded.
  Enabling top-surface processing no longer disables internal entries on
  unrelated flat walls: eligibility is checked against the actual paths.
- A raised brick path now retains its Brick wall preview role even when it
  also carries shrinkage compensation or a surface modifier. Flat realignment
  caps still use their ordinary role, because they are no longer raised bricks.
- A separate multi-path seam-splitting bug discarded brick, transition, inset
  and compensation metadata while retaining the raised Z. Splits now use the
  complete path metadata constructor; the regression failed before this fix
  and passes afterward. This repairs extrusion semantics as well as colour.

Exported-G-code regressions cover Classic/Arachne, both wall orders, zero and
nonzero seam gaps, straight/curved walls, unchanged non-connector wall material,
continuous positive extrusion through the seam, coincident internal start/end,
no repeated collinear wall extrusion, and no proper centerline crossings.
Scarf fallback/qualification and variable-height material-preservation tests
also pass. The combined seam/brick-preview filter passes 677,100 assertions in
8 cases. The brick-preview test failed before the label fix: all compensated
brick length was incorrectly classified under another role.

The subsequent same-height brick and per-path qualification changes pass
682,646 assertions in 20 seam/brick cases, plus the 12-assertion metadata-split
regression. This includes conservation of the original wall geometry across
changing material boundaries, counting reserved entry/finish pieces as wall
material while excluding the added connectors. The saved dual-colour Benchy
uses inner-only bricking and hybrid top-surface processing. At that checkpoint,
different-height pairs were still excluded; the subsequent height-aware
connector change replaces that restriction for qualified constant-height
courses. The 1:2 slope bound is a local qualification, not a proof of clearance
for every nozzle/toolhead geometry or physical dimensional accuracy.

The arm64-ninja app was rebuilt. These are geometry/export checks, not physical
print validation of the resulting seam's dimensional accuracy.

## Height-aware brick corner checkpoint

The new corner regression fails in all 16 inner-only combinations against the
previous implementation. With XYZ connectors, all 20 seam/brick cases pass
(1,178,092 assertions), including both wall orders, Classic/Arachne, straight
and curved walls, nonzero seam gaps, and hybrid top-surface mode.

The saved dual-colour Benchy was re-sliced to
`build/brick-seam-benchy.oGM1Bw/plate_1.gcode`. An independent modal G-code walk
finds 206 extruding seam entries, 201 with a Z change, 206 returns and 206
finishes coincident with the internal starting point within 0.003 mm. No
retraction, travel, or Z-only interruption occurs between those entries and
finishes. This audit confirms the paired paths, not whole-model printability;
the export still contains non-extruding seam fallbacks where qualification
does not produce a paired connection.
