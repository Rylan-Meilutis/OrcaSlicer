# September 18 print-quality and material-sync work

Reproduction: `/Users/rylan/Downloads/ksr_fdmtest_v4_0.2mm_PLA_Prusa CORE One HF_2h59m.gcode`.
The reported faults are physical sagging and dimensional error, not merely preview artifacts.
Non-planar work remains paused.

## Changes in this pass

- Extend the existing dissimilar-support slowdown to support courses deposited
  directly onto a different support material. Use actual selected filament slots,
  lower-course bead footprints and matching contact Z; exclude same-height side
  contacts. Cache footprints rather than reconstructing them for each path.
  Existing model-above-interface speed recovery remains supported. The legacy
  option key and disabled default remain compatible; the label is now
  "Slow dissimilar support contacts".
- Cap the intermediate-wall thickening multiplier at 100% in Brick mode, and
  do not change wall order for a now-capped multiplier. Scheduled brick volume
  is preserved; explicit flow reductions still work. Ordinary planar wall
  thickening remains unchanged.
- Prompt for a compatible material profile during OctoPrint sync when the
  material has no usable saved mapping. Stage choices until all prompts are
  accepted, then save global brand/material rules. Cancellation does not change
  the project or mappings. Existing roll-specific rules retain precedence.
- Retain provider entries without roll IDs for brand/material mapping. Saved
  material rules remain visible in the editor without a current inventory roll.
- Treat material capacity as a maximum, not a required count. OctoPrint physical
  multi-tool projects may retain fewer materials; normal count refresh must not
  repopulate all unused tools. Other vendors' physical-tool refresh is unchanged.

## Evidence and remaining work

Build succeeds for the app and both test executables. The targeted support,
brick and bridge-flow suite passes 17 cases / 8,593 assertions; the project
capacity regression passes 10 assertions for virtual and physical tool setups.
The brick-volume regression isolates seams/retractions, compares 100% and 110%
on identical brick geometry to 0.0001 mm filament, and separately allows 0.01%
between planar and shifted-course geometry. The old bit-identical cross-geometry
comparison failed by about 0.00985 mm filament even at 100% flow.
GUI compilation passes, but the new sync dialog has not been interactively
validated against a live OctoPrint provider. `/Applications` was not updated.

The supplied file sets inner-wall flow to 110%, arc flow to 112.5%, thick external
bridges on, and dissimilar-interface slowdown off. A modal G-code walk finds
1,189 extruding inner-to-outer seam transitions, 1,190 extruding returns, and 61
non-extruding transition fallbacks. Therefore inner seams are not universally
absent, but fallback qualification still needs analysis. Do not remove crossing
or clearance checks merely to force an L.

The existing bridge flow computes a circular cross-section. Regular external
bridge density already controls straight-line spacing. User additionally wants
separate bridge-overlap and sloped-wall-overlap controls to reduce sagging around
30-degree overhangs, accepting slower printing. Sloped-wall changes must retain
the outer model contour and update adjacent inner-wall/infill occupancy, rather
than just overextruding or moving the outside surface.

Physical short-arc sagging is NOT fixed by this pass. The original arc-overhang
reference explicitly warns that even 5 mm/s may be too fast for tiny arcs to cool:
https://github.com/stmcculloch/arc-overhang#6-suggested-print-settings
Cooling/time validation is separate from geometric support and crossing checks.

Print-time tool mapping remains UNIMPLEMENTED. OctoPrint currently uses the
generic send dialog; GCodeWriter::toolchange emits logical filament IDs. A real
mapping implementation must update generation and custom-template context,
temperatures, used-tool/calibration lists and spool metadata consistently. Do not
add a cosmetic mapping dialog or regex-rewrite only T commands. It must support
material/color suggestions, user overrides, cancellation, and re-export/reslice
when mapped nozzle geometry changes. No release or installation is authorized
by this pass; the build must not be described as completing all requested work.

## Follow-up: independent overlap and short-arc controls

- `bridge_line_overlap` (Quality > Bridging) replaces, rather than compounds,
  external bridge density when nonzero. Twenty percent corresponds to 125%
  density. Strand cross-section is unchanged, and arc paths retain their own
  spacing control.
- `overhang_wall_overlap` (Quality > Overhangs) reduces the inward step of extra
  overhang perimeter courses. It does not displace the outside contour or
  enlarge the bead. It requires Extra perimeters on overhangs and overhang
  detection; ordinary wall spacing is deliberately unchanged. The control is
  disabled when the arc selector owns these regions.
- `arc_overhang_min_path_time` (Quality > Arc overhangs) caps the speed of a
  complete short arc path, not each segment. The 0.1 mm/s lower speed bound is
  documented. No dwell, added extrusion, geometry change, or faster move is
  introduced. This is an optional tuning control, not verified physical sag
  elimination.
- All three default to zero for compatibility, participate in process-profile
  persistence, and invalidate the appropriate generation/export stage.
- Print/upload now uses the selected physical printer's host credentials, the
  same source as spool synchronization, instead of only the machine preset.
  This fixes a send-path discrepancy, **not** print-time tool remapping.

The more precise seam audit matches the exact `outer wall seam transition`
comment (excluding `move to first ...` travel comments). It finds 29 actual
non-extruding transitions in the Downloads KSR file, not 61: 23 lack an
immediately preceding primer and six follow a primer. The latter occur at
source lines 187896, 203659, 211731, 232896, 242509, and 353510. Existing
1,189 extruding connections are not evidence that these fallbacks are fixed.
No connector clearance checks have been relaxed.

Follow-up verification: the three new slicing/export cases pass 203 assertions,
including both Classic and Arachne, both arc roles, short and long paths, and
disabled controls. JSON persistence passes 10 assertions; the existing material
capacity case adds another 10. The initial synthetic arc test was corrected to
use the exporter's required fill collection and the sliced object's local
coordinates. The arc check requires actual preview extrusion moves and checks
their feedrate and total path length, so missing paths cannot pass vacuously.
All three changed GUI translation units compile. The full application relink
was paused for test debugging and has not been completed or installed.

The broader arc, brick, bridge-flow and dissimilar-support regression run passes
44,455 assertions in 52 cases. The Python G-code audit's 32 tests also pass.
The GUI's wall-order override now matches the brick flow cap: a saved multiplier
above 100% no longer forces inner/outer/inner ordering while Brick is selected.
Live OctoPrint interaction and physical cooling/adhesion remain unverified;
print-time tool remapping and remaining inner-seam fallbacks are still pending.
