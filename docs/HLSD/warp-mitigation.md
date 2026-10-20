# Warp and curl risk

Warp mitigation is conservative material screening, not a simulation of bed
adhesion or a probability of print failure. The slicer cannot measure plate
cleanliness, adhesive condition, drafts, chamber uniformity, or actual deformation.
Physical validation remains necessary.

The process setting `automatic_warp_mitigation` defaults to false for project and
preset compatibility. The Others page also offers a one-click application of the
same settings. Screening considers model material assignments, material painting,
feature assignments and modifiers rather than unused filament slots or support
materials alone. Susceptible matrix families include ABS, ASA, PC, polyamides,
PP and high-temperature engineering plastics. Filled variants remain candidates;
this does not imply they shrink as much as unfilled polymers.

The adhesion allowance is a bounded heuristic: 6 mm below a 40 mm XY bounding-box
diagonal, 8 mm below 100 mm, and 10 mm above that. These are not experimentally
calibrated physical thresholds. A continuous outer brim replaces painted or
mouse-ear-only coverage, preserving existing inner brims and any larger width.
The separation gap is capped at 0.1 mm. Painting data remains in the project.
Rafts and belt-printer adhesion are not replaced. Explicit object overrides
remain authoritative. The process-wide adhesion settings affect all objects
without such overrides; the largest susceptible object's span sets the allowance.
A wider brim may require additional clearance when arranging the plate.
Curled-perimeter slowdown and overhang-speed handling are
enabled, without changing nozzle/bed temperature, chamber targets or bridge
cooling. Automatic changes operate on the effective slicing configuration, not
the saved input preset.

Material-risk warnings remain visible in both Prepare and Preview. They distinguish
limited continuous adhesion from residual risk with adhesion measures in place;
enabling a brim, raft or automatic mitigation does not certify that a print is safe.
Revalidation replaces the warning and removes it when susceptible materials are no
longer used.

ValveBoard_case's ASA/ASA/TPU configuration with painted brim, zero general width
and a 0.2 mm gap motivates a regression fixture. Tests verify screening and the
continuous-brim policy, not a promise that the physical case cannot lift.

Material guidance:

- [Prusa ASA guide](https://help.prusa3d.com/article/asa_1809):
  temperature gradients and enclosure considerations.
- [Prusa enclosure guide](https://help.prusa3d.com/article/enclosure-guidepost_366332):
  materials susceptible to cooling shrinkage.
- [Prusa polycarbonate guide](https://help.prusa3d.com/article/polycarbonate-pc_165812):
  warping and increased bed-contact measures.
