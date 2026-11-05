# Non-planar selection and transitions

Experimental surface-following slicing remains quarantined: the normal UI does
not offer its modes, and print validation rejects imported enabled settings.
The internal planner and tests retain the implementation independently of that
availability gate. This is fixed-nozzle, three-axis surface following, not
rotary-axis or tilted-chunk slicing.

## Course ownership

Each candidate course owns its fill graph. Sparse courses project native sliced
source paths; required solid and finishing courses use the native Fill factory
on the source region's wall-bounded fill footprint, clipped to the patch.
Configured patterns, source path direction and extrusion roles are retained;
buried fill is not obtained by thinning or cloning the finishing raster.
Projection must retain the complete native bead footprint, including gap fill;
it cannot discard narrow features to make a coverage check pass. The required
footprint is recorded from the planar source subdivision before projecting Z.
Subdivision retains original vertices exactly and rounds only inserted samples
to the coordinate grid; explicit courses do not use lossy path simplification.
Only wall pieces retain stable
cross-course identities for wall morphing. Crossing, support and toolhead
clearance checks inspect every course independently, and a failed candidate
keeps the original planar geometry. Routing may reconnect exact existing wall
endpoints but does not invent new serpentine fill or bridge connectors.

For a complete heightfield roof, each supporting slice must match the selected
roof clipped at that slice's native plane. Only this verified case reuses the
last complete planar perimeter graph for every transition course. The interior
stair-cut edges of ordinary horizontal slices are not physical walls of a
surface-following course. Partial regions and continuing features do not enter
this replacement path. Each course has its own nominal owner and physical
height profile; local wall separation controls flow and must satisfy configured
minimum and maximum layer heights, predecessor support and toolhead clearance.
Wall mutations are transactional and restore the complete planar graph on
rejection. For this verified full-footprint case, the native fill graph replaces
the ordinary fill throughout each generated owner above the anchor. Clipping
only at the former horizontal nozzle plane would leave planar material above
the physical course. Partial patches retain local bead-at-height replacement;
the complete planar anchor is never removed by this owner rule.
Closed native wall courses start at an existing highest-Z vertex, preserving
their directed XYZ segments and segment-specific flow. This reduces avoidable
departure contact with the same loop; it does not exempt the route from
clearance checks. The non-planar writer respects the configured infill-first
ordering for ordinary surface courses instead of always putting retained walls
before their infill.

A sparse source is not relabeled as a required solid shell. Anchored
sparse-to-solid bridge planning remains required for cases where
the deposited source courses cannot support that transition. Rejecting these
cases is a safe fallback, not evidence that the requested roof was generated.

The course-ownership design follows the same separation of derived geometry
and print ordering used by Unlayered Slicer's declarative cut segments
(`CutPlaneSegmentation.hpp`, revision
`79a9bbd447613a52668e7891b9adafb9a937f3f8`). No implementation was copied.
That upstream mechanism slices tilted planar chunks, whereas this planner
follows curved surfaces. Its episode-based sampled hotend warning is not used
as a clearance certificate for same-axis surface-following paths.

## Selection

The region settings `nonplanar_selection_min_z` and
`nonplanar_selection_max_z` bound complete connected surface patches relative
to the transformed object's bottom. Zero maximum means unbounded. Their
defaults preserve automatic selection. They use the existing object/part
configuration mechanism, not a separate face-painting format.

Selection occurs after connected slope and narrow crown discovery. A patch
crossing a boundary is rejected whole, preserving the native planar paths.
The band selects the exposed surface, not the supporting transition courses;
those may extend below the selection. Angle, support, shell and toolhead
clearance validation still apply. Selection cannot enable a disabled method.

## Transitions

Non-planar travel uses machine-space departure, horizontal traverse and vertical
approach legs. Every leg is tested against the configured toolhead geometry or
fallback envelope and recorded deposited material. The planner first tries the
higher endpoint Z, then a height above the recorded material. It rejects a
blocked departure or approach even if the horizontal traverse is clear.

Non-planar jobs retain previously recorded occupancy across nominal layers:
a curved course can extend above several later layer planes. Ordinary
contouring retains its existing layer-local tracking behavior.

Recorded variable-Z beads retain their vertex heights. For detailed horizontal
toolhead slices, each travel/bead segment pair is tested in its two-parameter
domain: the height inequality clips a unit square, whose affine XY image is
checked against the toolhead footprint and bead radius. This avoids replacing
the downhill end by the segment's maximum height. It is a nominal swept-disc
bead envelope, not a softened-material contact model. Unknown or inconsistent
height profiles and the fallback cone retain conservative constant-top bounds;
there is no exemption for the last deposited bead.
Native clearance uses the writer's serialized XYZ precision, after applying
the exact plate offset, for both recorded deposition and travel endpoints.
Unrounded planner heights cannot be mixed with rounded emitted motion. This
qualification is limited to the identity Cartesian mapping; belt and remapped
kinematics require their own frame-aware clearance model and are rejected.

The writer emits the validated orthogonal route directly, without adding
unvalidated wiping, spiral lifts or perimeter-avoidance detours. Coincident
wall fragments retain the existing continuation exception. Non-planar
extrusions and planned transitions must remain within the printer's vertical
build limits; an unsafe transition stops export instead of silently flattening
the region. This does not certify custom G-code motion or all other travel
sites, such as tool changes.
