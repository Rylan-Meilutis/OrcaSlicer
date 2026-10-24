# Non-planar selection and transitions

Experimental surface-following slicing remains quarantined: the normal UI does
not offer its modes, and print validation rejects imported enabled settings.
The internal planner and tests retain the implementation independently of that
availability gate. This is fixed-nozzle, three-axis surface following, not
rotary-axis or tilted-chunk slicing.

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

The writer emits the validated orthogonal route directly, without adding
unvalidated wiping, spiral lifts or perimeter-avoidance detours. Coincident
wall fragments retain the existing continuation exception. Non-planar
extrusions and planned transitions must remain within the printer's vertical
build limits; an unsafe transition stops export instead of silently flattening
the region. This does not certify custom G-code motion or all other travel
sites, such as tool changes.
