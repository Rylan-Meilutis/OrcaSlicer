# Per-print filament bindings

In a process profile, select **Per print** for a role's filament instead of
choosing a fixed numbered slot. This is available for support, support interface,
outer walls, inner walls, sparse infill, internal solid infill, and top/bottom
surfaces. The wipe-tower selector retains its existing behavior.

An **Area bindings** section below the sidebar filament list then shows the
enabled roles. Choose a filament for each role in the current project. Automatic
uses the filament assigned to the most printable part instances (ties choose the
lowest slot). This is a part-count rule, not a painted-area or material-volume
estimate.

The process profile stores the Per print choice (`-1`), not the project's slot.
The project stores `project_filament_bindings` and `project_filament_roles`.
During slicing the role choices resolve to positive, one-based filament slots.
The project metadata restores the Per print choice when loading a saved 3MF.
Removing a slot renumbers bindings; unavailable slots fall back to Automatic.
Existing Default and fixed-slot choices are unchanged.

Regression coverage includes profile JSON persistence, project 3MF round-trip,
slot removal, automatic choice, and resolved wall/infill roles in a sliced model.
Cloud synchronization uses the normal profile JSON path; live cloud transport
has not been tested as part of these regressions.
