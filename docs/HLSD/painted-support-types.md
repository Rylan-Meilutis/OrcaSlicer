# Painted support types

Support painting can assign Normal or Tree to different contact regions of one object.
Mixed assignments use the shared column/branch planner: normal regions become fixed polygonal
columns, while tree regions produce branching circular nodes. They share collision
volumes, support layer heights, interfaces, area unioning, and toolpath generation.
Two independent support generators are not combined, because their layer grids and
collision decisions would not agree.

Painting never changes the object's support type or style. Unpainted contacts use the
configured default: Normal produces columns, Slim/Strong produces branches, and
Hybrid retains its area-based heuristic. Automatic/manual detection remains unchanged.
The painter's enable action only enables supports, without replacing their configuration.
Tree paint under an Organic default is an explicit local enforcer for that same Organic
planner. Mixing Normal paint with Organic remains unsupported and produces a specific
validation error rather than silently replacing Organic. Disabling support preserves
the annotation without generating it.

## Contact partition

Paint is projected onto the object's sliced contact surfaces. Normal and Tree regions
are also support enforcers, in both automatic and manual support modes, so the ordinary
small-overhang filter cannot remove explicitly painted demand. The final contact area
is partitioned geometrically, not classified by its centroid. Blockers take precedence.
If projections from different facets overlap, Normal takes precedence over Tree.
The existing enforcement margin carries the selected type so its fringe does not
create a different support type. Small Normal contacts bypass
the automatic large-area threshold; large Tree contacts bypass automatic conversion
to normal columns.

Typed paint applies to positive-area downward-facing support contacts. It does not
create a new contact footprint on a purely vertical face. Generic Enforce remains
available for the existing vertical-point support behavior. This feature does not
combine Organic branches with normal columns. Mixed normal/tree jobs share a layer grid
and collision planner; a Normal default retains normal columns and its material/pattern
settings, but does not promise identical classic Grid/Snug support geometry.

## Project compatibility

The support facet annotation retains its existing subdivision encoding and adds two
leaf states. Both 3MF readers and writers preserve the typed data in a versioned
extension. The legacy support-paint attribute contains the same subdivisions with
both typed states mapped to generic Enforce. Older readers therefore retain support
demand rather than misinterpreting a typed leaf. Unknown or malformed extended data
falls back to the legacy annotation. No typed annotation means the existing support
generation path is unchanged.
