# Rooting at vertical material interfaces

Rooting is an opt-in process/object setting in the Multimaterial tab. It assigns
buried volumes of the lower material to the material beginning above it. Five
overlapping circular branches (a central stem and four spreading branches)
converge into a narrower stem at the contact. The wider buried footprint provides
an undercut against vertical pull-out; it is not a guarantee of joint strength.
Flexible materials can still deform out of an anchor, and material combinations
need physical qualification.

The controls are root depth, branch diameter, spacing and protective skin. Old
profiles and projects default to rooting disabled; no existing keys change meaning.
The options serialize as normal process/object settings. Any change invalidates
slicing, since rooting changes material ownership before perimeter generation.

Only contacts between different filament regions within one multipart or painted
object qualify. Independently placed objects and generated supports do not.
Filaments marked as support material or soluble are excluded on either side of
the joint, even if assigned to a model part. Ordinary model filaments that are
also used by generated supports can still form roots at model-to-model contacts;
support and interface toolpaths themselves are never inputs to rooting.
Both regions must use one filament consistently for walls and fill, with nonzero
wall, top-shell and bottom-shell counts. Mixed feature-filament overrides are
rejected rather than producing roots with the wrong material in their core.
Original slices are captured before beam interlocking. Roots are applied after
beams, rejecting roots whose protective envelope intersects modified material.
Artificial beam interfaces therefore cannot recursively seed roots. Roots never
seed other roots during the same slice either.

All edits are below the original interface. The upper part, including its top
layers, is unchanged. A root's entire envelope must fit under the inset contact
patch and inside the lower material throughout its depth and the protective floor
below it. This keeps exposed ledges, side walls, cavities and bottom skin intact.
Roots are rejected whole rather than trimmed into disconnected anchors. Narrow
contacts and shallow bases can therefore receive no roots.

Depth is quantized inward to complete existing layers. Branches shift laterally
by no more than half their vertical rise, and branch diameter must be at least
twice the larger participating nozzle diameter. Spacing is increased if needed
to leave host material between roots. Material polygons are transferred, not
overlaid: the combined model volume remains unchanged and material regions stay
disjoint. Standard slicing, tool ordering and G-code export handle the new regions.

Rooting also enables internal interface-shell processing, without rewriting the
saved Interface shells option. This supplies top/bottom shell layers around the
buried material boundary using the configured shell counts. Use nonzero top and
bottom shell counts and appropriate walls/infill: the feature does not override
those settings or certify a structural fastener joint. It can add many tool changes.

Regression coverage in `test_multifilament.cpp` checks material-volume conservation,
disjointness, preserved exposed regions, wide buried anchors, vertical continuity,
rejection of shallow/narrow/cavity cases, and restoration when disabled. Geometric
checks do not replace pull-out testing or inspection of the final sliced toolpaths.
