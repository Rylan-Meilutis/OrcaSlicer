# Rooting at vertical material interfaces

Rooting is an opt-in process/object setting in the Multimaterial tab. It assigns
buried volumes of the lower material to the material beginning above it. The broad
attachment follows the inset upper model footprint, including concavities and
holes; it is not a cylindrical plug. Thick primary limbs divide recursively into
smaller lateral branches. The wider buried footprint provides
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
beams, reserving accepted root envelopes for native host material and roots.
Beams remain outside these envelopes; they cannot replace the root's protective
floor or interrupt its stem in the receiving material. Third-material regions
are never reclaimed.
Artificial beam interfaces therefore cannot recursively seed roots. Roots never
seed other roots during the same slice either.

Native upper geometry, including its top layers, is unchanged. Inside the
attachment footprint, beam-altered courses are restored to the native upper
material so beams cannot cut the root connection. Buried branches may extend
beyond the contact patch only where a full protective top skin remains.
Per-course clearance checks include host side walls, native material above and
below the branch, and third-material regions. This preserves ledges, holes,
cavities and bottom skin. A branch that does not fit is shortened or stopped;
it is never clipped into disconnected pieces. Narrow
contacts and bases too shallow for connected branches can therefore receive no roots.

Depth shortens to the available host material above a complete protective floor,
and is quantized inward to complete existing layers. The configured diameter is
the minimum tip diameter and must be at least twice the larger participating nozzle
diameter before a rounded terminal cap. Four primary limbs feed binary forks,
with at most three generations (16 terminal branches) per seed. Parent radii
scale with the square root of the number of daughter limbs rather than narrowing
all limbs into a tip-width stem. Small/shallow contacts use fewer generations.
Tapered sphere sweeps produce connected, rounded junctions and long horizontal
cross-sections, rather than circles that merely translate with layer height.
Horizontal reach is independent of a fixed vertical branch slope and targets
roughly twice the usable root depth, shortened as needed by the actual base.
All growth is deterministic and bounded: three generations, six length attempts
per branch, a minimum printable tip radius, and cancellation between courses.
Spacing separates attachment centers; same-material branches may merge into a
connected network. Material polygons are transferred, not
overlaid: the combined model volume remains unchanged and material regions stay
disjoint. Standard slicing, tool ordering and G-code export handle the new regions.

Rooting also enables internal interface-shell processing, without rewriting the
saved Interface shells option. This supplies top/bottom shell layers around the
buried material boundary using the configured shell counts. Use nonzero top and
bottom shell counts and appropriate walls/infill: the feature does not override
those settings or certify a structural fastener joint. It can add many tool changes.

Regression coverage in `test_multifilament.cpp` checks material-volume conservation,
disjointness, preserved exposed regions, model-shaped attachments, holes,
lateral branch spread and connections between adjacent courses,
rejection of shallow/narrow/cavity cases, and restoration when disabled. Geometric
checks do not replace pull-out testing or inspection of the final sliced toolpaths.

Preview represents roots as ordinary walls and fill in the owning filament, not
a separate line type. Use filament coloring and inspect layers below the original
material contact. Export tests check both G-code and parsed preview moves, with
beam interlocking enabled and disabled.
