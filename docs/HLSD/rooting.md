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
beams. Where a contiguous beam joint extends below the native contact, roots
start below its lowest course. A model-shaped collar joins the beam teeth from
below without replacing the alternating beam band. Root depth is the additional
depth below that band, limited by the base's remaining protective floor.
Beam additions are snapshotted before growing roots so one root cannot become
another root's beam attachment. Third-material regions are never reclaimed.
Artificial beam interfaces therefore cannot recursively seed roots. Roots never
seed other roots during the same slice either.

Native upper geometry, including its top layers, is unchanged when beams are
disabled. In combined mode the beam-only geometry above the root attachment
is preserved instead. Buried branches may extend
beyond the contact patch only where a full protective top skin remains.
Per-course clearance checks include host side walls, native material above and
below the branch, and third-material regions. This preserves ledges, holes,
cavities and bottom skin. A branch that does not fit is shortened or stopped;
it is never clipped into disconnected pieces. Narrow
contacts and bases too shallow for connected branches can therefore receive no roots.

Depth shortens to the available host material above a complete protective floor,
and is quantized inward to complete existing layers. The configured diameter is
the minimum tip width in XY and must be at least twice the larger participating nozzle
diameter before a rounded terminal cap. Shallow roots use vertically flattened
sweeps without reducing that XY width, and lateral reach has a width-based floor
so reducing depth does not collapse the root into a narrow vertical peg.
Four primary limbs retain their tips
while smaller side branches emerge at staggered positions along each limb.
Those branches recursively grow their own side branches, with at most three
generations (28 limbs) per seed. The trunk's fit is not tied to the sum of all
terminal cross-sectional areas, which would suppress branching at small contacts.
Radii taper to the printable minimum. Small/shallow contacts use fewer generations.
Tapered sphere sweeps produce connected, rounded junctions and long horizontal
cross-sections, rather than circles that merely translate with layer height.
Horizontal reach is independent of a fixed vertical branch slope and targets
roughly twice the usable root depth, shortened as needed by the actual base.
All growth is deterministic and bounded: three generations, six length attempts
per branch, a minimum printable tip radius, and cancellation between courses.
Spacing separates attachment centers. Before accepting each limb, the generator
checks the accumulated roots against each native host component. It shortens or
rejects a limb that encloses a disconnected host pocket or severs a printable
host web. The web check erodes the remaining base by half the larger of the
protective-skin width and two host-nozzle diameters. The model-shaped attachment
plate is exempt; the constraint applies to the buried branching network below it.
Same-material branches may merge only while preserving these host connections.
Material polygons are transferred, not
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
preserved beam courses with shallow roots beneath them,
rejection of insufficient-depth/narrow/cavity cases, and restoration when disabled. Geometric
checks do not replace pull-out testing or inspection of the final sliced toolpaths.

Preview represents roots as ordinary walls and fill in the owning filament, not
a separate line type. Use filament coloring and inspect layers below the original
material contact. Export tests check both G-code and parsed preview moves, with
beam interlocking enabled and disabled.
