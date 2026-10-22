# Rooting at material interfaces

Rooting is an opt-in process/object setting in the Multimaterial tab. It assigns
buried volumes of the lower material to the material beginning above it. The broad
attachment follows the inset upper model footprint, including concavities and
holes; it is not a cylindrical plug. Thick primary limbs divide recursively into
smaller lateral branches. The wider buried footprint provides
an undercut against vertical pull-out; it is not a guarantee of joint strength.
Flexible materials can still deform out of an anchor, and material combinations
need physical qualification.

Contacts require different nonempty filament material types (case-insensitive,
ignoring surrounding whitespace), not merely different filament slots or colors.
Equal material types and unknown types do not generate roots. This filtering is
per pair, including objects containing three or more materials. Conventional
prime towers are generated outside the model slicing pipeline; model-backed belt
purge towers are explicitly excluded. Objects without eligible material pairs do
not acquire rooting-only interface shells or root-failure warnings. Changing a
material type or an object's purge-tower role invalidates rooting geometry.

Failed Z contacts are aggregated into one warning per object listing affected
parts. A single affected part retains its direct selection link; a multi-part
summary links to the owning assembly. Successful contacts and excluded pairs
are not listed as failures.

XY rooting is an independent, default-off switch for sideways contacts. It grows
rounded, branching fans into the adjoining material, like roots entering a cliff.
The smaller native region on a course supplies the root material; the larger
region is the host (region order breaks equal-area ties). This is a geometric
heuristic, not a prediction of load direction. The shared depth setting limits
sideways penetration, while spacing separates both attachment sites and Z bands.
Two primary limbs produce recursive side branches, with three generations and
six shortening attempts per limb. Exposed side, top and bottom skins, third
materials and printable host webs constrain the complete swept geometry.
The swept trunk and limbs are trimmed at the protected envelope. The remaining
geometry must preserve host webs and connect to its owning material; floating
fragments reject the attachment. Connectivity may pass through adjacent
courses, allowing attachment to alternating beam teeth. Existing beam additions
are preserved; sideways growth starts at the outside of their attachment band.
At beam joints, the neck fits the beam tooth and widens behind the neighboring
teeth before branching. Growth follows the native material boundary normal,
not the diagonal tooth edges. Primary limbs must extend beyond short stubs;
otherwise that attachment is rejected. The swept tips stay within the selected
penetration depth.
Vertical and XY rooting can be enabled together.

Layer clearance and host-core preparation run in parallel, with each worker
owning a separate layer's output. XY clearance envelopes are cached once per
material pair and layer, shared read-only while intersecting neighboring skins,
then released before growing trees. Tree acceptance and material transfer stay
ordered because each accepted root changes the remaining host-web constraints.
Slicing progress identifies rooting separately from mesh intersection and beam
interlocking, reporting completed Z contact layers or XY material pairs. These
are work counters, not estimates of remaining time; contact complexity varies.

The controls are root depth, base diameter, tip diameter, spacing and protective skin. Old
profiles and projects default to rooting disabled; no existing keys change meaning.
The tip retains the `rooting_width` key and existing saved values. The independent
`rooting_base_width` defaults to 3 mm, while the tip defaults to 1.2 mm. A missing
base key uses that default; a base smaller than the tip uses the tip diameter.
Both vertical and XY branches interpolate from the base diameter toward the tip
across their branching generations, rather than repeatedly clamping a tapered
parent to one common width. The options serialize as normal process/object settings. Any change invalidates
slicing, since rooting changes material ownership before perimeter generation.

Only contacts between different filament regions within one multipart or painted
object qualify. Independently placed objects and generated supports do not.
Z rooting excludes root-owning components connected to the first native object layer:
bed-supported material, including inset first-layer lettering and its surrounding
body, does not need downward anchors. This does not exclude a bed-supported host
under a separate elevated region, or change XY rooting. Objects without eligible
Z contacts do not report a failed Z-root placement warning. Each disconnected native
Z contact records its own placement result, so successful roots elsewhere in an
assembly cannot hide a failure. Micron-wide clipping slivers are excluded using
the same 0.001 mm rounding tolerance as XY contacts; changes within 0.01 mm of
the previous course are treated as continuing material rather than new contacts.
Failed contacts are intersected with native part
slices to identify the affected volumes. Warnings carry stable volume IDs, and their
links re-resolve and select those parts in the assembly tree.
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
cavities and bottom skin. Vertical roots trim the obstructed side of a branch
to this protected interior rather than discarding its unobstructed side. Clipped
components and subsequent offshoots must remain connected upward to the tree;
floating fragments are discarded. Growth is shortened or stopped if the host
web or attachment cannot be preserved. Narrow
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
Narrow annular contacts allow a locally smaller attachment down to two nozzle
diameters, retaining the model-shaped collar and the full protective skin.
This attachment allowance does not reduce the requested terminal branch diameter;
only the side meeting a protected boundary is trimmed. Seed
search includes their inset contours instead of relying on a bounding-box center
that may lie inside a screw hole.
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
Host-web checks run independently across layers using TBB; a shared rejection
flag only short-circuits work, never changes acceptance. Candidates and limbs
are still accepted serially in their deterministic order. Geometry workers poll
cancellation but do not emit progress callbacks.
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
