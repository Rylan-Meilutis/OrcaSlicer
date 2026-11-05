# Spatial thermoelastic equilibrium

Deposition rasterization accepts an explicit nominal cross-section. The default
rectangular envelope preserves existing callers; a rounded rectangle uses a
central strip and semicircular sides of diameter equal to layer height, requiring
width at least height. Its interior section area is `w*h-(1-pi/4)*h*h`.
Neither option establishes the actual extruded shape from imported display
widths. In particular, bridge and scarf paths cannot silently inherit that
interpretation. Sloped paths retain vertical thickness, not a normal-to-path
section. Rounded end caps and overlapping paths form a geometric union.

G-code scenario results retain three distinct volume diagnostics: raw commanded
E volume where processor data is available (without firmware flow overrides),
nominal section area integrated along XY paths, and occupied raster volume.
Missing commanded areas are counted, not replaced with inferred amounts. These
quantities differ with end caps, overlaps, display-width inference and grid
resolution; neither density nor empty space is altered to force agreement.

`ThermoelasticLayers` provides a compact geometry-only snapshot for horizontal
rectangular-envelope beads. It unions nominal capsule footprints per material,
splits vertical intervals at every bead boundary, and retains polygon holes and
disconnected islands. Coordinates use the native 1 nm grid; round-cap accuracy
is an explicit arc tolerance. Long paths are processed in bounded polygon
batches without a safety offset. Segment, group, boundary-point and slab-work
budgets abort instead of returning a truncated body. Cancellation is checked
between geometry operations.

Different-material overlap is an error carrying its slab interval, area and
material slot. Such overlap cannot silently become extra mass, an arbitrary
material assignment or a welded interface. Sloped and rounded-section beads
are unsupported by this representation. These polygons are not an activated
thermal history, a finite-element mesh, or evidence of mechanical convergence;
they are not connected to automatic mitigation selection.

`section_deposition` evaluates a single horizontal cut without building a
bead-resolved volume mesh. Sloped rectangular envelopes retain their straight
cut planes and original endpoint caps. Explicit rounded-rectangle
beads use the section width `w-h+2*sqrt(d*(h-d))`, where `d` is depth below the
bead top, retaining the nominal corner voids. Sloped rounded beads interpolate
top Z along the straight path and retain vertical bead height. Bounded ellipse
arcs describe their section sides, with half-caps only at the original endpoints.
This is an explicit nominal shape, not an inference of variable scarf flow. Section
integration is checked against analytic straight-path and end-cap volumes.
These cuts preserve material separation and geometry budgets but are not a
homogenized constitutive model or a mechanical force result.

A cut may receive explicit native material ownership for the same plane.
Ownership resolves only intersections between different deposited envelopes;
nonconflicting beads, unfilled interiors, and existing holes stay unchanged.
Every conflict point must have exactly one claimed owner whose deposited
envelope covers that point. Missing or overlapping claims reject the cut.
The Boolean union is checked for conservation and resolved overlap area is
reported separately. Ownership masks do not become extra solid volume, a
blended modulus, or evidence of a welded interface. Without masks the default
remains rejection of cross-material overlap.

`ThermoelasticOwnership` captures masks from the current `LayerRegion::slices`
after slicing, perimeters, and infill have completed. This retains native
rooting/interlocking ownership instead of reverting to original CAD or raw
slices. The caller must hold the matching slicing generation immutable and
provide its coordinate placement and logical-filament-to-mechanical-material
mapping explicitly. All region extrusion roles, including bottom surfaces,
must resolve to one mechanical material. Physical tool numbers are not used
as material labels.

This adapter conservatively rejects purge-into-object/infill overrides, mixed
extrusion, spiral or height-modified paths, thick bridges, and combined infill:
these need per-extrusion or cross-layer provenance beyond a planar region
label. It does not derive support/brim ownership or attach current masks to
an unrelated imported G-code file. The masks describe intended native
ownership, not measured redistribution of overlapping extrudate.

`DepositionProvenance` provides a bounded, generation-local registry of numeric
object, instance, layer, region, logical-filament, height-transformation and
material-assignment IDs. The assignment distinguishes the writer's ordinary
nominal-role pass from per-instance wiping overrides after entity filtering.
It is internal export evidence, not a measurement of nozzle material purity;
unknown or wiping assignments cannot authorize paths when purging is enabled.
Its decoder requires a sealed registry and rejects stale, malformed, nested or
unfinished scopes and mismatched logical filaments. An independently captured
expected control-event transcript is copied before decoding and checked exactly,
with a four-MiB storage cap. Completeness is available only after successful
finalization; no-transcript decoding never certifies it. Live exports can opt in
to recording perimeter and infill regions after per-copy/tool assignment. The
serial producer owns a lazy registry; the output thread captures only provisional
IDs and never reads that mutable registry. After the pipeline joins, the sealed
registry and independent expected/observed transcripts must agree, and every
attributed move must match its source's final logical filament. Any failure
discards all provisional IDs. Published identities are not ownership masks.
Support, brim, tower and separate nonplanar-surface emission remain unowned.
Spiral and custom timelapse exports do not enable capture. Height-changing
features are conservatively labeled transformed, not planar.

`ThermoelasticSourceOwnership` copies bounded native masks during the matching
live export, after source verification and while the sliced `Print` remains
immutable. The result retains the snapshot alongside the exact registry; copy
and move delivery preserve both. Cartesian instance shifts are retained, while
writer/processor plate and nozzle offsets cancel. Belt and IMEX-shifted exports
do not publish this snapshot. Unsupported sources remain explicitly unavailable
instead of invalidating unrelated native sources. Snapshot resource limits or
unsupported input leave ownership unavailable without changing exported G-code.

Queries require the same registry, the actual logical filament, and the full
native bead-height interval for every supplied conflict contributor. Material
properties are mapped separately at query time. The masks are used only to
partition conflicting deposited envelopes; they cannot certify bead geometry,
bond strength, or bed adhesion. A transformed perimeter additionally requires
matching authenticated command dimensions, its entire design envelope within
the native layer interval, and a centerline contained in the native mask, except
for the separately authenticated brick-course mapping below.
Region-wide transformed provenance can also contain ordinary planar infill
when other paths use scarfing, Z-contouring or combined infill. Such infill is
admitted separately only for nominal roles with unmixed material. If settings
permit purging or infill/Z transformations, the source must explicitly belong
to the ordinary emission pass and each move must match authenticated command
dimensions. Its live source, material, endpoints, width, full native-layer
height and centerline containment must match. Actual combined-height or
variable-Z infill does not pass these native-plane checks. With ordinary
settings, the existing live planar-move evidence remains supported; optional
dimensions or a fabricated caller path never replace authenticated evidence.
Constant-height linear brick courses carry a generator-emitted transform kind
alongside their authenticated command dimensions. Their source layer identifies
the emission, not permission to extend that layer's mask vertically. Eligible
sources capture immutable references from the current and at most two contiguous
neighboring intervals, matched by object, instance, region and filament. Query
height selects the actual containing reference after live command checks; the
whole bead envelope must fit the captured interval window, and its XY centerline
must pass containment against the selected mask. Missing or mismatched references,
arcs, tapered commands and mixed lower-surface definitions reject. Shared references
are deduplicated under the existing point/work budgets. These are native design
references, not evidence of previously deposited material or physical bonding.
Command capture preserves public XYZ/E and the scheduled nominal path dimensions;
it does not infer physical thickness from transition-flow multipliers.
Legacy ContourZ paths carry a dedicated generator marker, distinct from the
shared variable-Z transport flag. Eligible ordinary walls and top skin retain
the independently known native layer bottom as an anchor; actual emitted top
coordinates determine their varying nominal heights. Capture requires the
original layer height and nominal plane and excludes other mixed transforms.
An unmarked variable-Z path cannot acquire that anchor merely from its role or
flow ratio. Native-ownership restrictions remain separate from this command
geometry evidence.
Positive-extrusion pure-Z brick commands retain their legacy preview travel type
and timing. A separate marker is published only after command-chain and source
authentication. The toolpath consumer independently checks that evidence and
retains their filament amount and elapsed time as unresolved deposition. A
vertical centerline supplies no supported XY bead footprint: it cannot become a
fabricated solid or silently disappear from the mechanical inventory.

An independent, explicitly selected overlap model assigns shared envelope space
to the first deposited command. It authenticates command chronology and live
geometry, then adds only previously unoccupied area to each material. Exact
Boolean checks preserve the occupied union at every prefix and disjoint material
interiors; same-material commands merge without artificial internal boundaries.
Each individual bead cut retains the sampler's quantized envelope and clipping,
but this policy composes those cuts exactly rather than rounding cross-bead
intersections through an integer union. Before its first checkpoint, this sampler
freezes bounded copies of the commands, geometry settings, clipping polygon and
ownership descriptor. Composition and source authentication use the same frozen
inputs; copying is charged to the existing sampler work budget. After complete
source matching, private per-command cut-admission records avoid reconstructing
present cuts merely to test whether they are empty. These records retain the
singleton construction and normalized boundary limits; a stricter resolver limit
falls back to reconstruction. Every nonempty source still requires live
authentication, and missing contributors still require independent cuts. Once
prefix and live-source authentication succeeds, a genuinely single-material request can
reuse its immutable exact sampled union as the sole owner, with zero competing
material area and the same conservative occupied-plus-owner storage bound. A
fresh private request identity binds the returned proof. Mixed-material requests
can reuse frozen material unions when each material occupies exactly one
chronological run and each union fits the streamed primitive's incoming limit.
Material returns or oversized unions fall back to individual command replay.
Both paths retain the complete occupied-union comparison against the sampled
union. The sampler accumulates each
material's cuts in bounded binary runs to avoid repeatedly overlaying a growing
prefix. All runs share the existing live-vertex limit; deterministic consolidation
can free redundant storage before admitting another cut. Every overlay is charged
to the same work budget. This changes union evaluation order, not chronological
ownership, and does not promise identical resource admission to a sequential union.
Sampler diagnostics separate exact-cut construction, material-union composition,
and other charged work without changing the total allowance. Local cursors
aggregate completed samples and retain any in-flight reservation in the labelled
other category. These are geometric-complexity counters, not elapsed-time
measurements or estimates of physical error.
For a sole material owner,
one exact union supplies both the owner and occupied-union handles; disjointness
is vacuous, while any previously recorded cross-material ambiguity is retained.
The live sampler reconstructs command cuts incrementally. Its streamed primitive
bounds each incoming cut, retained partition geometry, command count and
cumulative work separately; it does not store the complete batch of cuts.
The retained-batch API keeps its separate cumulative input-vertex limit.
The streamed path accumulates consecutive nonempty cuts of one material in
bounded binary runs before allocating their shared region. It consolidates
pending runs before requesting an unknown incoming cut when the maximum cut
would exceed the live allowance; the actual incoming payload is still checked.
Foreign ownership is constant within that run,
so this preserves first-deposition ownership and competing-material overlap.
Each requested prefix flushes its own final partial run. Live bounds include
the pending run and any next-material lookahead cut at stable streaming
boundaries. While certifying a consumed run, old operands and replacement
geometry are bounded operation scratch. Work and storage admission
can differ from the unbatched reference even when both produce identical sets.
Temporary Boolean arrangements remain subject to the exact operation limits;
the retained-vertex bound is not a claim about all simultaneous scratch storage.
Its sampler-bound proof records the deposited-order policy separately from
native-boundary assumptions. It neither models molten-material displacement nor
establishes bonding, adhesion, or a bound on stress error. Final-state material
priority cannot substitute for chronology: clipping an early prefix to a later
material's final ownership could otherwise erase already deposited material.
The existing explicit stationary-extrusion allowance remains a separate modeled
omission, with unchanged absolute/local/fraction bounds and a reported inventory.
The exact partition proof covers retained envelopes, not that omitted extrusion;
it cannot admit a moving vertical deposit through the stationary allowance.
Generated conventional brim commands have a separate auxiliary purpose and no
fabricated object source. Their typed command envelope, filament, live child
chain and generation must authenticate before deposited-order ownership can
include them. A role label or a source-zero move alone is not sufficient.

Prepared resolvers retain a bounded cache of successful transformed-source XY
containment predicates. Certificates are scoped to an immutable snapshot and
keyed by source, selected reference identity, rounded directed XY endpoints and
explicit boundary allowance.
All live command, material, span, dimension and query-height checks still run
before every lookup. Fixed direct-mapped storage is charged to the resolver's
workspace; collisions evict entries, and full keys establish hits. Lookups and
insertions consume work and check cancellation. Ownership sections, query-Z
validity and physical contact are never cached by these certificates.
Native centerline queries interpret the original oriented chains using the
slicer's nonzero fill rule. Simple valid polygons use the exact polygon union;
self-crossing chains use exact-intersection nonzero winding reconstruction,
with a separate 1024-input-vertex, 4096-output-vertex and cumulative-work bound.
Disjoint bounding-box components are processed independently. An exact edge-pair
preflight bounds intersection vertices before triangulation, rejecting excessive
crossings before they can allocate an unbounded arrangement.
For a complete local deposited-envelope pair, interface reconstruction first
clips the original oriented native chains to a surrounding rectangle in exact
rational arithmetic. The rectangle is padded beyond the full interface band;
its artificial edges cannot establish ownership or boundary-band coverage
inside the deposited region. All contours and holes participate in one nonzero
winding reconstruction. This removes remote boundary complexity without
rounding crop intersections, changing local fill semantics, or modifying the
source snapshot. Cropped chains retain the independent reconstruction limits.
Rational crossings stay rational in the coverage predicate; masks are neither
rounded nor rewritten. Oppositely oriented duplicate chains cancel globally,
including their zero-area traces. General simple-ring predicates remain strict.
Other transformed sources and purge-ambiguous regions remain ineligible.
The snapshot is not attached to imported G-code and does not
by itself enable automatic mechanical mitigation.

`ThermoelasticConflictOwnership` authenticates a single reported conflict against
the matching processor result. It reconstructs the full nominal timeline,
checks every spatial contributor (including unowned helper paths), and requires
the complete conflict mask to be covered by actual different-material
intersections. Unrelated brim paths do not invalidate a contact. Returned native
masks are clipped to that conflict, so they cannot implicitly authorize another
contact with unchecked contributors. Unresolved deposition, partial masks,
unauthenticated transformed contributors and exhausted budgets remain unavailable. This bridge
does not enable the study or Auto mode on unsupported geometry.
The explicit interface-sensitivity factory applies the native nonzero winding
interpretation only to authenticated native masks. Deposited envelopes and
claims remain subject to simple-ring validity checks. Opposite native chains
cannot create coverage by being repaired independently; their signed winding
is combined before any interface proof or rounding.

The prepared model resolver uses the identical observation-time prefix and role
selection as the section study. A bounded XYZ AABB tree identifies all possible
contributors; actual sections, not bounding boxes, determine ownership. The
section study always queries original nominal beads, including when approximate
integration partitions are selected. Its conflicts require unique native masks.
A separate resolver mode can preserve authenticated original unique material in
a caller-supplied expanded conflict, but does not authorize original void or
replace the study's original-geometry evaluation. Integer Boolean partitioning must
preserve the original occupied union exactly by default. Local studies may
explicitly allow a one-to-four native-grid-unit boundary band for integer
intersection representation. Every changed polygon must be contained in that
band; individual islands and holes must correspond and survive erosion and
dilation. Materials remain strictly disjoint. The reported coordinate allowance
is zero when a bounded exact comparison proves the original and represented
occupied sets identical; this does not clear previously accumulated source or
interface uncertainty. Input island-survival checks still apply. Bounded exact
changes are certified over their full rational difference, including slivers
absent from rounded Boolean output. Their added and removed areas are converted
to outward-rounded bounds only after exact unit conversion. This certificate
does not replace the material and occupied-set morphology checks.
For a changed representation, the reported coordinate allowance
includes three additional grid units for offset arc sag and Euclidean vertex
rounding through the offset, band difference, and containment operations.
Morphological contact tests reserve two further units for their own rounding.
Before modifying geometry in this opt-in mode, bounded CGAL regularized
predicates can disprove rounded Boolean overlap, missing-coverage, or union
difference diagnostics using the original integer polygon operands. Simple
contours, hole orientation and validity are checked before set operations;
quadratic work and local vertex bounds apply. Intersection emptiness permits
component-pair AABB pruning; difference and equality retain whole-set semantics.
Exact predicates emit no rounded replacement geometry and cannot excuse a
genuine positive-area ownership ambiguity.

An explicit analysis-only interface scenario may evaluate two opposite
assignments of a native-ambiguous two-material contact. It is separate from the
strict default. The exact arrangement derives the complete deposited conflict
and its native claims without intermediate integer clipping, requires complete
native coverage, and confines the ambiguous region within 7 nm of both exact
native material-union boundaries. Exclusive deposited material and uniquely
owned contact material retain their identities. Each full material pair is
constructed before one shared rational-to-integer rounding. Ordered edge
correspondence bounds each vertex displacement by sqrt(0.5) nm; linear boundary
homotopy confines the entire representation change to original-edge capsules
of radius 1 nm. Ring validity and exact disjointness remain mandatory.
Each separately constructed input operand is limited to 1,024 vertices.
Validation is performed per operand rather than on a concatenation that is
never used by the actual Boolean calculation. Validation, arrangement
construction and subsequent Boolean operations all consume the same cumulative
work allowance; this does not relax the boundary or ownership proof.
Construction reserves quadratic input work. Overlays of already planar
arrangements reserve cross-operand edge pairs (including both endpoints of
collinear overlaps) and input traversal, rather than charging intra-operand
intersection pairs again. Output arrangement vertex limits remain independent.
Empty-set identities avoid overlays, and complete-envelope studies reuse their
exclusive deposited regions across both interface variants. If the exact native
claims are disjoint, both variants are identical and reuse one complete
rounding/topology certificate instead of repeating its overlays. Identity
differences and repeated exclusive-region unions are also omitted in that case.
Overlay work accounting charges all cross-edge bounding-box tests and reserves
intersection endpoints only for pairs whose outward-rounded boxes overlap;
the actual Boolean operations still use exact geometry. Native coverage is
proved by subtracting each native mask from the local conflict, avoiding a
whole-mask union. Rejected uncovered regions report their exact-area estimate,
piece count and bounds, plus a bounded whole-region check against both native
boundary bands. Rejection diagnostics also check the deposited-envelope bands
to distinguish a bead-edge representation issue from a native-mask disagreement.
Without an explicit native-boundary sensitivity policy, neither small area nor
successful band coverage authorizes assigning a material outside the native
masks. If rational-to-integer conversion moves no vertices, the
existing exact arrangement supplies the topology certificate. Before rounding,
both materials split shared straight edges at all exactly incident interface
vertices, so a T-junction present on only one side cannot produce inconsistent
rounded boundary chains. This adds no area to either exact region. Moved vertices
still require component and hole correspondence checks. Correspondence first
rejects pairs with disjoint outward-rounded bounding boxes. A rational point
strictly interior to both exact regions certifies positive-area overlap;
failure to find such a witness falls back to the exact Boolean intersection.
These shortcuts do not infer empty overlap from failed sampling. If integer
union output fails component/hole matching, an exact certificate compares the
original unions and matches every component and its own holes one-to-one.
Equal aggregate topology counts are diagnostic only, not a certificate. All
boundary-band, surviving-feature and contact checks remain mandatory.
Local topology also verifies rounded positive-overlap diagnostics with an exact
predicate on the input polygons. A proved empty intersection can continue;
positive overlap, however small, still rejects. The predicate has independent
vertex and work bounds charged to the topology budget. Single validated
polygons bypass redundant material union operations. Rejections identify the
height interval, diagnostic area and material rather than reporting only a
generic overlap.
In explicit boundary-sensitivity scenarios the section sampler also checks
exact output disjointness, including when its rounded conflict diagnostic is empty, so a
sub-grid intersection cannot bypass reconciliation and enter mechanics.
The predicate retains its 1,024-vertex default. Explicit dense section queries
can use up to 2,048 vertices only within the same four-million quadratic-work
ceiling, reserved before constructing exact arrangements. Holes wholly outside
the other operand's bounds are excluded from intersection-predicate copies,
never from deposited geometry.
Strict ownership still requires exact native coverage. An explicit single-sided
native-boundary analysis option can admit an uncovered *deposited* conflict
component only when its entire area is certified inside one native 7 nm
boundary band and exact segment-distance tests exclude the other material's
band everywhere. Failed containment alone is never treated as separation.
Near both boundaries remains ambiguous under that single-sided option; the
separate exact-only dual-boundary sensitivity described below requires whole
coverage by both bands and explicitly different material assignments. Near
neither boundary is unsupported. Both single-sided
priority variants receive the same uniquely eligible material; no deposited
void is filled and the full deposited union is still checked independently.
The exact helper reports the admitted area, and G-code-history results identify
that this assumption was enabled. This is a bounded representation assumption,
disabled by default and forbidden with strict ownership, not measured material
redistribution or evidence of interface bonding.
The sampler's explicit per-section work allowance can reach eight million for
the interface calculation and this separate topology certificate; its normal
default remains two million, and the interface routine retains its four-million
ceiling. Prepared ownership still has a separate cumulative allowance across
all queries. The exact topology certificate bounds each input to 2,048 vertices,
construction to four million work units, and arrangements to 4,096 vertices.

Only the generation-verified prepared resolver can mint the opaque replacement
token, bound to one live sampler request and its complete unchanged source
prefix. The section operation installs those full regions without clipping or
priority subtraction, then always checks occupied-union representation and
component, hole and material morphology safeguards. Interlayer and bed-contact
qualification remains downstream. Reused tokens, mixed callback modes, more
than two materials, native coverage gaps not admitted by an explicit boundary
sensitivity policy, and exhausted exact-operation budgets
are refused. The two assignments require independent numerical convergence;
their difference is a sensitivity measure, not an exhaustive material-uncertainty
bound, calibrated adhesion prediction or permission to enable Auto mitigation.

Section queries use native 1 nm Z planes by default. An explicit half-native
sampling mode evaluates 0.5 nm planes, allowing an interior midpoint for a true
1 nm deposited interval. Bead endpoint and thickness quantization and all XY
polygon coordinates remain unchanged. Direct sections, sweep activation and
ownership plane identity use the same selected resolution; a prepared native
resolver cannot be reused with a different mode. This is quadrature support,
not evidence of sub-nanometre physical printing accuracy.
Local topology validates contour/hole winding and net positive area with exact
integer sums. Global-coordinate floating-point shoelace sums are not used for
this decision: translating a native sliver must not erase its area or reverse a
tiny hole. This validation does not repair or discard invalid polygons.

An explicit `LocalSectionLimits::partition_increment_m` can construct bounded
partition envelopes inside each streaming cursor. Zero retains exact endpoint
partitioning. This reuses the existing analysis-only envelope sampler and applies
consistently to reference, deposited-material, thermal and mechanical streams.
Original beads are unchanged and are evaluated at the selected section planes;
missed positive beads retain exact local endpoints. Invented bed contact is
rejected, and caller-supplied envelopes cannot be combined with this option.
Geometry metrics mark envelope use, and G-code history results report the chosen
increment. It is an unqualified quadrature approximation requiring independent
refinement of volume, connectivity and fields, not an Auto-mode default.

When integer Boolean containment itself leaves a residual, an exact bounded
integer test may instead triangulate a changed polygon using existing vertices
and certify each triangle inside one original boundary-segment capsule. The
bounded triangulator's output requires unique input vertices, positive-area
triangles, and exact directed-edge cancellation matching every
original contour/hole edge. Area equality alone is not its certificate. Convexity then
contains the whole triangle, not just its vertices. Bounded longest-edge
subdivision uses exact dyadic midpoints (24 fractional bits) and 256-bit squared
distances; exhausting depth or work refuses the geometry. Different nearest
segments for different triangle vertices do not constitute a certificate.
Exactly collinear triangles have zero area and may be omitted only after both
directed chains are split at the same original collinear vertices; positive
triangles must still cover the complete original boundary.
If a rounded diagnostic ring is self-intersecting and cannot obtain that
positive-triangle certificate, a conservative fan cover is used instead. Fan
excess may be pruned only by exact dyadic predicates proving that no boundary
touches or crosses the subtriangle and its constant nonzero-winding occupancy
is zero. Boundary touches remain ambiguous and require coverage or refinement.
Ordered normalization of residual
material intersections is permitted only after authenticated ownership; each
changed material must pass the same band and topology checks, and the final
material polygons must have strictly empty intersections.
Adjacent slab contacts must survive that allowance, and absent contacts must
remain absent under dilation. Artificial tile faces use interval endpoint
margins rather than removing the tile plane. Narrow unresolved contacts fail.
The history reports the maximum coordinate allowance and conservative
added/removed volume bounds: absolute exact fan-triangle areas integrated
through sampled slab thickness. Signed areas of crossing diagnostic rings are
not treated as occupancy measurements. These are representation diagnostics,
not bounds on physical displacement, force, or adhesion. Such studies are
explicitly approximate and still require geometry and basis convergence.

Native studies can use a noncopyable `SampledSectionConflict`, constructed only
by the section operation that produced the actual conflict. Its borrowed
segments, clip, plane, bed and arc approximation are valid only during that
synchronous callback. The prepared resolver checks every segment index against
its immutable original observation-prefix geometry and verifies spatial-index
completeness across the entire sampled tile/plane. Registry replacement, omitted
deposition, altered coordinates/materials, or mismatched section precision
invalidates this evidence. This avoids recomputing the same integer polygon
intersection in a different order merely to prove its provenance. It does not
skip native source eligibility, unique ownership, deposited-envelope checks or
occupied-union conservation. An externally constructed overlap diagnostic
cannot request this path and retains the ordinary geometric proof.

Source-geometry diagnostics inventory actual parsed horizontal, sloped and
unresolved moves per verified source, including observed Z and nominal height
ranges. They skip synthetic display vertices when selecting a segment start and
retain unowned extrusion separately. These bounded summaries help distinguish
active path transformations from conservative configuration flags; they do not
reclassify transformed sources or reconstruct scarf thickness.

Internal control lines are removed before both public-file output and processor
line numbering. The export stream buffers incomplete lines up to one MiB and
batches ordinary lines between control records; it does not write a second
instrumented G-code file. Source IDs are stored on move records so timing-stage
copies retain their association. Imported comments cannot activate tracking.
Preserving the markers alone is insufficient. Pressure equalization carries
private source boundaries separately from its numerical line window and reinserts
them after output rewriting; adaptive pressure advance ignores them in feedrate
lookahead. Spiral transitions can duplicate comments with transformed moves and
remain unsupported. Processors require byte-equivalence and scope-completeness
verification before provenance can pass through them. The transcript detects
removed or duplicated whole scopes, but cannot detect ordinary commands moving
across boundaries or exchanges of occurrences with identical source IDs.

`visit_deposition_sections` streams increasing cut planes using a height-interval
sweep. It validates all input beads, including those outside the requested
planes, and retains original bead order in each section. This produces the same
bounded polygon unions as independent cuts without repeatedly scanning the
whole toolpath. Only the active section is retained; plane count and intersecting
bead visits have separate hard budgets. Consumer exceptions and cancellation
propagate. Previously delivered sections are partial geometry if a later cut
fails, and must not qualify a simulation. This indexing does not replace
through-thickness quadrature convergence. An optional per-plane provider can
supply matching ownership masks, using the same conflict-only rules.

## Section-integrated equilibrium

The section solver explicitly applies AMD ordering before numerical LDLT.
An optional offline sparse-factor budget admits larger background grids only
after an elimination-tree reach counts the factor's structural entries in that
same ordering. The reach uses linear auxiliary storage and has a separate work
limit; excessive fill is rejected before numerical factor allocation. Explicit
zero matrix entries count as structural entries. Numerical factorization uses
natural ordering on the already ordered matrix, and the observed factor count
must match the preflight. Positive pivots, unscaled/scaled equilibrium and
componentwise backward error remain mandatory: a sparse pattern does not prove
that the mechanical problem is well-conditioned.

Normal grid and assembly defaults are unchanged. Explicit sparse-factor mode
allows at most 64000 occupied cells, 96000 nodes, 192000 solved coordinates and
19.2 million assembly triplets, with a maximum 64 million factor entries
(approximately 733 MiB for double values and 32-bit indices). Geometry,
integration, symbolic work, matrices and assembly storage retain independent
bounds. The factor budget alone is not a bound on total process memory. This
option is not enabled by Auto or by preview and is not a convergence certificate.
The optional streamed integration work ceiling is 3.2 billion operations in this
mode, versus 600 million otherwise; callers still select the actual budget.
Moment groups may be explicitly raised to 125,000 in sparse-factor mode; their
default stays 31,250. Cumulative input points may be explicitly raised to 64
million, while live polygon and component-label limits do not grow with that
CPU-work allowance. Cancellation remains checked during
integration. A larger ceiling does not reserve storage: callers must separately
request every increased budget, and symbolic factor preflight remains mandatory.

For refinement diagnostics, bed reactions can also be projected onto an explicit
fixed XY Q1 reporting grid. The reporting planes align with the source basis and
use integer-multiple spacings. Applying the transpose interpolation operator to
the nodal reactions preserves coarse-grid virtual work, total force and moment.
It distinguishes a fixed-scale restraint functional from a sum of positive
nodal forces whose measurement scale changes during refinement. Opposite signs
inside a reporting basis may cancel, so this quantity is not local traction,
adhesion capacity, or evidence that an unresolved tensile peak is absent. It
does not replace the existing raw-force convergence criterion or qualify a field.

The local section solver can use an explicit, graded vertical Q1 grid while
keeping regular X/Y spacing. The Z planes start at the bed, increase strictly
on the native nanometre grid, and enclose all occupied slabs. At most 128 planes
are accepted; a grid that would omit any solid is rejected. Geometry is clipped
and integrated on the original deposited domains, with each cell's actual height
used for basis derivatives, moments, node coordinates and sample recovery.
No material or support is added between planes. Empty plane lists preserve the
regular grid. This is an analysis-only option, not a change to print layers.

Explicit grading can concentrate vertical resolution near the bed without
paying for the same resolution throughout the full height. It does not certify
convergence: plane insertion and lateral refinement still need independent
checks. The regular-grid directional selector rejects explicit graded grids
rather than returning an ineffective Z-spacing recommendation. Uniform-plane
equivalence, nested energy bounds, affine stress recovery and whole-solid
coverage are regression-tested.

The local section solver requires every component to reach the bed by default.
An explicit diagnostic `allow_free_components` mode also solves unattached
components, without joining them to adjacent material or inventing bed bonds.
Uniform eigenstrain on a free component is solved analytically by its exact
affine contraction, even with heterogeneous elasticity. Its recovered nodal
equilibrium is checked, and no bed force is assigned. This avoids solving an
ill-conditioned tiny cut fragment when its stress-free solution is known.
Nonuniformly contracting free components receive six independent scalar displacement references,
selected by rank-revealing QR of its rigid-motion basis. Since the selected
reference matrix has rank six, any displacement field can satisfy these
references by adding only translation and rotation; no strain mode is removed.
Recovery checks equilibrium at the reference coordinates as well as the solved
degrees of freedom. A nonzero reference reaction rejects the result. Reference
forces never contribute to reported bed reactions.

Free-component displacements have an arbitrary rigid-body frame and are not
predictions of where a loose printed island would move. Node and sample labels
identify whether their component reaches the bed. Free-component counts remain
visible in diagnostics; this mode cannot qualify adhesion or make a disconnected
deposition printable. Additional singular modes and numerical failures still
reject the solve. Neither Auto nor the preview enables this diagnostic mode.

`ThermoelasticSectionTopology` normalizes a supplied piecewise-constant slab
solid. Equal-material footprints are unioned; cross-material intersections and
misaligned overlapping Z intervals reject. Occupied islands connect across
adjacent slab faces only through positive area. Each resulting solid component
must contact the bed. A private-constructor certificate carries these components
into `ThermoelasticSections`; a coarse bounding box cannot anchor an island
above a real hole. This certifies the supplied slab geometry only, not the
accuracy of the sections used to construct it.
For zero-uncertainty exact footprints, the same retained grouping and union-find
algorithm uses rational unions, intersections, containment and component
extraction. Mixed input promotes native regions rather than rounding exact ones.
Canonical boundary validation rejects pinches within a component, including
point-touching islands represented as one pinched contour. It never bonds
components through a point. Inherited boundary uncertainty
still requires local contact qualification and is not cleared by exact storage.

Authenticated exact interface publication can additionally retain an opaque
complete-partition proof. It certifies unchanged total occupancy relative to
the admitted sampled envelopes, not physical bead accuracy. Local topology
checks the shared proof, complete material membership, sampling band and exact
footprint equality before treating occupied-boundary representation uncertainty
as zero. Material-allocation sensitivity remains on emitted pieces. Splitting
pieces invalidates the group proof; it is not a per-material permission flag.
Reference activation retains it only when the complete activated partition is
identical to the reference partition. General material-specific clipping cannot
inherit it because changing ownership can change which regions are active.
Exact reference activation shares one bounded per-tile sampling schedule across
the reference and all deposited materials. The planner includes endpoints and
coverage requirements from both inputs, while each cursor evaluates only its
own original beads. This prevents independently chosen midplanes from changing
the intersection merely because materials were traversed separately. The
schedule retains only its latest tile's intervals, shares the index budget,
and preserves partial-deposition clipping and all partition-proof checks.

The streamed topology engine represents native and exact tile-face traces
separately. Exact Boolean and overlap operations preserve rational spans,
including spans smaller than one native coordinate unit. Positive-length face
contact and positive-area Z contact join the existing union-find; isolated
touching points do not. Exact traces use the same frontier expiry and recycling
as native traces, without rounded endpoint conversion or unbounded caching.

The thermal skin stream can select exact processing for its entire traversal,
promoting native bands as needed. It subtracts occupied horizontal coverage and
oppositely oriented vertical contact spans before integrating exposed area.
Tile seams retain opaque rational spans and native Z intervals; bounding boxes
are indexing aids only. Adjacent seam bands compact only after exact coverage
checks in both directions. The same completion, cancellation, work and storage
limits apply, and no exposed-skin result is complete until the matching topology
stream finishes. Optional normal bins change convection direction metadata,
not surface geometry or integrated area moments.

`ThermoelasticSections` uses a regular Q1 hexahedral displacement basis over
that solid. Each physical component has independent basis nodes. Cell matrices
integrate the clipped material polygons and their Z intervals using polynomial
moments, retaining voids rather than filling coarse cells or applying a density
multiplier to modulus. The caller supplies constitutive tensors and eigenstrain.
The 27 monomial integrals are accumulated by cell, physical component and
material before stiffness and thermal-load assembly. This exploits additivity
over disjoint slabs without merging constitutive laws or changing occupied
geometry. Each stored group follows a bounded integration-work charge; the
existing topology, input, cell and unknown limits still apply.
Exact cut-cell moments are integrated at the actual cell origin, then rebased
to the occupied-support conditioning anchor. Outward bounds select candidate
cells and conservative support extents, never substitute material geometry.
Exact field membership distinguishes strict interior from boundary points.
Released-bed cells retain exact clipped footprints under a separate aggregate
trace budget. Contact inspection evaluates the cell's bilinear displacement
on every rational contour and hole edge, including interior edge stationary
points. Outward minimum bounds detect penetration and outward maximum bounds
report lift; rounded vertices cannot create or erase contact. These checks
support prescribed release and conservative compressive re-contact, not a
material-specific adhesion or fracture law. Refinement normalization uses
exact occupied area rather than the area of its indexing enclosure.
The retained local-topology factory and incremental `LocalTopologyBuilder` use
one connectivity engine. It emits borrowed normalized material pieces with
stable provisional union-find labels, retains face traces and connectivity metadata, and releases
historical polygons after their last vertical-contact use. Only successful
finish publishes an opaque stream identity and final component mapping after
lateral-contact and bed-anchoring checks. Callback failures and reentry poison
the builder. Consumers must separately establish complete source traversal and
match the same stream identity before using accumulated moments; the mapping
alone is neither a mechanical result nor printer-path provenance. Existing
cumulative work and geometry limits remain the default for both adapters.
Uncertainty-sensitive vertical contacts are recorded as bounded label pairs and
never union components. After all certified vertical and lateral unions, a pair
is admissible only if its endpoints already share a root through independent
certified paths. Otherwise topology remains unresolved. At most 100,000 pairs
are retained, additionally bounded by the streaming label or retained pair-work
limit; diagnostics distinguish observed ambiguous and proved-redundant edges.
This certifies only the global component partition, not the ambiguous local
contact area, conductance, stiffness or stress. Nominal occupied geometry,
surface integration and refinement requirements remain unchanged.
The local island erosion check likewise requires one surviving connected solid,
not an unchanged number of voids in the eroded diagnostic shape. Holes can merge
or open while that solid remains connected. Original holes and coordinates
remain in all emitted pieces and integration; no eroded geometry is substituted.
Complete erasure, solid splitting and uncertainty-sensitive contacts still
reject. This does not relax the separate source or represented-geometry
component/hole correspondence certificates, nor certify stiffness or heat flow
through a thin web.
Labels are allocated only for islands without a proved contact to the previous
band; inherited labels can merge but are never recycled. Historical face records
are separate from current/previous geometric islands. Interior islands with no
tile-face trace need no face record. Contiguous Z runs coalesce only when their
retained exact traces, uncertainty and already-proved root agree, preserving all
later positive-area lateral contact tests. Delayed lateral joins still require
final label remapping. Diagnostic processed-island counts are distinct from
stored labels and face records; they are not persistent geometry identifiers.
The closed sampler selects validated lexicographic tile order. Incoming neighbor
faces join current negative faces online; positive faces remain only until their
target tile is processed or proved absent by the ordering. Consumed faces release
their intervals and empty record slots are recycled. Arbitrary-order callers
retain the original final lateral-join pass. Both modes require positive Z and
trace overlap, preserve uncertainty checks and reject incomplete traversal;
neither removes stable labels needed by final moment remapping.
An explicit streaming topology budget separately limits cumulative processed
points/work and live polygon payload (including the borrowed tile, copied
groups, current/previous vertical polygons and temporary Boolean outputs).
Island metadata and retained tile-face trace intervals have independent caps
checked before insertion. Polygon outputs are counted after bounded-input
Boolean calls and before retention; this is not a preallocation bound on
Clipper's workspace, vector capacity or process RSS. Historical payload is
released without resetting cumulative counters. Partial statistics do not
certify connectivity. Mechanics piece/mapping and sparse solve limits do not
change when only the explicit topology budget is selected; its slab limit is
separately configured and remains 10,000 by default.
The closed nominal-deposition wrapper consumes complete coalesced material
bands through checked begin-tile/band/end-tile events. Only the current sampler
partition and current/previous topology islands remain live. An open tile cannot
be finalized; missing, reordered or mismatched packets poison the builder.
The retained whole-tile adapter uses the same band processing and still groups
unordered input slabs before processing. Gaps never create vertical contact.
Native Boolean output may encode a point-attached void as a repeated vertex in
an outer contour. The local factory decomposes such walks into exact original
directed-edge cycles, using bounded fixed-width multiprecision predicates.
Positive simple cycles remain separate islands; negative cycles remain holes.
A hole may share exactly one original endpoint with its unique outer contour,
but crossings, overlapping edges, multiple contacts, nested holes and
ambiguous containment reject the geometry. The same private certification runs
after unions, so point contacts are never welded and voids are never filled.
Repeated vertices in a hole likewise split into exact clockwise directed cycles;
positive island cycles remain unsupported. Single-endpoint contacts between
holes are allowed only when their contact graph is a forest and each connected
void cluster has at most one attachment to its owning outer boundary. A closed
void cycle or a chain joining two outer points can disconnect solid and is
rejected rather than represented as one mechanically connected island.
No new coordinates or rounded Boolean repair are introduced by decomposition.
Closed tile and nominal-deposition solve wrappers feed this engine directly
into the same private moment accumulator and sparse solver as retained
topology. Sample membership is resolved while each borrowed polygon is alive;
only bounded provisional component IDs are retained. Complete source traversal
and matching final connectivity precede component remapping, physical cell/DOF
limits and field assembly. The nominal-deposition wrapper returns traversal
metrics only together with a successful single fixed-basis field. It does not
perform convergence studies. No runtime UI or Auto path enables the larger
offline streaming budgets. By default it preserves the retained solver's
shared work and cumulative point limits. An explicit streaming-only integration
budget separates total polygon/moment processing from sparse assembly: at most
600 million integration work units and 32 million processed input points,
while each live input/clip remains at most 200,000 points. Provisional moment
groups are capped at 31,250 outside sparse-factor mode and sample memberships
at two million. Explicit
streaming integration permits at most 250,000 processed pieces (500,000 with
explicit sparse-factor opt-in) and 250,000 provisional component-map entries.
The higher piece ceiling counts cumulative traversal, not retained geometry;
retained/default integration keeps its 10,000 limits.
Final physical components remain capped at 10,000 and all sparse cell/node
bounds apply after component remapping. Remapping may temporarily hold two bounded group
maps. Sparse cell/node/DOF/triplet limits, topology limits and sampler limits do
not increase implicitly. The section solver's
explicit offline ceilings are 4,000 occupied cells, 6,000 nodes, 12,000 solved
DOFs and 1.2 million assembly triplets; defaults remain 1,000 cells, 6,000 nodes,
6,000 DOFs and one million triplets. A dense triangular factor at the DOF ceiling
would need about 0.81 GiB for double values and 32-bit indices alone. Geometry,
assembly matrices and allocator overhead are additional, so this is not a total
memory estimate. The larger bounds are opt-in offline budgets, not automatic
refinement defaults. Independently explicit offline topology ceilings are
128 million processed points and four billion CPU-work units. The higher work
ceiling is opt-in for offline studies, not a larger memory allowance. Slabs, persistent
labels, retained face records and current-band islands each have a 250,000
ceiling; their defaults remain 10,000 slabs and 100,000 for the other counts.
The processed-point/work defaults remain 8 million/20 million. Its
200,000-point default live polygon budget and two-million trace bound do not
grow with cumulative processing. Optional partial statistics report processed work/points and peak
logical storage counts even on refusal; they are not partial geometry or a field.
Positive-area bed traces are fully bonded; reactions are generalized basis
forces, not a local traction field or measured adhesion. Sparse factorization
rejects singular or numerically unresolved systems without artificial pins or
a stiffness floor. Geometry, assembly, and unknown counts have hard budgets.
Positive finite stiffness diagonals provide congruence equilibration
`B = S K S`, `g = S f`, `u = S y`, with `Sii = 1/sqrt(Kii)`. This invertible
change of coordinates neither adds stiffness nor removes weak modes. The same
positive-pivot threshold applies to the scaled system. Both original and scaled
equilibrium residuals must pass, together with a row-wise backward-error check
against `abs(K)*abs(u)+abs(f)`, so a stiff component cannot mask errors in a
softer independent component. Reactions and energy use the original operators.

A coarse basis can restrict deformation between connected arms separated by
voids. Numerical equilibrium in that subspace is therefore not a mesh-converged
mechanical result. Independent basis and slab-geometry refinement, appropriate
thermal/material scenarios, and an interface law remain necessary for physical
qualification. Background nodes outside the occupied solid are basis
coefficients, not physical material or predicted surface locations.

`refine_deposited_sections` halves a nested background grid on identical
certified slab geometry. Elastic energy and integrated tensile generalized bed
reaction must stabilize over two successive comparisons with increased basis
size and increased occupied-cell resolution in every axis. A grid that adds
lateral nodes but remains one cell through the thickness cannot qualify from
unchanging scalar values. Roundoff tolerances scale with input eigenstrain energy and eigenstress
times bed area. Incomplete or resource-limited studies retain only scalar
history, not a qualifying field; caller cancellation propagates. The result is
explicitly basis convergence, not validation of deposition geometry or adhesion.
Affine patch, loaded tetrahedral-reference, translated-grid, separate-body,
singularity, and resource regressions check this contract.

`probe_section_refinement` compares a reference grid with independent X, Y,
and Z halvings. It retains scalar energy and generalized tensile reaction,
not four displacement fields. Eigenstrain-scaled roundoff is excluded when
ranking the largest relative change. All probes must complete before a finer
spacing is recommended; resource exhaustion and cancellation cannot silently
favor a cheaper direction. This diagnostic helps choose aspect-aware grids,
but does not replace basis convergence or qualify adhesion. Coarse Q1 bed
reactions can average tensile and compressive regions together, so energy can
stabilize before the positive tensile-reaction sum does.

`probe_local_deposition_refinement` applies the same selector to the closed
streaming pipeline. It repeats sampling on the identical original segments and
optional partition envelopes, retaining only one field at a time. Geometry,
topology, integration and assembly limits apply to each solve; a probe consumes
at most four such budgets. It returns no field or convergence qualification.
Ownership and progress callbacks must be deterministic across runs. Exceptions
from ownership, including deferred conflict resolution, propagate as cancellation
or errors rather than being mistaken for an internal budget refusal. A partially
completed probe retains diagnostics but cannot recommend a refinement direction.

Section studies can explicitly enable up to three initial-grid probes. Each
selects one axis for additional resolution, after which the unchanged full-axis
refinement test starts from the selected spacing. Diagnostic probes never
count as convergence comparisons. All probe solves consume the same cumulative
work budget and report progress; incomplete probes or insufficient subsequent
refinement return no field. Seed selection is off by default. A synthetic
contraction regression exercises this path at a stated 20% refinement tolerance
and compares its energy and generalized tensile reaction against a separate
finer-Z grid within 5%; this is a numerical check, not an ASA adhesion threshold.

`ThermoelasticSectionStudy` adds independent section-height refinement around
every bead top and bottom Z breakpoint. It preserves material ownership and
requires stable occupied volume, elastic energy, tensile generalized reaction,
and component count over two successive geometry comparisons. Horizontal
rectangular envelopes without ownership masks are piecewise constant between
their breakpoints and need only basis refinement. Rounded or sloped beads do
not use this shortcut. Cumulative section, polygon, visit, and work budgets
bound the study; incomplete runs publish no displacement field.
Ownership polygon work reserves one cut's maximum before execution and retains
only the consumed charges afterwards. This bounds each cut without reserving
all cuts' mutually exclusive worst cases at once. Primitive operations that do
not report consumed work still reserve their maxima.

`ThermoelasticGCodeSections` applies this study to parsed deposition, with an
explicit processor-filament-to-mechanical-material map and nominal bead shape.
An optional observation time truncates deposition without changing the source
G-code. Tower, custom, and unclassified extrusion are excluded; support and
adhesion paths remain. Shared nominal-deposition admission rejects unresolved
stationary extrusion by default and reports any explicitly bounded omission.
Section G-code histories also reject authenticated zero-thickness collapsed
envelopes by default. An explicit `CollapsedEnvelopeAllowance` can omit only
reauthenticated collapsed spans that satisfy the shared per-span, total-volume,
total-fraction and contiguous same-source neighbor limits. The identical option
is passed to fixed-reference construction, every chronological prefix and native
ownership preparation. The result separately reports omitted collapsed-span
count and commanded volume; it does not add that volume to neighboring beads or
manufacture minimum thickness. These are declared numerical omission limits,
not stress-error bounds, reconstructed physical shapes or adhesion guarantees.
Section histories accept a separate `exposure_times_s` schedule for changes to
the full boundary law. These events reuse the active assembly and its node
numbering, replace only boundary operators, and neither reconstruct geometry
nor add capacity or deposition heat. Coincident activation/exposure times use
one activation event; changes before the first material activation are reflected
in that activation's full law. Cooldown changes are supported. The combined
input schedule is bounded to 10,000 entries, with independent thermal stepping
and boundary-operator budgets. An empty exposure schedule preserves the
activation-only behavior. Fan command times alone do not establish the physical
heat-transfer coefficients supplied by the boundary law.
Processor timing selects geometry snapshots; it is not measured thermal
history. Numerical convergence of this nominal solid does not qualify the
unknown bed interface, hot-polymer properties, or actual airflow, and does not
by itself authorize Auto mitigation changes.

The local study can explicitly coarsen only its Z integration partition: bead
endpoint tops round upward and bottoms downward to choose enclosing slab
boundaries. The interval grid is halved at every geometry refinement. Every
sample still evaluates the unchanged original bead dimensions, path and material;
rounded cross-sections are not enlarged or reshaped. Positive beads without any
interior sampled plane and partitions introducing bed contact are rejected.
Reported boundary displacement measures the quadrature partition, not a moved
physical surface. The resulting midpoint-prism geometry remains approximate:
stable volume, component count, elastic energy and generalized tensile reaction
over successive refinements are numerical checks, not proof of unresolved
continuous contact topology or bed adhesion.

Rejected section overlaps retain a diagnostic mask of complete polygons,
bounded to 4096 points and explicitly marked when partial. The G-code adapter
can intersect that mask with source beads to identify contributing move indices,
G-code lines, final filament IDs, and extrusion roles. Bounds restrict polygon
checks and reports, including a per-filament cap so an early tessellated curve
does not hide the other material. Incomplete diagnostics are labeled; they do
not assign material ownership, average ASA and TPU, or enable Auto by themselves.

`ThermoelasticSolver` solves three-dimensional, small-strain elastic equilibrium
on a supplied conforming tetrahedral mesh. Each element has an isotropic Young's
modulus, Poisson ratio and prescribed eigenstrain. Positions and displacement use
metres, stress and modulus use pascals, reactions use newtons, and stored energy
uses joules. No material constants, thermal history or adhesion strength are
silently inferred by this solver.

The constitutive equation is `sigma = D (B u - epsilon_star)`. Engineering shear
strains occupy the last three components (xy, yz, xz). For isotropic thermal
loading the normal eigenstrains are `alpha * (T - T_reference)` and shear
eigenstrains are zero. Spatial variation is represented by separate element
eigenstrains, not by colouring a single whole-object risk value.

The assembled system is `K u = f`, with element contributions
`K_e = V B^T D B` and `f_e = V B^T D epsilon_star`. This is the standard
[weakly coupled thermoelastic formulation](https://comet-fenics.readthedocs.io/en/latest/demo/thermoelasticity/thermoelasticity.html).
Prescribed displacement degrees of freedom are eliminated exactly; there are no
arbitrary penalty springs. Sparse LDLT factorization is checked for positive,
adequately scaled pivots, and the free-degree equilibrium residual is checked.
Invalid geometry, nonfinite data, unrestrained components and failed solves throw
instead of producing a zero-stress result. Cancellation checkpoints surround
factorization and occur during assembly and recovery; factorization itself is
not interruptible.

The result contains nodal displacement, reactions on constrained degrees of
freedom, element stress and total elastic energy. Element stresses are deliberately
not smoothed across voids or material interfaces. Boundary reactions describe the
force needed to maintain the supplied restraint; they are not an adhesion limit
or a detachment probability. A fully bonded-bed calculation holds those nodes
fixed, so its zero bed displacement cannot be interpreted as proof against lift.

Linear tetrahedra are overly stiff in coarse bending meshes. Analytic tests cover
free contraction, fully constrained hydrostatic stress/energy, an affine shear
patch, rejected singular/invalid systems, and cantilever bending with refinement.
An imposed longitudinal strain gradient `epsilon_xx = g*z` also has the exact
stress-free displacement `u_x = g*x*z`, `u_z = -g*x*x/2` with a fixed `x=0` end.
The gradient-curl regression refines the tetrahedra in the length/thickness plane
and requires decreasing tip-displacement error, below 2% at its finest mesh.
This tests thermal curvature directly, rather than relying only on force-balance
residuals or a qualitative upward-displacement check.
Those tests verify the numerical core, not a printed part. A caller needs mesh
convergence, representative thermal histories, deposition/infill geometry,
boundary-condition scenarios and material uncertainty before using the results
to rank mitigations. Contact, fracture, viscoelastic relaxation, plasticity and
large deformation are not solved here. Near-incompressible material additionally
needs an appropriate formulation to avoid volumetric locking.

The existing `WarpMitigation` screen remains distinct from the equilibrium
calculation; its results must not be presented as spatial solver output.
`mechanically_improves` accepts an automatic candidate ranking only when every
matched scenario reduces both integrated tensile demand and elastic energy
beyond the two results' discretization tolerance bands. It rejects missing,
unconverged and incomparable results. It does not itself generate candidates,
authorize weaker adhesion, or justify an infill rotation with isotropic material
properties. Candidate generation and geometry-raster convergence are required
before this comparator can drive automatic toolpath changes.

## Periodic thermal transport closure

`ThermalHomogenization` solves three affine-temperature-gradient conduction
problems on explicit solid tetrahedra. It shares periodic node and boundary-trace
validation with the mechanical homogenizer. Each tetrahedron supplies a positive
definite world-frame conductivity tensor and positive volumetric heat capacity;
voids are insulating and do not acquire fictitious mass. Face-connected solids
have ideal thermal contact. No air transport, radiative exchange, contact
resistance or material parameters are inferred.

The caller supplies a conforming, non-overlapping tetrahedral mesh. Shared
admission rejects unused or coincident nodes, degenerate or duplicate elements,
over-incidence of faces, excess total volume and mismatched periodic traces;
it does not detect arbitrary tetrahedron intersections or interior hanging nodes.
The volume sum alone is not a proof that an arbitrary imported mesh is valid.

Periodic temperature fluctuations have one additive-temperature gauge per
connected periodic component. An isolated solid island therefore retains heat
capacity but has zero steady macroscopic transport. The conductivity tensor is
computed from volume-averaged flux response; capacity is averaged over the total
cell volume. Symmetry, passivity and microscopic/macroscopic energy consistency
are checked against the affine Voigt scale, which remains meaningful for
zero-transport directions. Admission is bounded to 2,048 nodes and 32,768
tetrahedra, with cancellation between assembly and solve stages.

Solid, series/parallel laminate and explicit-void analytical tests verify this
numerical closure. It is not a calibrated printed-material model or a transient
subcell equilibration model. Applying it to a coarse print grid requires a valid
representative periodic geometry and scale separation; local exterior surfaces,
bed contact, nonperiodic walls and deposition activation cannot simply be replaced
by a filled homogenized envelope.

## Geometry-preserving Galerkin reduction

`solve_reduced` accepts an explicit caller-provided interpolation from coarse
displacement unknowns to fine nodes. It assembles the restricted strain operator
`B P` on the original tetrahedra, preserving their occupied volume, material,
eigenstrain and void boundaries. It does not construct a coarse solid envelope,
automatically choose interpolation or invoke periodic homogenization.

Positive interpolation weights must reproduce constants and coordinates. Each
basis has connected positive support within one fine-mesh component, preventing
a basis from tying unbonded islands or opposite slit arms through a distant
connection with zero support. Constrained fine nodes require exact unit rows;
their prescribed displacement is eliminated exactly. Input sizes, interpolation
width, free reduced unknowns, assembly work and sparse triplets have hard limits.

The reduced residual measures `P^T r`, not fine equilibrium. `ReducedResult`
therefore separates recovered fine fields, the norm of residual forces at free
fine nodes and generalized constrained reactions. Generalized reactions include
weighted residual contributions from unconstrained fine nodes and are not a
local bed-traction distribution. The fine constrained forces alone are likewise
not qualified contact demands until comparison with a richer interpolation or
fine reference establishes reduction accuracy. A two-zone eigenstrain coupon
tests this distinction explicitly. Nested interpolation on a hollow printed
shell provides a fine-reference convergence test without changing geometry.

`occupied_line_basis` constructs a bounded interpolation along one selected
coordinate axis using actual tetrahedral edges of an occupied-grid mesh. Missing
edges split lines; coincident but distinct node IDs are not welded. Line endpoints,
stride anchors and explicitly pinned nodes remain independent unknowns. All bed
nodes are pinned by the deposition simulation, preserving its exact constraints.
This reduces unknowns without filling cavities or bridging slots. It does not
reduce the original tetrahedral integration or thermal-grid storage requirements.
The deposition simulation retains its 12,000 fine-node admission cap; this option
reduces factorization cost within that cap, not the admitted whole-part size.

An explicitly configured `SimulationLimits::reduction` enables nested line-basis
comparisons in `simulate_deposition_grid`. Relative displacement and elastic-energy
changes must satisfy caller-supplied tolerances on two consecutive comparisons;
an identity basis instead supplies the exact fine-discretization reference.
Insufficient comparisons or budget return an unqualified reduction state.
Reduced fields are returned separately from `BedResult`, never converted to local
bed traction. Basis stabilization is not spatial/thermal convergence, validation
of material or adhesion data, or automatic mitigation approval. The default path
continues to solve the unreduced bonded-bed scenario.

## Bonded-bed boundary scenario

`solve_bonded_bed` constrains all three displacement components on supplied
exterior tetrahedron faces in the bed plane. It rejects duplicate, internal,
off-plane and degenerate contact faces. Contact area comes from those faces,
not from a bounding rectangle that would incorrectly bridge holes or gaps.
Disconnected components without contact remain mechanically unsupported and
are rejected by the equilibrium solve.

The bed's downward nodal reactions are reported as tensile restraint demand in
newtons. Dividing by lumped nodal contact area provides a traction field in
pascals. Positive and negative normal reactions are not cancelled when computing
total tensile demand: a compressive center does not remove a corner's need for
adhesion. The field describes the fully bonded scenario only, not progressive
peeling. No adhesion capacity or detachment safety factor is assumed. Local
corner peaks can be mesh-dependent; they are not standalone acceptance criteria.

## Prescribed partial bed release

`solve_partially_bonded_bed` accepts the complete original footprint and an
explicit subset of triangle indices that remain bonded. Every footprint face
is validated as a distinct exterior tetrahedron face on the bed, including the
released faces. Nodes belonging to retained faces remain fixed in all three
directions; shared crack-front nodes therefore remain fixed too. Released faces
carry no applied traction. Components without enough remaining restraints are
unresolved, not stabilized using artificial pins.

The result separates retained restraint demand from released area and signed
footprint displacement. Positive displacement is lift; negative displacement is
penetration into the rigid plate. This traction-free release model does not
enforce unilateral compressive contact. Any penetration beyond a
coordinate-scaled floating-point roundoff threshold marks the scenario
`nonpenetrating=false`; such a result is **inadmissible as a bed-contact
solution**, rather than evidence of low lift risk. Linear footprint triangles
interpolate their nodal gap, so checking the footprint nodes also checks their
interiors under this displacement approximation.

The retained region is a supplied boundary condition, potentially taken from
an observed adhered region. It is not inferred from temperature, a guessed
adhesion strength or an image alone. This API neither predicts when peeling
begins nor evolves a peel front, and its lift still requires mesh convergence
and defensible constitutive and thermal inputs. It preserves the fully bonded
scenario when all footprint faces are retained. Releasing restraints reduces
the constrained elastic energy but does not imply that a practical mitigation
has improved adhesion.

`solve_partially_bonded_contact` adds frictionless unilateral contact to the
same prescribed released footprint. Its bounded active-set iteration adds a
normal restraint when a released node penetrates the bed, and removes that
restraint if its reaction would pull the node down. Retained bond nodes stay
fixed in all three directions; released contact nodes constrain only Z. A
result is returned only when the nodal gaps and normal reactions satisfy the
supplied absolute tolerances. Nonconvergence, an exhausted iteration budget or
unrestrained rigid-body modes throw without publishing a partial result.
Compressive contact reactions are reported separately from retained adhesive
tension. Contact does not restore adhesion or introduce tangential friction.
The model uses the original linearized footprint, not large-sliding contact or
a surface-to-surface finite-deformation formulation. It does not determine
which bonds fracture.

`bed_release_demands` evaluates strictly nested prescribed contact-loss regions
at one fixed eigenstrain/material state. With no applied mechanical loads, it
uses the decrease in minimized stored elastic energy divided by newly released
area as a finite-increment energy-release demand (J/m²). This is the resistance
an interface would need for that prescribed increment, not a measured PEI
fracture energy or a local crack-tip estimate. It does not combine states at
different temperatures or deposition stages in one energy difference. Increasing
energy after releasing constraints, bed penetration and unrestrained components
invalidate the study. There are at most 16 release steps and no mesh output on
disk. Mesh and release-increment convergence remain necessary; matching a
photographed peel boundary is not independent validation of an adhesion law.
Callers may explicitly supply `BedContactLimits` to use unilateral compressive
contact for every release state. Energy differences then compare equilibria
with the same rigid-bed inequality constraints and progressively fewer bonded
faces. The default remains traction-free release. An unconverged contact state
rejects the whole study, including previously completed release steps.

`evolve_bed_release` evaluates conditional finite-region debonding at one fixed
thermal/material state. The caller supplies disjoint candidate footprint
patches, an initial bonded footprint and an explicit fracture work per unit
area for each patch. Every trial re-solves unilateral contact. The chosen
increment maximizes the positive difference between the decrease in elastic
energy and the new fracture work, beyond absolute and relative decision
tolerances. Accepted release is irreversible; the next trials start from that
updated contact state. No mechanical external loads are included, so elastic
energy is the relevant stored potential energy in this comparison.

Completion means that no *individual supplied patch* has a favorable next
increment. It does not prove stability against combined releases, differently
shaped cracks or cohesive crack initiation. Patch size, mesh resolution,
constitutive uncertainty and interface fracture data require independent
qualification. Fracture energy is not inferred from filament name, bulk
tensile strength, brim width or a photograph. All contact iterations across
all trial states share a bounded solve budget; an incomplete study throws and
does not publish its intermediate sequence as a completed result. Completely
unrestrained components remain unresolved instead of acquiring artificial pins.

Continuation on the same mesh and patch partition uses the previous retained
face set as its initial bonded footprint. Entirely released patches retain
their identities but are skipped; a partially bonded patch is rejected rather
than silently transferring damage onto a different partition. Unloading can
close a gap through compressive contact, but does not restore adhesion or charge
fracture work again on reloading. This carries contact damage between explicit
states only: it does not supply viscoelastic memory, deposition onto the
deformed surface, or a converged observation schedule.

The deposition-grid simulation also accepts prescribed `BedReleaseEvent`s.
Each event names one stable input-grid bottom cell, an explicit release time
at or after its activation, and effective gap heat exchange. Events merge with
deposition and boundary-change times; they do not add mass or reset cell
temperature. From release onward, that bottom face uses the supplied gap
exchange instead of the global bed condition, including after subsequent
bed/chamber changes. Half-cell conduction resistance remains in series with
the surface exchange. Release is irreversible and duplicate releases are
rejected. These cell IDs are specific to the supplied grid, not portable to a
refined grid without an explicit geometric transfer.

The final mechanical snapshot removes released bottom faces from adhesion and
solves frictionless unilateral contact. It is returned as `partial_mechanical`,
not as a fully bonded result, so existing fully bonded convergence/ranking
code cannot silently consume it. Input-order voxel identity maps the active
thermal cells to mesh bed triangles even when deposition order differs.
Reduction with release histories is rejected. This path is a prescribed-damage
experiment, not a coupled predictor: it does not infer release times or change
gap heat transfer from the solved displacement, and the supplied gap exchange
does not automatically change on compressive re-contact.

`simulate_peeling_history` adds sampled feedback using explicit grid-cell
patches and interface fracture-work scenarios. At each observation it solves
the thermal snapshot, maps the same patch identities onto the active mesh,
and evaluates finite-patch release. Newly accepted releases become thermal
contact-loss events at that observation time for subsequent snapshots. The
thermal history is replayed with those events, preserving mass and deposited
energy rather than resetting temperatures at observation boundaries. New
upper cells may activate between observations; patch cells must already be
active at the first observation. Fully released patches never heal.

At most sixteen bounded thermal runs are made. All mechanical snapshot and
contact-trial solves share a global budget. Only scalar history, release events
and final fields are retained; failure rejects the study rather than publishing
partial success. Final mechanics lives in the release result, not a stale
fully bonded snapshot. Release is permitted only at supplied observation
times, so observation refinement is essential; these are conditional sampled
results, not established onset times. Thermal/grid/patch convergence, interface
uncertainty, gap-dependent heat transfer and deformed stress-free deposition
remain necessary before predictive Auto use. In particular, passing feedback
bookkeeping tests does not validate a ValveBoard adhesion margin.

`study_peeling_time` qualifies time discretization independently on a fixed
grid, patch partition and material/interface scenario. It uses three nested
observation schedules, inserting interval midpoints twice. At each schedule,
an independent thermal refinement halves the step cap and quarters the local
temperature-error tolerance. The cap also accounts for the largest interval
actually accepted across that history, so an inactive configured cap cannot
manufacture an identical nominal refinement. Adaptive thermal error control
is required.

Agreement includes final cell temperatures and nodal displacement, common-time lift, tensile force,
stored energy and released area, sampled force/energy/lift maxima, final heat
and boundary heat. Force and energy immediately before each release update
are included, so stress relieved by damage does not erase the pre-release
demand from sampled maxima. The released cell identities must agree, and their release
times must agree within an explicit absolute time tolerance. Equal detached
area at different locations or materially different release times is not a
qualified result. Kelvin temperatures use an absolute tolerance, not a relative
tolerance inflated by their temperature reference. A result is published only
after all three thermal comparisons and both observation refinements agree.
Six history runs share limits on total observations and mechanical solves;
budget exhaustion and cancellation propagate without partial publication.
Temporal qualification does not qualify spatial discretization, patch choice,
constitutive behavior or the actual bed interface, and is not an adhesion
safety verdict.

`study_peeling_space` applies the temporal study to three nested discretizations
of the same occupied voxel volume. Each refinement replaces a parent with eight
equal-volume children, retaining its material and activation state. Air gaps
remain gaps, including at negative world coordinates. Bottom release patches
inherit only their four bottom-face children, retaining physical area and
fracture work per area. The grouping of candidate release regions does not
change. This refines a prescribed voxel model, not the original G-code
rasterization or physical extrusion/activation history.

Both successive spatial comparisons must agree in sampled lift, force, energy,
released area, pre-release demand, heat and patch release identity/time. Final
child temperatures are volume-averaged onto their parent cells for comparison;
node indices from different meshes are never equated. This does not qualify
local crack-tip stress singularities. No field is published when either time
or spatial comparisons fail. All eighteen possible history runs share budgets
for observations and equilibrium solves; refinement checks cell capacity before
allocating children, and the existing mechanical node limit remains in force.
The qualified result includes its refined grid because cell IDs are not those
of the input grid. Constitutive/interface uncertainty and alternative patch
partitions still require separate evaluation.

`study_gcode_peeling` connects parsed deposition to that coupled study through
the same bounded preparation used by `simulate_gcode_deposition`. Preparation
clips deposition at the final observation, retains first-hit cell activation
provenance, and rejects unresolved extrusion instead of inventing missing
material geometry. The input must explicitly acknowledge modeled timing and
supply thermal, mechanical and interface properties. World-XY, half-open
release regions select bottom-cell centers, optionally by material; empty or
overlapping selections are rejected. Selected regions remain the same physical
voxel footprints during solver refinement. This connects the checked solver
to actual parser output, but does not turn initial rasterization assumptions
or a supplied interface scenario into measured print behavior.

Stationary positive extrusion remains a strict rejection by default. An
explicit `StationaryExtrusionAllowance` can admit only zero-displacement
extrusion immediately following a same-material, same-role bead at its exact
endpoint. It does not admit vertical paths or isolated dots. Absolute omitted
volume, its fraction of the total commanded line volume, and its cumulative
fraction of the preceding G-code command are bounded independently; arc
tessellation pieces are counted together for that command. Admission ceilings
are 0.001 cubic millimetres, one part per million globally, and one part per
thousand locally. These are approximation admission limits, not demonstrated
print-accuracy guarantees. Nominal geometry and thermal mass are unchanged,
and returned snapshot/peeling results report the stationary command count and
commanded volume separately. The nominal-volume discrepancy is not hidden by
changing density, scaling another bead, or calling the missing deposition
spatially resolved. Sensitivity to this explicit approximation remains part of
physical model qualification.

Attached brim area is not a dimensional-success criterion. The `BrimLift`
mechanism test retains the entire bed contact of a thin brim while releasing
the adjoining thicker case strip. It compares brim widths without changing the
case or its prescribed contraction gradient, and checks displacement against a
0.2 mm layer height. A connected, attached brim may coexist with case lift large
enough to interfere with subsequent layers. Halving prescribed strain must
halve displacement and quarter elastic energy in this linear model; that
scaling does not establish how a fan command changes actual locked-in strain.
This deliberately simplified strip is neither the ValveBoard geometry nor a
calibrated ASA/PEI interface. It does not predict when contact releases.

## Mesh refinement and spatial lookup

`mesh_voxels` splits occupied cells into six conforming tetrahedra and leaves
unoccupied cells empty. Nodes are joined only through shared faces. In
particular, an edge contact does not become a simulated weld just because the
two regions connect elsewhere in the body. Bed faces belong only to occupied
cells at the bottom of the grid. The caller supplies the grid, material data,
strain field and a node budget; exceeding the budget is an unresolved analysis,
not permission to merge gaps or return zero stress.

`solve_refined_voxels` subdivides those cells while preserving their volume,
material boundaries and piecewise-constant strain field. Integrated tensile
restraint and elastic energy must both stabilize over two successive
refinements. A single solve never qualifies as converged. This check concerns
discretization of the supplied model, not the accuracy of its physical inputs,
and does not establish convergence of singular corner-traction peaks.

`Thermoelastic::Field` uses a bounding-volume tree followed by tetrahedral
barycentric containment tests. It interpolates displacement within an element
and returns element stress without smoothing across material interfaces. It
does not extrapolate into voids. Maximum principal tensile stress is distinct
from von Mises stress: hydrostatic tension has a nonzero principal stress even
when von Mises stress is zero.

The preview renderer accepts per-move tensile stress in MPa. Missing samples
are gray, not zero, and the legend identifies unavailable analysis. Re-reading
postprocessed G-code invalidates samples even when the move count is unchanged.
The renderer does not itself construct a physical model or infer an adhesion
safety margin.
The configured field describes the final deposited geometry. Moving the preview
layer selector filters displayed paths; it does not recompute stress at that
point in the printing history. The legend identifies this as a final-state
static estimate.

## Deposited-path adapter

`voxelize_deposition` rasterizes deposited segments using their
reported width and height. The XY footprint is a capsule, the top elevation
is interpolated along the path, and vertical thickness remains constant;
this is a deposition-envelope approximation, not a resolved
bead cross-section or inter-bead weld model. Grid cells are selected by center
containment. The grid must resolve at least two cells across the width and
height. Infill voids are not replaced with the CAD solid envelope. Conflicting
material assignments in a cell, vertical-only paths, below-bed envelopes,
invalid dimensions and exceeded
resource budgets produce an unresolved analysis instead of an averaged material.
Sloped paths are not flattened to their ending height. Variable scarf thickness
and extrusion-volume-dependent cross sections are not inferred from the nominal
processor height, so this support does not establish physical accuracy for a
scarf joint. The section-based adapter separately accepts explicit linearly
varying endpoint design heights. Its ruled normal sections retain original
endpoint half-caps and use concavity-bounded refinement of rounded sides; this
is not a swept variable-radius disk union. Zero thickness is allowed at one
endpoint, but a wholly zero-thickness bead is rejected. The voxel adapter samples
the same explicit affine thickness at each cell's projected path fraction,
retaining the tapered empty space. Its minimum-resolution guard uses the maximum
endpoint thickness: no finite grid resolves the zero-height tip exactly, so
independent occupancy/grid refinement remains necessary.

Opt-in live exports capture scarf and constant-height perimeter command
dimensions before post-processing in a bounded private registry. Ordinary
perimeter commands in transformed source scopes retain their original XY-only
output; native-planar scopes need no redundant command records. Fitted XY arcs
retain their G2/G3 output and record the quantized center and direction. Child
chains check circular membership and angular order, including the analytic radial
drift bound of the legacy preview's truncated rotation. The explicitly commanded
last endpoint retains its own quantized position. Stationary/full-circle and
Z-contoured commands remain unsupported. Parent storage remains capped at 80 MiB;
the command-count ceiling accounts for the arc and lower-surface fields. Complete ordered child chains must match their
quantized parent commands before these dimensions are published to mechanical
consumers. Private controls are stripped from public G-code. Missing, reordered,
displaced or incomplete commands invalidate the captured dimensions. Imported
G-code has no such registry. Nominal design volume is integrated separately
from commanded E volume, including partial-time cuts; neither is a measurement
of the deposited bead or a qualification of interface bonding.

Leading scarf ramps whose Z and extrusion ratios coincide may also carry an
explicit nominal lower surface from the matching live native layer. This requires
the path height and nominal Z to match that layer; a base is not inferred from
an arbitrary tapered command. Authenticated endpoint heights use double precision
and are computed from each emitted child top minus this lower surface. The base
is transformed into the parsed coordinate frame without display-float rounding.
Anchored commands retain sparse precise endpoint Z spans, separately from the
float preview vertices. Consumers require the same sealed source and command
generation, matching move/source/command IDs, both original preview endpoint
projections, and command membership with identical derived heights. Missing or
inconsistent spans reject instead of falling back to rounded display Z. Ordinary
unanchored commands keep their existing coordinate representation.
This prevents command rounding or pressure-smoothing subdivisions from creating
false penetration below the native layer. It does not change extrusion volume or
public G-code. Negative envelopes reject. Exactly zero-height anchored spans
retain their authenticated identity, timing and positive commanded extrusion as
unresolved collapsed envelopes, not as zero-thickness beads or ordinary-height
fallback geometry. Spatial analysis rejects these spans; the stationary endpoint
allowance does not admit them. Bounded diagnostics report their volume and nearby
same-source bead candidates without authorizing omission. Unanchored return
ramps retain their separate flow-derived nominal envelope. The lower-surface
assumption is design geometry, not measured bead thickness or weld strength.

Transformed source ownership requires the command centerline to remain inside
its native material region. If integer line clipping reports an exterior
fragment, a bounded exact predicate partitions the segment at every rational
boundary crossing and checks every intervening interval against the closed
material union. Native oriented chains are cropped rationally to the segment
bounding box padded beyond the requested boundary allowance by two nanometres,
then resolved together with nonzero winding. Artificial crop edges cannot
certify boundary-band coverage. Input scanning is work-bounded; streaming the
first clip plane avoids allocating a full rational copy of remote chains, while
intermediate, retained-chain and arrangement limits remain unchanged. This is
local coverage, not a certificate of unrelated rings' validity. Only complete exact coverage can override the rounded clipping
result; material masks are not expanded and finite gaps remain exterior.
Rejected centerlines also report whether every exterior interval fits entirely
inside a4nm capsule around an original union boundary edge. This uses exact
endpoint distances and capsule convexity, not samples or rounded offsets.
The diagnostic band result is distinct from strict coverage. The general
predicate permits an explicit0..7nm sensitivity radius while retaining the
strict result separately. Source ownership defaults to strict coverage; an
explicit source-boundary sensitivity can admit only the exactly certified
near-boundary intervals, without changing source masks or command authentication.
The conflict resolver permits that option only alongside an explicitly
non-strict interface study. These are numerical sensitivity cases, not physical
material assignments or qualified curl predictions; ordinary Auto setup does
not enable the allowance. Original-void and deposited-partition checks remain
independent and unchanged.
Source verification has its own bounded work allowance, and successful work is
also charged against the prepared resolver's cumulative budget.

`analyze_gcode` obtains those segments from processor moves, excluding travel,
custom purge paths and the prime tower. Both static analysis entry points ignore
display-only acceleration subdivisions and seam-marker positions when tracking
physical endpoints, counting materials, and sampling stress. These markers
retain unavailable preview samples rather than becoming extra deposition.
The static adapter uses explicitly supplied material
and eigenstrain scenarios. After equilibrium refinement converges, it samples
principal tensile stress at the bead-center midpoint of each extrusion move;
missing samples remain unavailable. The sample is not the maximum over an
entire long segment. Samples are published together, never partially, and old
samples are cleared before starting or on cancellation.

Refining tetrahedra within an occupied grid checks equilibrium discretization,
not convergence of the rasterized deposition geometry. Comparing different
grid resolutions, physical loading scenarios and interface assumptions is still
necessary before using this adapter to rank automatic changes.

## Configured scenarios and bounded execution

`SectionIntegration` supplies shared occupied-volume monomial integration for
section mechanics and thermal elements. Polygon holes are subtracted, and Z
subintervals are integrated in the background cell's normalized coordinates.
The thermal element uses consistent trilinear heat capacity and a symmetric
positive-definite conductivity tensor. It does not replace empty volume with a
density-scaled solid. Element assembly must keep physical components separate;
moment integration itself is not a connectivity certificate.

`exact_section_volume_moments` provides a bounded rational alternative for
Boolean cut cells. It forms native nonzero unions, performs the requested
Boolean operation and cell clipping without rounding intersection vertices,
then integrates all 27 normalized monomials exactly using directed boundary
integrals and analytic Z integration. Only the final SI moments are converted
to doubles, accompanied by outward conversion intervals. Holes retain their
signed contribution. Input, arrangement and work limits bound each query.
This primitive does not certify contacts or replace the integer-boundary
growth and exposed-surface pipeline; its conversion intervals are not total
simulation-error bounds.

Area-only metadata and occupied-island volume use `exact_section_horizontal_area`
and `exact_section_volume`. These scalar queries preserve the same exact cell
clipping and Z interval semantics without constructing the other 26 moments.
An outward integer enclosure wholly inside the cell proves that XY clipping is
unnecessary; other cuts use the exact Boolean intersection. Directed rational
boundary area and integer Z thickness are combined before conversion to an SI
representative and outward interval. Tangent holes keep their signed area.
Full moment queries remain necessary for mechanical and thermal operators.

`ExactSection` is an immutable, shared, bounded rational set for composing
Boolean operations and extracting connected components without an integer
polygon round-trip. CGAL types remain private to the geometry implementation.
Its integer bounding box is an outward enclosure for indexing only; moment
integration uses the canonical rational set. Each operation reports work for
the caller's cumulative budget, and consumers must account for live retained
vertices. This representation is not itself a connectivity or thermal-exposure
certificate and does not enable automatic admission of an unsupported study.
Boolean identities involving an empty operand or the same immutable handle
return an existing handle or an empty set without rebuilding an arrangement.
They debit metadata work and retain budget, cancellation and edge-limit checks.
Distinct nonempty sets with equal edge counts and outward bounds may prove identity by comparing
their complete, sorted exact directed boundary edges, including holes and
disconnected components. This charged, bounded comparison does not round or
merge edges; differing subdivisions fall back to the fully budgeted overlay.
Consumers account for returned vertices even when their storage is shared.
For an actual overlay, a charged sweep of closed outward edge boxes counts
possible cross-operand contacts. Each operand is already planar; the input
edge count plus twice this candidate count bounds the split-edge event model.
The sweep selects X or Y using a charged linear estimate of projected edge
density. This changes search order only: the same closed two-dimensional box
predicate decides candidate pairs. Invalid estimates and ties use X; the
heuristic does not guarantee less work for every geometry.
The reservation includes a logarithmic sweep factor and input copying, with
output traversal charged separately. Preprocessing, including rejected active
box pairs, consumes work too. Touching and collinear boxes remain candidates;
zero candidates still require the exact overlay because one region may contain
the other. These are bounded geometric-complexity units, not a CPU-time or
rational-bit-complexity guarantee. Independent edge, memory and cancellation
limits remain necessary.
The selected exact interface-partition factory retains the existing native
ownership and interface-sensitivity checks, but publishes rational material
sets before integer rounding. Publication requires disjoint materials and
exact conservation of the complete deposited union, with an aggregate retained
vertex limit. This selected-output path constructs only its requested assignment
and uses event-sensitive overlay accounting; the legacy two-variant path still
constructs both assignments. Native source coverage or explicitly enabled
boundary sensitivity, and both interface-band checks, precede selection. The
inherited interface allowance remains explicit; zero new
representation allowance does not qualify the chosen material scenario as
physical truth. Existing integer consumers do not automatically use this path.

The exact-only `allow_dual_native_boundary` option is a separate, default-off
material-location sensitivity. Its candidate region is the intersection of
both complete deposited envelopes minus both native reference masks. Every
connected candidate must pass the existing whole-region exact certificate
against **both** native-union boundary bands at the unchanged 7 nm allowance.
Vertex proximity, small area, or failure to prove separation is insufficient.
The selected priority assigns that already deposited region to the first or
second material; final exact disjointness and deposited-union conservation are
still required. No modeled occupied region is added or removed. Exclusive
deposition and regions supported by only one boundary are not authorized by
this option; single-sided sensitivity remains independently controlled.

This option requires authenticated non-strict sampler ownership and retained
exact sections. Strict and rounded paths reject it. The exact result and its
occupied-partition metadata retain the admitted ambiguity area at the sampling
plane; that area is not a cumulative physical area or volume bound. History
results retain the enabled policy and selected priority. Both priority histories
need separate assessment before any future qualification, alongside temporal,
spatial and constitutive checks. Their difference is conditional numerical
sensitivity, not an exhaustive physical uncertainty bracket or evidence of
bonding. A history using only one priority is not sufficient qualification.

`SectionFootprint` stores either native polygons or an exact section, never both.
An exact empty set is authoritative, not a signal to use an integer fallback.
Legacy consumers request native geometry explicitly and reject exact handles;
this prevents partially integrated geometry from silently disappearing or being
rounded during assembly. Vertex counts and outward bounds can be queried without
changing the stored representation.

`ExactSectionBoundary` retains directed rational segments with material on the
left, including clockwise hole boundaries. Boundary intersection and subtraction
operate on positive-length collinear intervals; isolated point crossings are not
contacts. Callers explicitly choose whether right-hand spans must have opposed
orientation. Plane selection returns existing boundary spans on a native plane,
not an interior cross-section: tile-face extraction therefore requires prior
exact tile clipping. These operations establish nominal boundary coincidence,
not robustness to inherited interface uncertainty. Their storage and work are
bounded independently from the owning section's arrangement.

The horizontal-surface moment operation uses the same exact clipping and XY
integrals, with SI area and fixed normalized Z powers. Its caller removes
covered/contact regions, determines orientation, and assigns each exposed
face to one incident cell so shared boundary planes are not counted twice.
Vertical-face moments clip directed rational spans to the cell and integrate
the normalized polynomial factors exactly. Physical segment length can be
irrational: floating square-root estimates are bounded by exact squared
comparisons and tightened to adjacent doubles. An exact-rational seed refresh
and bounded bisection handle cancellation in the initial floating estimate.
One certified length interval per clipped span scales all nonnegative rational
polynomial coefficients, with outward summation. Material-side cell
ownership assigns cell-aligned spans to the incident cell containing material;
an explicit all-incident option leaves assignment to the caller. These area
moments do not supply directional airflow normals or decide which spans are
exposed. Covered material and tile interfaces must be removed first.
Retained thermal assembly promotes a mixed exact/native partition to exact
geometry once. It removes opposite-plane coverage from horizontal faces and
subtracts opposed boundary spans within each common Z interval from vertical
faces. Single-span handles supply bounded outward-normal metadata for directional
heat transfer without rounding endpoints. Cell-face ownership prevents duplicate
Robin area on grid boundaries. Added-volume assembly intentionally supplies no
exposure certificate; a deposited increment alone cannot determine which faces
of the whole body remain exposed. Streamed skin traversal has its own geometry
requirements and does not acquire exact support from this retained adapter.

Horizontal polygon faces and vertical exposed edge spans have separate area
moments for Robin heat exchange. The consumer must remove covered faces and
internal tile/material interfaces before applying convection. Reservoir
temperatures and heat-transfer coefficients are explicit inputs. Surface
integration represents the supplied piecewise-prismatic geometry; it does not
establish convergence of surface area for a curved bead represented by a
staircase.

The adaptive thermal-history solver accepts these eight-node operators and
assembles their action across shared active node IDs without constructing a
dense global matrix. Internal conduction must be symmetric positive
semidefinite with zero row sums and zero reservoir load. Boundary operators
carry their own load and heat-accounting channel. Inputs with duplicate or
inactive node IDs, invalid operators or exceeded budgets publish no field.
Iterative convergence requires both the global relative residual and every
equation's residual scaled by its local load and diagonal-temperature term.
The true residual is recomputed before accepting a step. A large heat capacity
elsewhere therefore cannot hide an unresolved temperature in a tiny cut cell.
Newly deposited nodal mass adds its heat and capacity to existing nodes using a
capacity-weighted temperature; it never resets the previous node temperature.
The nodal history uses row-sum mass lumping, not the consistent capacity matrix,
so spatial refinement must also assess this projection. Covered boundary
operators are removed when an exposure event replaces the current boundary set.
An explicit `LumpedDiagonal` quadrature option also lumps transverse gradient
weights into twelve nonnegative, axis-directed conductances per background cell.
The weights come from the consistent operator's action on each coordinate field;
the resulting operator must preserve that affine action. Only grid-aligned
diagonal conductivity is supported, and negative edge weights are rejected, not
clamped. Robin row-sum lumping preserves the reservoir load and constant-field
heat exchange. Occupied volume, connectivity, material properties and boundary
temperatures are unchanged, but mixed spatial modes differ: spatial refinement
is still required. A three-dimensional insulated cosine-mode regression checks
second-order spatial consistency on a fully occupied cube with diagonal
anisotropic conductivity. That check does not establish convergence for cut
cells, material interfaces or the coupled mechanical response. Assemblies and
G-code results report the selected quadrature;
the default remains consistent integration. This numerical option does not
qualify a material model or establish a bed-adhesion margin.
An explicit coalescence event can join previously separate thermal nodes at
physical contact. It supplies a dense old-to-new map before adding the new mass;
each old node contributes all of its volume, capacity, heat and source power to
one new node. Existing conduction operators are reindexed, internal self-links
are removed, and collapsing distinct corners of an element is rejected. The
incoming exposure uses the new numbering. This permits independent histories
before contact without introducing a fictitious high-conductance penalty link.
The geometry producer must certify coincident physical nodes at that event;
membership in the final printed component alone is insufficient. Coalescence
work is separately bounded, and observations use the latest node numbering.
Adaptive error control and timed observations use the same history engine as
the voxel network. Element-step work is budgeted separately from node-step work.

`assemble_section_thermal_volume` connects a retained `SectionTopology` to this
solver on a regular native-coordinate grid. It clips occupied polygons into
cells, accumulates moments separately by material and physical component, and
assigns shared node identities only within the certified component. Nodal heat
capacity and initial heat retain each material's contribution; disconnected
pieces in one background cell do not exchange heat through a shared node.
The assembly also constructs exposed horizontal-face moments by subtracting
opposing slab footprints. This excludes internal layer/material contacts while
preserving ledges, unsupported lower faces, holes and bed contact. Face groups,
input/clipped points, nodes, occupied cells and integration work are bounded.
Vertical exposure is assembled separately by subtracting collinear contact
intervals in native edge coordinates and Z. Partial-height neighbors cover only
their shared height; polygon holes retain their exposed walls. Parametric cell
clipping preserves oblique edge lengths without rounding intersection vertices.
A boundary on a grid plane belongs only to its material-side cell, so grid cuts
do not create artificial convection surfaces. Area moments and outward normals
are retained for explicit airflow boundary assignment.

This is a fixed-domain initial condition, not inferred print activation: all
supplied geometry enters at time zero at explicitly supplied material
temperatures. No boundary coefficient is assigned automatically. Consumers must
use both horizontal and vertical exposed-face lists when assigning whole-part
airflow; omitted faces must not silently become insulated. Evolving exposure
and batched deposition use the streamed history adapter described below;
actual machine timing is not inferred by this fixed-domain assembly.
This path is not yet enabled for full-case Auto mitigation;
analytical assembled tests do not qualify ValveBoard adhesion.

`assemble_local_deposition_thermal_volume` assembles nominal rounded and
variable-height deposition packets without retaining the complete polygon
history. It borrows each packet from the local sampler, integrates moments
under provisional component labels, and merges those moments only after the
same topology stream completes and supplies its final connectivity mapping.
Physical cell and material-group caps apply after that merge; provisional
groups, cumulative processed points, integration work and topology frontier
storage have independent bounds. Cancellation or incomplete connectivity
publishes no assembly. Disconnected strands retain separate thermal nodes even
when they occupy the same background cell.

The optional `include_exposure` setting also constructs complete skin through
`SectionThermalSkin`. Horizontal surfaces subtract the opposing material
partition in the adjacent Z band, retaining exposed undersides and gaps.
Bounds are computed once per band to skip disjoint horizontal footprints;
exactly equal adjacent footprints need no Boolean subtraction, even across
different materials. Neither check introduces a geometric tolerance.
Vertical material contacts are removed within the current band using one shared
exact supporting-line index. Owner and edge/Z intervals exclude self contacts
and disjoint spans; single-piece bands skip the index. The index is bounded to
250,000 edges independently of the per-piece patch budget. Bands with at most
16 pieces first exclude islands whose bounds touch no other piece, avoiding
index construction for isolated contours without rounding their boundaries.
At tile seams,
native edge-coordinate/Z rectangles are matched against the neighboring tile;
only uncovered remnants become external surfaces. A lexicographic frontier
releases missing-neighbor traces when traversal proves there is no neighbor.
Consecutive identical traces merge only for the same material/component owner,
never across gaps. Collinear subdivisions also merge when they have exactly the
same Z interval and owner; this does not infer connectivity from edge contact.
Both sorting passes consume the bounded work allowance. Compaction also runs before a full live frontier rejects a
new trace, rather than retaining redundant adjacent bands until tile completion.
At a full frontier, the source topology builder may supply its current
union-find representatives, in the original label namespace. This bounded,
work-charged snapshot merges only connections already proved by that builder;
it is not a connectivity certificate. Retained traces can then consolidate
under those representatives. Final source-matched certification is still
required, and failed or cyclic mappings invalidate the complete skin stream.
Pending neighboring traces are compacted before seam subtraction uses scratch
fragments. A single opposing rectangle that covers a whole trace cancels it
directly, and disjoint traces emit directly; neither case needs a scratch copy.
Partial overlap uses a depth-first fragment stack capped at the smaller of 4096
entries and the configured frontier count, independently of the retained frontier.
Finished fragments emit immediately instead of retaining
the entire subtraction result. Each fragment resumes after its generating cover;
the original subtraction and traversal order is preserved. Stack exhaustion
rejects the stream, never truncates it.
A bounded spatial tree indexes edge coordinate, height, plane,
and axis to skip unrelated rectangles; candidates retain their original subtraction
order. Native coordinates are exactly representable in the index. Index construction
and queries consume the same bounded work allowance, with a separate 64 MiB
construction-workspace ceiling per index.
Cancellation, invalid ordering, reentry or exhausted bounds
poison the surface stream; consumers discard its provisional callbacks unless
the same topology stream completes successfully.
Work-limit diagnostics report cumulative charges for each surface stage, so
the operation that encounters the limit is not mistaken for its dominant cost.
The default surface-work budget remains 20 million operations. Explicit offline
callers may request up to 2 billion cumulative operations independently of
the live geometry and seam-storage limits; cancellation remains active.

Surface area moments are mapped through the final connectivity certificate and
must belong to an occupied material group. Repeated vertical patches with the
same cell, material, component and normal combine into one moment group. No
geometric voids are filled by this aggregation; boundary coefficients use that
cell/normal discretization and require independent refinement. The streaming
assembly sets `exposure_complete` only when both surface orientations complete.
Vertical face integration separates nine XY line moments from three analytic
Z integrals. XY moments are computed once per intersected grid column and reused
only for that same patch's Z cells. No historical geometry cache or coordinate
approximation is involved. Material/component ownership remains in the destination
face groups. Factored Z power differences avoid subtractive cancellation in thin
bands. XY integration and Z products consume separate portions of the work budget.
XY monomials use midpoint-centered closed forms: odd parameter powers vanish,
leaving only the second and fourth parameter moments. Independent Gaussian
quadrature checks all 27 moments, including tiny edges near cell corners.
Polygon area moments use signed fan triangles with exact native-coordinate
determinants, avoiding subtraction of cell-scale edge integrals for tiny
translated faces. A Duffy map and three-by-three Gaussian quadrature integrate
all nine monomials exactly up to floating-point error. The signed fan supports
concave rings; contour and hole orientation is evaluated exactly once per ring.
Holes retain their negative area contribution. Analytic triangles, translated
anisotropic holed rectangles, and tiny translated convection faces check the
moments and boundary operators.
Volume-band integration factors the Z power differences and obtains thickness
from the integer endpoint difference, avoiding cancellation in one-nanometre
bands near a cell's upper plane.
With angular exposure bins, a bounded 4096-entry accumulator adds the nine XY
moments of faces sharing the same cell, material, component, normal and exact
Z interval. Each batch is extruded once into 27 moments before accumulation
in the final face group. A full buffer flushes; it does not drop faces or
change boundaries. Exact-normal mode integrates directly. Both XY additions
and flushed Z products consume work, and final group limits remain enforced.
The streaming thermal integrator defaults to20 million work units and allows
explicit offline studies up to2 billion cumulative units. This CPU allowance
does not change point, group, cell or node storage limits; cancellation remains
active throughout geometry integration and assembly.
Without `include_exposure`, the flag remains false, as it does for incremental
additions. Empty face lists must not silently mean whole-part insulation.
Comparisons against retained fixtures check area, heated-bed/ambient exchange,
nodal temperatures and stored heat, not actual printer airflow or adhesion.
These checks do not enable full-case Auto mitigation.

`qualify_section_growth` compares successive qualified prefix snapshots using
native-coordinate polygon differences split at overlapping Z boundaries. Every
previous occupied region must remain covered by the same material in one next
component; moved, removed, relabeled or split old material is rejected. It
returns the old-to-new component map and only the newly occupied prisms. Equal
volume alone is never sufficient evidence of monotone growth. Point, Boolean
work and retained-piece budgets bound this comparison.

`SectionThermalSequence` converts these retained snapshots into thermal events
at caller-supplied physical times. Constitutive thermal properties and the grid
remain fixed; explicit deposition temperatures apply only to new material.
Existing nodes are coalesced through the certified component map at identical
native grid positions, then new nodal mass and conduction are added. Exposure
is rebuilt from the complete next geometry, not from the new subset, so covered
faces stop receiving external exchange. The exposure callback is mandatory:
an empty result explicitly chooses insulation instead of inventing airflow.
The final assembly is reindexed to match the history and can feed the mechanical
load adapter after successful evolution. Rejected transitions do not commit
geometry, time or event count. This adapter retains two bounded snapshots during
comparison; it is not the full-case tile-streaming geometry/timing adapter.
Accumulated conduction elements and time integration retain the history budgets.

`LocalSectionCursor` is the resumable form of the same sampler used by retained
and packet traversal. It builds the segment-to-tile index once, advances one
lexicographic tile at a time, and maintains cumulative geometry/work budgets.
It borrows immutable source paths and partition envelopes. Exceptions or
reentry poison the cursor; partial metrics are diagnostics, not certification.
The retained wrapper preserves its smaller retained point/index ceilings.
Index candidates come from each bead's bounding rectangle, but an exact
segment-versus-expanded-box separating-axis test removes tiles the bead cannot
reach. Tangencies remain candidates; radius conversion and native coordinate
rounding are covered conservatively. This avoids filling a diagonal bead's
entire bounding rectangle with empty index entries. Candidate tests remain
charged even when rejected. The default tile cap is 4096; explicit offline
studies may select up to 16384 while retaining independent reference, point
and work ceilings.

`visit_local_section_growth` advances two such cursors in sorted tile order
and certifies each snapshot with its own incremental topology builder. Only
the two current tiles' polygons are retained. Native polygon subtraction over
overlapping Z bands proves material-preserving coverage and emits the added
subset. A Z-band index avoids scanning the whole tile history for each prism.
Old-to-new component links are resolved only after both topology streams finish;
provisional labels that later join do not cause false split failures. Missing,
moved or relabeled old geometry rejects even when aggregate volumes match.
Callbacks remain provisional until successful return, and separate bounds
cover tile points, output points/pieces, links and comparison work.

`assemble_local_thermal_growth` uses that same paired traversal to assemble the
complete next field/skin and only the added material. Both use the next stream's
final component mapping; an unchanged snapshot has an empty added assembly.
`SectionThermalStreamSequence` uses the same conservative event builder as the
retained reference. It maps old nodes at identical spatial positions through
certified joins, adds only new capacity and heat, and replaces exposure with
the complete next skin. It stores path inputs and the assembled field, not
polygon snapshots. Replayed ownership providers must capture immutable source
data. Failed updates leave the previous time, source and assembly unchanged.
Changing sampling planes can alter the represented old solid; the growth check
rejects that rather than accepting a volume-only approximation. Sampling and
activation schedules therefore still require independent convergence checks.

`simulate_section_gcode_history` is a read-only retained-prefix reference adapter
from processor G-code to this sequence and its final mechanical field. The
caller explicitly chooses rectangular nominal planar beads, acknowledges modeled
processor timing, maps processor filament IDs to constitutive materials, and
supplies deposited-polymer temperatures and exposure operators. Tower, custom
and unclassified paths are excluded; mixed roles, unresolved extrusion without
an admitted stationary allowance, unsupported shapes, material overlap and
non-monotone geometry reject. No source G-code is rewritten.

The strictly increasing activation schedule ends at the modeled print end.
Each interval's new occupied volume enters at its endpoint, so this is batched
activation, not a claim of continuously resolved motion or deposition timing.
Activation spacing must be refined independently of the thermal integrator's
time step and the geometric and mechanical discretizations. Exposure is held
constant between scheduled endpoints and through the optional cooldown; its
provider must include relevant boundary-change times in that schedule. Prefix
traversal work is cumulatively bounded in addition to per-snapshot geometry,
assembly, history and mechanics budgets. Failed thermal histories return no
mechanical field.

The optional `stream_tile_size_m` selects streamed prefix activation and final
mechanics on the same regular grid. It admits rounded/scarf sampling, accepts
explicit native ownership by processor time and sampled Z, and feeds the
successful evolved nodal field into the streamed mechanical solver. Thermal
assembly, skin, paired growth and mechanical integration retain independent
explicit budgets. Retained rectangular history remains the reference path.
`use_native_ownership` prepares matching live-generation evidence separately
for each activation prefix. Old/new resolvers retain their immutable source
geometry during growth comparison; final mechanics reuses the final resolver.
Imported files without matching provenance reject this option. Custom ownership
and automatic native preparation are mutually exclusive. Strict ownership is
the default; explicitly requested interface-sensitivity variants use fresh
sampler-bound partitions, not reusable masks or silent union-loss tolerances.
Parsed-source tests compare retained/streamed heat balance and exercise rounded
geometry through nonzero cooling and elastic response, tile-size invariance,
and adaptive integration. This is not full-size
ValveBoard resource/convergence qualification, an inferred fan coefficient,
or validation of actual adhesion and material laws; Auto is not promoted by
these reference tests.

Native ownership is not necessarily monotone over deposition time. A lone
nominal bead may initially occupy an interface region that a later foreign
bead causes the native partition to assign to another material. Such a history
is rejected: transferring its volume to the new material would discard or
reinterpret previously deposited mass and heat. Batching that contact into a
single activation avoids the transition numerically but is not temporal
convergence and must not be used as a qualification shortcut. A consistent
time-dependent deposition/interface model is required for those contacts.

`activate_reference_sections` provides a bounded retained reference model for
investigating that issue. It intersects cumulative deposited source regions
with the same material's fixed reference region, splits at every source and
reference Z boundary, and qualifies the resulting topology. Final reference
geometry alone never activates a region. Same-material repeated source regions
are unioned before activation and excluded-volume accounting; foreign nominal
overlap cannot retroactively replace previously active material. The caller
must keep the reference fixed, prove monotone growth between results, and prove
complete reference coverage at the final event. Work, live/output points and
partition counts are bounded, and cancellation returns no partial geometry.

Independent integer intersections can round a junction on only one side of a
material interface. When this produces conflicting active regions in a
two-material band, the clipping kernel reconstructs both intersections in
exact arithmetic and rounds their shared vertices and T-junctions together.
It certifies component/hole correspondence and disjointness after rounding;
it does not remove the conflict by assigning a material priority or activating
the entire reference. This bounded fallback shares the interface partition's
rounding implementation and charges its work to the activation budget.
Activation tries nine global axis-bias pairs from zero and plus/minus one eighth
nanometre, nearest rounding first, then four one-axis quarter-nanometre biases.
Every trial clears
its shared vertex map, uses the same mapping for both materials, and must prove
exact component/hole correspondence and disjointness. An exact squared-distance
check retains the one-nanometre vertex bound; native integer vertices do not
move. All attempts share the original four-million-work ceiling. Geometric
candidate failures may try the next phase; cancellation and resource failures
propagate immediately. Native interface partitioning keeps its original nearest
rounding policy.
Material-to-component assignment is checked before publishing contact records
or pieces for a band. If the integer occupied union excludes input material,
one bounded exact union of the original material chains is attempted. Its
output is accepted only when every boundary coordinate is exactly an integer;
fractional output is rejected, not snapped. Assignment and the usual component
checks are then repeated against that reconstructed union. Original material
pieces are preserved, and the retry shares the existing work and live-point
budgets.
The same reconstruction handles an independently clipped island that vanishes
under its contact-search uncertainty band, including single-material bands.
It recomputes from the original operands rather than deleting small polygons.
If the exact reconstruction retains that positive component, the ambiguous
topology remains rejected. Already-valid single source polygons bypass an
unnecessary union operation.
The additional integer representation moves vertices by at most one nanometre.
The slab's contact-search allowance is not a cumulative physical error bound.
Changing activation boundaries can still change rounded junction locations,
so independent snapshots must pass the same monotone-growth check; shared
rounding alone is not a certificate of chronological consistency.

This is an explicit fixed-interface approximation, not conservation of the
full commanded extrusion volume or a model of molten-bead displacement. It
reports nominal source volume excluded by the reference separately. Thermal
activation conserves heat and capacity within this chosen reference geometry;
it must not attribute the excluded nominal volume's heat to another material.
This retained helper requires independently qualified source topologies;
the streamed counterpart supplies the optional G-code path described below.
Neither is enabled in Auto. The reference tests cover
staged foreign-material activation, unchanged snapshots and heat accounting;
they do not qualify actual ASA/TPU interface deformation or adhesion.

`visit_reference_activation` applies the same clipping kernel to streamed
nominal paths. It advances the fixed reference and one sampler per deposited
material in common tile order. Resumable plane sweeps retain only a coalesced
height band and its pending successor per sampler. Reference activation merges
their height boundaries and clips one common band at a time, including
source-only bands whose excluded volume must still be counted. It does not
retain a tile's full vertical polygon history. Concurrent
samplers share tile and reference capacity by actual index occupancy, so sparse
materials do not reserve unused capacity. Index-work, section and segment-visit
allowances are divided between the concurrent samplers;
clipping work and processed/output points are cumulative across the traversal.
The reference and active geometry have separate completed connectivity
certificates. Sources do not receive artificial bed attachment. Piece and
packet callbacks are provisional until both streams finish successfully.

An explicit complete-deposited-reference option avoids resampling nominal source
unions when the entire deposited vector exactly matches the reference, including
source indices and endpoint dimensions. It freezes the bounded matching input
before callbacks and preserves the shared endpoint/partition-plane schedule.
Every forwarded band must carry an unchanged complete first-deposition partition
proof with zero boundary uncertainty. Its identical reference and active geometry
share one connectivity builder and its successfully completed certificate; packet
consumers still run, and copying the bounded final label mapping is charged.
Partial prefixes and mismatched inputs retain separate builders and ordinary
reference clipping. Excluded nominal volume remains explicit: the private proof
stores the sum of exact per-material differences between original and owned
unions, integrated over each band. This counts triple overlap twice and repeated
same-material coverage once, without cancellation from subtracting rounded areas.
Native and custom unproven packets cannot use this option.

First-deposition reference histories explicitly select exact per-command
composition for their activation source samplers, in both thermal and mechanical
assembly. Each source uses the same normalized, quantized individual cuts as the
fixed reference and unions them exactly before reference clipping; it does not
round cross-command intersections through a native polygon union. This
single-material geometry option is default-off, requires retained exact sections,
and grants no ownership or qualification certificate. Native reference sampling
is unchanged. These histories also explicitly select occupied-interior face
decomposition for exact topology. Material and occupied unions retain separate
point-touching faces and preserve tangent void boundaries without filling them.
Strictly separated polygon components are reconstructed locally after explicit
constructor-precondition checks and an independent occupied-face count. Weak or
tangent rings retain arrangement-face selection instead; the optimization does
not redefine connectivity. Both paths share the same work and retained-vertex
limits.
Positive-area vertical contacts and positive-length tile-face contacts remain
required; a point contact never supplies support. The default exact policy still
rejects nonmanifold boundary contacts, and native topology is unchanged.

`assemble_reference_thermal_volume` consumes that stream to assemble active
volume and complete exposed skin, including removal of tile seams and internal
material contacts. Thermal nodes are merged using the active stream's final
component mapping, and excluded nominal volume remains separate from deposited
heat. It is a snapshot assembly, not a temporal growth certificate.
Vertical contact extraction indexes exact native supporting lines rather than
comparing every pair of boundary edges. Collinear subdivisions share a key;
parallel separated edges do not, including one-nanometre air gaps.

Streamed thermal assembly optionally uses uniform angular quadrature for
vertical surface normals (`vertical_normal_bins`). Zero retains exact normals;
nonzero resolutions are multiples of four from 8 to 4096. Surface area moments,
material/component ownership and the occupied domain are unchanged. Only the
normal supplied to directional heat-transfer laws is rounded to the nearest
direction, with angular error at most pi/N. This preserves isotropic heat
exchange; a positive projected-airflow law has an integrated area error bounded
by total side area times 2 sin(pi/(2N)). General directional boundary laws require
their own angular convergence study. The resolution is retained in assembly and
G-code results, not silently inferred from a resource failure.
Tests compare streamed material volumes, heat and exposed area against
the retained reference and check partial activation from parsed G-code.

`ReferenceActivationCursor` makes the same traversal resumable at tile
boundaries. It owns per-material source paths, borrows an immutable reference,
preserves budgets across advances, and rejects incomplete finish or reentry.
`visit_reference_section_growth` pairs two such cursors and uses the ordinary
same-material coverage and component-link proof on their active packets. It
does not infer growth from aggregate volumes. `assemble_reference_thermal_growth`
integrates the complete next field and skin plus only that proven increment.
The thermal sequence optionally owns a fixed reference and immutable ownership
provider; failed appends leave its prior field, time and excluded-volume
diagnostic unchanged. Per-event ownership replacement is prohibited in this
mode.

`SectionGCodeScenario::use_fixed_interface_reference` explicitly selects this
analysis approximation. The adapter prepares the final reference once, fixes
native ownership at that timeline endpoint, charges its source traversal, and
activates each scheduled prefix against it. Thermal history uses the same
conservative event builder as the ordinary path. Mechanical recovery consumes
the identical activated reference geometry and evolved field through
`solve_reference_deposition_sections`, not the unmodified nominal bead union.
The result labels this model and reports final excluded nominal volume.
Existing behavior is unchanged when the option is disabled. Staged native
interface and parsed cooling tests exercise the full path; they do not qualify
its physical interface approximation, full-case resource/convergence behavior,
or actual bed adhesion. It remains outside automatic mitigation selection.

`section_thermal_load` transfers a successfully completed assembled thermal
history into the mechanical load callback. It retains eight-corner temperature
variation per occupied cell/material/component, converts it through the explicit
expansion tensor and stress-free temperature, and owns the resulting field.
Missing coverage, failed histories, different grid bounds or occupied volumes
are rejected. The consumer must use the same qualified topology: component IDs
are local to that snapshot and matching volume is not a geometry fingerprint.

The mechanical callback accepts either world-space affine strain or cell-local
trilinear corner strain. Trilinear strain is converted to eight polynomial
coefficients and integrated with the same occupied degree-two moments used for
stiffness. Initial strain energy, equivalent force, stress bounds and sampled
stress all retain cross terms; no affine temperature fit is introduced. General
trilinear loads bypass the free-body affine/quadratic analytic shortcuts and use
the equilibrium solve. Conversion is exact only when the expansion tensor is
constant over the cell's temperature range; an expansion-law breakpoint inside
that range requires subdivision and is rejected. A stress-free reference in a
different interval remains supported through integration of the expansion law.
These properties verify field transfer, not the unknown physical parameters or
the adequacy of a particular spatial or temporal resolution.

Transient heat histories retain signed heat exchange by boundary kind: ambient,
bed, part fan, auxiliary fan, dock fan, exhaust, and mixed exchange. Adaptive
integration commits only accepted half-step energies; rejected trials do not
contribute to this ledger. Concurrent surface reservoirs share one solid
half-cell conduction resistance and are reported as mixed exchange rather than
misattributed to ambient convection. These are heat-transfer accounting channels,
not causal stress contributions or fan-PWM calibration. Temporal peeling
comparisons require each channel to stabilize independently, so cancellation
between bed heating and airflow cooling cannot conceal a net-energy error.

`analyze_configured_gcode` constructs four no-relaxation scenarios from the
actual deposited paths: the published ASA CTE endpoints, each with either
uniform chamber temperature or a linear bed-to-chamber temperature gradient
over each face-connected component's height. A taller disconnected object
therefore does not change a shorter object's thermal load. The stress-free temperature is capped by Tg and the
configured nozzle temperature. These prescribed temperatures are scenario
inputs, not a transient heat-transfer solution. Chamber and bed targets are not
measurements of deposited-polymer temperature. The MPa preview is the maximum
principal tensile stress across these scenarios, not a failure probability.

ASA CTE, Tg and room-temperature modulus come from the representative material
table documented in [warp mitigation](warp-mitigation.md). Poisson ratio 0.35 is
the explicit modelling assumption in the
[ASA constitutive study](https://www.mdpi.com/1996-1944/19/3/470), not a measured
value for an arbitrary spool. The no-relaxation, room-temperature-modulus
scenario is deliberately stiff, but is not a guaranteed worst-case envelope.
Other materials, filled grades and mixed ASA/TPU paths are unavailable rather
than assigned these constants. Interfacial adhesion and TPU's nonlinear,
near-incompressible response require different constitutive data and treatment.

For explicit studies, callers may supply a `MaterialScenario` for each filament
slot, including its exact configured type, provenance, Young's modulus, Poisson
ratio, CTE endpoints, stress-free temperature and a chamber-temperature validity
limit. Missing or mismatched slots remain unsupported; supplied invalid
properties are rejected. This enables heterogeneous linear-elastic scenarios,
not an automatic TPU material model. It does not supply nonlinear behavior,
temperature-dependent stiffness, interfacial debonding or protection against
near-incompressible locking. The four cases pair all lower CTE endpoints and all
upper endpoints with the two temperature fields; they are correlated scenarios,
not an exhaustive combination or guaranteed mixed-material worst case.
Candidate comparisons require identical material, provenance and temperature
scenario identities; a changed constitutive assumption cannot masquerade as an
improved toolpath.

Every scenario must converge before any stress samples are published. Missing
material data, unset temperatures, chamber temperatures at or above the Tg
range, invalid geometry, resource limits and lack
of convergence have distinct result states. Cancellation propagates, with old
stress samples already invalidated. Default execution allows 12,000 mesh nodes
and 250,000 deposited segments; callers cannot raise the node limit above
60,000 through this entry point. It writes no mesh, cache or G-code copy to
disk. Exceeding a limit must not cause automatic coarse voxelization that
closes infill gaps or artificial bonding between components.
Raster traversal also has a two-million-visit work cap, including empty and
already occupied cells, so overlapping paths cannot create unbounded work
without exceeding the occupancy limit. A refinement stopped by its resource
budget is distinguished from one that uses its allowed refinements without
converging. ASA color IDs are canonicalized only when their material type and
all scenario temperature inputs match; incompatible overlapping cells still
fail instead of being averaged.

A full case at bead-resolved spacing can exceed this budget. A scalable
continuum/shell alternative would need section-wise homogenization validated
against resolved paths, geometry and discretization convergence, material
interfaces and a temperature history. Neither the coarse raster nor a
converged equilibrium for a smaller unrelated specimen establishes the
case's resistance to corner lift.

## Bounded periodic stiffness homogenization

`homogenize_periodic` computes a small representative cell's effective 6×6
elastic stiffness from its conforming tetrahedral solid geometry, including
explicit empty space. Opposing boundary nodes and triangle traces must match
exactly. It rejects coincident nodes instead of guessing a weld. A cell with
unrestrained islands or a macroscopic mechanism is unavailable, not stabilized
by artificial supports. Input eigenstrains must be zero for this stiffness-only
calculation.

The displacement is decomposed into an imposed affine macroscopic strain and
periodic fluctuations. `solve_affine_periodic` eliminates shared fluctuation
degrees of freedom exactly, assembling `Tᵀ K T` and `Tᵀ(f − K u_affine)`; it fixes
one fluctuation class's translation and introduces no penalty stiffness. Six
engineering-strain cases give volume-averaged stresses divided by the full
cell volume, including its voids. Microscopic strain energy must agree with
macroscopic work (Hill–Mandel consistency). Effective stiffness must be symmetric
within numerical tolerance and positive definite. This follows the
[periodic elasticity formulation](https://bleyerj.github.io/comet-fenicsx/tours/homogenization/periodic_elasticity/periodic_elasticity.html).

Default cells are limited to 1,024 nodes; the absolute limit is 2,048 nodes and
32,768 tetrahedra, with a conservative dense-fill/triplet workspace preflight
of 512 MiB. Each of the six solves is discarded after accumulating its column.
Cancellation occurs between cases and in assembly/recovery, but not inside a
sparse factorization. No mesh or large intermediate is written to disk.

Tests recover homogeneous isotropic elasticity, series/parallel laminate
stiffness, porous-cell volume and energy consistency. Separate perforated-body
tests enforce exact free thermal contraction and covariance under joint
rotation of geometry and restraint. These checks do not validate any printed
material constants, weld areas, nonlinear TPU behavior, transient thermal
history or a whole-case reduced model. In particular, periodic stiffness does
not manufacture a directional thermal expansion coefficient from raster angle;
uniform isotropic eigenstrain remains a stress-free affine contraction of a
single-material free skeleton.

`homogenize_periodic_eigenstrain` closes an explicitly supplied per-element
eigenstrain load on the same periodic geometry. It snapshots bounded inputs
before callbacks, computes stiffness with zero loads, and performs a seventh
solve at zero macroscopic strain with the supplied loads. The volume-averaged
stress gives the effective stress offset; the recovered energy is normalized
by total cell volume, including voids. For engineering strain `E`, the reduced
energy density is `0.5 E^T C E + E^T stress_offset + energy_offset`.
The effective free strain is `-C^-1 stress_offset`. The remaining energy at
that strain is retained: incompatible constituent strains can leave internal
stress even when average stress vanishes. Only subtraction-level roundoff may
be clamped to zero; inconsistent equilibrium or negative residual energy is
rejected. This closure supplies neither a temperature law nor a certificate
that a region of a printed part is periodic. Wall, interface, bed-contact and
transient-history applicability still require separate treatment.

`solve_hybrid_periodic` applies these closures to explicitly selected elements
of a conforming macro mesh. It preserves the supplied nodes, connectivity,
restraints and all unassigned materials. Closure consistency and bounded input
sizes are checked before solving; callbacks cannot mutate the frozen inputs.
The returned macro equilibrium excludes irreducible microscopic energy, which
is reported separately as the sum of residual energy density times assigned
tetrahedral volume. Their sum is the total energy. This adapter does not infer
periodic regions or add bonds: shared nodes already imply perfect bonding, and
the caller must establish which interior regions admit the approximation.

## Explicit anisotropic elasticity

Elements may supply a symmetric positive-definite 6-by-6 elasticity matrix in
place of Young's modulus and Poisson's ratio. The legacy isotropic path is
unchanged when this optional matrix is absent. Matrix ordering is xx, yy, zz,
xy, yz, xz, mapping engineering strain to physical stress. `rotate_elasticity`
rotates the tensor using a proper orthogonal local-to-world orientation, with
the engineering-shear factors included. Eigenstrain remains an independently
prescribed world-frame quantity; rotating stiffness does not invent anisotropic
thermal expansion. Voxel construction retains its explicit isotropic inputs;
callers supplying homogenized matrices must assign them to the resulting
elements deliberately.

## Thermal strain coupling

The configured-analysis airflow summary retains decoded events only when the
entire supplied stream is read successfully. A byte, line, or event budget
failure retains diagnostic counts but publishes no partial timeline. Event
line IDs start at the supplied stream position; callers must establish that
this was the complete matching G-code before attaching processor timing.
Retention alone does not incorporate airflow into the static analysis.

`ThermalHistory::time_airflow_events` aligns decoded fan states with
modeled G-code line completion times. It requires complete processor timing,
includes travel and attributed waits, and excludes display-only vertices.
Events sharing a line retain their order and follow all motion records on that
line, including arc subdivisions. The bounded, cancellable adapter preserves
unknown fan states and automatic-controller flags. The caller must establish
that events and move timing come from the same final G-code; coincident filenames
or line counts are not source authentication. These timestamps are not measured
firmware timing, and PWM is not converted into a heat-transfer coefficient.
Spatial convection and material assumptions remain explicit scenario inputs.

`apply_timed_airflow` attaches this trace to a caller-supplied spatial boundary
law. Fan changes merge into the separate exposure schedule, leaving deposition
activations and material properties untouched. Equal-time commands take effect
in their original order; the last state remains through cooldown. Before the
first command the state is unknown. Invalid schedules and cancellation leave
the scenario unchanged. Matching total duration is a consistency check, not
proof that the trace came from the same G-code. This adapter supplies no
PWM-to-convection calibration or automatic material assumptions.

`compare_section_cooling` runs two labeled boundary-law scenarios sequentially
on the same source and common section-history inputs. Only the exposure law
and its change times differ. It checks that occupied and omitted material,
deposition heat, activation counts and modeled duration agree before returning
compact final-time energy and tensile-demand diagnostics. Common callbacks must
be deterministic; matching scalar inventories is not a geometry fingerprint.
Each run retains its own work limits and only one full field is held at a time.
Failures publish no partial pair. This comparison does not establish peak print
stress, numerical refinement, an Auto ranking, or an adhesion safety margin.

`study_gcode_deposition` independently compares time-step refinement at each
spatial grid, then requires two successive spatial refinements to stabilize
integrated tensile demand, elastic energy, deposited/stored/boundary heat and
bed-contact area. The same explicit physical scenario and observation time are
used throughout. Positive local thermal error control is mandatory; a successful
linear solve alone does not qualify the field. Absolute tolerances carry their
physical units, alongside the declared relative tolerance. Budget exhaustion,
missing full mechanical output or nonconvergence publishes no qualified field.
At most four grids and two time-integration settings per grid are evaluated,
retaining scalar comparison metrics instead of all intermediate meshes.
The finer time-step cap is at most half the largest interval actually accepted
by the coarse run, since deposition events can make nominal caps inactive.

This numerical study concerns one modeled observation. It does not establish
peak demand over the print, convergence of local crack-tip stress, actual
adhesion strength, or validity of the supplied material and airflow inputs.
Success records observed stability over the tested resolutions, not a rigorous
error bound; features narrower than the finest raster can remain unresolved.
Reduced mechanical results remain separate and cannot bypass this qualification
by being relabeled as full bed-contact output.

`study_gcode_history` applies that qualification independently to an explicit,
strictly increasing sequence of modeled observation times. It retains scalar
force and energy samples, not a growing collection of meshes. Separate indices
identify the largest sampled tensile demand and largest sampled elastic energy;
they need not occur at the same time. A final low-stress state does not erase an
earlier maximum. Every observation must qualify before any history or peak is
published. Failure or exhaustion of the shared simulation-run budget returns
diagnostic progress only. The maximum is 16 observations and 128 simulation runs.
These are sampled maxima, not an upper bound between samples: temporal sampling
of mechanical demand remains distinct from thermal time-integration convergence.
The caller still needs matched scenario inputs and history sampling before
using these demands to rank automatic candidates.

`ThermomechanicalSimulation` couples an occupied deposition grid to finite-volume
thermal history and a final mechanical snapshot. Every occupied cell has one
explicit activation time and initial polymer temperature. Duplicate coordinates
are rejected, so repeated paths do not duplicate deposited mass. New cells are
appended at activation; only active face neighbors exchange heat. The shared-face
conductance uses the two half-cell thermal resistances in series, plus an
explicit inter-material resistance. Voids and disconnected islands are retained.

Each activation or boundary-scenario event rebuilds the exposed faces, removing
surfaces buried by new deposition. Bed faces use an explicit bed-contact heat
exchange coefficient; other exposed faces use declared ambient, part-fan,
auxiliary, dock or exhaust exchanges. These coefficients are physical scenario
inputs, not inferred from PWM. Exchanges optionally select signed axis-normal
faces and an inclusive finite world-space box. A face is selected by its center;
this mesh-resolved approximation does not compute fractional coverage of faces
cut by the box. Defaults select all exposed non-bed faces. Explicit rear,
part-fan or dock zones can change with boundary events, without changing deposited
mass, occupied geometry or mechanical bed footprint. Selector evaluations count
toward the bounded schedule work. This exposed-face approximation does not
equate face and cell-center temperatures: conductance is
`area / (half_cell_length / conductivity + 1 / h)`, with zero `h` insulated.
The same series resistance applies at the bed. This approximation does not
apply the solid resistance independently to overlapping exchange mechanisms:
on each face, supplied coefficients are summed and reservoir temperatures are
weighted by those coefficients before adding the one half-cell resistance.
The thermal energy result is aggregate, not attribution to individual fans.
Callers must avoid double-counting overlapping empirical convection laws.
This approximation does not
solve airflow, occlusion or chamber circulation: spatial coefficients and
reservoir temperatures remain prescribed inputs, not a predicted fan field. Snapshot
eigenstrain is derived from final cell temperatures and explicit expansion data,
then assigned to six tetrahedra per active voxel. Optional anisotropic stiffness
is assigned in the same order. Supplied stiffness is already world-frame;
`local_to_world` rotates expansion only, so local stiffness must be explicitly
rotated with `rotate_elasticity` before assignment. Mechanical material contacts are perfectly bonded;
the thermal contact resistance does not imply a mechanical interface law.

Schedule entries and exposure reconstruction work are bounded, as are thermal
cells and mechanical nodes. Thermal failure returns no mechanical field;
unrestrained mechanical islands, unsupported expansion domains and oversized
meshes fail rather than being silently merged or coarsened. This coupling does
not yet evolve mechanical stress during deposition or relax thermal strain.
Its mechanical reference is the original undeformed grid, not stress-free new
layers deposited onto an already-deformed body. Rebuilding complete exposure at
each unique cell activation can produce quadratic exposure-reconstruction work;
the entry budget bounds that work. Events are streamed into the thermal solver,
so previous boundary snapshots do not accumulate in memory. Grouped
activation requires an explicitly justified temporal discretization, not silent
batching of toolpath times.
Large bead-resolved cases remain subject to the mechanical node limit. Activation
times and initial polymer temperatures must be supplied by an independently
supported deposition mapping; nozzle setpoints are not automatically equated
with deposited polymer temperatures.

`ThermoelasticThermalStrain` integrates caller-supplied piecewise-constant
symmetric expansion tensors over an explicit temperature domain. Its result is
world-frame engineering eigenstrain, including rotated shear components. The
stress-free temperature is an input, not inferred from nozzle temperature or
glass transition. Missing intervals and extrapolation are rejected. A synthetic
coupling regression passes temperatures from `ThermalHistory::step`, with
separate explicit part-fan and exhaust heat-exchange boundaries, into thermal
strain and then constrained anisotropic equilibrium. This verifies the data
path and signs; it neither calibrates PWM to convection nor represents a full
deposition history. Cooling-induced strain must not be multiplied by an
arbitrary relaxation factor: coupled viscoelastic equilibrium requires the
history-dependent tangent and internal-stress right-hand side.

## Isothermal scalar relaxation history

`ThermoelasticRelaxation` supplies an explicit generalized-Maxwell material-point
kernel. It is a constitutive building block, not a viscoelastic FE coupling or an
automatic TPU material assignment. The caller supplies instantaneous modulus,
reference temperature and up to 64 nonnegative Prony fractions with positive
relaxation times. Fractions sum to at most one; the remainder is the equilibrium
spring. Input strain excludes thermal/free strain.

Each update integrates a linear strain ramp exactly: a branch retains
`exp(-dt/tau)` of its previous stress and adds
`E_branch * (tau/dt) * (1-exp(-dt/tau)) * delta_strain`. Zero duration represents
an instantaneous strain step. Short steps use `expm1` to avoid cancellation.
Invalid updates leave the history unchanged. The supplied temperature must equal
the reference temperature; no unmeasured time-temperature shift is invented.

Analytical step, ramp, hold, unloading and time-partition tests verify the kernel.
A separate regression uses the eight-term Filaflex 82A spectrum measured at
20 degrees Celsius in [Polymers 2021, 13, 2365, Table 7](https://pmc.ncbi.nlm.nih.gov/articles/PMC8309542/).
It verifies the published normalized relaxation, not the properties of the
user's TPU, hot deposition, or bed adhesion. Coupled stress evolution still
requires appropriate tensor shear/bulk constitutive data, temperature histories
and assembly of the history-dependent FE equations.

## Prescribed contact loss in section mechanics

Section mechanics accepts explicit rectangular bed-release regions in world
coordinates. Their edges must align with the XY background cell faces. The
rectangles stay fixed when refining the basis; they remove bed restraint, not
material or stiffness. The union selects whole bed-cell traces. Their actual
occupied areas are integrated from the clipped polygons, including holes;
overlapping rectangles never double-count released area. Nodes shared with a
retained trace remain fixed, including the crack boundary.

This supplies a traction-free partial-detachment scenario at a prescribed
material/eigenstrain state. Comparing stored elastic energy against the bonded
solution at that same state gives finite-area release work, not a singular
corner stress. Dividing by released area gives an average driving demand, not
the unknown ASA/PEI interface resistance or a prediction of crack initiation.
Both the basis and release-region discretization require independent checks.
No automatic release propagation or Auto qualification follows from this API.

Previously attached components must retain restraint. A region that removes
all restraint from a component is rejected, even when analysis of pre-existing
free components is enabled. The original topology certificate is still checked.
Penetration of the occupied footprint rejects a traction-free solve. Release
studies retain bounded polygonal bed traces, not the full volume geometry:
200,000 trace points by default, with an explicit ceiling of two million. These
traces use the same clipped polygons as the area integrals. A bilinear Q1 scalar
field is harmonic and attains its extrema on the boundary; each polygon edge
has a quadratic restriction. Endpoints and interior stationary points therefore
give the exact footprint extrema up to floating-point arithmetic. Testing only
polygon vertices can miss an edge-interior extremum. Testing background cell
corners instead can incorrectly reject penetration outside deposited material.

The actual polygonal-footprint maximum lift is reported separately from the
larger cell-corner upper bound. Neither is a measured lift height or independently
converged displacement. Empty release lists preserve the fully bonded model.

An explicitly selected `SectionCellContact` study permits frictionless re-contact
in a conservative Q1 subspace. A bounded active-set solve constrains only normal
displacement at contacting bed-cell corners. Tangential displacement remains
free; positive compressive reactions are reported separately from bonded-bed
reactions. Complementarity is checked with displacement/force-scaled numerical
tolerances. Cycles, exhausted iteration or cumulative work budgets, nonfinite
reactions and unresolved equilibrium return no field. Contact adds neither
adhesive healing nor a penalty spring or artificial stiffness floor.

For cut cells, corner constraints also restrict the polynomial extension outside
occupied material. This is a more restrictive admissible displacement space
than contact only on the actual footprint, so it can overestimate stored energy
and underestimate release work. It is not a conservative prediction of safety
against peeling. The field is explicitly marked as a cell-contact bound and
requires independent basis/geometry refinement; its multipliers are generalized
constraint forces, not local pressure on the plate. Contact does not qualify
the field for Auto or establish the actual interface's fracture resistance.

### Spatially varying prescribed shrinkage

Section mechanics also accepts a spatial load provider on the certified occupied
cell/material groups. Queries include world bounds, occupied volume, material and
the current physical-component identity. Component IDs belong to this geometry
snapshot, not a persistent thermal timeline. The provider supplies either affine
world-space strain or trilinear corner strain for the requested cell; it replaces prescribed strain
but never stiffness, material ownership, voids or bed contact. Missing coverage
must fail rather than fall back to unloaded material. Responses are cached once
per group and used consistently in assembly, energy and stress recovery. Free-body
shortcuts compare the supplied affine spatial loads, not stale material defaults;
general trilinear loads use the equilibrium solve.
Directional probes retain the provider and propagate its exceptions.

The transient heat driver can expose requested observation times during one
streamed history. It splits integration at those times and applies all coincident
deposition/boundary events before observation. Active cells are borrowed only for
the synchronous callback; no complete sequence of thermal fields is retained.
Callbacks may drive section solves through the spatial load provider. All derived
snapshots remain provisional until the complete thermal history succeeds, and
each mechanical/thermal discretization still needs independent convergence checks.
An adaptive transient two-cell heat-flow benchmark checks this handoff against
the exact piecewise axial contraction, including energy and unchanged volume.

These interfaces do not construct the full-case section thermal network, its
activation sequence, heat-transfer coefficients or adhesive fracture law. They
do not enable Auto qualification or substitute a fitted gradient for unavailable
thermal coverage.

Section materials may prescribe an affine world-coordinate engineering strain,
`epsilon*(x) = epsilon0 + gradient*x`, with coordinates in metres. This represents
an explicit spatial loading scenario, including a through-thickness shrinkage
gradient; it does not infer a thermal history from bed/nozzle targets. A zero
gradient retains the uniform-strain model and its stress-free free-component
shortcut. Attached bodies and incompatible material interfaces use the Q1
equilibrium solve.

The strain-displacement polynomial is multilinear. Its product with this linear
strain field, and the quadratic eigenstrain energy, fit the existing degree-two
cut-polygon moments. The implementation therefore integrates these loads and
their energy without cell-centre sampling or new geometry approximation. Point
stress recovery subtracts the prescribed strain at the requested world position.
Refinement roundoff scales include integrated prescribed-strain energy and an
eigenstress bound over each occupied background cell, not just the constant term
at world zero. Gradient assembly is charged against the existing work budget.

Analytical checks include a vertical linear eigenstrain whose exact displacement
is quadratic: Q1 recovers its nodal values and its energy error is
`E*g^2*A*sum(dz^3)/24`. Graded Z cells, translated bed coordinates and floating
bodies exercise the coordinate and constraint handling. A released slender strip
with linear through-thickness axial shrinkage is checked against its natural
curvature. These checks verify numerical mechanics, not the actual ASA temperature
gradient, stress relaxation, or ASA-to-PEI fracture resistance.

`linear_thermal_eigenstrain` converts an explicit affine temperature field on a
bounded spatial domain to this strain representation. It integrates the supplied
temperature-dependent expansion law from the declared stress-free temperature
and rotates tensor shear into world-frame engineering strain. The reference
temperature may span expansion-law intervals, but the current spatial temperature
range must have constant expansion; otherwise domain subdivision is required.
There is no extrapolation or silent affine fit across property changes. A steady
thermal-slab benchmark feeds its solved temperature gradient through this
conversion and checks section curvature against the continuum solution. This is
a coupling check, not an implicit heat-history model for full-size section runs.
An attached/released coupon also checks linear lift and quadratic release-energy
scaling under proportional thermal loads, with unchanged geometry and contact
conditions. Neither benchmark supplies a fracture-resistance law.

A genuinely free connected component with the same linear prescribed strain
throughout admits an exact compatible quadratic displacement, even if its
elasticity changes between materials. The identity
`u_i,jk = epsilon_ij,k + epsilon_ik,j - epsilon_jk,i` constructs it; its elastic
strain is identically zero. These components bypass cut-cell factorization,
contribute zero elastic energy/reaction, and are sampled analytically rather
than interpolating their nodal values with Q1. They are counted separately from
uniform affine free components. A finite-difference strain-recovery test exercises
all spatial and engineering-shear coefficients. Any material-side difference in
the prescribed strain or gradient disables this analytic path and retains the
equilibrium solve with rigid-body reference coordinates. No artificial stiffness
or bed attachment is introduced to stabilize detached fragments.

For a detached component contained entirely in one background cell, Q1 basis
coordinates use the component's occupied bounding box. This is an invertible
axis-affine change of the same trilinear polynomial space, not a coarser model
or a removal of strain modes. Integrated moments accumulate near occupied
support before being rebased; thermal strain is evaluated at the new basis
corners and stress samples use the same coordinates. Bed-connected and
multi-cell components retain the shared background basis. Free-component nodes
therefore need not lie on the background lattice; they carry no bed reactions.
The same pivot and equilibrium checks still apply. This conditioning treatment
does not resolve topology ambiguity or qualify material/adhesion assumptions.

The sparse solve also checks row-wise relative backward error, in addition to
global and equilibrated residuals. Up to three residual-correction solves may
reuse the existing LDLT factor when the row-wise check fails. Corrections target
the original assembled equations through the same congruence scaling; they do
not add stiffness or relax the acceptance threshold. Work and factor-storage
ceilings remain enforced. Exhausted corrections or nonfinite arithmetic return
no field. The reported correction count is cumulative across contact iterations,
while the reported backward error belongs to the final equilibrium. This is
iterative refinement of a linear solve, distinct from mesh refinement and from
physical validation (see [LAPACK Working Note 165](https://www.netlib.org/lapack/lawnspdf/lawn165.pdf)).
