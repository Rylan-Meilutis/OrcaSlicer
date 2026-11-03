# G-code deposition geometry export

The **Export geometry** action is available in the print/export
button's action selector and in **File → Export**. Its save dialog selects STL,
STEP or OBJ. It reconstructs the current
plate's generated toolpaths, not its original CAD volumes. Preview filament
visibility and layer-slider filters do not restrict the export. Infill, walls,
supports, skirts, brims and the prime tower share one STL; disconnected printed
islands remain disconnected shells in that file.

The fast STL option streams separate closed beads into one binary STL. It keeps
the nominal inter-line gaps, but overlapping beads are not fused and material
identities are not retained. This avoids whole-print volumetric reconstruction;
receiving simulation software must handle the overlaps explicitly.

The fused STL option is a union of nominal rounded extrusion beads. The user chooses
0.1 mm or 0.05 mm spatial sampling; gaps smaller than that sampling may disappear.
There is no silent coarsening when a resource limit is reached. A bead is the
Minkowski sum of its three-dimensional segment, a horizontal disc of radius
`(width - height) / 2`, and a sphere of radius `height / 2`, shifted downward by
half its height. Sloped and scarf paths retain their Z motion; a bounded convex
distance minimization handles these paths without flattening them. Horizontal
paths use a direct projection. Invalid bead dimensions and insufficient sampling
resolution are rejected. Extruder-only commands cannot specify a spatial purge shape, and
positive-extrusion wiping or other nondeposition moves may have no bead dimensions.
These commands are not reconstructed; omission of moving extrusion commands or
stationary extrusion vertices requires confirmation. Ordinary retraction recovery
does not trigger that prompt.
Preview processing can infer or clamp bead dimensions and transform startup
coordinates. In particular, nominal width/height metadata can remain constant
while scarf extrusion flow varies: this reconstruction does not infer the
actual varying cross-section or conserve commanded volume in those regions.
The selected sampling spacing is not a guaranteed surface-error tolerance.
Reconstruction uses bounded tiles on one global sampling lattice. Overlapping
tile halos provide the same field at tile boundaries; whole mesh faces have a
single owning tile. Collinear zero-area mesh faces are resolved by exact edge
subdivision across the halo, retaining that ownership and checking opposite-oriented edge incidence before
removing those faces. Unresolved topology stops the export rather than dropping
faces and leaving holes. Exact coplanar horizontal interior faces are reduced without
moving curved surfaces or tile boundaries. This reduces redundant triangles,
not the spatial sampling resolution, and does not fill holes to reduce file size.
This is neither an exact physical print nor the original
solid model: actual extrusion spreading, shrinkage and inter-layer bonding are
not measured by G-code. STL also discards material identities. External mechanical
analysis must establish its own material regions, bond assumptions, boundary
conditions and mesh-resolution convergence.

The GUI takes a bounded, owned snapshot of extrusion segments, independent of
the mutable preview result. Reconstruction and file writing run as a cancellable
UI job with progress updates. Memory, segment, sampling-work and output limits
bound the operation. The GUI allows a 1 GiB reconstruction workspace, up to
4 GiB of geometry output and 50 billion work visits. These are ceilings, not estimates
of total application memory or guarantees that every plate fits. Detailed
exports can take minutes and produce files that are too large for some external
tools. Cancellation is checked during tile traversal, field
sampling, planar reduction and file writing. Each tile's OpenVDB meshing call
is not interruptible internally, so cancellation can wait for that bounded call
to return. The exporter never holds the entire output mesh in memory.

Output is written to a private temporary file beside the selected destination.
Tiles are streamed into one seekable binary STL; its final triangle count is
patched into the header only after all tiles complete. The application checks
free space before writing, during progress updates and at timed checkpoints,
retaining 30 GiB of headroom plus the remaining output allowance. The preflight
reserves the entire output budget; subsequent worker-side checks subtract the
private file's current size and fail closed if that size cannot be verified.
Only successful,
uncanceled completion replaces the destination. Partial temporary output is
removed on cancellation or failure; the existing destination is preserved.

## Separate-bead STL, STEP and OBJ

STEP and OBJ export each extrusion segment or joined path as its own closed, faceted
nominal bead, including the same roles as STL and independent of preview filters.
Arcs have already been tessellated by the G-code processor, so one arc command
can produce multiple objects. Names retain the first segment sequence, G-code ID,
processor extruder ID and color ID; these identifiers are not material property
definitions or a guarantee of a physical-tool mapping.

The endpoint cross-sections sample the horizontal-disc-plus-sphere envelope at
8 or 16 angular intervals for OBJ, and 16 or 32 for separate-bead STL.
Compact OBJ is the default format in the save dialog; its eight-sample option
has coarser curved surfaces and is explicitly selected by the user. OBJ coordinates
are rounded to six decimal places in millimetres (at most 0.5 nanometres per axis),
with trailing zeroes omitted. Relative face indices reduce text size without changing
connectivity. STEP retains its existing numeric precision.
A bounded cache stores the convex cross-section
polytope for each exact width/height pair. Its forward/backward facets and swept
silhouette construct the segment hull, including slopes, without a new 3D hull
calculation for every line. Dimensions and angular detail are not quantized.
Construction uses local coordinates before translating to the
G-code location. It does not intentionally fill intervening gaps or merge distinct paths.
Adjacent beads may overlap; their volumes must not simply be summed as the
physical printed volume. Simulation tools must establish interfaces or unions
explicitly. Angular sampling is not an absolute surface-error tolerance, and
neither format captures measured flow, thermal deformation or bonding.

STEP uses millimetre units and a separate product definition for each bead.
Compact STEP stores AP242 tessellated solids with indexed coordinates and
triangles; receiving software must support tessellated STEP. The alternative
CAD-solids option stores AP214 planar-face faceted B-reps. Both use 16 angular
samples. Adjacent exactly coplanar triangles share a polygonal CAD face, avoiding
duplicate planes and placements without moving the boundary. Patches with multiple
boundary loops or ambiguous incidence retain their triangles. Real-valued text
retains 13 significant digits without redundant trailing zeroes and includes the
decimal point required by STEP even for integer-valued coordinates. Formatting
uses a thread-local C numeric locale and restores the caller's locale.
B-reps can still exceed the output budget on detailed prints; export stops
rather than silently reducing fidelity.
These representations are not analytic cylinders or editable original CAD features.
OBJ uses millimetre coordinates and named
objects with separately indexed, closed triangle meshes. Neither stores a
whole-print CAD model in memory: a bounded segment mesh and text block are
written at a time. Byte limits are checked before each block, and cancellation
is checked between segments. This prevents the exporter from accumulating
millions of CAD objects in memory, but a receiving CAD application may still
need substantial memory to load the resulting assembly. The GUI uses the same
temporary-file publication and disk headroom checks as STL.

### Continuous-path compaction

Exact collinear planar segments with identical width, height, extruder and color
are merged without changing the nominal occupied envelope. Optional path
compaction offers 0.005 mm and 0.001 mm tolerances, or unchanged paths. A bounded
Ramer–Douglas–Peucker pass limits the two-sided centreline displacement. It does
not cross layer, material, dimension or discontinuity boundaries. Gentle bends
can then share cross-section rings instead of repeating rounded caps at every
G-code sample. Miters add at most the selected tolerance again; cross-section
faceting error is additional. Sharp turns and detected distant self-contact fall
back to separate beads. Each joined mesh is limited to 128 segments. Chunk ends
remain capped, and separate meshes can still overlap.

This is not a global Boolean union or conformal material-interface reconstruction.
Gaps comparable to the approximation tolerance require inspection. The retained
air spaces and material identifiers do not establish bonding, dielectric properties
or mechanical constitutive data. In particular, CADFEKO's STL mesh import defaults
to PEC surfaces, not polymer dielectric regions; its media must be configured
explicitly. STEP CAD-solid regression checks use OpenCascade import and solid
validity checks; they do not establish compatibility with a particular Feko release.
