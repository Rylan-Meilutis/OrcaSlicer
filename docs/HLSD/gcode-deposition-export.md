# G-code deposition STL export

The **Export printed extrusion as STL** action is available in the print/export
button's action selector and in **File → Export**. It reconstructs the current
plate's generated toolpaths, not its original CAD volumes. Preview filament
visibility and layer-slider filters do not restrict the export. Infill, walls,
supports, skirts, brims and the prime tower share one STL; disconnected printed
islands remain disconnected shells in that file.

The geometry is a union of nominal rounded extrusion beads. The user chooses
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
32 GiB of STL output and 50 billion work visits. These are ceilings, not estimates
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
