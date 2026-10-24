# Spatial thermoelastic equilibrium

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

`voxelize_deposition` rasterizes horizontal deposited segments using their
reported width and height. The footprint is a capsule and the vertical section
is rectangular; this is a deposition-envelope approximation, not a resolved
bead cross-section or inter-bead weld model. Grid cells are selected by center
containment. The grid must resolve at least two cells across the width and
height. Infill voids are not replaced with the CAD solid envelope. Conflicting
material assignments in a cell, nonplanar paths, invalid dimensions and exceeded
resource budgets produce an unresolved analysis instead of an averaged material.

`analyze_gcode` obtains those segments from processor moves, excluding travel,
custom purge paths and the prime tower. It uses explicitly supplied material
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
