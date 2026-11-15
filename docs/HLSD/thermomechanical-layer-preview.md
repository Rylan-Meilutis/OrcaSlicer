# Mechanical fields at a selected print layer

`preview_mechanical_layer` evaluates an explicit thermal/material/boundary
scenario at the end of the selected processor layer's last model-deposition
command. Layer identity comes from the processor, not a comparison against Z:
sloped paths and different extrusion heights need not lie on one plane. The
observation uses modeled command timing, not measured printer timing. Subsequent
deposition and final cooldown do not contribute to that observation.

The bridge calls the transient deposition study and requires its independent
time-step and spatial comparisons to qualify integrated force, energy, heat and
contact area. An unsupported, incomplete or resource-limited result contains no
sample map. Cancellation propagates, including cancellation expressed with the
same exception type used for resource limits. A completed linear solve by itself
does not authorize publication.

For a qualified scenario, each included physical deposition move is sampled at
its nominal bead midpoint. Samples contain the element stress tensor, maximum
principal stress, von Mises stress and displacement in SI units. They are not
maxima along long moves. Missing spatial coverage remains explicitly unsampled;
there is no interpolation across holes or extrapolation to the nearest solid.
Bed samples are separate: node positions, tributary contact areas, actual
reaction vectors and area-normalized tensile restraint demand. Stress is never
multiplied by an arbitrary area to produce a purported curl force.

Integrated discretization qualification does not establish pointwise stress or
sharp-edge traction convergence. The fully bonded model measures required bed
restraint, not available adhesion, crack onset or a lift guarantee. Its mechanical
state is an elastic equilibrium snapshot driven by the simulated temperature
field; it does not model stress-free deposition onto an already-deformed body or
complete viscoelastic process history.

Scenario inputs are supplied explicitly. Filament names, fan PWM and heater
targets do not supply missing constitutive properties, heat-transfer coefficients
or a bed fracture law. The transient backend retains its fine-grid resource
limits and only accepts explicit, authenticated variable-thickness deposition; the nominal
section solver is not substituted as though it provided thermal history.
Consequently a large mixed ASA/TPU case may remain unavailable even when its
nominal scarf geometry can be sectioned successfully. The returned standalone
arrays are associated with the source move indices and selected observation;
consumers must invalidate them when either the source or scenario changes.
# Prescribed-strain section preview

`preview_section_layer` is a separate static scenario path. It selects the
model geometry deposited by the end of a processor layer, then applies explicit
material elasticity tensors and eigenstrains through the bounded section study.
It does not calculate temperature history or infer strain from fan commands.

Requested move midpoints are sampled while the certified slab topology and Q1
cell basis are alive. Material membership and polygon holes are checked before
interpolation. Displacement follows Q1 shape functions; stress is computed from
their displacement gradient minus the specified eigenstrain. Points outside the
sampled material are counted as unavailable, never assigned zero stress.

Only selected-layer move samples are retained. The geometry/basis study must
satisfy its integrated convergence criteria before any field is published;
this does not independently establish convergence of local peak stresses.
The integrated tensile reaction is a sum of generalized constrained-basis
forces. It is not a local peel traction, adhesion capacity or prediction of
detachment. No per-face traction is fabricated from these reactions.

## Tile-local section quadrature

An explicit `SectionStudySettings::local_tile_size_m` selects spatially local
sectioning. Each XY tile indexes only nominal beads whose conservative envelope
bounds meet it. Its Z partition retains those beads' endpoint top and bottom
coordinates, including linearly varying scarf thickness. Midplane footprints
are clipped to the tile before cross-material conflict reconciliation. Empty
tiles and empty intervals contribute no material; internal infill voids and
polygon holes are not filled. Index work, references, active-bead visits,
retained polygon points and mechanical assembly have independent hard bounds.

The synchronous tile visitor shares the retaining sampler's implementation but
releases each complete tile after its consumer returns. It returns scalar
metrics only after the complete traversal succeeds. Cumulative geometry and
work limits still apply; a consumer sees unqualified geometry, not a topology
certificate or a partial mechanical result. This allows consumers to reduce
geometry incrementally without requiring all tile polygons to remain alive.

Adjacent sampled intervals share one stored prism only when their complete
material partitions have identical native-coordinate contours, holes, material
identities and boundary uncertainty. Every sampling plane is still evaluated.
This compacts the sampled solid without bridging empty intervals or independently
merging one material across a change in another. The retained prisms may be
taller than the sampling interval; their volume and Q1 polynomial integrals are
unchanged. Comparisons consume the bounded indexing-work budget.

Neighboring tiles need not share Z partitions. The topology certificate joins
occupied islands through positive-area horizontal faces, or positive-length
tile-boundary traces with positive Z overlap. Edge/point-only contacts do not
anchor disconnected material. Every component needs actual positive-area bed
contact, and independent components have independent mechanical degrees of
freedom. Generalized Q1 integration uses the clipped prisms rather than a solid
tile or a density-scaled modulus.

The local path still requires independent basis convergence and two successive
geometry refinements with stable volume, elastic energy, integrated tensile
reaction and component count. Horizontal rectangular envelopes without ownership
reconciliation are exactly piecewise constant between their retained boundaries.
Numerical qualification applies only to the supplied nominal representation and
prescribed strain scenario, not the printed bead shape, material bond or adhesion.
The same fine-geometry/assembly caps remain enforced; selecting local quadrature
does not permit unbounded memory or silently return an unconverged coarse field.

`initial_z_quantization_m` is a separate opt-in integration-partition
approximation. Endpoint tops round upward and bottoms downward only to choose
sampling-interval boundaries on a bed-relative grid; the increment halves
between refinement levels. Cross-sections always evaluate the original beads,
including their original rounded profile and zero-height tips. The report
retains the grid increment and maximum partition-boundary shift. A positive
bead missed by the coarse midplanes adds its original local endpoints to the
partition. An internal artificial cut one native unit from a retained exact
endpoint is removed so that it does not create an interval with no integer
midplane. Retained endpoints and the enclosing extremes are never removed.
Coverage is rechecked after repartitioning, under the same work and section
limits; unresolved native-scale intervals are rejected rather than discarded.
The explicit half-native plane mode instead retains these cuts and samples a
one-nanometre interval at its half-nanometre interior. Only query-plane resolution
changes; bead boundaries, dimensions and XY geometry retain their native grid.
A partition envelope that introduces bed contact is rejected. Printing G-code
and material identities are unchanged.

Sub-grid gaps and contacts still require independent resolution, so successful
refinement has the distinct `ApproximateGeometryAndBasisConverged` state. It
requires two successive stable geometry comparisons and separately converged
mechanical bases, even for horizontal rectangular input. Stability is numerical
evidence for the specified quadrature family, not a physical error bound or proof
that all features below the finest increment are resolved.

## Authenticated material conflicts

`GCodeSectionScenario::use_native_ownership` prepares original nominal geometry
once for the same modeled observation prefix, excluding the same tower/custom
roles and applying the same explicit stationary-extrusion admission. The worker
borrows an immutable source result and verifies its registry/snapshot identity on
each query. A bounded existing AABB tree indexes original XYZ bead envelopes;
construction and queries remain cancellable. Total candidate and polygon work
is bounded across queries, and candidate order retains original deposition order.

Ownership is deferred until sectioning has assembled the complete actual
conflict, including tile clipping. Every contributing original bead is examined;
later deposition and unrelated unowned material do not authorize or veto that
contact. Original multi-material overlap uses the authenticated native masks.
Every contributing source still needs authentication. Conflict area in original void is rejected, never
assigned a plausible-looking owner. Native layer, source, material and transformed-path restrictions
remain enforced. Resolved masks partition conflicts only and must leave the
deposited union unchanged. Imported G-code lacking matching live provenance
cannot use this evidence path.

## Preview publication

The GUI's section-analysis request snapshots the matching processor result under
its mutex before submitting a cancellable UI-worker job. No plate or viewer
pointer is used by the numerical worker. Finalization checks the slice identity
and filename before replacing the scalar map; replacing or resetting the loaded
slice invalidates the previous field. Scenario selection remains explicit.
An invalid sample in a matching result clears the previous field and publishes
a failure state; it cannot silently leave a previous layer's colors visible.

Rendered vertices retain their processor move identity, including synthetic
path-start vertices. A result therefore updates both ends of a rendered path
without confusing G-code line numbers with processor moves. Missing samples are
cleared, not carried forward from another layer. Scalar replacement changes CPU
data only and invalidates the visibility-keyed range cache; the next render
updates the color range and GPU colors in its normal
OpenGL context without reconstructing toolpath geometry. The section legend
distinguishes prescribed static contraction from transient cooling history and
does not label generalized bed reaction as local peeling traction.
