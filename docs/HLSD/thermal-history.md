# Transient thermal history

`ThermalHistory::step` advances the temperatures of an explicitly supplied active-cell network. A cell has volume, volumetric heat capacity, initial temperature in kelvin and source power. Symmetric connections supply thermal conductance in W/K. Boundary connections exchange heat with an explicitly specified reservoir temperature.

The finite-volume energy balance uses backward Euler:

`C_i (T_i,new - T_i,old) / dt = Q_i + sum_j G_ij (T_j,new - T_i,new) + sum_b G_ib (T_b - T_i,new)`.

Positive cell heat capacities and nonnegative conductances yield a symmetric positive-definite system. A matrix-free conjugate-gradient solve uses the diagonal as a preconditioner. Input size and iteration limits bound work and memory; cancellation propagates. The true residual is checked after iteration, and unsuccessful solves publish no temperature field. Reported stored, source and boundary energies provide an independent discrete conservation check.

The caller owns deposition activation: only deposited cells enter the network. Adding a cell contributes its actual mass and initial energy, not an ambient-temperature fictitious solid. Occupied contacts and exposed boundaries must be updated with activation. Temperature-dependent properties, contact changes, radiation linearization and substepping require caller-controlled re-evaluation; this core solves one fixed-coefficient step and does not silently invent these inputs.

Boundary labels distinguish bed, ambient, part-fan, auxiliary-fan, dock-fan and exhaust effects. Labels do not supply a calibration or change the equations. Their conductances must represent effective heat transfer for the actual exposed area and airflow. Fan PWM is not a heat-transfer coefficient. Multiple boundary contributions must not double-count the same natural/forced convection; directional airflow, heated air and dock shielding require explicit scenarios or measured data. The core contains no assumed material properties or fan curves.

An enclosure-air reservoir can be represented explicitly with its heat capacity,
connections to part/bed surfaces, and an exhaust heat exchange with room air.
This distinguishes exhaust heat loss from circulation inside the enclosure.
It is a lumped-zone approximation, not CFD: rear-fan jets, inlet drafts and
spatial temperature variation require multiple zones or separately bounded
surface coefficients. A controller target is not an imposed uniform air
temperature. Exhaust flow, controller behavior and room temperature remain
explicit inputs; no automatic mapping from the target or fan PWM is supplied.

Regression tests cover internal energy conservation, source/boundary balance, activation energy, analytical two-cell exchange, first-order time convergence to lumped cooling, invalid inputs, budgets and cancellation. These are numerical verification, not experimental validation of a printer or material. Published transient-temperature benchmarks and bounded sensitivity studies remain necessary before using this model to choose print settings. Thermal histories alone do not establish adhesion strength, nozzle-collision safety or a guaranteed absence of corner lift.

## Fan command interpretation

`AirflowTimeline` retains changes in commanded part, auxiliary, dock and exhaust
fan PWM using an explicit printer-specific channel map and installed-firmware
policy (`GenericM106`, `RmeLegacy` or `RmeDockCooldown`). Aliases must be declared;
there is no universal assumption that P1/P2/P3 denotes a particular fan.
An unassigned fan command marks the history incomplete rather than silently
discarding another source of airflow. Initial fan states are unknown.

Legacy RME firmware drives P6 directly: tool parking does not change that request,
and M979 is flagged as unsupported rather than interpreted as automation. Its
M106 handler requires S/A/R; regulator-only N/G commands do not imply full PWM.
The optional dock-cooldown policy starts with unknown control mode. M979 selects
that mode; parking/tool changes and manual requests during automatic cooldown
cannot establish actual dock PWM. Explicit manual mode survives tool changes.
Both RME policies leave rear exhaust output unresolved, including after M141 or
M191: chamber feedback, filtration floors and safety overrides can change it.
A temperature target is not a measured temperature or airflow. These events
require a runtime controller model or uncertainty
scenario. PWM is never converted directly into a heat-transfer coefficient.
The command trace is an input utility, not yet an automatically calibrated
boundary-condition model or a deposition scheduler. Combined fan effects on the
same surface must be supplied as a defensible effective boundary condition, not
double-counted by independently adding overlapping air streams.

Configured static G-code mechanics accepts an optional matching source stream,
explicit fan-channel map and installed-firmware policy. The bounded scan consumes
from the stream's current position without reopening a filename or guessing
firmware from a printer name. It reports final command state, unresolved roles
during E-bearing movement, and incomplete scans caused by byte/event/line limits
or stream errors. Cancellation propagates. The summary survives a mechanical
material or geometry failure, so unsupported mechanics does not conceal missing
airflow inputs. Without this explicit source the summary is `NotProvided`.
These command observations do **not** enter the current prescribed static
temperature scenarios or their force ranking; the summary explicitly records
that distinction. They are not measured airflow or a calibrated convection law.

`ThermalDeposition` advances explicit physical-time events with the thermal core.
Events append newly deposited mass at its initial temperature and connect only
already active cells. Each event replaces the exposed boundary list, allowing a
covered surface to stop convecting. Steps end at event times, and a maximum time
step bounds temporal discretization. Event, step, active-cell and cumulative
cell-step limits bound storage and work; cancellation propagates and failed
histories publish no temperature field. Deposited sensible energy (using a
Kelvin reference), boundary heat and volumetric source heat are accounted for
separately against final stored energy.

Optional local time-error control compares one backward-Euler step with two
half steps and accepts the half-step temperatures only when the maximum cell
temperature difference meets the supplied kelvin tolerance. Rejected trials
commit neither temperatures nor heat accounting. All trial solves, including
rejections, count against the step and cell-work budgets. Event boundaries are
never crossed. Results distinguish accepted intervals, solver calls, rejected
trials, and the maximum accepted local error estimate. A zero tolerance retains
the fixed-step API. This is a local truncation-error estimate, not a global error
bound or experimental validation; temporal refinement must still be checked
for derived mechanical quantities. Linear-solver residual convergence alone
does not establish time-discretization accuracy.

Events can be pulled from a bounded producer, keeping only the current exposed
boundary list instead of retaining a copy for every activation. The spatial
adapter uses this streaming path. Exposure-reconstruction work still has a
cumulative budget; streaming reduces storage, not the number of solves needed
for finely separated activation times.

This scheduler does not infer deposition time, activation geometry, material
properties, surface exposure or effective convection from G-code. Its caller
must supply these consistently; an existing cell cannot be silently replaced
with fresh hot mass. Synthetic tests verify activation and energy accounting,
not a printed-material calibration or an adhesion prediction.

`ThermomechanicalGCode` connects processor moves to this scenario evaluator.
It requires explicit acknowledgement of modeled move timing and explicit
deposited-polymer temperatures; nozzle setpoints are never substituted. The
processor duration includes some lumped delays, so physical deposition phase
timing remains unresolved even when cumulative durations close to its total.
Display-only acceleration vertices are omitted without shortening source paths.
An optional observation time stops the simulation within the modeled timeline;
only cells activated by that time enter its thermal and mechanical body. This
supports snapshots during deposition rather than treating the final part as
present from the start. It does not remove the processor's timing uncertainty.

Cell activation uses the nominal bead footprint's first-hit segment and path
fraction. The nominal vertical cross-section is not measured bead shape;
same-material overlaps add neither duplicate mass nor reheating. Explicit
surface-exchange scenarios are uniform over the exposed faces they address,
not a resolved fan-flow field. The final elastic solve uses reference geometry
and does not reconstruct a stress-free growing body or establish available bed
adhesion. These assumptions accompany results and do not authorize automatic
mitigation selection by themselves. A hidden real-export test audits extraction
and modeled-time closure only; it is not full-case physical validation.
