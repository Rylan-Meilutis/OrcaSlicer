# Thermal mechanical reference data

Material-history calculations require properties tied to material grade,
manufacturing process, temperature, loading rate and measurement method. The
reference datasets below support numerical verification and explicitly named
parameter studies. They are not automatic defaults for an unidentified filament.
In particular, neither an upper measured modulus nor the absence of relaxation
defines a universal worst-case ASA corner-lift scenario.

## ASA

[Temperature-Dependent Plastic Behavior of ASA: Johnson–Cook Plasticity Model
Calibration and FEM Validation](https://pmc.ncbi.nlm.nih.gov/articles/PMC12897946/)
measures injection-molded BASF Luran S778T at a tensile strain rate of 0.01/s.
Reported Young's moduli are 2487±196 MPa at −10°C, 2206±103 MPa at 23°C and
2102±156 MPa at 65°C. The paper assumes Poisson's ratio 0.35 rather than measuring
it. These data do not establish printed orthotropy, a relaxation spectrum, or
constitutive behavior near the glass transition.

The [Stratasys ASA manufacturer datasheet](https://www.stratasys.com/siteassets/materials/materials-catalog/fdm-materials/asa/mds_fdm_asa_0921a.pdf)
reports printed-specimen mean coefficients of thermal expansion of 69.38
µm/(m·K) in XY and 63.55 µm/(m·K) in XZ/ZX, measured over −50 to 90°C by ASTM E831.
Its DSC glass-transition inflection is 103.55°C. These are reference-grade
measurements, not confidence bounds for all commercial ASA. Extrapolating the
CTE or glassy stiffness to a 110°C bed-contact region is not justified by this
dataset; bed setpoint also does not establish polymer temperature.

The [2025 Stratasys ASA datasheet](https://www.stratasys.com/siteassets/materials/materials-catalog/fdm-materials/asa/mds_fdm_asa_0925a.pdf)
also supplies matched natural-ASA conductivity and diffusivity measurements:

| Temperature, °C | Conductivity, W/(m·K) | Diffusivity, mm²/s |
| ---: | ---: | ---: |
| 0 | 0.1685 | 0.108 |
| 30 | 0.1642 | 0.096 |
| 60 | 0.1622 | 0.087 |
| 90 | 0.1563 | 0.077 |

The named reference registry derives volumetric heat capacity as conductivity
divided by diffusivity at each matching temperature. It interpolates within the
table only. Neither this transport dataset nor the Luran modulus dataset is a
measured constitutive law for the other's material grade.

[Struž et al., 2024](https://doi.org/10.3390/polym16213018) report dynamic
mechanical measurements on printed Prusa Galaxy Black ASA. Their Table 9 gives
storage modulus from 1210 to 976 MPa over 20–40°C; Figure 21 places the storage
modulus transition inflection at 105.7°C. The specimens use 100% rectilinear
infill with alternating raster directions, 0.3 mm layers and a 260°C nozzle.
The three-point bending test uses 3 Hz loading and a 2°C/min heating ramp.
These measurements concern a printed grade rather than injection-molded resin,
but storage modulus is frequency-dependent. They are not substituted for static
Young's modulus or an hours-long relaxation law in the reference registry.
In particular, the plotted softening curve does not establish the stress-free
temperature or residual strain of a newly deposited road.

## Named unfilled TPU references

The [UltiMaker TPU 95A TDS, version 3.010](https://ase.au.dk/fileadmin/www.ase.au.dk/Filer/Laboratorier_og_vaerksteder/MakerLab/Ultimaker-TDS_TPU_95A_v3.010_ENG.pdf)
reports printed tensile modulus 26 MPa, specific gravity 1.22 and CTE
100 µm/(m·K). The CTE temperature interval is unspecified. These constants are
retained as reference information, not silently extended into temperature laws.
Its −24°C glass transition is not an appropriate stress-free solidification
temperature for a deposition model.

[Delbart et al., 2024](https://doi.org/10.1016/j.compositesa.2024.108381)
model unfilled PolyFlex 95A separately from conductive CB-PLA. Their Table 3
uses initial modulus 8.4 MPa, conductivity 0.15 W/(m·K), specific heat
1210 J/(kg·K), and CTE 110 µm/(m·K). These are published model assumptions,
including a nonlinear mechanical response, not a complete measured
temperature-dependent TPU calibration. The study explicitly disregards
printing residual stresses. The registry therefore does not present its
constants as a validated deposition law or extrapolate them automatically.

`ThermalMaterialReferences` uses explicit reference IDs, never generic filament
names. Per-property evidence distinguishes measurements, derived values and
published assumptions. Missing density, Poisson ratio, reference/deposition
temperature, relaxation, interface resistance and boundary exchange remain
missing. Scenario authors must state how they complete such inputs. Min/max
outputs across selected scenarios describe those scenarios only: they are not
population confidence bounds or guaranteed bounds on corner lift. Thermal
gradients, stiffness contrasts and contact changes need not vary monotonically
with individual input parameters.

## TPU relaxation verification

[Viscoelastic Behaviour of Flexible Thermoplastic Polyurethane Additively
Manufactured Parts](https://pmc.ncbi.nlm.nih.gov/articles/PMC8309542/) measures
Filaflex 82A at 20±0.5°C. Table 7 provides the normalized Prony spectrum used by
the scalar relaxation regression:

| Relaxation time, seconds | Modulus fraction |
| ---: | ---: |
| 0.054063 | 0.14292 |
| 0.14614 | 0.021163 |
| 0.39504 | 0.029420 |
| 1.0678 | 0.035212 |
| 2.8866 | 0.013225 |
| 7.8028 | 0.036500 |
| 21.092 | 0.00080567 |
| 57.015 | 0.063713 |

The convention is `E(t)/E(0) = 1 - sum(e_i * (1 - exp(-t/tau_i)))`.
The reported average `E(0)/E(10 seconds)` is 1.390±0.018: approximately 72% of
initial modulus remains, not 72% lost. Agreement verifies implementation of this
constitutive equation. It does not validate another TPU grade or thermal history.

[A Hybrid Experimental-Numerical Framework for Identifying Viscoelastic
Parameters of 3D-Printed Polyurethane Samples](https://www.techscience.com/cmc/v86n3/65492/html)
provides a separate REC TPU Easy Flex dataset, printed at 100% infill and tested
at 23°C. Its compression-rate-fitted standard-linear-solid model has equilibrium
modulus 31.1 MPa, relaxing spring 16.7 MPa and viscosity 109.5 MPa·s. The
creep-derived fit instead uses 25.9 MPa, 10.1 MPa and 1378 MPa·s. These distinct
fits illustrate the importance of loading protocol; their coefficients must not
be combined into an invented material calibration.

Dynamic storage modulus depends on oscillation frequency and is not generally
the instantaneous or long-time static modulus. Likewise, Shore hardness does
not uniquely identify TPU chemistry or its temperature-dependent relaxation.
An isothermal reference spectrum needs independently supported shift factors
or temperature-specific data before use at other temperatures. The scalar
kernel therefore rejects updates away from its explicit reference temperature.

## Bed-interface quantities

[Print Bed Adhesion in Fused Filament Fabrication of ABS, PETG and PA12](https://doi.org/10.1002/masy.70120)
investigates these polymers on several surfaces including PEI, but tests shear
adhesion at room temperature and 1 mm/min. Its 5.3 MPa ABS maximum uses PVA.
That value is neither an ASA-on-PEI hot peel threshold nor a fracture energy.
Surface-energy measurements likewise do not provide a mixed-mode cohesive law.

The mechanical boundary model distinguishes prescribed bonded contact from
prescribed release. Finite-increment stored-energy reduction per released area
is an interface driving demand, not interface resistance. Predicting propagation
requires appropriate normal and shear traction–separation behavior, fracture
energies and mixed-mode dependence at the relevant temperature, rate and surface
condition. No resistance is inferred from a material name or a bed setpoint.

Unknown polymer temperatures, relaxation, surface condition and adhesion can
trade off against one another when fitting a lift observation. One photograph or
a near-detachment height cannot identify them independently. Scenario comparison
and numerical convergence therefore remain separate from physical validation of
corner-lift onset and propagation.
