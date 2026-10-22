# Thermal deformation and bed adhesion

The reduced-order model evaluates connected cross sections through the first
layer of each printable, transformed model instance. Disconnected feet are not
joined across their empty bounding box. Each contact is represented by a strip
whose length is twice its greatest contour distance from its centroid. The solid
base thickness is the lesser of object height and configured bottom-shell
thickness (including the bottom-layer count). Holes remain holes in the section;
their stress concentrations and spatially varying plate stiffness are not solved.
The model is a screening calculation, not finite-element analysis.

For applicable amorphous polymers, solid thermal strain is alpha * delta_T.
The stress-free reference is bounded by glass transition and deposition
temperature, rather than treating the entire nozzle-to-room drop as elastic
shrinkage. For a linear bed-to-chamber gradient through thickness h, free strip
curvature is alpha * delta_T / h; centre-to-edge bow is curvature * length^2 / 8.
Uniform cooling separately produces free dimensional contraction strain * length.
Fully restrained uniaxial elastic energy per area is E * strain^2 * h / 2,
with thickness converted to metres for J/m².

For semicrystalline resins, the stress-free reference is bounded by the published
melting endpoint instead of Tg: crystals can transmit stress above Tg. Extending
solid CTE and room-temperature modulus to this reference is a no-relaxation
screening scenario, not a measured cooling curve. Separately, mass conservation
gives a full amorphous-to-crystalline isotropic contraction bound of
`1 - cbrt(rho_amorphous / rho_crystalline)`. The displayed interval runs from zero
to that bound times contact span; it does not assume the printed part becomes
fully crystalline. This term is not added to thermal contraction or curvature:
phase fractions and through-thickness crystallization kinetics are not known.
Uniform phase shrinkage cannot establish curl. The interval therefore does not
alone trigger automatic mitigation. These calculations cover thermal screening
and phase-volume uncertainty, not a solved crystallization-induced warp field.

These are limiting scenarios, not confidence intervals or observed deformation.
Relaxation and restraint can reduce deformation to zero. Chamber and bed targets
are not measured polymer temperatures. With no chamber target the ambient
scenario is 20 C, explicitly identified in the UI. Drafts, transient reheating,
anisotropy, infill compliance and temperature-dependent modulus are not resolved.
The small-deflection approximation becomes qualitative at large bow.
Multi-material objects use the envelope of assigned model materials across their
contact sections, not an invented homogenized composite modulus. Modifiers and
feature assignments are included conservatively; spatial material interfaces are
not mechanically simulated. Supports are not treated as model material.
Unknown data never implies safe geometry.

## Published data

Published references below supply solid CTE ranges (1/K), glass transition
ranges (C) and bulk room-temperature Young's modulus (Pa):

| Material | CTE | Glass transition | Maximum modulus |
| --- | --- | --- | --- |
| ABS | 80–100e-6 | 95–105 | 3.0e9 |
| ASA | 85–105e-6 | 95–105 | 2.9e9 |
| PC | 75–80e-6 | 140–150 | 2.4e9 |
| HIPS | 80–90e-6 | 101.1 | 2.0e9 |
| PETG | 51e-6 | 80 | 2.1e9 |
| PMMA | 90–110e-6 | 45–115 | 3.3e9 |
| PLA | 80–85e-6 | 45–65 | 2.8e9 |
| PP | 130–180e-6 | -20–20 | 1.8e9 |
| PA6 | 80–90e-6 | 45–80 | 2.8e9 |
| PA11 | 85–120e-6 | 40–55 | 1.4e9 |
| PA12 | 120–140e-6 | 40–50 | 1.4e9 |
| POM | 160–180e-6 | -85–-75 | 3.2e9 |
| PPS | 50–70e-6 | 85–100 | 3.7e9 |
| PEEK | 50–70e-6 | 145–155 | 3.7e9 |

- [ABS](https://analyzing-testing.netzsch.com/en-US/polymers-netzsch-com/commodity-thermoplastics/abs-acrylonitrile-butadiene-styrene-copolymer)
- [ASA](https://analyzing-testing.netzsch.com/en-AU/polymers-netzsch-com/commodity-thermoplastics/asa-acrylonitrile-styrene-acrylate-copolymer)
- [PC](https://analyzing-testing.netzsch.com/en-US/polymers-netzsch-com/engineering-thermoplastics/pc-polycarbonate)
- [PMMA](https://analyzing-testing.netzsch.com/en-US/polymers-netzsch-com/engineering-thermoplastics/pmma-polymethylmethacrylate/): Tg depends on tacticity; the upper endpoint is the screening scenario.
- HIPS: [Versalis E05](https://versalis.eni.com/assets/documents/versalis/en/portfolio/polimeri-e-intermedi/stirenici/schede-tecniche/gpps-hips/TDS_ENG_HIPS_E05.pdf) supplies 90e-6/K CTE and 1.9 GPa modulus; [INEOS 495N](https://www.ineos-styrolution.com/Product/HIPS_Styrolution-PS-495N_SKU200400070622_lang_es_ES.html) supplies 80e-6/K and 2 GPa. The 101.1 C transition comes from the unfilled HIPS filament in [Sieradzka et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8623337/), not Vicat or HDT. Combining different grades is an explicit representative approximation.
- PETG: [Eastman Eastar 6763](https://productcatalog.eastman.com/tds/ProdDatasheet.aspx?pn=Eastar+6763+copolyester&product=71040786) supplies 51e-6/K CTE, 80 C Tg and 2.1 GPa injection-molded tensile modulus. CTE is measured over -30 to 40 C; carrying it to Tg is a constant-CTE screening extrapolation, not a measured high-temperature curve. These are representative grade values, not universal PETG bounds.

- PLA: [ORNL thermal-property measurements](https://www.osti.gov/servlets/purl/1463519) for solid CTE; [NETZSCH PLA](https://analyzing-testing.netzsch.com/es/polimeros-netzsch-com/commodity-thermoplastics/pla-polylactide-polylactic-acid/) for Tg, modulus and melting endpoint. Its cooling example also illustrates why full crystallization must not be assumed.
- [PP](https://analyzing-testing.netzsch.com/id/polimer-netzsch-com/commodity-thermoplastics/pp-isotactic-polypropylene), [PA6](https://analyzing-testing.netzsch.com/en-US/polymers-netzsch-com/engineering-thermoplastics/pa6-polyamide-6), [PA11](https://analyzing-testing.netzsch.com/es/polimeros-netzsch-com/engineering-thermoplastics/pa11-polyamide-11), [PA12](https://analyzing-testing.netzsch.com/nl/polymeren-netzsch-com/engineering-thermoplastics/pa12-polyamide-12/) reference solid properties and transitions.
- [POM homopolymer](https://analyzing-testing.netzsch.com/en-US/polymers-netzsch-com/engineering-thermoplastics/pom-homo-polyoxymethylene-homopolymer), [PPS](https://analyzing-testing.netzsch.com/da/polymerer-netzsch-com/high-temperature-resistant-thermoplastics/pps-polyphenylenesulfide/), [PEEK](https://analyzing-testing.netzsch.com/ja/polymers-netzsch-com/high-temperature-resistant-thermoplastics/peek-polyetheretherketone) reference solid properties and transitions.

| Resin | Melting endpoint C | Amorphous / crystalline density g/cm³ | Phase-density reference |
| --- | --- | --- | --- |
| PLA | 160 | 1.248 / 1.290 | [PLA morphology study](https://www.sciencedirect.com/science/article/abs/pii/S0032386116305948) |
| PP | 165 | 0.855 / 0.936 | [Isotactic PP study](https://research.chalmers.se/publication/549830/file/549830_Fulltext.pdf) |
| PA6 | 235 | 1.08 / 1.23 | [Verdu, polymer properties](https://sam.ensam.eu/bitstream/handle/10985/18689/PIMM_OMSC_2011_VERDU.pdf?isAllowed=y&sequence=3) |
| PA11 | 190 | 1.01 / 1.18 | Same reference; polymorph-dependent endpoint |
| PA12 | 180 | 0.99 / 1.085 | [In-situ PA12 deformation study](https://www.mdpi.com/2073-4360/11/6/1001) |
| POM | 190 | 1.215 / 1.492 | [POM feedstock study](https://www.mdpi.com/1996-1944/13/7/1535) |
| PPS | 290 | 1.32 / 1.43 | [Toray resin properties](https://www.plastics.toray/de/technical/torelina/tec_029.html) |
| PEEK | 345 | 1.265 / 1.400 | [PEEK morphology reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC2040108/) |

Upper endpoints are used for the displayed scenario. These representative bulk
ranges are not guaranteed bounds for every spool. `PA`/`Nylon` uses the independent
envelope of PA6/11/12 properties, flagged as a family scenario, not identified
nylon chemistry. Named nylon grades use their own references. Aromatic nylons,
copolymers, foamed grades, unknown materials and proprietary blends without a
matching reference remain unassessed rather than silently treated as ABS.

CF/GF variants of the listed resins receive a separately labelled **resin-only
screening scenario**. It takes no credit for reinforcement reducing contraction,
does not assign a fictitious fiber fraction, and does not report the matrix stored
energy as composite energy. Fiber alignment, loading, moisture and stiffness can
change actual deformation; the resin scenario is not a guaranteed directional
bound or a composite constitutive model. See [NETZSCH's anisotropic-filler
measurements](https://analyzing-testing.netzsch.com/en-AU/blog/2021/why-the-effect-of-anisotropic-fillers-on-thermal-expansion-is-process-dependent/).
Case, spaces and CF/GF separator spelling are normalized, but arbitrary additive
suffixes are not dropped. No spool-density inference of crystallinity is used.
Published molding shrinkage is not added to CTE strain as isolated crystallization
strain. Arc-cooling heat capacity is not a substitute for these mechanical data.

## Decisions and mitigation

A contact whose estimated free bow exceeds its configured layer height prompts
an adhesion precaution. This resolution-based decision is distinct from bed
failure prediction. Small ABS parts do not trigger merely because they are ABS.
Dimensional contraction is reported separately; uniform contraction is not curl.
Unsupported-edge curl remains the existing toolpath-stage analysis, not an
inference from material name or the base-bow model.

Bed detachment remains **unknown**: configuration does not specify a measured
peel/fracture law for the actual polymer, plate, adhesive and surface condition.
Stored elastic energy is not compared to fabricated adhesion energy or to
published shear strength (different units and failure mode).

Automatic mitigation defaults off. When enabled, or explicitly applied in Others,
elevated base-bow potential selects a continuous outer brim and overhang/curled-
perimeter slowdown. Suggested brim width is a bounded policy of twice estimated
bow, limited to 3–10 mm; it is not a solved adhesive width. Larger existing brims
and inner coverage remain. Separation gap is capped at 0.1 mm. Rafts, belt-printer
adhesion and explicit object overrides are preserved. Temperatures, fan settings
and material strength are not changed. Automatic changes affect the effective
slicing configuration, not saved presets.

An unmitigated elevated estimate produces a warning. Adequate continuous brim or
raft produces an informational note, not a warning or certification of safety.
The note only confirms mitigation and asks the user to validate changes before
printing. The Others page opens the detailed assessment in a modal dialog, with
an explicit action to apply recommendations; opening or closing it changes no settings.
Revalidation clears stale messages. Tests cover dimensional consistency, thermal
and geometric scaling, disconnected contacts, small parts, unknown materials and
warning/note behavior. They do not establish real-world failure probabilities.

## Cosmetic top relief

`cosmetic_top_relief` is an independent, experimental, opt-in process/region
setting in Localized shrinkage mitigation. It is never enabled by the thermal
screening model. Its 0.2 mm gap and 10 mm section defaults describe toolpath
geometry, not calibrated mechanical performance. It does not alter chamber,
bed, nozzle or fan targets.

After infill grouping, exposed top infill is inset and divided by X/Y grooves.
Walls, buried solid infill and bottom surfaces are not cut. Treatment requires
at least two configured top shell layers and actual non-bridge solid fill or
wall coverage of the same region immediately underneath. Small surfaces and unsupported skins
are left unchanged. Region identity is part of top-fill grouping when enabled,
so a neighboring region cannot inherit the option. Classified lower-layer fill
and completed perimeters are read without touching concurrently generated infill.
Ironing is skipped for opted-in regions so it cannot reconnect grooves.

The sections remain vertically bonded to the structural layer below: this is
not a detached skin and does not eliminate stress transfer. Visible grooves,
surface-strength and sealing tradeoffs require user validation. Preserving the
underlying toolpaths does not establish identical part strength or demonstrate
that a particular lid will stop warping. The disabled default preserves existing
projects and toolpaths; changing these options invalidates infill and downstream
operations without regenerating walls.

## Buried zigzag relief

`buried_top_relief` is a separate opt-in mode using the same gap and section-size
controls. It overrides cosmetic top relief when both are enabled. Only internal
solid infill connected through solid layers to an exposed top within the selected
depth is cut. At least four configured top shell layers are required. Bridges,
bottom surfaces, sparse infill and perimeters are unchanged. A protective rim is
retained around walls and holes. Small areas that cannot fit sections are skipped.

The channels form an object-aligned zigzag grid, with adjustable tooth depth and
pitch. The grid stays aligned across buried layers. Channel geometry accounts for
extrusion width, and is also removed from gap-fill coverage. The exposed top and
the two solid layers immediately beneath it remain uncut; ironing remains
available. This continuous three-layer cap limits surface imprinting from the
slots below. Relief does not add solid layers or cut thin roofs that cannot
accommodate the cap. Near stepped roofs, the protected nearer top takes priority
over a higher eligible top.

The default 0.8 mm relief depth, 1 mm tooth depth, 4 mm tooth pitch and 1 mm rim
are geometric test values, not calibrated shrinkage compensation. Zigzags do not
establish preserved strength, and the continuous rim and cap still transfer
stress. The cap may show depressions over channels. Neither relief mode is
automatically enabled by risk screening; print validation is required.
