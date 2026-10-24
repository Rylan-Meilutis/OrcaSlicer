# OctoPrint tool dispatch

For non-Bambu printers with independently selectable tools, project filament
slots describe materials, not fixed physical tool numbers. Ordinary slicing
and preview use the logical slots. Print/Send, G-code export and sliced 3MF
export prompt for a one-to-one assignment to physical tools. Single-extruder
multi-material printers and Bambu dispatch retain their existing behavior.

Independent-tool machines expose add/remove material controls regardless of
their print host. Deleting a material does not change nozzle count, diameters
or offsets, and preset refresh does not repopulate the palette to hardware
capacity. OctoPrint supplies automatic loaded-spool suggestions; other
connections use explicit manual assignments with an inventory warning. This
path maps independent T-tools, not sub-nozzles within a single toolhead.

## Different nozzle sizes by feature

Keep a separate logical filament slot for each nozzle/material combination;
slots may use exactly the same filament profile and color. Set the desired
feature's filament option to Per print, then select its slot under Filament
Bindings. For example, bind outer walls to slot 1, inner walls to slot 2 and
infill to slot 3; at dispatch map those slots to loaded 0.2, 0.4 and 0.6 mm tools.
Support base and interface have independent bindings too.

The mapping dialog shows configured nozzle diameters alongside loaded material
and color, and lists default feature bindings beneath each project material.
Diameters come from the machine profile, not from a guess based on spool data.
The private job remaps feature assignments before slicing. Standard region and
support flow calculations therefore use the destination nozzle. Auto widths
scale with nozzle size; explicit widths are deliberately preserved. The shared
layer height must be valid for all tools in use. Tests check actual exported
tool IDs and extrusion widths, as well as support base/interface assignments.

The dialog opens before the upload destination or export file picker and reads
the selected physical printer's loaded spool inventory afresh. Only materials
in the selected print's used-extruder set are shown (including support and
explicit layer tool changes), not the whole project palette. Unused entries
receive private, distinct filler assignments without generating tool use. It
suggests distinct tools by material and color, with nozzle-size preference when
the project has a manual nozzle assignment from sync; the user can change every
assignment. Other hosts without automatic inventory permit explicitly verified
manual assignment. OctoPrint requires a live or last-known loaded-tool snapshot.
Known material mismatches require confirmation. Each dispatch prompts again,
so changing loaded spools cannot silently reuse an old assignment. Multi-plate
jobs are dispatched separately. Virtual mixed-material slots are not supported
by this one-material/one-tool mapping.

The dialog presents project material swatches alongside loaded-tool selectors,
with live per-row messages for matching materials, incompatible materials,
unloaded tools and duplicate assignments. These messages do not rely on color
alone. Filament swatches retain their actual colors in dark mode. The material
list scrolls within a bounded viewport; the confirmation controls stay outside
the scrolling area.

Empty physical slots display as Unloaded and cannot be confirmed after a
live or cached inventory read. The upstream mixed-filament add controls are hidden
for independent-tool OctoPrint dispatch, while existing mixed entries remain
editable for project compatibility. Sync from machine first presents a checklist
of loaded tools using their original physical numbers, configured diameters,
names and colors. The
checked count determines the palette size (at least one and no more than printer
capacity). Unloaded tools are omitted. Cancelling this selection changes neither
the project nor saved profile mappings. Sync compacts only the selected materials
into the logical palette rather than preserving holes for unloaded tools. For
an existing model it asks before replacing the palette; removed tail entries
use the existing filament deletion/remapping path and fall back to material 1.
An empty or failed inventory does not erase the project's palette.
Printer-settings callbacks, including those invoked while loading a preset,
do not add or remove project filaments to match nozzle count or color capacity.
Sync retains selected tools in the existing manual `filament_map`, preserving
nozzle/material pairs even when only non-adjacent tools are selected. Duplicate
profiles or colors remain separate slots for different nozzle sizes. Nozzle
diameters stay in the printer profile; sync does not overwrite that geometry.
Dispatch refreshes inventory and confirms assignments rather than treating this
preference as permission to print on an unchecked tool.

## Loaded-spool snapshots

Print/Send refreshes OctoPrint assignments before asking for the upload location,
independently of the manual Sync and metadata-embedding options. Independent
tools use these assignments in the mapping dialog; single-nozzle jobs validate
their used slot indices against the snapshot without replacing project profiles.
The confirmed snapshot is reused for metadata rather than queried again after
destination selection.

Successful reads (including manual Sync) are saved immediately in AppConfig's
`octoprint_loaded_spools` section. Versioned snapshots retain physical indices,
unloaded slots, spool metadata and a timestamp. Cache keys hash the host URL and
provider endpoint, so a different printer/provider does not inherit assignments
and credentials are not stored in cache keys. A successful empty assignment
replaces the old snapshot; an inventory-only or malformed response does not.

A failed refresh can use the persisted snapshot after restarting Orca. The
mapping dialog warns that these are last-known materials and shows the snapshot
timestamp; single-nozzle dispatch requires acknowledgment of that warning.
Known-unloaded tools remain invalid even when offline. Without usable live or
cached assignments, dispatch stops with an error rather than inventing a tool
mapping. The entire provider-probing sequence has a 15-second request budget.
Manual Sync still requires a live reading, so it cannot silently replace a
project palette with stale inventory. `[OctoPrintSpoolCache]` tests cover snapshot
round trips, empty assignments, corrupt data and host/provider isolation.

`PresetBundle::tool_mapped_config` composes an independent job from the current
presets and plate overrides. It remaps colors, material profiles, purge tables,
resolved area bindings, model/volume/range assignments, paint and layer-slider
tool changes. An implicit object material is made explicit before remapping.
Destination nozzle geometry stays with the printer tool. Profile variants are
resolved through the usual `Print::apply` pipeline, not by copying flattened
temperature arrays or rewriting finished G-code. Literal custom G-code is not
rewritten: templates should use the slicer's tool placeholders when they mean
the currently selected material, and literal tool numbers when they mean a
specific physical tool.

`BackgroundSlicingProcess` owns the mapped `Print` until dispatch completes or
is cancelled. An identity mapping with an identical complete effective
configuration reuses the sliced G-code through a separate working copy. Unused
physical-slot material padding is excluded from that comparison. Changed plain
G-code assignments can reuse motion when their effective material rows, nozzle
properties, feature roles, purge tables and sequences are equivalent. Unknown
configuration differences require regeneration.

The generator captures custom-template inputs, runtime overrides, output state,
globals and random-generator state. Immutable parser settings are shared between
snapshots; retained data is capped at 16 MiB and 2048 expansions. Marked template
spans and a fingerprint of the final generated file establish provenance. Dispatch
requires a fresh slice from this generator to obtain that in-memory certificate;
an older exported G-code file alone cannot establish script provenance. Dispatch
first performs a read-only preflight: it reproduces each original expansion, then
replays it with physical tool IDs and relocated arrays. Literal custom parameters
are never globally replaced. Changed macro motion or writer/global state rejects
reuse. Replay rebuilds effective retraction settings from the destination nozzle
defaults and material overrides, preserving the parser's extrusion-state outputs.
The exported derived nozzle-map metadata describes the physical slots even though
the destination Print has not been sliced. Replay preserves the source index of implicit legacy machine-limit vectors
(normal/silent modes), independently of the relocated physical tool. Explicit
machine-limit indices and ordinary generation keep their existing semantics.
The motion check covers recognized positioning, cleaner/parking, coordinate,
extrusion, acceleration and driver-control commands; it is not a complete
firmware-state interpreter. Firmware calibration commands such as RME's `G427`
and `M976` are evaluated from their captured templates with the mapped tool lists.
Outside those spans an explicit native-command dialect identifies tool
operands; an unknown command rejects reuse rather than guessing what its parameters
mean. The copied motion remains byte-identical. Physical configuration, filament
lists and usage statistics are updated separately.

Within a template span, dispatch replaces only lines that changed during replay,
not the entire rendered template. Original and replayed line counts must match;
every changed line must occur exactly as often in the sealed final span as in the
original expansion. Replacements follow original occurrence order in a single
pass, so a newly written tool ID is never remapped again. Missing, rewritten or
additional matching lines reject reuse. All other final bytes are preserved,
including relocated fan commands, split moves, timing updates and computed mesh
bounds inserted by the generator's post-processing.

The preflight edit plan is capped at 32 MiB. A bounded streaming pass writes one
private destination and verifies the source fingerprint again before the output
can be published. Disk-space checks reserve system headroom and account for the
additional spool-metadata rewrite. IMEX, SEMM, prime towers, nozzle preheating, logical nozzle-group
templates, non-Marlin dialects, external processors and changed-tool sliced-3MF
packaging use regeneration. Rejected preflight asks for Regenerate or Cancel on
the main thread; it does not silently start a second slice. Different nozzle
diameters receive a specific explanation before regeneration.

Regeneration produces a separate temporary
G-code file; the editable project, ordinary preview and its memory-mapped G-code are
not replaced. Cancellation reaches both prints, and cancelled jobs are checked
before upload enqueueing. Post-processing receives the mapped configuration;
filename statistics come from the mapped print when regenerated and the original
print when reusing its motion. Connection credentials are
not copied into the exported job configuration. Progress identifies mapped-job
slicing, G-code generation, spool metadata, packaging and copying. Script and
plugin names are reported when those processors actually run.

Physical-index spool metadata is embedded when spool-name embedding is enabled.
Sliced 3MF packaging uses the mapped print and its G-code processor result,
including used physical tool IDs and material metadata. It does not package
the old logical preview file for changed mappings. Identity reuse packages the
unchanged print and its existing processor result. Upload receives its own temporary copy so GUI
cleanup cannot remove a queued job's input.

Regression coverage uses `[ToolMapping]` in the core and FFF suites. It checks
project preservation, painting and role remapping, invalid assignments,
T0/T7 commands, used-tool masks, temperature ordering, purge matrices and the
G-code/filament metadata inside a sliced 3MF.
The opt-in FFF tag `[.ToolMappingProject]` reads `ORCA_TOOL_MAPPING_PROJECT` and
checks an actual three-material/eight-tool project's applied configuration without
slicing or exporting the model.
Setting `ORCA_TOOL_MAPPING_SLICE=1` additionally creates one fresh logical slice,
certifies its real printer scripts, and verifies that physical remapping retains
the motion fingerprint and slicing timestamp without processing the mapped print.

The mapping editor separates brand/material rules from exact-roll overrides.
Profile choices are material-filtered and searchable; the current selection
stays pinned when filtered out, so searching cannot silently replace a rule.
The first-sync profile prompt also offers a searchable compatible-profile list.

The Objects-list instance context menu exposes Change filament. Assigning a
filament to only some copies first separates those copies using the existing
instance-to-object path, under the same undo snapshot. Object geometry, painting,
placement and settings are copied; sibling copies retain their assignments.
This changes the object's base filament, not an individual painted-color map.

Manual GUI verification required: select one of two instances, change its
filament, verify the sibling is unchanged, then undo and verify both instances
are linked again. In the mapping editor, search, change tabs, clear search and
cancel without changing saved rules; repeat with Save and reopen.
