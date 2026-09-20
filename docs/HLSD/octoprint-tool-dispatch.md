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
suggests distinct tools by material and color; the user can change every
assignment. Missing inventory permits explicitly verified manual assignment.
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
successful inventory read. The upstream mixed-filament add controls are hidden
for independent-tool OctoPrint dispatch, while existing mixed entries remain
editable for project compatibility. Sync from machine compacts loaded materials
into the logical palette rather than preserving holes for unloaded tools. For
an existing model it asks before replacing the palette; removed tail entries
use the existing filament deletion/remapping path and fall back to material 1.
An empty or failed inventory does not erase the project's palette.

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
is cancelled. It re-slices the job and generates a separate temporary G-code
file; the editable project, ordinary preview and its memory-mapped G-code are
not replaced. Cancellation reaches both prints, and cancelled jobs are checked
before upload enqueueing. Post-processing receives the mapped configuration;
filename statistics come from the mapped print. Connection credentials are
not copied into the exported job configuration.

Physical-index spool metadata is embedded when spool-name embedding is enabled.
Sliced 3MF packaging uses the mapped print and its G-code processor result,
including used physical tool IDs and material metadata. It does not package
the old logical preview file. Upload receives its own temporary copy so GUI
cleanup cannot remove a queued job's input.

Regression coverage uses `[ToolMapping]` in the core and FFF suites. It checks
project preservation, painting and role remapping, invalid assignments,
T0/T7 commands, used-tool masks, temperature ordering, purge matrices and the
G-code/filament metadata inside a sliced 3MF.

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
