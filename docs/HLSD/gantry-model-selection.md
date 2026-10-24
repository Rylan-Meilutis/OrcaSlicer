# Gantry model selection and sequential spacing

The printer's Gantry model selector offers automatic printer matching and the
bundled entries in `resources/data/printer_gantries/geometries.json`. Explicit
built-in selections store `builtin:<filename>` and load the same curated slices
as automatic matching. They do not reinterpret visualization meshes.

Import gantry model accepts a nozzle-centred STL and validates its collision
bands before changing the profile. Saving a printer profile copies an external
mesh into `<profile>.assets/gantry-<unique-id>.stl` beside the profile JSON.
Subsequent saves reuse that file. Save As creates an independent asset for the
new profile. The original imported file is never moved or overwritten. Built-in
models do not require copies. Legacy custom JSON geometry remains readable,
but selecting or importing a model clears that old override; the raw JSON
editor is no longer exposed.

The selector persists enum values (`builtin:<filename>`), not its translated
display labels. Selecting a model clears legacy JSON before the change is
published. Gantry keys and printer notes invalidate the plater clearance cache;
Arrange reads the edited printer preset so changes also apply before slicing.
Each catalog entry supplies its own curated collision bands; an imported STL
supplies measured bands. Neither path copies INDX clearance constants.

The stored custom-model path is absolute. Moving profiles to another computer
therefore requires copying the assets and reselecting the model; cloud asset
transfer is not implemented by this local-file mechanism.

Sequential arrangement and validation use separate X/Y nozzle-to-head reaches
from the collision model, with half the resulting rectangular envelope on each
object. The legacy radius is a fallback only and is hidden for modeled printers;
its stored value is retained for compatibility. Skirts and explicit arrange
spacing remain additional allowances. The INDX model therefore reserves 28 mm
sideways and 90 mm fore/aft, rather than spreading the rear reach sideways.
This is still conservative in print order: it reserves the larger reach on
both sides. Rod/lid constraints remain separate. Modeled arrangement preserves
existing object rotations, because rotating an expanded envelope would also
rotate the machine geometry. Unmodeled printers retain the existing rotation
and circular-clearance behavior.

Project print order offers Use profile, By layer, and By object, applied to all
existing plates. Overrides use the established per-plate configuration stored
in 3MF projects, without modifying the process preset. Projects whose plates
have different overrides display Custom per plate. New plates retain the
existing plate-creation behavior; the control does not alter printer defaults.
