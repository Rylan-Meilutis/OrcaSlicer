# Preview nozzle/toolhead visibility

The expanded G-code legend has a **Show nozzle/toolhead** checkbox beneath
the view selector. It controls the existing model, including the full gantry
when that asset is selected. The application preference
`preview_show_tool_model` defaults to visible when absent; it is not a printer,
process, or project setting. Position readouts remain visible.

Manual GUI verification checklist (requires an interactive preview):

1. Open sliced G-code and move the toolpath slider to a visible location.
2. Uncheck the option: the model disappears, but paths and the position readout
   remain. Check it again: the model returns at the current tool position.
3. Repeat with the full gantry model selected, and in Filament and Line Type
   views. Folding the legend must not reset the choice.
4. Change layers, reload G-code, then restart the app normally: the visibility
   preference should remain unchanged.
5. With a fresh app configuration, the model should remain visible by default.

This control changes rendering only; no slicing or G-code data is modified.
