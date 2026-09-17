# INDX / RME custom G-code templates

These templates preserve the custom Cygnas X1 INDX startup and tool-change settings.
They require the corresponding RME firmware commands (including `M976` and `M870`)
and a calibrated INDX cleaner. They are not generic Prusa printer defaults.

Startup requires RME firmware supporting the explicit physical-tool list in
`G427 T... R2 P3`. The list comes from `is_extruder_used`: a job using tools 0
and 7 emits `G427 T0,7 R2 P3`. This template assumes Orca tool IDs correspond
directly to physical INDX tool IDs; do not use that list unchanged with remapping.
Older firmware without this parameter may still calibrate all tools.

- `indx_rme_start.gcode`: printer preset **Machine start G-code**.
- `indx_rme_toolchange.gcode`: printer preset **Change filament G-code**.

These files do not overwrite saved or cloud printer presets. Copy each template
into its corresponding preset field, save the preset, and re-slice the project.

Bucket priming replaces a bed purge line. After purging, each sequence retracts,
waits for motion completion, wipes with `G12 S1`, waits again, ejects the new
pellet with `G12 S30`, and exits with `G12 S91`. Exit alone is not a wipe.
The initial leftover-pellet ejection also follows a quick wipe.

Startup homes XY once, before vent control and tool pickup, then homes Z after
tool pickup. Do not add another XY or all-axis home later in the template.

Auto-PA (`M976`) runs after homing and before tool-offset calibration (`G427`).
The bed heats during PA; the template waits for its requested temperature
before offset calibration and restores the initial print tool only after both
calibration stages. No startup command turns off the bed. Optional chamber
warm-up still raises the bed target and then restores the printing target;
its `M104 S0` turns off the **nozzle**, not the bed, and is retained.
Heater shutdown at the end of the print must remain in the end G-code.

Do not insert `[filament_start_gcode]` into the change-filament template: Orca
already evaluates and emits the selected filament's start G-code after the tool
change. Raw substitution duplicates it and can leave expressions in the output.

The `[INDX]` tests exercise startup with single, nonconsecutive, and all tools and tool changes with the
wipe tower enabled/disabled. Physical pellet release and nozzle priming still
need verification on the printer. Existing exported G-code must be re-sliced.
