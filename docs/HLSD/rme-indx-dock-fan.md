# RME INDX dock-fan cooldown protocol

This is an Orca-to-firmware contract. It requires a matching firmware implementation; the inspected RME firmware does not implement it. Install that implementation before using the enabled RME INDX preset. An older firmware may report or ignore the unknown command and will not provide automatic cooldown.

## Commands

`M979 S1 C180` enables firmware-managed dock-fan cooling with a minimum duration of 180 seconds after the latest potentially hot tool parking event. `C` is an integer duration in seconds, valid from 30 through 1800. The command is emitted once before machine start G-code, not repeatedly at layer changes. `M979 S0` returns to manual management; it must not itself force the fan off or bypass a pending thermal safety requirement. Other parameter combinations must be rejected without changing state.

M979 is a dedicated proposed RME command, absent from the inspected firmware command table. It is not an extension of M106: existing `M106 P6 R` does not provide automatic dock cooling and must not be used as a substitute. The S and C operands contain no tool identifiers, so physical-tool remapping leaves the command unchanged.

## Firmware responsibilities

- Activate the dock fan **before** depositing a potentially hot nozzle into its dock. Include successful parking through P0, implicit tool changes, startup operations and recovery paths; do not rely on slicer-generated toolchange comments.
- Restart the deadline on every potentially hot parking event. A recently parked tool must not inherit another tool's expired deadline. Use monotonic, wrap-safe real wall-clock time, not estimated G-code duration or layer count.
- Stop the fan only after both the minimum duration and a validated thermal safety condition have been met. Three minutes is not proof that a nozzle is cold. Unknown temperature, sensor failure or unavailable dock-temperature measurement must retain cooling unless a separately validated conservative cooling model establishes safety.
- Once the dock is safe, keep the fan off throughout a long gap, and restart it before the next hot parking event. The controller runs asynchronously and must not stall printing or inject a dwell.
- Continue pending cooldown after the final tool change, print completion, cancellation and idle. Print-end cleanup must not unconditionally clear the fan. While automatic cooling has pending safety work, legacy `M106 P6 S0` and print-end fan-stop calls must not defeat it. Fault/emergency handling and independent hardware protection retain priority.
- Define startup/reboot behavior conservatively for unknown hot dock contents. Re-enabling automatic mode must not erase an outstanding safety deadline. Manual mode selection must not be treated as an immediate safe-to-stop assertion.
- Make state transitions testable: cold idle, hot park, repeated parks, time rollover, long printing gaps, last park, end/cancel, manual selection, invalid arguments and unavailable thermal data. Bench validation of the cooling condition remains a firmware/hardware responsibility.

## Orca settings

The printer capability `support_rme_indx_dock_fan_cooldown` defaults to false globally. When false, Orca emits no M979 commands. With the capability enabled, `rme_indx_dock_fan_cooldown` selects S1 or S0, and `rme_indx_dock_fan_cooldown_time` supplies the minimum duration (default 180 seconds). These controls are in Printer Settings → Basic information → Accessory.

The RME INDX system preset enables the capability and automatic mode by default; stock and other printer presets remain unchanged. Existing customized presets and saved projects can override inherited values and must be checked explicitly. Changing these options regenerates G-code, not slicing geometry. No slicer thermal estimate is used to turn the fan off, and no printer commands are sent directly by editing the settings.
