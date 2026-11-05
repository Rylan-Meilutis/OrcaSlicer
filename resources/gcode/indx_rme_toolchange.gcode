{local speed_tc = min(travel_speed, 350.0) * 60}

  {if enable_prime_tower == false}
    ; Change Tool[previous_extruder] -> Tool[next_extruder] (layer [layer_num])

    {if is_nil(retract_length_toolchange[current_extruder])}
      G1 E-[retraction_length[current_extruder]] F{retraction_speed[current_extruder] * 60}
    {endif}

    G1 F{speed_tc}
    M204 S7000

    {if not (layer_num < 0 and next_extruder == initial_tool)}
      G27 W3 Z{z_hop[current_extruder]} P2 R{retract_toolchange} V{retraction_speed[current_extruder]} A{travel_slope[current_extruder]}

      P0 S1 L0 D0
    {endif}

    T{next_extruder} S1 L0 D0
    ; Orca emits the evaluated filament-start G-code after this tool change.

    {local target_temp = nozzle_temperature[next_extruder]}

    {if layer_z <= (initial_layer_print_height + 0.0001)}
      {target_temp = nozzle_temperature_initial_layer[next_extruder]}
    {endif}

    M104 S{target_temp}

    G12 S90 ; enter cleaner
    G12 S1 ; wipe before clearing a leftover pellet
    M400 ; finish wiping before ejecting
    G12 S30 ; eject leftover pellet

    M106 S{255 / 100 * fan_max_speed[next_extruder]}
    M109 S{target_temp}

    M906 P1 ; Set extruder current
    M83

    {if is_nil(retract_length_toolchange[next_extruder])}
      G1 E{retract_toolchange + retraction_length[next_extruder]} F{20 * 60}
    {else}
      G1 E{retract_toolchange + retract_length_toolchange[next_extruder]} F{20 * 60}
    {endif}

    M400

    M201 E{default_acceleration * line_width * layer_height / (3.14 * filament_diameter[next_extruder] * filament_diameter[next_extruder] / 4)}

    ;EXCLUDE_E_START
    G1 E3.4 F{0.6 * 60 * filament_max_volumetric_speed[next_extruder] / 3.14 / filament_diameter[next_extruder]}
    ;EXCLUDE_E_END

    ;FLUSH_START
    G91 ; use relative coordinates
    M83 ; extruder relative mode

    G1 E{(0.4 * filament_minimal_purge_on_wipe_tower[current_extruder]) / (3.14 * filament_diameter[next_extruder] * filament_diameter[next_extruder] / 4)} F{0.8 * 60 *
    filament_max_volumetric_speed[next_extruder] / 3.14 / filament_diameter[next_extruder]}

    G1 Y-1.5 E{(0.3 * filament_minimal_purge_on_wipe_tower[current_extruder]) / (3.14 * filament_diameter[next_extruder] * filament_diameter[next_extruder] / 4)} F{0.8 * 60 *
    filament_max_volumetric_speed[next_extruder] / 3.14 / filament_diameter[next_extruder]}

    G1 Y1.5 E{(0.3 * filament_minimal_purge_on_wipe_tower[current_extruder]) / (3.14 * filament_diameter[next_extruder] * filament_diameter[next_extruder] / 4)} F{min(250,
    60*filament_max_volumetric_speed[next_extruder]/3.14/filament_diameter[next_extruder]/filament_diameter[next_extruder]*4)}

    G90 ; use absolute coordinates
    M83 ; extruder relative mode
    ;FLUSH_END

    M201 E[machine_max_acceleration_e]
    M906 P0 ; Revert extruder current
    M400

    G1 E-{1.5 * retraction_length[next_extruder]} F{retraction_speed[next_extruder] * 60}
    {e_retracted[next_extruder] = 1.5 * retraction_length[next_extruder]}
    M400

    G1 F{speed_tc}

    G12 S1 ; wipe the nozzle before ejecting the pellet
    M400 ; finish wiping before ejecting
    G12 S30 ; eject the newly purged pellet
    G12 S91 ; Exit cleaning station

    G0 Z{layer_z + 1} ; Lift

    {if layer_z <= (initial_layer_print_height + (close_fan_the_first_x_layers[next_extruder] - 1) * layer_height + 0.0001)}
      M106 S0
    {endif}

    G1 F{speed_tc}

  {else}

    {if previous_extruder != -1}
      ; wipe tower enabled
      ; Change Tool[previous_extruder] -> Tool[next_extruder] (layer [layer_num])

      {if is_nil(retract_length_toolchange[current_extruder])}
        G1 E-{retraction_length[current_extruder]} F{retraction_speed[current_extruder] * 60}
      {endif}

      G1 F{speed_tc}
      M204 S7000
      M104 S0

      G27 W3 Z{z_hop[current_extruder]} P2 R{retract_toolchange} V{retraction_speed[current_extruder]} A{travel_slope[current_extruder]}

      P0 S1 L0 D0

      T{next_extruder} S1 L0 D0

      {local target_temp = nozzle_temperature[next_extruder]}

      {if layer_num == 0}
        {target_temp = nozzle_temperature_initial_layer[next_extruder]}
      {endif}

      M104 S{target_temp}

      G12 S90 ; enter cleaner

      G750 Y85 F21000 A
      G750 Y93 F21000 A
      G750 Y82 F21000 A
      G750 Y98.5 F21000 A
      G750 Y75 F21000 A
      G750 Y98.5 F21000 A

      M109 S{target_temp}

      {if is_nil(retract_length_toolchange[next_extruder])}
        G1 E{retract_toolchange + retraction_length[next_extruder]} F{20 * 60}
      {else}
        G1 E{retract_toolchange + retract_length_toolchange[next_extruder]} F{20 * 60}
      {endif}

      ;EXCLUDE_E_START
      G1 E1.2 F{0.6 * 60 * filament_max_volumetric_speed[next_extruder] / 3.14 / filament_diameter[next_extruder]}
      ;EXCLUDE_E_END

      G1 E-{1.5 * retraction_length[next_extruder]} F{retraction_speed[next_extruder] * 60}
      {e_retracted[next_extruder] = 1.5 * retraction_length[next_extruder]}
      M400
      G12 S1 ; wipe the nozzle before ejecting the pellet
      M400 ; finish wiping before ejecting
      G12 S30 ; eject the newly purged pellet
      G12 S91 ; exit cleaner

      G0 Z{layer_z + 0.5} ; Lift
      G1 F{speed_tc}
    {endif}

  {endif}
