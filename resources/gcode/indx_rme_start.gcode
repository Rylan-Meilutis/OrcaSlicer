; CORE One INDX start G-code for OrcaSlicer + RME firmware
; Based on the current Prusa INDX startup with RME CORE One enhancements

M75 ; start print timer / serial print mode

{
global retract_toolchange = 8;
}

; number of used tools
{global used_tools = 0}
{if is_extruder_used[0]}{used_tools = used_tools + 1}{endif}
{if is_extruder_used[1]}{used_tools = used_tools + 1}{endif}
{if is_extruder_used[2]}{used_tools = used_tools + 1}{endif}
{if is_extruder_used[3]}{used_tools = used_tools + 1}{endif}
{if is_extruder_used[4]}{used_tools = used_tools + 1}{endif}
{if is_extruder_used[5]}{used_tools = used_tools + 1}{endif}
{if is_extruder_used[6]}{used_tools = used_tools + 1}{endif}
{if is_extruder_used[7]}{used_tools = used_tools + 1}{endif}

{local probe_temp = "" + ((filament_notes[initial_tool]=~/.*MBL160.*/) ? 160 : (filament_notes[initial_tool]=~/.*HT_MBL10.*/) ? (nozzle_temperature_initial_layer[initial_tool] - 10) : (filament_type[0] == "PC" or filament_type[0] == "PA") ? (nozzle_temperature_initial_layer[0] - 45) : (filament_type[initial_tool] == "FLEX") ? 170 : (filament_type[initial_tool]=~/.*PET.*/) ? 170 : 170) + "";}

M17 ; enable steppers
M151 W255 D0 T500 ; initialize printer lighting/status LEDs
M862.3 P "COREONEINDX" ; printer model check
M862.5 P2 ; g-code level check
M862.6 P"Input shaper" ; FW feature check
M862.6 P"INDX lock" ; FW feature check
M115 U6.6.0+15528

M591 S0 ; disable stuck detection

M555 X{(min(print_bed_max[0], first_layer_print_min[0] + 32) - 32)} Y{(max(0, first_layer_print_min[1]) - 4)} W{((min(print_bed_max[0], max(first_layer_print_min[0] + 32, first_layer_print_max[0])))) - ((min(print_bed_max[0], first_layer_print_min[0] + 32) - 32))} H{((first_layer_print_max[1])) - ((max(0, first_layer_print_min[1]) - 4))}

; inform about nozzle diameter
{if (is_extruder_used[0])}M862.1 T0 P{nozzle_diameter[0]}{endif}
{if (is_extruder_used[1])}M862.1 T1 P{nozzle_diameter[1]}{endif}
{if (is_extruder_used[2])}M862.1 T2 P{nozzle_diameter[2]}{endif}
{if (is_extruder_used[3])}M862.1 T3 P{nozzle_diameter[3]}{endif}
{if (is_extruder_used[4])}M862.1 T4 P{nozzle_diameter[4]}{endif}
{if (is_extruder_used[5])}M862.1 T5 P{nozzle_diameter[5]}{endif}
{if (is_extruder_used[6])}M862.1 T6 P{nozzle_diameter[6]}{endif}
{if (is_extruder_used[7])}M862.1 T7 P{nozzle_diameter[7]}{endif}

G90 ; use absolute coordinates
M83 ; extruder relative mode

; Home XY before vent control or tool pickup (once only)
G28 XY

; --- RME chamber vent control ---
{if chamber_temperature[initial_tool] > 35}
M870 C ; close top vent for heated-chamber materials
{else}
M870 O ; open top vent for low-temperature materials
{endif}

M140 S[first_layer_bed_temperature] ; set bed temp

{if used_tools > 1 and (is_extruder_used[0] and filament_type[0]=="PLA") or (is_extruder_used[1] and filament_type[1]=="PLA") or (is_extruder_used[2] and filament_type[2]=="PLA") or (is_extruder_used[3] and filament_type[3]=="PLA") or (is_extruder_used[4] and filament_type[4]=="PLA") or (is_extruder_used[5] and filament_type[5]=="PLA") or (is_extruder_used[6] and filament_type[6]=="PLA") or (is_extruder_used[7] and filament_type[7]=="PLA")}
M141 S20 ; set nominal chamber temp
{elsif chamber_minimal_temperature[initial_tool]==0}
M141 S{chamber_temperature[initial_tool]} ; set nominal chamber temp
{endif}

; Pick tool for Z homing
{if is_extruder_used[0]}T0 S1 L2 D0
{elsif is_extruder_used[1]}T1 S1 L2 D0
{elsif is_extruder_used[2]}T2 S1 L2 D0
{elsif is_extruder_used[3]}T3 S1 L2 D0
{elsif is_extruder_used[4]}T4 S1 L2 D0
{elsif is_extruder_used[5]}T5 S1 L2 D0
{elsif is_extruder_used[6]}T6 S1 L2 D0
{elsif is_extruder_used[7]}T7 S1 L2 D0
{endif}

M104 S120

; Home Z
G28 Z

G0 Z40 F10000
; try picking tools used in print
G1 F10000

{if (is_extruder_used[0] and filament_type[0] != "FLEX")}M574 S0 V35 T{nozzle_temperature[0]} F{(filament_max_volumetric_speed[0] > 0 ? min(outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[0]/2) * (filament_diameter[0]/2)), filament_max_volumetric_speed[0] / (3.14159 * (filament_diameter[0]/2) * (filament_diameter[0]/2))) : outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[0]/2) * (filament_diameter[0]/2)))}{endif}
{if (is_extruder_used[1] and filament_type[1] != "FLEX")}M574 S1 V35 T{nozzle_temperature[1]} F{(filament_max_volumetric_speed[1] > 0 ? min(outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[1]/2) * (filament_diameter[1]/2)), filament_max_volumetric_speed[1] / (3.14159 * (filament_diameter[1]/2) * (filament_diameter[1]/2))) : outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[1]/2) * (filament_diameter[1]/2)))}{endif}
{if (is_extruder_used[2] and filament_type[2] != "FLEX")}M574 S2 V35 T{nozzle_temperature[2]} F{(filament_max_volumetric_speed[2] > 0 ? min(outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[2]/2) * (filament_diameter[2]/2)), filament_max_volumetric_speed[2] / (3.14159 * (filament_diameter[2]/2) * (filament_diameter[2]/2))) : outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[2]/2) * (filament_diameter[2]/2)))}{endif}
{if (is_extruder_used[3] and filament_type[3] != "FLEX")}M574 S3 V35 T{nozzle_temperature[3]} F{(filament_max_volumetric_speed[3] > 0 ? min(outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[3]/2) * (filament_diameter[3]/2)), filament_max_volumetric_speed[3] / (3.14159 * (filament_diameter[3]/2) * (filament_diameter[3]/2))) : outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[3]/2) * (filament_diameter[3]/2)))}{endif}
{if (is_extruder_used[4] and filament_type[4] != "FLEX")}M574 S4 V35 T{nozzle_temperature[4]} F{(filament_max_volumetric_speed[4] > 0 ? min(outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[4]/2) * (filament_diameter[4]/2)), filament_max_volumetric_speed[4] / (3.14159 * (filament_diameter[4]/2) * (filament_diameter[4]/2))) : outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[4]/2) * (filament_diameter[4]/2)))}{endif}
{if (is_extruder_used[5] and filament_type[5] != "FLEX")}M574 S5 V35 T{nozzle_temperature[5]} F{(filament_max_volumetric_speed[5] > 0 ? min(outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[5]/2) * (filament_diameter[5]/2)), filament_max_volumetric_speed[5] / (3.14159 * (filament_diameter[5]/2) * (filament_diameter[5]/2))) : outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[5]/2) * (filament_diameter[5]/2)))}{endif}
{if (is_extruder_used[6] and filament_type[6] != "FLEX")}M574 S6 V35 T{nozzle_temperature[6]} F{(filament_max_volumetric_speed[6] > 0 ? min(outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[6]/2) * (filament_diameter[6]/2)), filament_max_volumetric_speed[6] / (3.14159 * (filament_diameter[6]/2) * (filament_diameter[6]/2))) : outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[6]/2) * (filament_diameter[6]/2)))}{endif}
{if (is_extruder_used[7] and filament_type[7] != "FLEX")}M574 S7 V35 T{nozzle_temperature[7]} F{(filament_max_volumetric_speed[7] > 0 ? min(outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[7]/2) * (filament_diameter[7]/2)), filament_max_volumetric_speed[7] / (3.14159 * (filament_diameter[7]/2) * (filament_diameter[7]/2))) : outer_wall_speed * (layer_height * (outer_wall_line_width - layer_height) + 3.14159 * (layer_height/2) * (layer_height/2)) / (3.14159 * (filament_diameter[7]/2) * (filament_diameter[7]/2)))}{endif}

; --- RME automatic pressure-advance calibration ---
; Run after native INDX homing and tool pickup, before G427 offset calibration.
; Format: physical INDX tool : logical Orca tool : filament type : temperature.
; Only tools actually used by the sliced job are included.
M976 A {if is_extruder_used[0]}0:0:{filament_type[0]}:{nozzle_temperature_initial_layer[0]}{endif}{if is_extruder_used[1]}{if is_extruder_used[0]},{endif}1:1:{filament_type[1]}:{nozzle_temperature_initial_layer[1]}{endif}{if is_extruder_used[2]}{if is_extruder_used[0] or is_extruder_used[1]},{endif}2:2:{filament_type[2]}:{nozzle_temperature_initial_layer[2]}{endif}{if is_extruder_used[3]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2]},{endif}3:3:{filament_type[3]}:{nozzle_temperature_initial_layer[3]}{endif}{if is_extruder_used[4]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2] or is_extruder_used[3]},{endif}4:4:{filament_type[4]}:{nozzle_temperature_initial_layer[4]}{endif}{if is_extruder_used[5]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2] or is_extruder_used[3] or is_extruder_used[4]},{endif}5:5:{filament_type[5]}:{nozzle_temperature_initial_layer[5]}{endif}{if is_extruder_used[6]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2] or is_extruder_used[3] or is_extruder_used[4] or is_extruder_used[5]},{endif}6:6:{filament_type[6]}:{nozzle_temperature_initial_layer[6]}{endif}{if is_extruder_used[7]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2] or is_extruder_used[3] or is_extruder_used[4] or is_extruder_used[5] or is_extruder_used[6]},{endif}7:7:{filament_type[7]}:{nozzle_temperature_initial_layer[7]}{endif}


; Keep the bed heating throughout PA; offsets require the requested bed temperature.
M190 R[first_layer_bed_temperature]
G427 T{if is_extruder_used[0]}0{endif}{if is_extruder_used[1]}{if is_extruder_used[0]},{endif}1{endif}{if is_extruder_used[2]}{if is_extruder_used[0] or is_extruder_used[1]},{endif}2{endif}{if is_extruder_used[3]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2]},{endif}3{endif}{if is_extruder_used[4]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2] or is_extruder_used[3]},{endif}4{endif}{if is_extruder_used[5]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2] or is_extruder_used[3] or is_extruder_used[4]},{endif}5{endif}{if is_extruder_used[6]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2] or is_extruder_used[3] or is_extruder_used[4] or is_extruder_used[5]},{endif}6{endif}{if is_extruder_used[7]}{if is_extruder_used[0] or is_extruder_used[1] or is_extruder_used[2] or is_extruder_used[3] or is_extruder_used[4] or is_extruder_used[5] or is_extruder_used[6]},{endif}7{endif} R2 P3 ; Calibrate only the used physical tools

; Restore the print tool once, after both calibration stages.
T{initial_tool} S1 L0 D0

{if chamber_minimal_temperature[initial_tool]!=0 and ! (used_tools > 1 and (is_extruder_used[0] and filament_type[0]=="PLA") or (is_extruder_used[1] and filament_type[1]=="PLA") or (is_extruder_used[2] and filament_type[2]=="PLA") or (is_extruder_used[3] and filament_type[3]=="PLA") or (is_extruder_used[4] and filament_type[4]=="PLA") or (is_extruder_used[5] and filament_type[5]=="PLA") or (is_extruder_used[6] and filament_type[6]=="PLA") or (is_extruder_used[7] and filament_type[7]=="PLA"))}
; Min chamber temp section
M104 S0
M140 S115 ; set bed temp for chamber heating
G1 Z10 F720 ; set bed position
G1 X242 Y0 F4800 ; set print head position
M191 S{chamber_minimal_temperature[initial_tool]} ; wait for minimal chamber temp
M141 S{chamber_temperature[initial_tool]} ; set nominal chamber temp
M140 S[first_layer_bed_temperature] ; set bed temp
{endif}

{if first_layer_bed_temperature[initial_tool]<=60}M106 S70{endif}
G0 Z40 F10000
M104 S{if is_nil(idle_temperature[initial_tool])}100{else}{idle_temperature[initial_tool]}{endif}
M190 R[first_layer_bed_temperature] ; wait for bed temp
M107

G29 G ; absorb heat

M109 S{probe_temp} ; wait for temp

;
; MBL
;
M84 E ; turn off E motor
G29 P1 ; invalidate mbl & probe print area
;G29 P1 X150 Y0 W100 H20 C ; probe near purge place
G29 P3.2 ; interpolate mbl probes
G29 P3.13 ; extrapolate mbl outside probe area
G29 A ; activate mbl

G0 Z10 ; add Z clearance

{if used_tools == 1 or (used_tools > 1 and enable_prime_tower == true)}
    M104 S{nozzle_temperature_initial_layer[initial_tool]}
    G12 S90 ; enter cleaner
    G12 S1 ; wipe before clearing a leftover pellet
    M400 ; finish wiping before ejecting
    G12 S30 ; eject leftover pellet
    M106 S{255/100*fan_max_speed[initial_tool]}
    M109 S{nozzle_temperature_initial_layer[initial_tool]}
    M906 P1 ; Set extruder current
    M83
    G1 E{retract_toolchange+retraction_length[initial_tool]} F{deretraction_speed[initial_tool] * 60}
    M400
    M201 E{default_acceleration*line_width*layer_height/(3.14 * filament_diameter[initial_tool] * filament_diameter[initial_tool] / 4)}
    G1 E{0.8*filament_minimal_purge_on_wipe_tower[initial_tool] / (3.14 * filament_diameter[initial_tool] * filament_diameter[initial_tool] / 4)} F{0.8*60*filament_max_volumetric_speed[initial_tool]/3.14/filament_diameter[initial_tool]}
    M400
    G1 E{0.2*filament_minimal_purge_on_wipe_tower[initial_tool] / (3.14 * filament_diameter[initial_tool] * filament_diameter[initial_tool] / 4)} F{min(250, 60*filament_max_volumetric_speed[initial_tool]/3.14/filament_diameter[initial_tool]/filament_diameter[initial_tool]*4)}
    M400
    M201 E[machine_max_acceleration_e]
    G4 P20
    G1 E-{1.5*retraction_length[initial_tool]} F{retraction_speed[initial_tool] * 60}
    {e_retracted[initial_tool] = 1.5*retraction_length[initial_tool]}
    G4 P20
    M906 P0 ; Revert extruder current
    G1 F{min(travel_speed, 350.0) * 60}
    M400 ; finish retraction before wiping
    G12 S1 ; wipe the nozzle before ejecting the pellet
    M400 ; finish wiping before ejecting
    G12 S30 ; eject the newly purged pellet
    G12 S91 ; Exit cleaning station
{endif}

M569 S0 E ; set spreadcycle mode for extruder
M591 S0 ; disable stuck detection
G92 E0 ; reset extruder position
M221 S100 ; set flow to 100%
