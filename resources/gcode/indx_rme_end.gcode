G90
{if layer_z < max_print_height}
G1 Z{min(max_layer_z + z_offset + 10, max_print_height)} F720
{endif}

; turn off extruder heaters
{if is_extruder_used[0]}M104 T0 S0{endif}
{if is_extruder_used[1]}M104 T1 S0{endif}
{if is_extruder_used[2]}M104 T2 S0{endif}
{if is_extruder_used[3]}M104 T3 S0{endif}
{if is_extruder_used[4]}M104 T4 S0{endif}
{if is_extruder_used[5]}M104 T5 S0{endif}
{if is_extruder_used[6]}M104 T6 S0{endif}
{if is_extruder_used[7]}M104 T7 S0{endif}

M140 S0 ; turn off heatbed
M141 S0 ; disable chamber control
M107 ; turn off fan
M106 P6 S0 ; turn off dock fan
P0 S1 ; park tool before vent motion
M400 ; finish parking before closing the vent
M870 C ; close top vent after print

G1 Y0  F10200
G1 X242 Y205 F10200 ; park
M400
G1 Z{min(max_layer_z + z_offset + 200, max_print_height)} F720
M400 ; wait
M572 S0 ; reset PA
M221 S100 ; reset flow percentage
M84 X Y E ; disable motors
M151 W0 D0 T500 ; turn off printer lighting/status LEDs
M77 ; stop print timer
; max_layer_z = [max_layer_z]
