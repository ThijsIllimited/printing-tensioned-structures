G90 ; use absolute coordinates
M83 ; extruder relative mode
M104 S[first_layer_temperature] ; set extruder temp
M140 S[first_layer_bed_temperature] ; set bed temp
M190 S[first_layer_bed_temperature] ; wait for bed temp
M109 S[first_layer_temperature] ; wait for extruder temp
;M106 S100 ; Turn on fan
M900 K0.20 ; set K-factor
G28 W ; home all without mesh bed level
G80 ; mesh bed leveling
G1 Z0.3 F720
G1 Y-3 F1000 ; go outside print area
G92 E0
G1 X60 E6 F1000 ; intro line
G1 X100 E6 F1000 ; intro line
G92 E0
G0 Z1.200 E0.0 F30000 ; move to point

