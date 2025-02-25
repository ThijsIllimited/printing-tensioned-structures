G90 ; use absolute coordinates
M83 ; extruder relative mode
M104 S[first_layer_temperature] ; set extruder temp
M140 S[first_layer_bed_temperature] ; set bed temp
M190 S[first_layer_bed_temperature] ; wait for bed temp
M109 S[first_layer_temperature] ; wait for extruder temp
;M106 S100 ; Turn on fan
M900 K[K-factor] ; set K-factor
G28 W ; home all without mesh bed level
G80 ; mesh bed leveling