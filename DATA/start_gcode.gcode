G90 ; use absolute coordinates
M83 ; extruder relative mode

M104 S[first_layer_temperature] ; set extruder temp
M140 S[first_layer_bed_temperature] ; set bed temp
M190 S[first_layer_bed_temperature] ; wait for bed temp

M900 K[K-factor] ; set K-factor

G28 W ; home all without bed mesh level

; Prusa XL-specific pre-print and leveling routine
G85 ; un-tare the load cell (if it was tared previously)
G0 X10 Y-3 Z0.2 F10000 ; move a bit away from the bed
G0 Z10 F1000 ; raise nozzle for cleaning
M109 S[first_layer_temperature] ; wait for extruder temp

; Nozzle cleaning routine
G0 X60 Y-3 Z0.2 F10000 ; move to nozzle cleaning position
G1 E15 F100 ; extrude 15mm
G1 X120 E15 F1000 ; move to nozzle cleaning
G0 X60 Z0.2 F10000 ; retract and move
G1 E15 F100 ; extrude 15mm
G1 X120 E15 F1000 ; move to nozzle cleaning

G0 Z10 F1000 ; raise nozzle for clearing
G28 W ; home without mesh bed leveling

G80 ; mesh bed leveling

M420 S1 Z0.2 ; activate the mesh bed leveling with a fade height of 0.2mm