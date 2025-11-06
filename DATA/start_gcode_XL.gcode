; external perimeters extrusion width = 0.45mm
; perimeters extrusion width = 0.45mm
; infill extrusion width = 0.45mm
; solid infill extrusion width = 0.45mm
; top infill extrusion width = 0.42mm
; first layer extrusion width = 0.50mm

M73 P0 R239
M73 Q0 S252
M201 X7000 Y7000 Z200 E2500 ; sets maximum accelerations, mm/sec^2
M203 X400 Y400 Z12 E100 ; sets maximum feedrates, mm / sec
M204 P4000 R1200 T5000 ; sets acceleration (P, T) and retract acceleration (R), mm/sec^2
M205 X8.00 Y8.00 Z2.00 E10.00 ; sets the jerk limits, mm/sec
M205 S0 T0 ; sets the minimum extruding and travel feed rate, mm/sec

M486 S0
M486 AFDM-flat-frame.STL
M486 S-1

;TYPE:Custom
M17 ; enable steppers
M862.3 P "XL" ; printer model check
M862.5 P2 ; g-code level check
M862.6 P"Input shaper" ; FW feature check
M115 U6.2.5+8912
G90 ; use absolute coordinates
M83 ; extruder relative mode
; set print area
M555 X65.5816 Y68.3251 W205.061 H205.061
; inform about nozzle diameter
M862.1 T0 P0.4 A0 F0

; turn off unused heaters
M104 T1 S0

M217 Z2 ; set toolchange z hop to 2mm, or zhop variable from slicer if higher
; set bed and extruder temp for MBL
M140 S[first_layer_bed_temperature] ; set bed temp
G0 Z5 ; add Z clearance
M104 T0 S[idle_temperature] ; set hotend temp to a lower temp for probing
M109 T0 S[idle_temperature] ; wait for hotend to reach idle temp

; Home XY
G28 XY
; try picking tools used in print
G1 F24000

; select tool that will be used to home & MBL
T0 S1 L0 D0
; home Z with MBL tool
M84 E ; turn off E motor
G28 Z
G0 Z5 ; add Z clearance

M104 T0 S[idle_temperature] ; set idle temp
M190 S[first_layer_bed_temperature] ; wait for bed temp

;G29 G ; absorb heat

; MBL
M84 E ; turn off E motor
G29 P1 ; invalidate mbl & probe print area
;G29 P1 X30 Y0 W50 H20 C ; probe near purge place
G29 P3.2 ; interpolate mbl probes
G29 P3.13 ; extrapolate mbl outside probe area
G29 A ; activate mbl
G1 Z10 F720 ; move away in Z
G1 F24000
P0 S1 L1 D0; park the tool
; set extruder temp
M104 T0 S[first_layer_temperature]
M109 T0 S[first_layer_temperature] ; wait for hotend to reach final temp

;
; purge initial tool
;
G1 F24000
P0 S1 L2 D0; park the tool
M109 T0 S[first_layer_temperature]
T0 S1 L0 D0; pick the tool
G92 E0 ; reset extruder position

G0 X30 Y-7 Z10 F24000 ; move close to the sheet's edge
G0 E30 X40 Z0.2 F170 ; purge while moving towards the sheet
G0 X70 E9 F800 ; continue purging and wipe the nozzle
G0 X73 Z0.05 F8000 ; wipe, move close to the bed
G0 X76 Z0.2 F8000 ; wipe, move quickly away from the bed
G1 E-1.5 F2400 ; retract

G92 E0 ; reset extruder position
M900 K[K-factor] ; set K-factor