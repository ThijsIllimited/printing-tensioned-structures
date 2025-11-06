
G1 Z10 F720 ; Move print head up
G1 X0 Y200 F3600 ; park
G4 ; wait
M221 S100 ; reset flow
M900 K0 ; reset LA
M104 S0 ; turn off temperature
M140 S0 ; turn off heatbed
M107 ; turn off fan
M84 ; disable motors