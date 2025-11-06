
import numpy as np
import os
# from src.g_code_generation_copy.gcode_generator import G_code_generator

def make_purge_line(printing_params, gen, point0, point1, ef_extra = 1.0, comment=''):
    """"
    ""Make a purge line at the start of the print."
    """
    g_code = '\n'
    g_code += gen.move_to_point(point0[0:2], point0[2] + printing_params['nozzle_lift'], comment='Move to start point')
    g_code += gen.move_to_point(point0[0:2], point0[2], comment='Lower Nozzle')
    g_code += gen.unretract()
    g_code += gen._print_line(
            point0=point0,
            point1=point1,
            move_to_start=False, # move to start point without extruding
            extrude_factor=printing_params['extrude_factor']*ef_extra,
            comment=comment)
    g_code += gen.retract()
    g_code += gen.wipe(2 * np.pi)
    g_code += '\n'
    return g_code

def add_path_to_gcode(cor_list, gen, printing_params, comment = ' '):
    g_code = f'\n;{comment}\n'
    # cor_list   = np.array(cor_list)
    # Move the coordinates to the center of the bed and add the layer height
    # cor_list[:,0] += bed_width/2
    # cor_list[:,1] += bed_height/2 - 1
    cor_list[:,2] += printing_params['layer_height']
    g_code += '\n'
    # Move the coordinates to the start point
    g_code += gen.move_to_point(cor_list[0][0:2], cor_list[0][2] + printing_params['nozzle_lift'], comment='Move to start point')
    # Lower the nozzle
    g_code += gen.move_to_point(cor_list[0][0:2], cor_list[0][2], comment='Lower Nozzle')
    # Unretract the filament
    g_code += gen.unretract()
    # Print the path
    for point0, point1 in zip(cor_list[:-1], cor_list[1:]):
        g_code += gen._print_line(
            point0=point0,
            point1=point1,
            move_to_start=False, # move to start point without extruding
            extrude_factor = printing_params['extrude_factor'],
            comment=comment)
    # Retract the filament
    g_code += gen.retract()
    # Wipe the nozzle
    g_code += gen.wipe_from_last_points(g_code)    
    # Raise the nozzle
    g_code += gen.move_to_point(point1[0:2], point1[2] + printing_params['nozzle_lift'], speed_factor=0.5, comment='Raise Nozzle')
    return g_code

def replace_brackets(line, temperature_settings):
    """"Replace the brackets and the brackets contents with the corresponding variable from the dictionary"""
    for key, value in temperature_settings.items():
        key_b = '[' + key + ']'
        if key_b in line:
            line = line.replace(key_b, str(value))
    return line

def save_gcode(path, g_code):
    with open(path, "w") as g_code_file:
        g_code_file.write(g_code)