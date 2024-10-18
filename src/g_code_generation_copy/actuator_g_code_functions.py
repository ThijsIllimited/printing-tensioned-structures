
# import sys
# sys.path.insert(0, r"C:\Users\tibor\Documents\GitHub\Python-modules")

import numpy as np
from tool_changer_functions import tool_change, take_photo

# TODO:
# uredi kako se uporabi speed in extrude factor v funkcijah

class Actuator_g_code_generator():
    """
    Class for actuator g_code generation.

    Params:
    g_code_generators = {
        gen_electrode,
        gen_dielectric,
        gen_electrode_iron,
        gen_dielectric_iron
    }
    regions = {
        skirt_1_rect,
        skirt_2_rect,
        electrode_1_reg,
        electrode_2_reg,
        electrode_1_per_rect,
        electrode_2_per_rect,
        contact_1_reg,
        contact_2_reg,
        insulation_1_reg,
        insulation_2_reg,
        dielectric_reg
    }
    printer_settings = {
        cooling,
        prime_macro,
        tools
    }
    """
    def __init__(self, g_code_generators, regions, printer_settings, layer_cam=False):
        
        self.g_code_generators = g_code_generators
        self.regions = regions
        self.printer_settings = printer_settings

        # g_code generator objects
        self.gen_electrode = g_code_generators['gen_electrode']
        self.gen_dielectric = g_code_generators['gen_dielectric']
        self.gen_electrode_iron = g_code_generators['gen_electrode_iron']
        self.gen_dielectric_iron = g_code_generators['gen_dielectric_iron']

        # printer settings
        self.tools = printer_settings['tools'] # dict: {'passive': 'T0', 'electrode': 'T1', 'dielectric': 'T3', etc.}
        self.cooling = printer_settings['cooling']
        self.prime_macro = printer_settings['prime_macro']

        self.layer_cam = layer_cam
        
    def print_perimeter(self, region_key, z):
        # checking material to define generator
        if self.regions[region_key]['material'] == 'electrode':
            generator = self.gen_electrode
        elif self.regions[region_key]['material'] == 'dielectric':
            generator = self.gen_dielectric
        g_code = generator.print_rectangular_perimeter(
            rectangle = self.regions[region_key]['coordinates'],
            z = z,
            start = self.regions[region_key]['start_pos'],
            speed_factor = self.regions[region_key]['speed_factor'],
            extrude_factor = self.regions[region_key]['extrude_factor'],
            comment = self.regions[region_key]['heading']
        )
        return g_code

    def print_surface(self, region_key, z):
        # checking material to define generator
        if self.regions[region_key]['material'] == 'electrode':
            generator = self.gen_electrode
        elif self.regions[region_key]['material'] == 'dielectric':
            generator = self.gen_dielectric
        # generating g_code
        g_code = generator.print_surface(
            surface = self.regions[region_key]['coordinates'],
            z = z,
            infill_angle = self.regions[region_key]['infill_angle'],
            start = self.regions[region_key]['start_pos'],
            perimeter = self.regions[region_key]['perimeter'],
            overlap_factor = self.regions[region_key]['overlap_factor'],
            speed_factor = self.regions[region_key]['speed_factor'],
            extrude_factor = self.regions[region_key]['extrude_factor'],
            comment = self.regions[region_key]['heading']
        )
        return g_code

    def iron_surface(self, region_key, z):
        # checking material to define generator
        if self.regions[region_key]['material'] == 'electrode':
            generator = self.gen_electrode_iron
        elif self.regions[region_key]['material'] == 'dielectric':
            generator = self.gen_dielectric_iron
        
        # # defining ironing angle (perpendicular to infill angle)
        # infill_angle = self.regions[region_key]['infill_angle']
        # if infill_angle == 0:
        #     ironing_angle = 90
        # else:
        #     ironing_angle = 0
        
        # generating g_code
        g_code = generator.print_surface(
            surface = self.regions[region_key]['coordinates_ironing'],
            z = z,
            infill_angle = self.regions[region_key]['infill_angle_ironing'],
            start = self.regions[region_key]['start_pos_ironing'],
            perimeter = False,
            overlap_factor = 1,
            speed_factor = 1,
            extrude_factor = 1,
            comment = self.regions[region_key]['heading'] + ' - ironing'
        )
        return g_code

    def actuator_electrode_1(self, z, speed_factor=1, extrude_factor=1): # speed and extrude factor are not used
        """Generates g-code for actuator layer with electrode 1.
        Sequence:
            - skirt_2
            - electrode_1
            - (ironing)
            - contact_2
            - tool_change
            - skirt_1
            - electrode_1_perimeter
            - insulation_2
        Args:
            z (float): layer height
            speed_factor (int, optional): Speed factor for whole layer. Defaults to 1.
            extrude_factor (int, optional): Extrude factor for whole layer. Defaults to 1.
        Returns:
            g_code_dict: g_code for each tool
        """
        # tools definition:
        T_el = self.tools['electrode']
        T_de = self.tools['dielectric']
        
        g_code_dict = {
            T_el: '',
            T_de: ''
        }
        
        # tool for electrodes
        g_code = f'; electrode 1 - start: z = {z} mm\n'
        # skirt 2
        g_code += self.print_perimeter(region_key='skirt_2', z=z)
        # electrode 1
        g_code += self.print_surface(region_key='electrode_1', z=z)
        if self.regions['electrode_1']['ironing']:
            # electrode 1 - ironing
            g_code += 'M98 P"brush.g"\n'
            g_code += 'M98 P"brush.g"\n\n'
            g_code += self.iron_surface(region_key='electrode_1', z=z)
            g_code += 'M98 P"brush.g"\n'
            g_code += 'M98 P"brush.g"\n\n'
            
        # start position for electrode region - alternating y start position
        x_start_ele = self.regions['electrode_1']['start_pos'][0]
        y_start_ele = ['y0', 'y1']
        start_positions_ele = [[x_start_ele, y_start_ele[0]], [x_start_ele, y_start_ele[1]]]
        last_start_pos_ele = self.regions['electrode_1']['start_pos']
        if last_start_pos_ele == start_positions_ele[0]:
            self.regions['electrode_1']['start_pos'] = start_positions_ele[1]
        else:
            self.regions['electrode_1']['start_pos'] = start_positions_ele[0]

        # contact 2
        g_code += self.print_surface(region_key='contact_2', z=z)
        g_code_dict[T_el] = g_code
        
        # tool for dielectric
        g_code = f'; insulation 2 - start: z = {z} mm\n'
        # skirt 1
        g_code += self.print_perimeter(region_key='skirt_1', z=z)
        # dielectric electrode 1 perimeter
        g_code += self.print_perimeter(region_key='electrode_1_per', z=z)
        g_code += 'M98 P"brush.g"\n\n'
        # insulation 2
        g_code += self.print_surface(region_key='insulation_2', z=z)
        g_code_dict[T_de] = g_code
        
        # start position for insulation region - alternating y start position
        x_start_ins = self.regions['insulation_2']['start_pos'][0]
        y_start_ins = ['y1', 'y0']
        start_positions_ins = [[x_start_ins, y_start_ins[0]], [x_start_ins, y_start_ins[1]]]
        last_start_pos_ins = self.regions['insulation_2']['start_pos']
        if last_start_pos_ins == start_positions_ins[0]:
            self.regions['insulation_2']['start_pos'] = start_positions_ins[1]
        else:
            self.regions['insulation_2']['start_pos'] = start_positions_ins[0]  
        
        return g_code_dict


    def actuator_electrode_2(self, z, speed_factor=1, extrude_factor=1):
        """Generates g-code for actuator layer with electrode 2.
        Sequence:
            - skirt_2
            - electrode_2
            - (ironing)
            - contact_1
            - tool_change
            - skirt_1
            - electrode_2_perimeter
            - insulation_1
        Args:
            z (float): layer height
            speed_factor (int, optional): Speed factor for whole layer. Defaults to 1.
            extrude_factor (int, optional): Extrude factor for whole layer. Defaults to 1.
        Returns:
            g_code_dict: g_code for each tool
        """
        # tools definition:
        T_el = self.tools['electrode']
        T_de = self.tools['dielectric']

        g_code_dict = {
            T_el: '',
            T_de: ''
        }

        # tool for electrodes
        g_code = f'; electrode 2 - start: z = {z} mm\n'
        # skirt 2
        g_code += self.print_perimeter(region_key='skirt_2', z=z)
        # electrode 2
        g_code += self.print_surface(region_key='electrode_2', z=z)
        if self.regions['electrode_2']['ironing']:
            # electrode 2 - ironing
            g_code += 'M98 P"brush.g"\n'
            g_code += 'M98 P"brush.g"\n\n'
            g_code += self.iron_surface(region_key='electrode_2', z=z)
            g_code += 'M98 P"brush.g"\n'
            g_code += 'M98 P"brush.g"\n\n' 
        
        # start position for electrode region - alternating y start position
        x_start_ele = self.regions['electrode_2']['start_pos'][0]
        y_start_ele = ['y1', 'y0']
        start_positions_ele = [[x_start_ele, y_start_ele[0]], [x_start_ele, y_start_ele[1]]]
        last_start_pos_ele = self.regions['electrode_2']['start_pos']
        if last_start_pos_ele == start_positions_ele[0]:
            self.regions['electrode_2']['start_pos'] = start_positions_ele[1]
        else:
            self.regions['electrode_2']['start_pos'] = start_positions_ele[0]        
        
        # contact 1
        g_code += self.print_surface(region_key='contact_1', z=z)
        g_code_dict[T_el] = g_code

        # tool for dielectric
        g_code = f'; insulation 1 - start: z = {z} mm\n'
        # skirt 1
        g_code += self.print_perimeter(region_key='skirt_1', z=z)
        # dielectric electrode 1 perimeter
        g_code += self.print_perimeter(region_key='electrode_2_per', z=z)
        g_code += 'M98 P"brush.g"\n\n'
        # insulation 1
        g_code += self.print_surface(region_key='insulation_1', z=z)
        g_code_dict[T_de] = g_code
        
        # start position for insulation region - alternating y start position
        x_start_ins = self.regions['insulation_1']['start_pos'][0]
        y_start_ins = ['y1', 'y0']
        start_positions_ins = [[x_start_ins, y_start_ins[0]], [x_start_ins, y_start_ins[1]]]
        last_start_pos_ins = self.regions['insulation_1']['start_pos']
        if last_start_pos_ins == start_positions_ins[0]:
            self.regions['insulation_1']['start_pos'] = start_positions_ins[1]
        else:
            self.regions['insulation_1']['start_pos'] = start_positions_ins[0]         
        
        return g_code_dict


    def actuator_dielectric(self, z, speed_factor=1, extrude_factor=1):
        """Generates g-code for actuator layer with dielectric.
        Sequence:
            - skirt_1
            - dielectric
            - (ironing)
            - tool_change
            - skirt_2
            - contact_1
            - contact_2
        Args:
            z (float): layer height
            speed_factor (int, optional): Speed factor for whole layer. Defaults to 1.
            extrude_factor (int, optional): Extrude factor for whole layer. Defaults to 1.
        Returns:
            g_code_dict: g_code for each tool
        """
        # tools definition:
        T_el = self.tools['electrode']
        T_de = self.tools['dielectric']
        
        g_code_dict = {
            T_el: '',
            T_de: ''
        }
        
        # tool for dielectric
        g_code = f'; dielectric - start: z = {z} mm\n'
        # skirt 2
        g_code += self.print_perimeter(region_key='skirt_1', z=z)
        # dielectric
        g_code += self.print_surface(region_key='dielectric', z=z)
        if self.regions['dielectric']['ironing']:
            # electrode 2 - ironing
            g_code += 'M98 P"brush.g"\n'
            g_code += 'M98 P"brush.g"\n\n'
            g_code += self.iron_surface(region_key='dielectric', z=z)
            g_code += 'M98 P"brush.g"\n'
            g_code += 'M98 P"brush.g"\n\n'
        g_code_dict[T_de] = g_code
        
        # start position for dielectric region - alternating y start position
        start_positions = [['x0', 'y0'], ['x1', 'y1']]
        last_start_pos = self.regions['dielectric']['start_pos']
        if last_start_pos == start_positions[0]:
            self.regions['dielectric']['start_pos'] = start_positions[1]
            self.regions['dielectric']['start_pos_ironing'] = start_positions[1]
        else:
            self.regions['dielectric']['start_pos'] = start_positions[0]
            self.regions['dielectric']['start_pos_ironing'] = start_positions[0]

        # tool for electrodes
        g_code = f'; contacts - start: z = {z} mm\n'
        # skirt 1
        g_code += self.print_perimeter(region_key='skirt_2', z=z)
        # contact 1
        g_code += self.print_surface(region_key='contact_1', z=z)
        # contact 2
        g_code += self.print_surface(region_key='contact_2', z=z)
        g_code_dict[T_el] = g_code
        
        return g_code_dict


    def generate_actuator_g_code(self, start_z, N, end_with_electrode=False):
        """Generates g-code for whole SDEA.

        Args:
            start_z (float): Height of base layer. 0 if printing directly onto print bed
            N (int): Number of active layers (active layer = electrode 1 + dielectric + electrode 2).
            end_with_electrode (bool, optional): Top layer is electrode. Defaults to False (top layer is dielectric).

        Returns:
            g_code_dict (dict)
            z_last (float): z height of the printed actuator
        """
        z = start_z # current layer height
        last_electrode = 2
        
        # tools definition:
        T_el = self.tools['electrode']
        T_de = self.tools['dielectric']
        
        g_code_dict = {}
        
        for i in range(N):
            
            if i == 0:
                speed_factor = 0.7
                extrude_factor = 1.0
            else:
                speed_factor = 1.0
                extrude_factor = 1.0
            
            # electrode - alternating 1, 2
            z += self.gen_electrode.layer_height
            if last_electrode == 1:
                g_code_electrode = self.actuator_electrode_2(z, speed_factor=speed_factor, extrude_factor=extrude_factor)
                last_electrode = 2
            elif last_electrode == 2:
                g_code_electrode = self.actuator_electrode_1(z, speed_factor=speed_factor, extrude_factor=extrude_factor)
                last_electrode = 1
            g_code_dict[round(z, 2)] = g_code_electrode # g_code for each tool for this layer
            #print(f'Electrode {last_electrode} at {z:2.2f} mm.')

            
            # dielectric
            z += self.gen_dielectric.layer_height
            g_code_dielectric = self.actuator_dielectric(z, speed_factor=speed_factor, extrude_factor=extrude_factor)
            g_code_dict[round(z, 2)] = g_code_dielectric # g_code for each tool for this layer
            #print(f'Dielectric at {z:2.2f} mm.')
        
        if end_with_electrode:
            z += self.gen_electrode.layer_height
            if last_electrode == 1:
                g_code_electrode = self.actuator_electrode_2(z, speed_factor=0.7*speed_factor, extrude_factor=1.05*extrude_factor) # 5% increased extrusion for top layer
                last_electrode = 2
            elif last_electrode == 2:
                g_code_electrode = self.actuator_electrode_1(z, speed_factor=0.7*speed_factor, extrude_factor=1.05*extrude_factor) # 5% increased extrusion for top layer
                last_electrode = 1
            g_code_dict[round(z, 2)] = g_code_electrode # g_code for each tool for this layer

        z_last = z
        
        return g_code_dict, z_last

    
#############################################################################
# Funkcija za izračun lokacij posamenih regij/površin SDEA
# Primer definiranih dimenzij:
# dimensions = {
#     'reference_point': [200, 100],
#     'active_layer_dims': [15, 5],
#     'contact_width': 4,
#     'insulation_width': 1,
#     'support_offset': [10, 5],
#     'spacer_dims': [8, 3],
#     'skirt_1_distance': 1,
#     'skirt_2_distance': 1,
#     'trace_widths': {
#         'dielectric': 0.2,
#         'electrode': 0.5,
#         'passive': 0.5
#     }
# }

def calculate_region_locations(dimensions):
    """
    x, y -> actual locations on print bed
    Lx, Ly -> lengths in x and y axis
    Ox, Oy -> offsets in x and y axis
    W -> trace widths
    """
    
    # unpacking params:
    x_ref, y_ref = dimensions['reference_point']
    Lx_ac, Ly_ac = np.asarray(dimensions['active_layer_dims'])
    Lx_con = dimensions['contact_width']
    Lx_ins = dimensions['insulation_width']
    Ox_sup, Oy_sup = dimensions['support_offset']
    Lx_spac, Ly_spac = dimensions['spacer_dims']
    O_sk1 = dimensions['skirt_1_distance']
    O_sk2 = dimensions['skirt_2_distance']
    
    W_el = dimensions['trace_widths']['electrode']
    W_de = dimensions['trace_widths']['dielectric']
    
    # dimesions calculations
    # lengths in X
    Lx_de = Lx_ac + 2 * Lx_ins + 2 * W_de
    Lx_el = Lx_con + Lx_ins + Lx_ac
    Lx_dea = 2 * Lx_con + Lx_de
    Lx_sup = Lx_dea + 2 * Ox_sup
    Lx_sk1 = Lx_sup + 2 * O_sk1
    Lx_sk2 = Lx_sup + 2 * O_sk2
    # lengths in Y
    Ly_de = Ly_ac + 2 * W_de
    Ly_el = Ly_ac
    Ly_con = Ly_de
    Ly_ins = Ly_de
    Ly_dea = Ly_de
    Ly_sup = Ly_dea + 2 * Oy_sup
    Ly_sk1 = Ly_sup + 2 * O_sk1
    Ly_sk2 = Ly_sup + 2 * O_sk2

    region_dims = {
        'dielectric': [Lx_de, Ly_de],
        'electrode': [Lx_el, Ly_el],
        'dea': [Lx_dea, Ly_dea],
        'support': [Lx_sup, Ly_sup],
        'skirt_1': [Lx_sk1, Ly_sk1],
        'skirt_2': [Lx_sk2, Ly_sk2]
    }
    
    # region locations
    # left vertices in X
    x_sk2_l = x_ref
    x_sk1_l = x_sk2_l + (O_sk2 - O_sk1)
    x_sup_l = x_sk1_l + O_sk1
    x_dea_l = x_sup_l + Ox_sup
    
    x_el1_l = x_dea_l + W_de
    x_el1_per_l = x_el1_l - W_de / 2
    x_con1_l = x_dea_l
    x_ins1_l = x_dea_l + Lx_con
    x_de_l = x_dea_l + Lx_con
    
    x_el2_l = x_dea_l + Lx_con + Lx_ins + W_de
    x_el2_per_l = x_el2_l - W_de / 2
    x_con2_l = x_dea_l + Lx_con + Lx_de
    x_ins2_l = x_con2_l - Lx_ins

    # right vertices in X
    x_sk1_r = x_sk1_l + Lx_sk1
    x_sk2_r = x_sk2_l + Lx_sk2
    x_sup_r = x_sup_l + Lx_sup
    x_dea_r = x_dea_l + Lx_dea
    
    x_el1_r = x_el1_l + Lx_el
    x_el1_per_r = x_el1_r + W_de / 2
    x_con1_r = x_con1_l + Lx_con
    x_ins1_r = x_ins1_l + Lx_ins
    x_de_r = x_de_l + Lx_de
    
    x_el2_r = x_el2_l + Lx_el
    x_el2_per_r = x_el2_r + W_de / 2
    x_con2_r = x_con2_l + Lx_con
    x_ins2_r = x_ins2_l + Lx_ins
    
    # left vertices in Y
    y_sk2_l = y_ref
    y_sk1_l = y_sk2_l + (O_sk2 - O_sk1)
    y_sup_l = y_sk1_l + O_sk1
    y_dea_l = y_sup_l + Oy_sup
    
    y_el1_l = y_dea_l + W_de
    y_el1_per_l = y_el1_l - W_de / 2
    y_con1_l = y_dea_l
    y_ins1_l = y_dea_l
    y_de_l = y_dea_l
    
    y_el2_l = y_dea_l + W_de
    y_el2_per_l = y_el2_l - W_de / 2
    y_con2_l = y_dea_l
    y_ins2_l = y_dea_l
    
    # right vertices in Y
    y_sk1_r = y_sk1_l + Ly_sk1
    y_sk2_r = y_sk2_l + Ly_sk2
    y_sup_r = y_sup_l + Ly_sup
    y_dea_r = y_dea_l + Ly_dea
    
    y_el1_r = y_el1_l + Ly_el
    y_el1_per_r = y_el1_r + W_de / 2
    y_con1_r = y_con1_l + Ly_con
    y_ins1_r = y_ins1_l + Ly_ins
    y_de_r = y_de_l + Ly_de
    
    y_el2_r = y_el2_l + Ly_el
    y_el2_per_r = y_el2_r + W_de / 2
    y_con2_r = y_con2_l + Ly_con
    y_ins2_r = y_ins2_l + Ly_ins
    
    # spacer
    x_center = x_sk2_l + Lx_sk2 / 2
    y_center = y_sk2_l + Ly_sk2 / 2
    x_spac_l = x_center - Lx_spac / 2
    x_spac_r = x_center + Lx_spac / 2
    y_spac_l = y_center - Ly_spac / 2
    y_spac_r = y_center + Ly_spac / 2
    
    region_locations = {
        'skirt_1': [[x_sk1_l, y_sk1_l],[x_sk1_r, y_sk1_r]],
        'skirt_2': [[x_sk2_l, y_sk2_l],[x_sk2_r, y_sk2_r]],
        'support': [[x_sup_l, y_sup_l],[x_sup_r, y_sup_r]],
        'spacer': [[x_spac_l, y_spac_l],[x_spac_r, y_spac_r]],
        'electrode_1': [[x_el1_l, y_el1_l],[x_el1_r, y_el1_r]],
        'electrode_2': [[x_el2_l, y_el2_l],[x_el2_r, y_el2_r]],
        'electrode_1_ironing': [[x_el1_l - 0.3, y_el1_l - 0.3],[x_el1_r + 0.3, y_el1_r + 0.3]],
        'electrode_2_ironing': [[x_el2_l - 0.3, y_el2_l - 0.3],[x_el2_r + 0.3, y_el2_r + 0.3]],
        'electrode_1_perimeter': [[x_el1_per_l, y_el1_per_l],[x_el1_per_r, y_el1_per_r]],
        'electrode_2_perimeter': [[x_el2_per_l, y_el2_per_l],[x_el2_per_r, y_el2_per_r]],
        'contact_1': [[x_con1_l, y_con1_l],[x_con1_r, y_con1_r]],
        'contact_2': [[x_con2_l, y_con2_l],[x_con2_r, y_con2_r]],
        'insulation_1': [[x_ins1_l, y_ins1_l],[x_ins1_r, y_ins1_r]],
        'insulation_2': [[x_ins2_l, y_ins2_l],[x_ins2_r, y_ins2_r]],
        'dielectric': [[x_de_l, y_de_l],[x_de_r, y_de_r]]
    }

    return region_dims, region_locations