
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
        skirt_dielectric,
        skirt_electrode,
        electrode_left,
        electrode_right,
        electrode_left_ironing,
        electrode_right_ironing,
        contact_left,
        contact_right,
        insulation_left,
        insulation_left,
        dielectric,
        dielectric_ironing
    }
    printer_settings = {
        cooling,
        prime_macro,
        tools
    }
    """
    def __init__(self, g_code_generators, regions, printer_settings, layer_cam=False):
        
        self.gen = g_code_generators    # g-code generation objects for each material
        self.regions = regions          # dict of regions
        self.printer_settings = printer_settings

        # printer settings
        # self.tools = printer_settings['tools'] # dict: {'passive': 'T0', 'electrode': 'T1', 'dielectric': 'T3', etc.}
        # self.cooling = printer_settings['cooling']
        # self.prime_macro = printer_settings['prime_macro']

        # self.layer_cam = layer_cam
        
    # def print_perimeter(self, region_key, z):
    #     # checking material to define generator
    #     if self.regions[region_key]['material'] == 'electrode':
    #         generator = self.gen_electrode
    #     elif self.regions[region_key]['material'] == 'dielectric':
    #         generator = self.gen_dielectric
    #     g_code = generator.print_rectangular_perimeter(
    #         rectangle = self.regions[region_key]['coordinates'],
    #         z = z,
    #         start = self.regions[region_key]['start_pos'],
    #         speed_factor = self.regions[region_key]['speed_factor'],
    #         extrude_factor = self.regions[region_key]['extrude_factor'],
    #         comment = self.regions[region_key]['heading']
    #     )
    #     return g_code

    # def print_surface(self, region_key, z):
    #     # checking material to define generator
    #     if self.regions[region_key]['material'] == 'electrode':
    #         generator = self.gen_electrode
    #     elif self.regions[region_key]['material'] == 'dielectric':
    #         generator = self.gen_dielectric
    #     # generating g_code
    #     g_code = generator.print_surface(
    #         surface = self.regions[region_key]['coordinates'],
    #         z = z,
    #         infill_angle = self.regions[region_key]['infill_angle'],
    #         start = self.regions[region_key]['start_pos'],
    #         perimeter = self.regions[region_key]['perimeter'],
    #         overlap_factor = self.regions[region_key]['overlap_factor'],
    #         speed_factor = self.regions[region_key]['speed_factor'],
    #         extrude_factor = self.regions[region_key]['extrude_factor'],
    #         comment = self.regions[region_key]['heading']
    #     )
    #     return g_code

    # def iron_surface(self, region_key, z):
    #     # checking material to define generator
    #     if self.regions[region_key]['material'] == 'electrode':
    #         generator = self.gen_electrode_iron
    #     elif self.regions[region_key]['material'] == 'dielectric':
    #         generator = self.gen_dielectric_iron
        
    #     # # defining ironing angle (perpendicular to infill angle)
    #     # infill_angle = self.regions[region_key]['infill_angle']
    #     # if infill_angle == 0:
    #     #     ironing_angle = 90
    #     # else:
    #     #     ironing_angle = 0
        
    #     # generating g_code
    #     g_code = generator.print_surface(
    #         surface = self.regions[region_key]['coordinates_ironing'],
    #         z = z,
    #         infill_angle = self.regions[region_key]['infill_angle_ironing'],
    #         start = self.regions[region_key]['start_pos_ironing'],
    #         perimeter = False,
    #         overlap_factor = 1,
    #         speed_factor = 1,
    #         extrude_factor = 1,
    #         comment = self.regions[region_key]['heading'] + ' - ironing'
    #     )
    #     return g_code


    def print_electrode_left_layer(self, z):
        """Generates g-code for actuator layer with electrode 1.
        Sequence:
            - skirt_electrode
            - electrode_left
            - contact_right
            - ironing
            - tool_change
            - skirt_dielectric
            - insulation_right
        Args:
            z (float): layer height
        Returns:
            g_code_dict: g_code for each tool
        
        TODO: alternating start positions
        
        """

        g_code = ''

        # print skirt_electrode
        g_code += self.gen.electrode.print_region(self.regions['skirt_electrode'], z_height=z)

        # print electrode_left
        g_code += self.gen.electrode.print_region(self.regions['electrode_left'], z_height=z)

        # print contact_left
        g_code += self.gen.electrode.print_region(self.regions['contact_right'], z_height=z)

        # brush before ironing
        g_code += 'M98 P\"brush_conductive.g\"\n'

        # iron the left electrode
        g_code += self.gen.electrode_ironing.print_region(self.regions['electrode_left_ironing'], z_height=z)

        # brush after ironing
        g_code += 'M98 P\"brush_conductive.g\"\n'

        # tool_change: electrode to dielectric
        g_code += tool_change(
            current_material='electrode',
            next_material='dielectric',
            printer_settings=self.printer_settings)

        # print skirt_dielectric
        g_code += self.gen.dielectric.print_region(self.regions['skirt_dielectric'], z_height=z)

        # print insulation_right
        g_code += self.gen.dielectric.print_region(self.regions['insulation_right'], z_height=z)

        return g_code


    def print_electrode_right_layer(self, z):
        """Generates g-code for actuator layer with electrode 1.
        Sequence:
            - skirt_electrode
            - electrode_right
            - contact_left
            - ironing
            - tool_change
            - skirt_dielectric
            - insulation_left
        Args:
            z (float): layer height
        Returns:
            g_code_dict: g_code for each tool
        
        TODO: alternating start positions
        
        """

        g_code = ''

        # print skirt_electrode
        g_code += self.gen.electrode.print_region(self.regions['skirt_electrode'], z_height=z)

        # print electrode_right
        g_code += self.gen.electrode.print_region(self.regions['electrode_right'], z_height=z)

        # print contact_left
        g_code += self.gen.electrode.print_region(self.regions['contact_left'], z_height=z)

        # brush before ironing
        g_code += 'M98 P\"brush_conductive.g\"\n'

        # iron the left electrode
        g_code += self.gen.electrode_ironing.print_region(self.regions['electrode_right_ironing'], z_height=z)

        # brush after ironing
        g_code += 'M98 P\"brush_conductive.g\"\n'

        # tool_change: electrode to dielectric
        g_code += tool_change(
            current_material='electrode',
            next_material='dielectric',
            printer_settings=self.printer_settings)

        # print skirt_dielectric
        g_code += self.gen.dielectric.print_region(self.regions['skirt_dielectric'], z_height=z)

        # print insulation_right
        g_code += self.gen.dielectric.print_region(self.regions['insulation_left'], z_height=z)

        return g_code


    def print_dielectric_layer(self, z):
        """Generates g-code for actuator layer with electrode 1.
        Sequence:
            - skirt_dielectric
            - dielectric
            - ironing
            - tool_change
            - skirt_electrode
            - contact_left
            - contact_right
        Args:
            z (float): layer height
        Returns:
            g_code_dict: g_code for each tool
        
        TODO: alternating start positions
        
        """

        g_code = ''

        # print skirt_dielectric
        g_code += self.gen.dielectric.print_region(self.regions['skirt_dielectric'], z_height=z)

        # print dielectric
        g_code += self.gen.dielectric.print_region(self.regions['dielectric'], z_height=z)

        # brush before ironing
        g_code += 'M98 P\"brush.g\"\n'

        # iron the dielectric
        g_code += self.gen.dielectric_ironing.print_region(self.regions['dielectric_ironing'], z_height=z)

        # brush after ironing
        g_code += 'M98 P\"brush.g\"\n'

        # tool_change: dielectric to electrode      
        g_code += tool_change(
            current_material='dielectric',
            next_material='electrode',
            printer_settings=self.printer_settings)
        
        # print skirt_electrode
        g_code += self.gen.electrode.print_region(self.regions['skirt_electrode'], z_height=z)

        # print contact_left
        g_code += self.gen.electrode.print_region(self.regions['contact_left'], z_height=z)

        # print contact_right
        g_code += self.gen.electrode.print_region(self.regions['contact_right'], z_height=z)

        return g_code


    def generate_actuator_g_code(self, start_z, N):
        """
        Generates g-code for whole SDEA.

        Args:
            start_z (float): Height of base layer. 0 if printing directly onto print bed
            N (int): Number of active layers (active layer = electrode 1 + dielectric + electrode 2).
            end_with_electrode (bool, optional): Top layer is electrode. Defaults to False (top layer is dielectric).

        Returns:
            g_code_dict (dict)
            z_last (float): z height of the printed actuator
        """
        z = start_z # height ot the support surface (print bed or previusly printed surface)
        last_electrode = 'right'
        
        g_code_dict = {}
        
        # loop over active layers
        for i in range(N):

            
            # electrode - alternating left/right
            z += self.gen.electrode.layer_height # lift printing for electrode layer height
            z = round(z, 2)
            if last_electrode == 'left':
                g_code_electrode = self.print_electrode_right_layer(z)
                last_electrode = 'right'
            elif last_electrode == 'right':
                g_code_electrode = self.print_electrode_left_layer(z)
                last_electrode = 'left'
            g_code_dict[z] = g_code_electrode
            
            # dielectric
            z += self.gen.dielectric.layer_height # lift printing for dielectric layer height
            z = round(z, 2)
            g_code_dielectric = self.print_dielectric_layer(z)
            g_code_dict[z] = g_code_dielectric

        z_last = z
        
        return g_code_dict, z_last

