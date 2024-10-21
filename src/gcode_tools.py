import numpy as np
from src.GcodeGenerator.gcode_generator import G_code_generator

class GCode(G_code_generator):    
    def _print_line(self, point0, point1, move_to_start=True,
                    extrude_factor=1, speed_factor=1, 
                    comment=None):
        """
        Generates a g-code for a line without wipe/retraction/unretraction. 
        Enables printing in 3 axis (printing over air).

        Params:
        point0          ... list: [x, y, z] in mm
        point1          ... list: [x, y, z] in mm
        extrude_factor  ... float: extrusion multiplier
        speed_factor    ... float: feed_rate = speed_factor * self.move_feedrate
        comment         ... string: c_code comment at end of line
        """

        if comment == None:
            comment = 'single line'
        
        x0, y0, z0 = point0
        x1, y1, z1 = point1

        point0 = np.array(point0)
        point1 = np.array(point1)
        line_length = np.sqrt(np.sum((point1-point0)**2))
        extrude_length = self.calculate_extrusion_length(line_length)

        g_code = ''
        if move_to_start:
            # 1) move to print point
            g_code += self.move_to_printing_point([x0, y0, z0])
        # 3) print line
        g_code += f'G1 '
        g_code += f'X{x1:.3f} '
        g_code += f'Y{y1:.3f} '
        g_code += f'Z{z1:.3f} '
        g_code += f'E{extrude_length * extrude_factor:.5f} '
        g_code += f'F{self.print_feedrate * speed_factor:.0f} '
        g_code += f'; {comment}\n'
        self.nozzle_locations.append(point1) # adding point1 to object history

        return g_code

    def xyz_from_gcode_line(self, line):
        """
        Extracts the x, y, z coordinates from a gcode line.
        """
        x = float(line[line.find('X')+1:line.find(' ', line.find('X'))])
        y = float(line[line.find('Y')+1:line.find(' ', line.find('Y'))])
        z = float(line[line.find('Z')+1:line.find(' ', line.find('Z'))])
        return np.array([x, y, z])

    def wipe_from_last_points(self, gcode):
        """
        Wipes accross the previous points.
        first part of the code Keeps extracting the previous nozzle locations from the gcode until the total wiped distance exceeds the wipe distance set in the printing parameters.
        Second part wipes back and forth accross the extracted points.
        """
        split_lines = gcode.splitlines()
        last_line = split_lines[-2]
        current_xyz = self.xyz_from_gcode_line(last_line)
        xyz_list = [current_xyz]
        total_len = 0
        k=-3
        while total_len <= self.wipe_len:
            previous_line = split_lines[k]
            previous_xyz = self.xyz_from_gcode_line(previous_line)
            xyz_list.append(previous_xyz)
            total_len += np.linalg.norm(current_xyz - previous_xyz)
            current_xyz = previous_xyz
            k -= 1
        
        # Wipe back and forth
        xyz_list = np.array(xyz_list)
        g_code = ''
        for xyz in xyz_list[1:]:
            g_code += f'G1 '
            g_code += f'X{xyz[0]:.3f} '
            g_code += f'Y{xyz[1]:.3f} '
            g_code += f'F{self.wipe_feedrate:.0f} '
            g_code += '; wipe back\n'
        for xyz in np.flip(xyz_list[:-1], axis=0):
            g_code += f'G1 '
            g_code += f'X{xyz[0]:.3f} '
            g_code += f'Y{xyz[1]:.3f} '
            g_code += f'F{self.wipe_feedrate:.0f} '
            g_code += '; wipe forth\n'
        return g_code