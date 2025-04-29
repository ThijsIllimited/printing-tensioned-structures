# from compas.utilities import geometric_key
from compas.geometry import midpoint_point_point
# from compas.utilities import i_to_rgb

from compas.datastructures import Network

import json

from ast import literal_eval

import src.src.uw_spider


# class FD_Network(Network):

#     def __init__(self):
#         super(Network, self).__init__()
#         self.default_edge_attributes.update({'q':{},
#                                              'f':{},
#                                              'cables':{},
#                                              })


class FD_Network(object):

    def __init__(self):
        self.vertices = {}
        self.vertices_scaled = {}
        self.edges = {}
        self.q = {}
        self.f = {}
        self.gkey_vk = {}
        self.gkey_ek = {}
        self.cables = {}
        self.fixed = []
        self.loads = {}
        self.l0 = {}
        self.l1 = {}
        self.ls1 = {}
        self.eps = {}
        self.elem_radials_center = {}
        self.elem_radials_catch2 = {}
        self.elem_gap_in = {}
        self.elem_gap_out = {}
        self.elem_Uturn_rad = {}
        self.elem_Uturn_spir = {}
        self.elem_spir_cent = {}
        self.elem_spir_catch = {}
        self.Spi_UT_elem2 = {}
        self.elem_frame_ring = {}
        self.elem_frame_sec = {}
        self.gap_size = {}
        self.gkey_spiral = {}
        self.gkey_center_n = {}
        self.description = {}
        self.R  = {}
        self.th = {}
        self.edges_arc = {}
        self.Radials_keys = {},
        self.UTurn_keys = {},
        self.Spiral_center_keys = {},
        self.Spiral_catching_keys = {},
        self.Frame_ring_keys = {},
        self.Frame_sec_keys= {},
        self.xyz_vec= {},

    def add_cable(self, pts, q):
        eks = []
        for i in range(len(pts) - 1):
            u = pts[i]
            v = pts[i + 1]
            mpt = midpoint_point_point(u, v)
            # ek = self.gkey_ek[geometric_key(mpt)]
            self.q[ek] = q
            eks.append(ek)
        self.cables[self.num_cables()] = eks

    def update_cable_q(self, cable, q):
        for ek in self.cables[cable]:
            self.q[ek] = q

    def num_vertices(self):
        return len(self.vertices)

    def num_edges(self):
        return len(self.edges)

    def num_cables(self):
        return len(self.cables)

    def vertex_list(self):
        vks = sorted(self.vertices.keys(), key=int)
        return [self.vertex_xyz(vk) for vk in vks]
    
    def vertex_scaled_list(self):
        vks = sorted(self.vertices_scaled.keys(), key=int)
        return [self.vertex_scaled_xyz(vk) for vk in vks]

    def edges_list(self):
        eks = sorted(self.edges.keys(), key=int)
        return [self.edges[ek] for ek in eks]

    def q_list(self):
        eks = sorted(self.q.keys(), key=int)
        return [self.q[ek] for ek in eks]

    def load_list(self):
        nks = sorted(self.loads.keys(), key=int)
        return [self.loads[nk] for nk in nks]

    @classmethod
    def from_lines(cls, lines, fixed=[]):
        fdnet = cls()
        for l in lines:
            uv = []
            for i in range(2):
                # gk = geometric_key(l[i])
                if gk not in fdnet.gkey_vk:
                    vk = fdnet.num_vertices()
                    fdnet.gkey_vk[gk] = vk
                else:
                    vk = fdnet.gkey_vk[gk]
                fdnet.vertices[vk] = {'x':l[i][0], 'y':l[i][1], 'z':l[i][2]}
                uv.append(vk)
            ek = fdnet.num_edges()
            fdnet.edges[ek] = uv
            # gke = geometric_key(midpoint_point_point(l[0], l[1]))
            fdnet.gkey_ek[gke] = ek
            fdnet.q[ek] = 1

        # fdnet.fixed = [fdnet.gkey_vk[geometric_key(pt)] for pt in fixed]

        fdnet.loads = {vk: [0, 0, 0] for vk in fdnet.vertices.keys()}

        return fdnet

    def vertex_xyz(self, key):
        return self.vertices[key]['x'], self.vertices[key]['y'], self.vertices[key]['z']

    def vertex_scaled_xyz(self, key):
        return self.vertices_scaled[key]['x'], self.vertices_scaled[key]['y'], self.vertices_scaled[key]['z']
    
    def rhino_index_to_key(self, elem):    
        return int(self.num_edges() - elem - 1)

    @property
    def data(self):
        data = {'vertices'  : {},
                'vertices_scaled'  : {},
                'edges'     : {},
                'q'         : {},
                'f'         : {},
                'gkey_vk'   : {},
                'gkey_ek'   : {},
                'cables'    : {},
                'fixed'     : self.fixed,
                'loads'     : {},
                'l0'        : {},
                'l1'        : {},
                'ls1'        : {},
                'eps'        : {},
                'elem_radials_center'   : {},
                'elem_radials_catch2'   : {},
                'elem_gap_in'           : {},
                'elem_gap_out'          : {},
                'elem_Uturn_rad'        : {},
                'elem_Uturn_spir'       : {},
                'elem_spir_cent'        : {},
                'elem_spir_catch'       : {},
                'Spi_UT_elem2'          : {},
                'elem_frame_ring'       : {},
                'elem_frame_sec'        : {},
                'gap_size'              : {},
                'gkey_spiral'           : {},
                'gkey_center_n'         : {},
                'description'           : {},
                'R'                     : {},
                'th'                    : {},
                'edges_arc'             : {},
                'Radials_keys'          : {},
                'Spiral_center_keys'    : {},
                'Spiral_catching_keys'  : {},
                'UTurn_keys'            : {},
                'Frame_ring_keys'       : {},
                'Frame_sec_keys'        : {},
                'xyz_vec'               : {},
                }

        for key in self.vertices:
            data['vertices'][repr(key)] = self.vertices[key]
        
        for key in self.vertices_scaled:
            data['vertices_scaled'][repr(key)] = self.vertices_scaled[key]

        for key in self.edges:
            data['edges'][repr(key)] = self.edges[key]

        for key in self.q:
            data['q'][repr(key)] = self.q[key]

        for key in self.f:
            data['f'][repr(key)] = self.f[key]

        for key in self.gkey_vk:
            data['gkey_vk'][repr(key)] = self.gkey_vk[key]

        for key in self.gkey_ek:
            data['gkey_ek'][repr(key)] = self.gkey_ek[key]

        for key in self.cables:
            data['cables'][repr(key)] = self.cables[key]

        for key in self.loads:
            data['loads'][repr(key)] = self.loads[key]

        for key in self.l0:
            data['l0'][repr(key)] = self.l0[key]

        for key in self.l1:
            data['l1'][repr(key)] = self.l1[key]

        for key in self.ls1:
            data['ls1'][repr(key)] = self.ls1[key]

        for key in self.eps:
            data['eps'][repr(key)] = self.eps[key]

        for key in self.elem_radials_center:
            data['elem_radials_center'][repr(key)] = self.elem_radials_center[key]
        
        for key in self.elem_radials_catch2:
            data['elem_radials_catch2'][repr(key)] = self.elem_radials_catch2[key]

        for key in self.elem_gap_in:
            data['elem_gap_in'][repr(key)] = self.elem_gap_in[key]

        for key in self.elem_gap_out:
            data['elem_gap_out'][repr(key)] = self.elem_gap_out[key]

        for key in self.elem_Uturn_rad:
            data['elem_Uturn_rad'][repr(key)] = self.elem_Uturn_rad[key]
        
        for key in self.elem_Uturn_spir:
            data['elem_Uturn_spir'][repr(key)] = self.elem_Uturn_spir[key]

        for key in self.elem_spir_cent:
            data['elem_spir_cent'][repr(key)] = self.elem_spir_cent[key]

        for key in self.elem_spir_catch:
            data['elem_spir_catch'][repr(key)] = self.elem_spir_catch[key]

        for key in self.Spi_UT_elem2:
            data['Spi_UT_elem2'][repr(key)] = self.Spi_UT_elem2[key]

        for key in self.elem_frame_ring:
            data['elem_frame_ring'][repr(key)] = self.elem_frame_ring[key]

        for key in self.elem_frame_sec:
            data['elem_frame_sec'][repr(key)] = self.elem_frame_sec[key]

        for key in self.gap_size:
            data['gap_size'][repr(key)] = self.gap_size[key]

        for key in self.gkey_spiral:
            data['gkey_spiral'][repr(key)] = self.gkey_spiral[key]
        
        for key in self.gkey_center_n:
            data['gkey_center_n'][repr(key)] = self.gkey_center_n[key]

        for key in self.description:
            data['description'][repr(key)] = self.description[key]
        
        for key in self.R:
            data['R'][repr(key)] = self.R[key]
        
        for key in self.th:
            data['th'][repr(key)] = self.th[key]
        
        for key in self.edges_arc:
            data['edges_arc'][repr(key)] = self.edges_arc[key]
        
        for key in self.Radials_keys:
            data['Radials_keys'][repr(key)] = self.Radials_keys[key]
        
        for key in self.Spiral_center_keys:
            data['Spiral_center_keys'][repr(key)] = self.Spiral_center_keys[key]
        
        for key in self.Spiral_catching_keys:
            data['Spiral_catching_keys'][repr(key)] = self.Spiral_catching_keys[key]

        for key in self.UTurn_keys:
            data['UTurn_keys'][repr(key)] = self.UTurn_keys[key]
        
        for key in self.Frame_ring_keys:
            data['Frame_ring_keys'][repr(key)] = self.Frame_ring_keys[key]

        for key in self.Frame_sec_keys:
            data['Frame_sec_keys'][repr(key)] = self.Frame_sec_keys[key]

        for key in self.xyz_vec:
            data['xyz_vec'][repr(key)] = self.xyz_vec[key]

        return data

    @data.setter
    def data(self, data):
        vertices        = data.get('vertices') or {}
        vertices_scaled = data.get('vertices_scaled') or {}
        edges           = data.get('edges') or {}
        q               = data.get('q') or {}
        f               = data.get('f') or {}
        gkey_vk         = data.get('gkey_vk') or {}
        gkey_ek         = data.get('gkey_ek') or {}
        cables          = data.get('cables') or {}
        self.fixed      = data.get('fixed') or []
        loads           = data.get('loads') or []
        l0              = data.get('l0') or {}
        l1              = data.get('l1') or {}
        ls1              = data.get('ls1') or {}
        eps              = data.get('eps') or {}
        elem_radials_center = data.get('elem_radials_center') or {}
        elem_radials_catch2 = data.get('elem_radials_catch2') or {}
        elem_gap_in         = data.get('elem_gap_in') or {}
        elem_gap_out        = data.get('elem_gap_out') or {}
        elem_Uturn_rad      = data.get('elem_Uturn_rad') or {}
        elem_Uturn_spir      = data.get('elem_Uturn_spir') or {}
        elem_spir_cent      = data.get('elem_spir_cent') or {}
        elem_spir_catch     = data.get('elem_spir_catch') or {}
        Spi_UT_elem2        = data.get('Spi_UT_elem2') or {}
        elem_frame_ring     = data.get('elem_frame_ring') or {}
        elem_frame_sec      = data.get('elem_frame_sec') or {}
        gap_size            = data.get('gap_size') or {}
        gkey_spiral         = data.get('gkey_spiral') or {}
        gkey_center_n       = data.get('gkey_center_n') or {}
        description         = data.get('description') or {}
        R                   = data.get('R') or {}
        th                  = data.get('th') or {}
        edges_arc           = data.get('edges_arc') or {}
        Radials_keys        = data.get('Radials_keys') or {}
        Spiral_center_keys  = data.get('Spiral_center_keys') or {}
        Spiral_catching_keys= data.get('Spiral_catching_keys') or {}
        UTurn_keys          = data.get('UTurn_keys') or {}
        Frame_ring_keys     = data.get('Frame_ring_keys') or {}
        Frame_sec_keys      = data.get('Frame_sec_keys') or {}
        xyz_vec             = data.get('xyz_vec') or {}



        self.vertices = {}
        for key in vertices:
            self.vertices[literal_eval(key)] = vertices[key]

        self.vertices_scaled = {}
        for key in vertices_scaled:
            self.vertices_scaled[literal_eval(key)] = vertices_scaled[key]

        self.edges = {}
        for key in edges:
            self.edges[literal_eval(key)] = edges[key]

        self.q = {}
        for key in q:
            self.q[literal_eval(key)] = q[key]

        self.f = {}
        for key in f:
            self.f[literal_eval(key)] = f[key]

        self.gkey_vk = {}
        for key in gkey_vk:
            self.gkey_vk[literal_eval(key)] = gkey_vk[key]

        self.gkey_ek = {}
        for key in gkey_ek:
            self.gkey_ek[literal_eval(key)] = gkey_ek[key]

        self.cables = {}
        for key in cables:
            self.cables[literal_eval(key)] = cables[key]

        self.loads = {}
        for key in loads:
            self.loads[literal_eval(key)] = loads[key]
        
        self.l0 = {}
        for key in l0:
            self.l0[literal_eval(key)] = l0[key]

        self.l1 = {}
        for key in l1:
            self.l1[literal_eval(key)] = l1[key]
        
        self.ls1 = {}
        for key in ls1:
            self.ls1[literal_eval(key)] = ls1[key]

        self.eps = {}
        for key in eps:
            self.eps[literal_eval(key)] = eps[key]

        self.elem_radials_center = {}
        for key in elem_radials_center:
            self.elem_radials_center[literal_eval(key)] = elem_radials_center[key]
        
        self.elem_radials_catch2 = {}
        for key in elem_radials_catch2:
            self.elem_radials_catch2[literal_eval(key)] = elem_radials_catch2[key]
        
        self.elem_gap_in = {}
        for key in elem_gap_in:
            self.elem_gap_in[literal_eval(key)] = elem_gap_in[key]
        
        self.elem_gap_out = {}
        for key in elem_gap_out:
            self.elem_gap_out[literal_eval(key)] = elem_gap_out[key]
        
        self.elem_Uturn_rad = {}
        for key in elem_Uturn_rad:
            self.elem_Uturn_rad[literal_eval(key)] = elem_Uturn_rad[key]
        
        self.elem_Uturn_spir = {}
        for key in elem_Uturn_spir:
            self.elem_Uturn_spir[literal_eval(key)] = elem_Uturn_spir[key]
        
        self.elem_spir_cent = {}
        for key in elem_spir_cent:
            self.elem_spir_cent[literal_eval(key)] = elem_spir_cent[key]
        
        self.elem_spir_catch = {}
        for key in elem_spir_catch:
            self.elem_spir_catch[literal_eval(key)] = elem_spir_catch[key]
        
        self.Spi_UT_elem2 = {}
        for key in Spi_UT_elem2:
            self.Spi_UT_elem2[literal_eval(key)] = Spi_UT_elem2[key]
        
        self.elem_frame_ring = {}
        for key in elem_frame_ring:
            self.elem_frame_ring[literal_eval(key)] = elem_frame_ring[key]
        
        self.elem_frame_sec = {}
        for key in elem_frame_sec:
            self.elem_frame_sec[literal_eval(key)] = elem_frame_sec[key]
        
        self.gap_size = {}
        for key in gap_size:
            self.gap_size[literal_eval(key)] = gap_size[key]
        
        self.gkey_spiral = {}
        for key in gkey_spiral:
            self.gkey_spiral[literal_eval(key)] = gkey_spiral[key]

        self.gkey_center_n = {}
        for key in gkey_center_n:
            self.gkey_center_n[literal_eval(key)] = gkey_center_n[key]

        self.description = {}
        for key in description:
            self.description[literal_eval(key)] = description[key]
        
        self.R = {}
        for key in R:
            self.R[literal_eval(key)] = R[key]
        
        self.th = {}
        for key in th:
            self.th[literal_eval(key)] = th[key]
        
        self.edges_arc = {}
        for key in edges_arc:
            self.edges_arc[literal_eval(key)] = edges_arc[key]

        self.Radials_keys = {}
        for key in Radials_keys:
            self.Radials_keys[literal_eval(key)] = Radials_keys[key]
        
        self.Spiral_center_keys = {}
        for key in Spiral_center_keys:
            self.Spiral_center_keys[literal_eval(key)] = Spiral_center_keys[key]
        
        self.Spiral_catching_keys = {}
        for key in Spiral_catching_keys:
            self.Spiral_catching_keys[literal_eval(key)] = Spiral_catching_keys[key]
        
        self.UTurn_keys = {}
        for key in UTurn_keys:
            self.UTurn_keys[literal_eval(key)] = UTurn_keys[key]

        self.Frame_ring_keys = {}
        for key in Frame_ring_keys:
            self.Frame_ring_keys[literal_eval(key)] = Frame_ring_keys[key]
        
        self.Frame_sec_keys = {}
        for key in Frame_sec_keys:
            self.Frame_sec_keys[literal_eval(key)] = Frame_sec_keys[key]
        
        self.xyz_vec = {}
        for key in xyz_vec:
            self.xyz_vec[literal_eval(key)] = xyz_vec[key]

    def to_json(self, filepath):
        with open(filepath, 'w+') as fp:
            json.dump(self.data, fp)

    def to_lines(self):
        lines = []
        for u,v in self.edges_list():
            line = self.vertex_xyz(u), self.vertex_xyz(v)
            lines.append(line)
        return lines

    @classmethod
    def from_json(cls, filepath):
        with open(filepath, 'r') as fp:
            data = json.load(fp)
        room = cls()
        room.data = data
        return room

    # PLOTTING      

    def rhino_plot(self, color=False, elables=False, scaled = False, plot_stress = False):
        import rhinoscriptsyntax as rs

        if color:
            f = [self.f[i] for i in self.f]
            if plot_stress:
                f = [self.eps[i] for i in self.eps]
            minf = min(f)
            maxf = max(f)

        for ek in self.edges:
            if scaled:
                u = self.vertex_scaled_xyz(self.edges[ek][0])
                v = self.vertex_scaled_xyz(self.edges[ek][1])
            else:    
                u = self.vertex_xyz(self.edges[ek][0])
                v = self.vertex_xyz(self.edges[ek][1])
            if u != v:
                l = rs.AddLine(u, v)
            if elables:
                rs.AddTextDot(str(ek), rs.CurveMidPoint(l))
            if color:
                c = (self.f[ek] - minf) / (maxf - minf)
                if plot_stress:
                    c = (self.eps[ek] - minf) / (maxf - minf)
                c = i_to_rgb(c, False)
                rs.ObjectColor(l, c)
            rs.ObjectName(l, str(self.q[ek]))

        for vk in self.fixed:
            rs.AddPoint(self.vertex_xyz(vk))

    def plotly_net_plot(self, color=False, elables=False, scaled=False, plot_stress=False):
        import plotly.graph_objects as go
        line_marker = dict(width=1.5)
        x, y, z, c = [], [], [], []
        text_positions = []
        edge_labels = []

        # Determine color scale if needed
        if color:
            f = [self.f[i] for i in self.f] if not plot_stress else [self.eps[i] for i in self.eps]
            minf, maxf = min(f), max(f)

        for ek in self.edges:
            u = self.vertex_scaled_xyz(self.edges[ek][0]) if scaled else self.vertex_xyz(self.edges[ek][0])
            v = self.vertex_scaled_xyz(self.edges[ek][1]) if scaled else self.vertex_xyz(self.edges[ek][1])
            if u != v:
                x.extend([u[0], v[0], None])
                y.extend([u[1], v[1], None])
                z.extend([u[2], v[2], None])
                if color:
                    c_val = (self.f[ek] - minf) / (maxf - minf) if not plot_stress else (self.eps[ek] - minf) / (maxf - minf)
                    c.append(c_val)
                if elables:
                    mid_point = [(u[i] + v[i]) / 2 for i in range(3)]
                    text_positions.append(mid_point)
                    edge_labels.append(str(ek))

        # Create lines
        lines = [go.Scatter3d(x=x, y=y, z=z, mode='lines', line=line_marker)]

        # Create text labels for edges if needed
        texts = [go.Scatter3d(x=[pos[0]], y=[pos[1]], z=[pos[2]], text=[label], mode='text') for pos, label in zip(text_positions, edge_labels)]

        # Setup layout
        layout = go.Layout(title='Rhino Plotly Plot', scene=dict(aspectmode='data'))

        # Combine data and layout into a figure
        fig = go.Figure(data=lines, layout=layout)

        # Display the figure
        fig.show()

    def plotly_plot(self, scale = False, arcs = False, check_order = False, junction_vertices = None):
        import plotly.graph_objects as go
        line_marker = dict(color='rgb(0,0,0)', width=1.5)
        x, y, z = [], [],  []
        if arcs:
            for key, points in sorted(self.xyz_vec.items()):
                for k in range(len(points)-1):
                    x.extend([points[k][0], points[k+1][0], [None]])
                    y.extend([points[k][1], points[k+1][1], [None]])
                    z.extend([points[k][2], points[k+1][2], [None]])
        else:
            if scale:
                vertices = self.vertex_scaled_list()
            else:
                vertices = self.vertex_list()
            edges = [[vertices[u], vertices[v]] for u, v in self.edges_list()]
            for u, v in edges:
                x.extend([u[0], v[0], [None]])
                y.extend([u[1], v[1], [None]])
                z.extend([u[2], v[2], [None]])
        if check_order:
            x = [item for item in x if item != [None]]
            y = [item for item in y if item != [None]]
            z = [item for item in z if item != [None]]
        lines = [go.Scatter3d(x=x, y=y, z=z, mode='lines', line=line_marker)]
        
        if junction_vertices:
            markers = [go.Scatter3d(x=[v[0] for v in junction_vertices], y=[v[1] for v in junction_vertices], z=[0.2*0.9]*len(junction_vertices), mode='markers', marker=dict(symbol = 'circle', color='red', size=3))]

        layout = go.Layout(title='fd plot',
                          scene=dict(aspectmode='data',
                                    xaxis=dict(
                                               gridcolor='rgb(255, 255, 255)',
                                               zerolinecolor='rgb(255, 255, 255)',
                                               showbackground=False,
                                               backgroundcolor='rgb(230, 230,230)',
                                               tickfont=dict(color='rgb(255,255,255)'),
                                               title=''),
                                    yaxis=dict(
                                               gridcolor='rgb(255, 255, 255)',
                                               zerolinecolor='rgb(255, 255, 255)',
                                               showbackground=False,
                                               backgroundcolor='rgb(230, 230,230)',
                                               tickfont=dict(color='rgb(255,255,255)'),
                                               title=''),
                                    zaxis=dict(
                                               gridcolor='rgb(255, 255, 255)',
                                               zerolinecolor='rgb(255, 255, 255)',
                                               showbackground=False,
                                               backgroundcolor='rgb(230, 230,230)',
                                               tickfont=dict(color='rgb(255,255,255)'),
                                               title='')
                                    ),
                          showlegend=False,
                            )
        if junction_vertices:
            fig = go.Figure(data=lines +  markers, layout=layout)
        else:
            fig = go.Figure(data=lines, layout=layout)
        name = 'eye = (x:0., y:0., z:10)'
        camera = dict(eye=dict(x=0., y=0., z=10), up=dict(x=0, y=1., z=0),)
        fig.update_layout(scene_camera=camera, title=name)
        fig.show()


if __name__ == '__main__':
    import os

    filepath = os.path.join(uw_spider.DATA, 'networks', 'Web_spirals_r56_s19_ecc_temp20.json')
    net = FD_Network.from_json(filepath)
    net.plotly_plot()
