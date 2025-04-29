import os
# from compas.utilities import geometric_key, reverse_geometric_key
from compas.geometry import midpoint_point_point
# from scipy.sparse import diags
# from scipy.sparse.linalg import spsolve
# from numpy import ones

from compas.rpc import Proxy
# from compas.utilities import i_to_rgb

from src.src.uw_spider.fd.fd_network import FD_Network

def fd_from_rhino_layers(lines, fixed, qs, loads, color=False):
    import rhinoscriptsyntax as rs
    lines = rs.ObjectsByLayer(lines)
    lines = [[rs.CurveStartPoint(l), rs.CurveEndPoint(l)] for l in lines]

    vdict = {}
    for l in lines:
        vdict[geometric_key(l[0])] = l[0]
        vdict[geometric_key(l[1])] = l[1]

    vi_dict = {k:i for i, k in enumerate(vdict)}

    vertices = []
    for k in vdict:
        vertices.append([vdict[k][0], vdict[k][1],vdict[k][2]])

    edges = [[vi_dict[geometric_key(l[0])], vi_dict[geometric_key(l[1])]] for l in lines]

    fixed = [rs.PointCoordinates(p) for p in rs.ObjectsByLayer(fixed)]
    fixed = [vi_dict[geometric_key(f)] for f in fixed]

    edict = {}
    for i, e in enumerate(edges):
        mpt = midpoint_point_point(vertices[e[0]], vertices[e[1]])
        edict[geometric_key(mpt)] = i

    qlist = [1] * len(edges)

    pls = rs.ObjectsByLayer(qs)
    for pl in pls:
        pts = rs.PolylineVertices(pl)
        q = rs.ObjectName(pl)
        for i in range(len(pts) - 1):
            mpt = midpoint_point_point(pts[i], pts[i + 1])
            qlist[edict[geometric_key(mpt)]] = q

    loads_ = [[0,0,0]] * len(vertices)
    for pt in rs.ObjectsByLayer(loads):
        l = float(rs.ObjectName(pt))
        loads_[vi_dict[geometric_key(rs.PointCoordinates(pt))]] = [0, 0, l]

    proxy = Proxy()
    numerical = Proxy('compas.numerical')
    xyz, q, f, l, r = numerical.fd_numpy(vertices, edges, fixed, qlist, loads_)
    f = [force[0] for force in f]
    minf = min(f)
    maxf = max(f)

    for i, ek in enumerate(edges):
        l = rs.AddLine(xyz[ek[0]], xyz[ek[1]])
        rs.ObjectName(l, f[i])
        if color:
            c = (f[i] - minf) / (maxf - minf)
            c = i_to_rgb(c, False)
            rs.ObjectColor(l, c)

    return xyz, q, f, l, r

def fd_from_rhino_network(network, color=True):
    vertices = network.vertex_list()
    edges = network.edges_list()
    fixed = network.fixed
    qlist = network.q_list()
    loads_ = [network.loads] * len(vertices)

    numerical = Proxy('compas.numerical')
    xyz, q, f, _, _ = numerical.fd_numpy(vertices, edges, fixed, qlist, loads_)

    for i, p in enumerate(xyz):
        network.vertices[i] = {'x': p[0], 'y': p[1], 'z': p[2]}

    for i in range(len(q)):
        network.q[i] = q[i][0]
        network.f[i] = f[i][0]

def fd_from_rhino_layers_net(lines_layer, fixed_layer=None, qs_layer=None, loads_layer=None, color=False, elables=False):
    import rhinoscriptsyntax as rs
    lines = rs.ObjectsByLayer(lines_layer)
    lines = [[rs.CurveStartPoint(l), rs.CurveEndPoint(l)] for l in lines]

    if fixed_layer:
        fixed = [rs.PointCoordinates(p) for p in rs.ObjectsByLayer(fixed_layer)]

    network = FD_Network.from_lines(lines, fixed)

    if qs_layer:
        pls = rs.ObjectsByLayer(qs_layer)
        for pl in pls:
            pts = rs.PolylineVertices(pl)
            q = rs.ObjectName(pl)
            for i in range(len(pts) - 1):
                mpt = midpoint_point_point(pts[i], pts[i + 1])
                ek = network.gkey_ek[geometric_key(mpt)]
                network.q[ek] = q

    if loads_layer:
        for pt in rs.ObjectsByLayer(loads_layer):
            l = float(rs.ObjectName(pt))
            xyz = rs.PointCoordinates(pt)
            nk = network.gkey_vk[geometric_key(xyz)]
            network.loads[nk] = [0, 0, l]

    vertices = network.vertex_list()
    edges = network.edges_list()
    fixed = network.fixed
    qlist = network.q_list()
    loads_ = network.load_list()

    numerical = Proxy('compas.numerical')
    xyz, q, f, _, _ = numerical.fd_numpy(vertices, edges, fixed, qlist, loads_)

    for i, p in enumerate(xyz):
        network.vertices[i] = {'x': p[0], 'y': p[1], 'z': p[2]}

    for i in range(len(q)):
        network.q[i] = q[i][0]
        network.f[i] = f[i][0]

    network.rhino_plot(color=color, elables=elables)
    return network

def fd_from_rhino_layers_net2(lines_layer, fixed_layer=None, qs_layer=None, loads_layer=None, color=False, elables=False, elem_radials_center = None, elem_radials_catch2 = None, elem_gap_in = None, elem_gap_out = None, elem_Uturn_rad = None, elem_spir_cent = None, elem_spir_catch = None, Spi_UT_elem2 = None, elem_frame_ring = None, elem_frame_sec = None, description = None):
    import rhinoscriptsyntax as rs
    from compas.geometry import  vector_average, distance_point_point
    lines = rs.ObjectsByLayer(lines_layer)
    lines = [[rs.CurveStartPoint(l), rs.CurveEndPoint(l)] for l in lines]

    if fixed_layer:
        fixed = [rs.PointCoordinates(p) for p in rs.ObjectsByLayer(fixed_layer)]

    network = FD_Network.from_lines(lines, fixed)

    if qs_layer:
        pls = rs.ObjectsByLayer(qs_layer)
        for pl in pls:
            pts = rs.PolylineVertices(pl)
            q = rs.ObjectName(pl)
            for i in range(len(pts) - 1):
                mpt = midpoint_point_point(pts[i], pts[i + 1])
                ek = network.gkey_ek[geometric_key(mpt)]
                network.q[ek] = q

    if loads_layer:
        for pt in rs.ObjectsByLayer(loads_layer):
            l = reverse_geometric_key(rs.ObjectName(pt))
            # l = float(rs.ObjectName(pt))
            xyz = rs.PointCoordinates(pt)
            nk = network.gkey_vk[geometric_key(xyz)]
            network.loads[nk] = [l[0], l[1], l[2]]
    network.Radials_keys            = {'name' : 'Radials_keys'}
    network.Spiral_center_keys      = {'name' : 'Spiral_center'}
    network.Spiral_catching_keys    = {'name' : 'Spiral_catching'}
    network.Frame_ring_keys         = {'name' : 'Frame_ring'}
    network.Frame_sec_keys          = {'name' : 'Frame_sec'}
    network.UTurn_keys              = {'name' : 'UTurn_keys'}
    
    if elem_radials_center:
        for radial in range(elem_radials_center['n_radials']):
            network.elem_radials_center[radial]   = elem_radials_center[str(radial)+'_n']
            network.elem_gap_in[radial]           = elem_gap_in[str(radial)+'_n']
            network.elem_radials_catch2[radial]   = elem_radials_catch2[str(radial)+'_n']
            network.elem_gap_out[radial]          = elem_gap_out[str(radial)+'_n']
            # if elem_Uturn_rad:
            #     network.elem_Uturn_rad[radial] = elem_Uturn_rad[str(radial)+'_n']
            #     Full_rad = elem_radials_center[str(radial)+'_n'] + elem_gap_in[str(radial)+'_n'] + elem_radials_catch2[str(radial)+'_n'] + elem_Uturn_rad[str(radial)+'_n'] + elem_gap_out[str(radial)+'_n']
            # else:
            Full_rad = elem_radials_center[str(radial)+'_n'] + elem_gap_in[str(radial)+'_n'] + elem_gap_out[str(radial)+'_n'] + elem_radials_catch2[str(radial)+'_n']
            key_list = []
            for index in Full_rad:
                key_list.append(network.rhino_index_to_key(index))
            network.Radials_keys[radial] = key_list

        for edge_id, edge in enumerate(elem_spir_cent['line_n']):
            network.elem_spir_cent[edge_id] = elem_spir_cent['line_n'][edge_id]
            network.Spiral_center_keys[edge_id] = network.rhino_index_to_key(edge)

        for edge_id, edge in enumerate(elem_spir_catch['line_n']):
            network.elem_spir_catch[edge_id] = elem_spir_catch['line_n'][edge_id]
            network.Spiral_catching_keys[edge_id] = network.rhino_index_to_key(edge)
            # key = network.rhino_index_to_key(edge)
            # network.Spiral_catching_keys[key] = elem_spir_catch['line_n'][edge_id]

        for edge_id, edge in enumerate(Spi_UT_elem2['line_n']):
            network.Spi_UT_elem2[edge_id] = Spi_UT_elem2['line_n'][edge_id]
        for turn in range(Spi_UT_elem2['turns'] ):
            key_list = []
            for edge in Spi_UT_elem2[str(turn) + '_n']:
                key_list.append(network.rhino_index_to_key(edge))
            network.UTurn_keys[turn] = key_list

        for edge_id, edge in enumerate(elem_frame_ring['line_n']):
            network.elem_frame_ring[edge_id] = elem_frame_ring['line_n'][edge_id]
            network.Frame_ring_keys[edge_id] = network.rhino_index_to_key(edge)

        for edge_id, edge in enumerate(elem_frame_sec['line_n']):
            network.elem_frame_sec[edge_id] = elem_frame_sec['line_n'][edge_id]
            network.Frame_sec_keys[edge_id] = network.rhino_index_to_key(edge)
            # key = network.rhino_index_to_key(edge)
            # network.Frame_sec_keys[key] = elem_frame_sec['line_n'][edge_id]
    if description:
        network.description = description   
    vertices = network.vertex_list()
    edges = network.edges_list()
    fixed = network.fixed
    qlist = network.q_list()
    loads_ = network.load_list()

    numerical = Proxy('compas.numerical')

    xyz, q, f, _, _ = numerical.fd_numpy(vertices, edges, fixed, qlist, loads_)

    for i, p in enumerate(xyz):
        network.vertices[i] = {'x': p[0], 'y': p[1], 'z': p[2]}
        

    for i in range(len(q)):
        network.q[i]        = q[i][0]
        network.f[i]        = f[i][0]
    
    for key, [vert0, vert1] in sorted(network.edges.items()):
        network.l1[key]       = distance_point_point(network.vertex_xyz(vert0), network.vertex_xyz(vert1))

    if elem_gap_in:
        gap_size = []
        for edge in elem_gap_in['radial_n']:
            vertices = network.edges[network.num_edges() - elem_gap_in[str(edge)+'_n'][0]-1]
            gap_size.append(distance_point_point(network.vertex_xyz(vertices[0]), network.vertex_xyz(vertices[1])))
        network.gap_size['gap_size_vec'] = gap_size
        network.gap_size['gap_size_mean'] = vector_average(gap_size)
        
    network.rhino_plot(color=color, elables=elables)
    network.xyz_vec = {}
    return network

def fd_network(network, color=True):
    from compas.numerical import fd_numpy

    vertices = network.vertex_list()
    edges = network.edges_list()
    fixed = network.fixed
    qlist = network.q_list()
    loads_ = [network.loads] * len(vertices)

    xyz, q, f, _, _ = fd_numpy(vertices, edges, fixed, qlist, loads_)

    for i, p in enumerate(xyz):
        network.vertices[i] = {'x': p[0], 'y': p[1], 'z': p[2]}

    for i in range(len(q)):
        network.q[i] = q[i][0]
        network.f[i] = f[i][0]

def fd_materialization(q, l, ea, **kwargs):
    from scipy.sparse.linalg import spsolve
    from scipy.sparse import diags
    from numpy import ones
    """Parameters
    ---------
    q : list
        Force density of edges.
    f : array
        Forces in the edges.
    l : array
        Lengths of the edges
    ea : array
        Axial stiffness times the cross sectional area element

    Returns
    -------
    l0 : array
        Initial lengths of the edges
    """
    v = len(q)
    Q = diags([q.flatten()], [0])
    L = diags([l.flatten()], [0])
    EA = diags([ea.flatten()], [0])
    I = diags(ones(v))
    A = spsolve(EA, Q.dot(L)) + I
    l0 = spsolve(A,l)

    return l0

if __name__ == '__main__':
    import os
    import rhinoscriptsyntax as rs
    import uw_spider

    rs.DeleteObjects(rs.ObjectsByLayer('Default'))
    for i  in range(30): print('')

    fd_from_rhino_layers_net2('lines', 'fixed', 'qs', 'loads', color=True)

    # rs.DeleteObjects(rs.ObjectsByLayer('Default'))
    # rs.CurrentLayer('Default')
    #
    # network = fd_from_rhino_layers_net('lines', 'fixed')
    #
    # filepath = os.path.join(uw_spider.DATA, 'networks', 'test_02.json')
    # network.to_json(filepath)
