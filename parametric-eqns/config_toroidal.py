import numpy as np

# Toroidal propeller parameter configuration

# Domains
s_domain = [0, 1]
t_domain = [0, 2]

# CAD resolutions
s_resolution_cad = 25
t_resolution_cad = 60

# Global scaling / hub / blades
global_scale = 7.5
hub_radius   = 5
hub_length   = 18.0
num_blades   = 3

# Airfoil params
m         = 0.04
p         = 0.4
thickness = 0.75

# Centerline control points / blade vector
loc_ctrl_point2 = [-0.25, 14.0, 35.50000000000001]
loc_ctrl_point3 = [-5, 0.75, 30]
blade_vector    = [-12.500000000000028, 16.25]  # [circumferential_offset, z_offset]

# Angle of attack polynomial coefficients
a_AoA = -0.01592653589793258
b_AoA = -0.01592653589793258
c_AoA = -0.01592653589793258
d_AoA = 3.134073464102066
e_AoA = -0.01592653589793258

# Scaling polynomials (x / y)
a_scX = 1
b_scX = 0
c_scX = -2
d_scX = 6.100000000000001
e_scX = 5.710000000000001

a_scY = 0
b_scY = 0
c_scY = -1
d_scY = 0
e_scY = 2.01

# Flags
apply_thickness_normal = False