import kiui
from kiui.mesh import Mesh

mesh = Mesh.load("todoial_propeller.obj", resize=False, clean=True)

mesh.clean_mesh(mesh.v, mesh.f)