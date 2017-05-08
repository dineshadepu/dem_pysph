"""100 spheres falling under gravity

Check basic equations working. Body force and integrator
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

# PySPH base and carray imports
from pysph.base.utils import get_particle_array_dem
from pysph.base.kernels import CubicSpline

from pysph.solver.solver import Solver
from pysph.sph.integrator import EPECIntegrator
from pysph.sph.integrator_step import DEMStep

from pysph.sph.equation import Group
from pysph.sph.rigid_body import (BodyForce)
from pysph.sph.molecular_dynamics import (LinearSpringForceParticleParticle,
                                          MakeForcesZero)
from pysph.solver.application import Application


def add_properties(pa, *props):
    for prop in props:
        pa.add_property(name=prop)


class FluidStructureInteration(Application):
    def initialize(self):
        self.dx = 1

    def create_particles(self):
        x = np.asarray([-2, 2])
        y = np.asarray([3, 3])
        u = np.asarray([1, -1])
        m_inverse = np.ones_like(x) * 1
        m = np.ones_like(x) * 1
        R = np.ones_like(x) * 1
        h = np.ones_like(x) * 3
        sand = get_particle_array_dem(x=x, y=y, m=m, m_inverse=m_inverse, R=R,
                                      u=u, h=h, name="sand")

        x = np.asarray([0])
        y = np.asarray([0])
        m_inverse = np.ones_like(x) * 1
        m = np.ones_like(x) * 1
        R = np.ones_like(x) * 1
        h = np.ones_like(x) * 3
        wall = get_particle_array_dem(x=x, y=y, m=m, m_inverse=m_inverse, R=R,
                                      h=h, name="wall")
        return [sand, wall]

    def create_solver(self):
        kernel = CubicSpline(dim=2)

        integrator = EPECIntegrator(sand=DEMStep())

        dt = 1e-4
        print("DT: %s" % dt)
        tf = 2
        solver = Solver(kernel=kernel, dim=2, integrator=integrator, dt=dt,
                        tf=tf, adaptive_timestep=False)

        return solver

    def create_equations(self):
        equations = [
            Group(equations=[
                # BodyForce(dest='sand', sources=None, gy=-9.81),
                LinearSpringForceParticleParticle(dest='sand',
                                                  sources=['wall', 'sand']),
                MakeForcesZero(dest='sand', sources=None)
            ]),
        ]
        return equations

    # def pre_step(self, solver):
    #     solver.dump_output()


if __name__ == '__main__':
    app = FluidStructureInteration()
    app.run()
    # x, y = create_fluid()
    # xc, yc, indices = create_cube()
    # xt, yt = create_boundary(1 * 1e-3)
    # plt.scatter(x, y)
    # plt.scatter(xc, yc)
    # plt.scatter(xt, yt)
    # plt.axes().set_aspect('equal', 'datalim')
    # plt.show()
    # xt, yt = create_boundary(1 * 1e-3)
    # xc, yc, indices = create_cube()
    # xf, yf = create_fluid_with_solid_cube()
    # plt.scatter(xt, yt)
    # plt.scatter(xc, yc)
    # plt.scatter(xf, yf)
    # plt.axes().set_aspect('equal', 'datalim')
    # plt.show()

#  LocalWords:  SummationDensityShepardFilter
