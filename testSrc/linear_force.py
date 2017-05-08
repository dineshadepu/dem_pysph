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
from pysph.sph.molecular_dynamics import LinearSpringForceParticleParticle
from pysph.solver.application import Application


def add_properties(pa, *props):
    for prop in props:
        pa.add_property(name=prop)


class FluidStructureInteration(Application):
    def initialize(self):
        self.dx = 1

    def create_particles(self):
        x = np.linspace(0, 1, 10)
        y = np.linspace(2, 3, 10)
        r = (x[1] - x[0]) / 2.
        x, y = np.meshgrid(x, y)
        x = x.ravel()
        y = y.ravel()
        m_inverse = np.ones_like(x) * 1
        m = np.ones_like(x) * 1
        R = np.ones_like(x) * r
        h = np.ones_like(x) * 4 * r
        sand = get_particle_array_dem(x=x, y=y, m=m, m_inverse=m_inverse, R=R,
                                      h=h, name="sand")

        x = np.linspace(0, 1, 10)
        y = np.linspace(-1, 1, 1)
        x, y = np.meshgrid(x, y)
        x = x.ravel()
        y = y.ravel()
        m_inverse = np.ones_like(x) * 1
        m = np.ones_like(x) * 1
        R = np.ones_like(x) * r
        wall = get_particle_array_dem(x=x, y=y, m=m, m_inverse=m_inverse, R=R,
                                      name="wall")
        return [sand, wall]

    def create_solver(self):
        kernel = CubicSpline(dim=2)

        integrator = EPECIntegrator(sand=DEMStep())

        dt = 1e-3
        print("DT: %s" % dt)
        tf = 2
        solver = Solver(kernel=kernel, dim=2, integrator=integrator, dt=dt,
                        tf=tf, adaptive_timestep=False)

        return solver

    def create_equations(self):
        equations = [
            Group(equations=[
                BodyForce(dest='sand', sources=None, gy=-9.81),
                LinearSpringForceParticleParticle(dest='sand',
                                                  sources=['wall', 'sand'])
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
