"""100 spheres falling under gravity

check basic equations working. body force and integrator
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

# pysph base and carray imports
from pysph.base.utils import get_particle_array_dem
from pysph.base.kernels import cubicspline

from pysph.solver.solver import solver
from pysph.sph.integrator import epecintegrator
from pysph.sph.integrator_step import demstep

from pysph.sph.equation import group
from pysph.sph.rigid_body import (
    bodyforce
     )
from pysph.solver.application import application


def add_properties(pa, *props):
    for prop in props:
        pa.add_property(name=prop)


class fluidstructureinteration(application):
    def initialize(self):
        self.dx = 1

    def create_particles(self):
        x = np.linspace(0, 1, 10)
        y = np.linspace(0, 1, 10)
        x, y = np.meshgrid(x, y)
        x, y = x.ravel(), y.ravel()
        m_inverse = np.ones_like(x) * 1
        m = np.ones_like(x) * 1
        sand = get_particle_array_dem(x=x, y=y, m=m, m_inverse=m_inverse,
                                      name="sand")

        return [sand]

    def create_solver(self):
        kernel = cubicspline(dim=2)

        integrator = epecintegrator(sand=demstep())
        # integrator = epecintegrator(fluid=wcsphstep())

        dt = 1e-3
        print("dt: %s" % dt)
        tf = 2
        solver = solver(kernel=kernel, dim=2, integrator=integrator, dt=dt,
                        tf=tf, adaptive_timestep=false)

        return solver

    def create_equations(self):
        equations = [
            group(equations=[
                bodyforce(dest='sand', sources=None, gy=-9.81),
            ]),
        ]
        return equations

    # def pre_step(self, solver):
    #     solver.dump_output()


if __name__ == '__main__':
    app = fluidstructureinteration()
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

#  localwords:  summationdensityshepardfilter
