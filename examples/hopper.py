"""100 spheres falling inside hopper

Check the complete molecular dynamics code
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
from pysph.sph.molecular_dynamics import (LinearSpringForceParticleParticle,
                                          MakeForcesZero, BodyForce)
from pysph.solver.application import Application


def add_properties(pa, *props):
    for prop in props:
        pa.add_property(name=prop)


def create_hopper(r):
    d = 2 * r
    x_start = 0.2
    x_final = x_start
    y_start = 0
    y_final = y_start
    theta = 60 * np.pi / 180.
    x = []
    y = []
    while x_final < 2:
        x.append(x_final)
        y.append(y_final)
        x_final = x_final + d * np.cos(theta)
        y_final = y_final + d * np.sin(theta)
    x_l = np.asarray(x)
    y_l = np.asarray(y)

    x_r = -np.asarray(x)
    y_r = np.asarray(y)

    x, y = np.concatenate([x_l, x_r]), np.concatenate([y_l, y_r])
    return x, y


class FluidStructureInteration(Application):
    def initialize(self):
        self.dx = 1

    def create_particles(self):
        x = np.linspace(-0.5, 0.5, 10)
        y = np.linspace(0.77, 1.77, 10)
        r = (x[1] - x[0]) / 2.
        x, y = np.meshgrid(x, y)
        x, y = x.ravel(), y.ravel()
        R = np.ones_like(x) * r
        _m = np.pi * 2*r * 2*r
        m = np.ones_like(x) * _m
        m_inverse = np.ones_like(x) * 1. / _m
        _I = 2. / 5. * _m * r**2
        I_inverse = np.ones_like(x) * 1. / _I
        h = np.ones_like(x) * r
        sand = get_particle_array_dem(x=x, y=y, m=m, m_inverse=m_inverse, R=R,
                                      h=h, I_inverse=I_inverse, name="sand")

        x, y = create_hopper(r)
        m = np.ones_like(x) * _m
        m_inverse = np.ones_like(x) * 1. / _m
        R = np.ones_like(x) * r
        h = np.ones_like(x) * r
        _I = 2. / 5. * _m * r**2
        I_inverse = np.ones_like(x) * 1. / _I
        wall = get_particle_array_dem(x=x, y=y, m=m, m_inverse=m_inverse, R=R,
                                      h=h, I_inverse=I_inverse, name="wall")
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
                BodyForce(dest='sand', sources=None, gy=-9.81),
                LinearSpringForceParticleParticle(
                    dest='sand', sources=['sand', 'wall'], k=1e4,
                    ln_e=abs(np.log(0.8)), m_eff=0.5, mu=.5),
                # MakeForcesZero(dest='sand', sources=None)
            ]),
        ]
        return equations

    # def pre_step(self, solver):
    #     solver.dump_output()


if __name__ == '__main__':
    app = FluidStructureInteration()
    app.run()
    # x, y = create_hopper(0.1)
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
