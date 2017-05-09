from pysph.base.utils import get_particle_array
from pysph.sph.equation import Equation
import numpy as np


def get_particle_array_dem(constants=None, **props):
    """Return a particle array for the DEM formulation.

    These are the additional properties to get_particle_array

        ['x0 ', 'y0', 'z0', 'u0', 'v0', 'w0', 'fx', 'fy', 'fz', 'R']

    Parameters
    ----------
    constants : dict
        Dictionary of constants

    Other Parameters
    ----------------
    props : dict
        Additional keywords passed are set as the property arrays.

    See Also
    --------
    get_particle_array

    """

    dem_props = [
        'x0', 'y0', 'z0', 'u0', 'v0', 'w0', 'wx ', 'wy', 'wz', 'wx0', 'wy0',
        'wz0', 'fx', 'fy', 'fz', 'R', 'm_inverse'
    ]

    pa = get_particle_array(constants=constants, additional_props=dem_props,
                            **props)

    # default property arrays to save out.
    # pa.set_output_arrays([
    #     'x', 'y', 'z', 'u', 'v', 'w', 'fx ', 'fy', 'fz', 'm', 'pid', 'gid',
    #     'tag', 'p'
    # ])
    pa.set_output_arrays([
        'x', 'y', 'z', 'u', 'v', 'w', 'wx', 'wy', 'wz', 'm', 'p', 'pid', 'tag',
        'gid', 'fx', 'fy', 'fz'
    ])
    # pa.set_output_arrays([
    #     'x', 'y', 'z', 'u', 'v', 'w', 'm', 'pid', 'gid',
    #     'tag', 'p'
    # ])

    return pa


class BodyForce(Equation):
    def __init__(self, dest, sources, gx=0.0, gy=0.0, gz=0.0):
        self.gx = gx
        self.gy = gy
        self.gz = gz
        super(BodyForce, self).__init__(dest, sources)

    def initialize(self, d_idx, d_m, d_fx, d_fy, d_fz):
        d_fx[d_idx] = d_m[d_idx] * self.gx
        d_fy[d_idx] = d_m[d_idx] * self.gy
        d_fz[d_idx] = d_m[d_idx] * self.gz


class LinearSpringForceParticleParticle(Equation):
    """Documentation for LinearSpringForce

    """

    def __init__(self, dest, sources, k=1e4, ln_e=1.0, m_eff=0.5):
        super(LinearSpringForceParticleParticle, self).__init__(dest, sources)
        self.k = k

        ln_e2 = ln_e * ln_e
        _tmp = np.sqrt(ln_e2 + np.pi * np.pi)
        self.eta = 2 * np.sqrt(m_eff * self.k) * ln_e / _tmp

    def loop(self, d_idx, d_m, s_idx, d_fx, d_fy, d_fz, VIJ, XIJ, RIJ, d_R,
             s_R):
        overlap = 0

        if RIJ > 0:
            overlap = d_R[d_idx] + s_R[s_idx] - RIJ

        if overlap > 0:
            # basic variables: normal vector
            _rij = 1. / RIJ
            nx = -XIJ[0] * _rij
            ny = -XIJ[1] * _rij
            nz = -XIJ[2] * _rij

            # conservative force
            # d_fx[d_idx] += 1e2 * overlap * nx
            # d_fy[d_idx] += 1e2 * overlap * ny
            # d_fz[d_idx] += 1e2 * overlap * nz

            # damping force
            # v_n = VIJ[0] * nx + VIJ[1] * ny + VIJ[2] * nz
            # d_fx[d_idx] += -2 * v_n * nx
            # d_fy[d_idx] += -2 * v_n * ny
            # d_fz[d_idx] += -2 * v_n * nz

            # conservative damping force at one in single equation
            v_n = (VIJ[0] * nx + VIJ[1] * ny + VIJ[2] * nz)
            k_multi_overlap = self.k * overlap
            eta_multi_normal_velocity = self.eta * v_n
            d_fx[d_idx] += (-k_multi_overlap - eta_multi_normal_velocity) * nx
            d_fy[d_idx] += (-k_multi_overlap - eta_multi_normal_velocity) * ny
            d_fz[d_idx] += (-k_multi_overlap - eta_multi_normal_velocity) * nz


class MakeForcesZero(Equation):
    """Documentation for LinearSpringForce

    """

    def loop(self, d_idx, d_fx, d_fy, d_fz):
        d_fx[d_idx] = 0
        d_fy[d_idx] = 0
        d_fz[d_idx] = 0
