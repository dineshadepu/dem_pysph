from pysph.base.utils import get_particle_array
from pysph.sph.equation import Equation
import numpy as np
from pysph.sph.integrator_step import IntegratorStep


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
        'x0', 'y0', 'z0', 'u0', 'v0', 'w0', 'tang_x', 'tang_y', 'tang_z',
        'tang_x0', 'tang_y0', 'tang_z0', 'vt_x', 'vt_y', 'vt_z', 'wx', 'wy',
        'wz', 'wx0', 'wy0', 'wz0', 'fx', 'fy', 'fz', 'torX', 'torY', 'torZ',
        'R', 'm_inverse', 'I_inverse'
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
        'gid', 'fx', 'fy', 'fz', 'torX', 'torY', 'torZ', 'I_inverse'
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

    def initialize(self, d_idx, d_m, d_fx, d_fy, d_fz, d_torX, d_torY, d_torZ):
        d_fx[d_idx] = d_m[d_idx] * self.gx
        d_fy[d_idx] = d_m[d_idx] * self.gy
        d_fz[d_idx] = d_m[d_idx] * self.gz

        d_torX[d_idx] = 0
        d_torY[d_idx] = 0
        d_torZ[d_idx] = 0


class LinearSpringForceParticleParticle(Equation):
    """Documentation for LinearSpringForce

    """

    def __init__(self, dest, sources, k=1e4, ln_e=1.0, m_eff=0.5, mu=0.5):
        super(LinearSpringForceParticleParticle, self).__init__(dest, sources)
        self.k = k

        ln_e2 = ln_e * ln_e
        _tmp = np.sqrt(ln_e2 + np.pi * np.pi)
        self.eta = 2 * np.sqrt(m_eff * self.k) * ln_e / _tmp

        self.kt = 2. / 7. * k
        self.etaT = 0.5 * self.eta

        self.mu = mu

    def loop(self, d_idx, d_m, d_wx, d_wy, d_wz, d_fx, d_fy, d_fz, d_torX,
             d_torY, d_torZ, d_tang_x, d_tang_y, d_tang_z, d_vt_x, d_vt_y,
             d_vt_z, VIJ, XIJ, RIJ, d_R, s_idx, s_R, s_wx, s_wy, s_wz):
        overlap = 0

        if RIJ > 0:
            overlap = d_R[d_idx] + s_R[s_idx] - RIJ

        if overlap > 0:
            # basic variables: normal vector
            _rij = 1. / RIJ
            nx = -XIJ[0] * _rij
            ny = -XIJ[1] * _rij
            nz = -XIJ[2] * _rij

            # radius multiplied by the angular velocity
            R_wijx = (d_R[d_idx] * d_wx[d_idx] + s_R[s_idx] * s_wx[s_idx])
            R_wijy = (d_R[d_idx] * d_wy[d_idx] + s_R[s_idx] * s_wy[s_idx])
            R_wijz = (d_R[d_idx] * d_wz[d_idx] + s_R[s_idx] * s_wz[s_idx])

            # angular relative velocity
            w_ij_x = R_wijy * nz - R_wijz * ny
            w_ij_y = -R_wijx * nz + R_wijz * nx
            w_ij_z = R_wijx * ny - R_wijy * nx

            # add angular relative velocity to linear relative velocity
            VIJ[0] += w_ij_x
            VIJ[1] += w_ij_y
            VIJ[2] += w_ij_z

            # normal force magnitude
            v_n = (VIJ[0] * nx + VIJ[1] * ny + VIJ[2] * nz)

            # tangential velocity
            vt_x = VIJ[0] - v_n * nx
            vt_y = VIJ[1] - v_n * ny
            vt_z = VIJ[2] - v_n * nz

            # tangential velocity magnitude
            vt_magn = pow(vt_x * vt_x + vt_y * vt_y + vt_z * vt_z, 1. / 2.)

            # tangential unit vectors
            tx = 0
            ty = 0
            tz = 0
            if vt_magn > 0:
                tx = vt_x / vt_magn
                ty = vt_y / vt_magn
                tz = vt_z / vt_magn

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
            # normal force calculation
            k_multi_overlap = self.k * overlap
            eta_multi_normal_velocity = self.eta * v_n
            fn_x = (-k_multi_overlap - eta_multi_normal_velocity) * nx
            fn_y = (-k_multi_overlap - eta_multi_normal_velocity) * ny
            fn_z = (-k_multi_overlap - eta_multi_normal_velocity) * nz
            fn_magn = pow(fn_x * fn_x + fn_y * fn_y + fn_z * fn_z, 1. / 2.)

            # tangential force addition
            ft_x = -(self.kt * d_tang_x[d_idx]) - self.etaT * vt_x
            ft_y = -(self.kt * d_tang_y[d_idx]) - self.etaT * vt_y
            ft_z = -(self.kt * d_tang_z[d_idx]) - self.etaT * vt_z
            ft_magn = pow(ft_x * ft_x + ft_y * ft_y + ft_z * ft_z, 1. / 2.)

            # check for Coulomb criterion
            f_colo = self.mu * fn_magn

            if ft_magn > self.mu * fn_magn:
                ft_x = -f_colo * tx
                ft_y = -f_colo * ty
                ft_z = -f_colo * tz

            # add normal and tangential force to global force
            d_fx[d_idx] += fn_x + ft_x
            d_fy[d_idx] += fn_y + ft_y
            d_fx[d_idx] += fn_z + ft_z

            # torque calculation
            tor_x = (ny * ft_z - nz * ft_y) * d_R[d_idx]
            tor_y = (-nx * ft_z + nz * ft_x) * d_R[d_idx]
            tor_z = (nx * ft_y - ny * ft_x) * d_R[d_idx]

            # add to global torque
            d_torX[d_idx] += tor_x
            d_torY[d_idx] += tor_y
            d_torZ[d_idx] += tor_z

            # post calculation works
            # assign tangential velocity to the particle
            d_vt_x[d_idx] = vt_x
            d_vt_y[d_idx] = vt_y
            d_vt_z[d_idx] = vt_z

        else:
            d_tang_x[d_idx] = 0
            d_tang_y[d_idx] = 0
            d_tang_z[d_idx] = 0


class MakeForcesZero(Equation):
    """Documentation for LinearSpringForce

    """

    def loop(self, d_idx, d_fx, d_fy, d_fz, d_torX, d_torY, d_torZ):
        d_fx[d_idx] = 0
        d_fy[d_idx] = 0
        d_fz[d_idx] = 0

        d_torX[d_idx] = 0
        d_torY[d_idx] = 0
        d_torZ[d_idx] = 0


class DEMStep(IntegratorStep):
    """RK2 step for integrating particles over time

    Use this integrator for DEM formulations. In the predictor step,
    the particles are advanced to `t + dt/2`. The particles are then
    advanced with the new force computed at this position.

    This integrator can be used in PEC or EPEC mode.

    The same integrator can be used for other problems. Like for
    example solid mechanics (see SolidMechStep)

    """

    def initialize(self, d_idx, d_x0, d_y0, d_z0, d_x, d_y, d_z, d_u0, d_v0,
                   d_w0, d_tang_x, d_tang_y, d_tang_z, d_tang_x0, d_tang_y0,
                   d_tang_z0, d_vt_x, d_vt_y, d_vt_z, d_u, d_v, d_w, d_wx,
                   d_wy, d_wz, d_wx0, d_wy0, d_wz0):
        d_x0[d_idx] = d_x[d_idx]
        d_y0[d_idx] = d_y[d_idx]
        d_z0[d_idx] = d_z[d_idx]

        d_u0[d_idx] = d_u[d_idx]
        d_v0[d_idx] = d_v[d_idx]
        d_w0[d_idx] = d_w[d_idx]

        d_wx0[d_idx] = d_wx[d_idx]
        d_wy0[d_idx] = d_wy[d_idx]
        d_wz0[d_idx] = d_wz[d_idx]

        d_tang_x0[d_idx] = d_tang_x[d_idx]
        d_tang_y0[d_idx] = d_tang_y[d_idx]
        d_tang_z0[d_idx] = d_tang_z[d_idx]

    def stage1(self, d_idx, d_x0, d_y0, d_z0, d_x, d_y, d_z, d_u0, d_v0, d_w0,
               d_u, d_v, d_w, d_tang_x, d_tang_y, d_tang_z, d_tang_x0,
               d_tang_y0, d_tang_z0, d_vt_x, d_vt_y, d_vt_z, d_fx, d_fy, d_fz,
               d_torX, d_torY, d_torZ, d_wx, d_wy, d_wz, d_wx0, d_wy0, d_wz0,
               d_m_inverse, d_I_inverse, dt):
        dtb2 = 0.5 * dt
        d_x[d_idx] = d_x0[d_idx] + dtb2 * d_u[d_idx]
        d_y[d_idx] = d_y0[d_idx] + dtb2 * d_v[d_idx]
        d_z[d_idx] = d_z0[d_idx] + dtb2 * d_w[d_idx]

        d_u[d_idx] = d_u0[d_idx] + dtb2 * d_fx[d_idx] * d_m_inverse[d_idx]
        d_v[d_idx] = d_v0[d_idx] + dtb2 * d_fy[d_idx] * d_m_inverse[d_idx]
        d_w[d_idx] = d_w0[d_idx] + dtb2 * d_fz[d_idx] * d_m_inverse[d_idx]

        d_wx[d_idx] = d_wx0[d_idx] + dtb2 * d_torX[d_idx] * d_I_inverse[d_idx]
        d_wy[d_idx] = d_wy0[d_idx] + dtb2 * d_torY[d_idx] * d_I_inverse[d_idx]
        d_wz[d_idx] = d_wz0[d_idx] + dtb2 * d_torZ[d_idx] * d_I_inverse[d_idx]

        d_tang_x[d_idx] = d_tang_x0[d_idx] + dtb2 * d_vt_x[d_idx]
        d_tang_y[d_idx] = d_tang_y0[d_idx] + dtb2 * d_vt_y[d_idx]
        d_tang_z[d_idx] = d_tang_z0[d_idx] + dtb2 * d_vt_z[d_idx]

    def stage2(self, d_idx, d_x0, d_y0, d_z0, d_x, d_y, d_z, d_u0, d_v0, d_w0,
               d_u, d_v, d_w, d_tang_x, d_tang_y, d_tang_z, d_tang_x0,
               d_tang_y0, d_tang_z0, d_vt_x, d_vt_y, d_vt_z, d_fx, d_fy, d_fz,
               d_torX, d_torY, d_torZ, d_wx, d_wy, d_wz, d_wx0, d_wy0, d_wz0,
               d_m_inverse, d_I_inverse, dt):
        d_x[d_idx] = d_x0[d_idx] + dt * d_u[d_idx]
        d_y[d_idx] = d_y0[d_idx] + dt * d_v[d_idx]
        d_z[d_idx] = d_z0[d_idx] + dt * d_w[d_idx]

        d_u[d_idx] = d_u0[d_idx] + dt * d_fx[d_idx] * d_m_inverse[d_idx]
        d_v[d_idx] = d_v0[d_idx] + dt * d_fy[d_idx] * d_m_inverse[d_idx]
        d_w[d_idx] = d_w0[d_idx] + dt * d_fz[d_idx] * d_m_inverse[d_idx]

        d_wx[d_idx] = d_wx0[d_idx] + dt * d_torX[d_idx] * d_I_inverse[d_idx]
        d_wy[d_idx] = d_wy0[d_idx] + dt * d_torY[d_idx] * d_I_inverse[d_idx]
        d_wz[d_idx] = d_wz0[d_idx] + dt * d_torZ[d_idx] * d_I_inverse[d_idx]

        d_tang_x[d_idx] = d_tang_x0[d_idx] + dt * d_vt_x[d_idx]
        d_tang_y[d_idx] = d_tang_y0[d_idx] + dt * d_vt_y[d_idx]
        d_tang_z[d_idx] = d_tang_z0[d_idx] + dt * d_vt_z[d_idx]
