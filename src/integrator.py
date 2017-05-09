class IntegratorStep(object):
    """Subclass this and implement the methods ``initialize``, ``stage1`` etc.
    Use the same conventions as the equations.
    """

    def __repr__(self):
        return '%s()' % (self.__class__.__name__)


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
