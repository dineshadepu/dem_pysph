from pysph.base.utils import get_particle_array
from pysph.sph.equation import Equation


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

    dem_props = ['x0 ', 'y0', 'z0', 'u0', 'v0', 'w0', 'fx', 'fy', 'fz', 'R']

    pa = get_particle_array(constants=constants, additional_props=dem_props,
                            **props)

    # default property arrays to save out.
    pa.set_output_arrays([
        'x', 'y', 'z', 'u', 'v', 'w', 'fx ', 'fy', 'fz', 'm', 'pid', 'gid',
        'tag', 'p'
    ])
    pa.set_output_arrays([
        'x', 'y', 'z', 'u', 'v', 'w', 'm', 'pid', 'gid',
        'tag', 'p'
    ])

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
