import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import pyproj
from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata


def rho_to_vert_geo(lonr, latr, lonp, latp):
    Mm, Lm = lonr.shape
    lon = np.zeros((Mm+1,Lm+1))
    lat = np.zeros((Mm+1,Lm+1))

    lon[1:-1, 1:-1] = lonp[:,:]
    lat[1:-1, 1:-1] = latp[:,:]

    #North edge
    lon[Mm,0:-2] = lonr[Mm-1,0:-1] - ( lonp[Mm-2,:] - lonr[Mm-1,0:-1] )
    lon[Mm,-2:] = lonr[Mm-1,-2:] - ( lonp[Mm-2,-2:] - lonr[Mm-1,-2:] )
    lat[Mm,0:-2] = latr[Mm-1,0:-1] - ( latp[Mm-2,:] - latr[Mm-1,0:-1] )
    lat[Mm,-2:] = latr[Mm-1,-2:] - ( latp[Mm-2,-2:] - latr[Mm-1,-2:] )

    #South edge
    lon[0,0:-2] = lonr[0,0:-1] - ( lonp[0,:] - lonr[0,0:-1] )
    lon[0,-2:] = lonr[0,-2:] - ( lonp[0,-2:] - lonr[0,-2:] )
    lat[0,0:-2] = latr[0,0:-1] - ( latp[0,:] - latr[0,0:-1] )
    lat[0,-2:] = latr[0,-2:] - ( latp[0,-2:] - latr[0,-2:] )

    #East edge
    lon[0:-2,Lm] = lonr[0:-1,Lm-1] - ( lonp[:,Lm-2] - lonr[0:-1,Lm-1] )
    lon[-2:,Lm] = lonr[-2:,Lm-1] - ( lonp[-2:,Lm-2] - lonr[-2:,Lm-1] )
    lat[0:-2,Lm] = latr[0:-1,Lm-1] - ( latp[:,Lm-2] - latr[0:-1,Lm-1] )
    lat[-2:,Lm] = latr[-2:,Lm-1] - ( latp[-2:,Lm-2] - latr[-2:,Lm-1] )

    #West edge
    lon[0:-2,0] = lonr[0:-1,0] - ( lonp[:,0] - lonr[0:-1,0] )
    lon[-2:,0] = lonr[-2:,0] - ( lonp[-2:,0] - lonr[-2:,0] )
    lat[0:-2,0] = latr[0:-1,0] - ( latp[:,0] - latr[0:-1,0] )
    lat[-2:,0] = latr[-2:,0] - ( latp[-2:,0] - latr[-2:,0] )

    return lon, lat

def get_ROMS_hgrid(gridid):
    """
    hgrid = get_ROMS_hgrid(gridid)
    Load ROMS horizontal grid object
    """

    nc = Dataset(gridid)

    #Check for cartesian or geographical grid
    spherical = nc.variables['spherical'][0]

    #if it is type byte, then convert to string
    try:
      spherical=spherical.decode('utf8')
    except:
      print('Assuming spherical is integer',spherical, type(spherical))

    #Get horizontal grid
    if ((spherical == 0) or (spherical == 'F')):
        #cartesian grid
        print('Load cartesian grid from file')
        if 'x_vert' in list(nc.variables.keys()) and 'y_vert' in list(nc.variables.keys()):
            x_vert = nc.variables['x_vert'][:]
            y_vert = nc.variables['y_vert'][:]
        elif 'x_rho' in list(nc.variables.keys()) and 'y_rho' in list(nc.variables.keys()) \
                 and 'pm' in list(nc.variables.keys()) and 'pn' in list(nc.variables.keys()):
            x_rho = nc.variables['x_rho'][:]
            y_rho = nc.variables['y_rho'][:]
            pm = nc.variables['pm'][:]
            pn = nc.variables['pn'][:]
            try: angle = nc.variables['angle'][:]
            except: angle = np.zeros(x_rho.shape)
            #compute verts from rho point, pm, pn, angle
            x_vert, y_vert = rho_to_vert(x_rho, y_rho, pm, pn, angle)
        else:
            raise ValueError('NetCDF file must contain x_vert and y_vert \
                     or x_rho, y_rho, pm, pn and angle for a cartesian grid')

        if 'x_rho' in list(nc.variables.keys()) and 'y_rho' in list(nc.variables.keys()) and \
             'x_u' in list(nc.variables.keys()) and 'y_u' in list(nc.variables.keys()) and \
             'x_v' in list(nc.variables.keys()) and 'y_v' in list(nc.variables.keys()) and \
             'x_psi' in list(nc.variables.keys()) and 'y_psi' in list(nc.variables.keys()):
            x_rho = nc.variables['x_rho'][:]
            y_rho = nc.variables['y_rho'][:]
            x_u = nc.variables['x_u'][:]
            y_u = nc.variables['y_u'][:]
            x_v = nc.variables['x_v'][:]
            y_v = nc.variables['y_v'][:]
            x_psi = nc.variables['x_psi'][:]
            y_psi = nc.variables['y_psi'][:]
        else:
            x_rho = None
            y_rho = None
            x_u = None
            y_u = None
            x_v = None
            y_v = None
            x_psi = None
            y_psi = None

        if 'pm' in list(nc.variables.keys()) and 'pn' in list(nc.variables.keys()):
            pm = nc.variables['pm'][:]
            dx = 1. / pm
            pn = nc.variables['pn'][:]
            dy = 1. / pn
        else:
            dx = None
            dy = None

        if 'dndx' in list(nc.variables.keys()) and 'dmde' in list(nc.variables.keys()):
            dndx = nc.variables['dndx'][:]
            dmde = nc.variables['dmde'][:]
        else:
            dndx = None
            dmde = None

        if 'angle' in list(nc.variables.keys()):
            angle = nc.variables['angle'][:]
        else:
            angle = None

        #Get cartesian grid
        hgrd = CGrid(x_vert, y_vert, x_rho=x_rho, y_rho=y_rho, \
                     x_u=x_u, y_u=y_u, x_v=x_v, y_v=y_v, \
                     x_psi=x_psi, y_psi=y_psi, dx=dx, dy=dy, \
                     dndx=dndx, dmde=dmde, angle_rho=angle)

        #load the mask
        try:
            hgrd.mask_rho = np.array(nc.variables['mask_rho'][:])
        except:
            hgrd.mask_rho = np.ones(hgrd.x_rho.shape)

    else:
        #geographical grid
        print('Load geographical grid from file')
        proj = Basemap(projection='merc', resolution=None, lat_0=0, lon_0=0)
        if 'lon_vert' in list(nc.variables.keys()) and 'lat_vert' in list(nc.variables.keys()):
            lon_vert = nc.variables['lon_vert'][:]
            lat_vert = nc.variables['lat_vert'][:]
        elif 'lon_rho' in list(nc.variables.keys()) and 'lat_rho' in list(nc.variables.keys()) \
                and 'lon_psi' in list(nc.variables.keys()) and 'lat_psi' in list(nc.variables.keys()):
            lon_rho = nc.variables['lon_rho'][:]
            lat_rho = nc.variables['lat_rho'][:]
            lon_psi = nc.variables['lon_psi'][:]
            lat_psi = nc.variables['lat_psi'][:]
            #compute verts from rho and psi point
            lon_vert, lat_vert = rho_to_vert_geo(lon_rho, lat_rho, lon_psi, lat_psi)
        else:
            raise ValueError('NetCDF file must contain lon_vert and lat_vert \
                  or lon_rho, lat_rho, lon_psi, lat_psi for a geographical grid')

        if 'lon_rho' in list(nc.variables.keys()) and 'lat_rho' in list(nc.variables.keys()) and \
              'lon_u' in list(nc.variables.keys()) and 'lat_u' in list(nc.variables.keys()) and \
              'lon_v' in list(nc.variables.keys()) and 'lat_v' in list(nc.variables.keys()) and \
              'lon_psi' in list(nc.variables.keys()) and 'lat_psi' in list(nc.variables.keys()):
            lon_rho = nc.variables['lon_rho'][:]
            lat_rho = nc.variables['lat_rho'][:]
            lon_u = nc.variables['lon_u'][:]
            lat_u = nc.variables['lat_u'][:]
            lon_v = nc.variables['lon_v'][:]
            lat_v = nc.variables['lat_v'][:]
            lon_psi = nc.variables['lon_psi'][:]
            lat_psi = nc.variables['lat_psi'][:]
        else:
            lon_rho = None
            lat_rho = None
            lon_u = None
            lat_u = None
            lon_v = None
            lat_v = None
            lon_psi = None
            lat_psi = None

        if 'pm' in list(nc.variables.keys()) and 'pn' in list(nc.variables.keys()):
            pm = nc.variables['pm'][:]
            dx = 1. / pm
            pn = nc.variables['pn'][:]
            dy = 1. / pn
        else:
            dx = None
            dy = None

        if 'dndx' in list(nc.variables.keys()) and 'dmde' in list(nc.variables.keys()):
            dndx = nc.variables['dndx'][:]
            dmde = nc.variables['dmde'][:]
        else:
            dndx = None
            dmde = None

        if 'angle' in list(nc.variables.keys()):
            angle = nc.variables['angle'][:]
        else:
            angle = None

        #Get geographical grid
        hgrd = CGrid_geo(lon_vert, lat_vert, proj, \
                         lon_rho=lon_rho, lat_rho=lat_rho, \
                         lon_u=lon_u, lat_u=lat_u, lon_v=lon_v, lat_v=lat_v, \
                         lon_psi=lon_psi, lat_psi=lat_psi, dx=dx, dy=dy, \
                         dndx=dndx, dmde=dmde, angle_rho=angle)

        #load the mask
        try:
            hgrd.mask_rho = np.array(nc.variables['mask_rho'][:])
        except:
            hgrd.mask_rho = np.ones(hgrd.lat_rho.shape)

    return hgrd

def get_ROMS_vgrid(gridid, zeta=None):
    """
    vgrid = get_ROMS_vgrid(gridid)
    Load ROMS vertical grid object. vgrid is a s_coordinate or
    a z_coordinate object, depending on gridid.grdtype.
    vgrid.z_r and vgrid.z_w (vgrid.z for a z_coordinate object)
    can be indexed in order to retreive the actual depths. The
    free surface time serie zeta can be provided as an optional
    argument. Note that the values of zeta are not calculated
    until z is indexed, so a netCDF variable for zeta may be passed,
    even if the file is large, as only the values that are required
    will be retrieved from the file.
    """

    nc = Dataset(gridid)

    #Get vertical grid
    try:
        h = nc.variables['h'][:]
    except:
        raise ValueError('NetCDF file must contain the bathymetry h')

    try:
        hraw = nc.variables['hraw'][:]
    except:
        hraw = None

    gridinfo_type = 'roms'
    if gridinfo_type == 'roms':
        Vtrans = 4
        theta_b = 0.4
        theta_s = 5.0
        Tcline = 100
        N = 40
        if Vtrans == 1:
            vgrid = s_coordinate(h, theta_b, theta_s, Tcline, N, hraw=hraw, zeta=zeta)
        elif Vtrans == 2:
            vgrid = s_coordinate_2(h, theta_b, theta_s, Tcline, N, hraw=hraw, zeta=zeta)
        elif Vtrans == 4:
            vgrid = s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw, zeta=zeta)
        elif Vtrans == 5:
            vgrid = s_coordinate_5(h, theta_b, theta_s, Tcline, N, hraw=hraw, zeta=zeta)
        else:
            raise Warning('Unknown vertical transformation Vtrans')

    elif  gridinfo_type == 'z':
        N = gridinfo.N
        depth = gridinfo.depth
        vgrid = z_coordinate(h, depth, N)

    else:
        raise ValueError('Unknown grid type')

    return vgrid


def get_ROMS_grid(gridid, zeta=None, hist_file=None,grid_file=None):
    #ROMS_gridinfo above.
    hgrd = get_ROMS_hgrid(gridid)
    vgrid = get_ROMS_vgrid(gridid, zeta=zeta)

    name='nenhum'
    #Get ROMS grid
    return ROMS_Grid(name, hgrd, vgrid)

class CGrid(object):
    """
    Curvilinear Arakawa C-Grid
    The basis for the CGrid class are two arrays defining the verticies of the
    grid in Cartesian (for geographic coordinates, see CGrid_geo). An optional
    mask may be defined on the cell centers. Other Arakawa C-grid properties,
    such as the locations of the cell centers (rho-points), cell edges (u and
    v velocity points), cell widths (dx and dy) and other metrics (angle,
    dmde, and dndx) are all calculated internally from the vertex points.
    Input vertex arrays may be either type np.array or np.ma.MaskedArray. If
    masked arrays are used, the mask will be a combination of the specified
    mask (if given) and the masked locations.
    EXAMPLES:
    --------
    >>> x, y = mgrid[0.0:7.0, 0.0:8.0]
    >>> x = np.ma.masked_where( (x<3) & (y<3), x)
    >>> y = np.ma.MaskedArray(y, x.mask)
    >>> grd = pyroms.grid.CGrid(x, y)
    >>> print(grd.x_rho)
    [[-- -- -- 0.5 0.5 0.5 0.5]
     [-- -- -- 1.5 1.5 1.5 1.5]
     [-- -- -- 2.5 2.5 2.5 2.5]
     [3.5 3.5 3.5 3.5 3.5 3.5 3.5]
     [4.5 4.5 4.5 4.5 4.5 4.5 4.5]
     [5.5 5.5 5.5 5.5 5.5 5.5 5.5]]
    >>> print(grd.mask)
    [[ 0.  0.  0.  1.  1.  1.  1.]
     [ 0.  0.  0.  1.  1.  1.  1.]
     [ 0.  0.  0.  1.  1.  1.  1.]
     [ 1.  1.  1.  1.  1.  1.  1.]
     [ 1.  1.  1.  1.  1.  1.  1.]
     [ 1.  1.  1.  1.  1.  1.  1.]]
    """

    def __init__(self, x_vert, y_vert, x_rho=None, y_rho=None, x_u=None, y_u=None, x_v=None, y_v=None, \
                    x_psi=None, y_psi=None, dx=None, dy=None, dndx=None, dmde=None, angle_rho=None):

        assert np.ndim(x_vert)==2 and np.ndim(y_vert)==2 and np.shape(x_vert)==np.shape(y_vert), \
            'x and y must be 2D arrays of the same size.'

        if np.any(np.isnan(x_vert)) or np.any(np.isnan(y_vert)):
            x_vert = np.ma.masked_where( (isnan(x_vert)) | (isnan(y_vert)) , x_vert)
            y_vert = np.ma.masked_where( (isnan(x_vert)) | (isnan(y_vert)) , y_vert)

        self.x_vert = x_vert
        self.y_vert = y_vert

        self.f = None
        self.spherical = 'F'

        mask_shape = tuple([n-1 for n in self.x_vert.shape])
        self.mask_rho = np.ones(mask_shape, dtype='d')

        # If maskedarray is given for verticies, modify the mask such that
        # non-existant grid points are masked.  A cell requires all four
        # verticies to be defined as a water point.
        if isinstance(self.x_vert, np.ma.MaskedArray):
            mask = (self.x_vert.mask[:-1,:-1] | self.x_vert.mask[1:,:-1] | \
                    self.x_vert.mask[:-1,1:] | self.x_vert.mask[1:,1:])
            self.mask_rho = np.asarray(~(~np.bool_(self.mask_rho) | mask), dtype='d')

        if isinstance(self.y_vert, np.ma.MaskedArray):
            mask = (self.y_vert.mask[:-1,:-1] | self.y_vert.mask[1:,:-1] | \
                    self.y_vert.mask[:-1,1:] | self.y_vert.mask[1:,1:])
            self.mask_rho = np.asarray(~(~np.bool_(self.mask_rho) | mask), dtype='d')

        if x_rho is None or y_rho is None or x_u is None or y_u is None or \
            x_v is None or y_v is None or x_psi is None or y_psi is None:
            self._calculate_subgrids()
        else:
            self.x_rho = x_rho
            self.y_rho = y_rho
            self.x_u = x_u
            self.y_u = y_u
            self.x_v = x_v
            self.y_v = y_v
            self.x_psi = x_psi
            self.y_psi = y_psi

        if dx is None or dy is None:
            self._calculate_metrics()
        else:
            self.dx = dx
            self.dy = dy

        self.xl = np.maximum(self.dx[0,:].sum(), self.dx[-1,:].sum())
        self.el = np.maximum(self.dy[:,0].sum(), self.dy[:,-1].sum())

        if dndx is None or dmde is None:
            self._calculate_derivative_metrics()
        else:
            self.dndx = dndx
            self.dmde = dmde

        if angle_rho is None:
            self._calculate_angle_rho()
        else:
            self.angle_rho = angle_rho

        self._calculate_angle()


    def _calculate_subgrids(self):
        self.x_rho = 0.25*(self.x_vert[1:,1:]+self.x_vert[1:,:-1]+ \
                           self.x_vert[:-1,1:]+self.x_vert[:-1,:-1])
        self.y_rho = 0.25*(self.y_vert[1:,1:]+self.y_vert[1:,:-1]+ \
                           self.y_vert[:-1,1:]+self.y_vert[:-1,:-1])
        self.x_u = 0.5*(self.x_vert[:-1,1:-1] + self.x_vert[1:,1:-1])
        self.y_u = 0.5*(self.y_vert[:-1,1:-1] + self.y_vert[1:,1:-1])
        self.x_v = 0.5*(self.x_vert[1:-1,:-1] + self.x_vert[1:-1,1:])
        self.y_v = 0.5*(self.y_vert[1:-1,:-1] + self.y_vert[1:-1,1:])
        self.x_psi = self.x_vert[1:-1,1:-1]
        self.y_psi = self.y_vert[1:-1,1:-1]

    def _calculate_metrics(self):
        'Calculates pm, pn, dndx, dmde from x_vert and y_vert'
        x_temp = 0.5*(self.x_vert[1:,:]+self.x_vert[:-1,:])
        y_temp = 0.5*(self.y_vert[1:,:]+self.y_vert[:-1,:])
        self.dx = np.sqrt(np.diff(x_temp, axis=1)**2 + np.diff(y_temp, axis=1)**2)
        x_temp = 0.5*(self.x_vert[:,1:]+self.x_vert[:,:-1])
        y_temp = 0.5*(self.y_vert[:,1:]+self.y_vert[:,:-1])
        self.dy = np.sqrt(np.diff(x_temp, axis=0)**2 + np.diff(y_temp, axis=0)**2)

    def _calculate_derivative_metrics(self):
        if isinstance(self.dy, np.ma.MaskedArray):
            self.dndx = np.ma.zeros(self.x_rho.shape, dtype='d')
        else:
            self.dndx = np.zeros(self.x_rho.shape, dtype='d')

        if isinstance(self.dx, np.ma.MaskedArray):
            self.dmde = np.ma.zeros(self.x_rho.shape, dtype='d')
        else:
            self.dmde = np.zeros(self.x_rho.shape, dtype='d')

        self.dndx[1:-1,1:-1] = 0.5*(self.dy[1:-1,2:] - self.dy[1:-1,:-2])
        self.dmde[1:-1,1:-1] = 0.5*(self.dx[2:,1:-1] - self.dx[:-2,1:-1])

    def _calculate_angle(self):
        if isinstance(self.x_vert, np.ma.MaskedArray) or \
           isinstance(self.y_vert, np.ma.MaskedArray):
            self.angle = np.ma.zeros(self.x_vert.shape, dtype='d')
        else:
            self.angle = np.zeros(self.x_vert.shape, dtype='d')

        angle_ud = np.arctan2(np.diff(self.y_vert, axis=1), np.diff(self.x_vert, axis=1))
        angle_lr = np.arctan2(np.diff(self.y_vert, axis=0), np.diff(self.x_vert, axis=0)) - np.pi/2.0
        # domain center
        self.angle[1:-1,1:-1] = 0.25*(angle_ud[1:-1,1:]+angle_ud[1:-1,:-1]\
                                     +angle_lr[1:,1:-1]+angle_lr[:-1,1:-1])
        # edges
        self.angle[0,1:-1] = (1.0/3.0)*(angle_lr[0,1:-1]+angle_ud[0,1:]+angle_ud[0,:-1])
        self.angle[-1,1:-1] = (1.0/3.0)*(angle_lr[-1,1:-1]+angle_ud[-1,1:]+angle_ud[-1,:-1])
        self.angle[1:-1,0] = (1.0/3.0)*(angle_ud[1:-1,0]+angle_lr[1:,0]+angle_lr[:-1,0])
        self.angle[1:-1,-1] = (1.0/3.0)*(angle_ud[1:-1,-1]+angle_lr[1:,-1]+angle_lr[:-1,-1])
        #conrers
        self.angle[0,0] = 0.5*(angle_lr[0,0]+angle_ud[0,0])
        self.angle[0,-1] = 0.5*(angle_lr[0,-1]+angle_ud[0,-1])
        self.angle[-1,0] = 0.5*(angle_lr[-1,0]+angle_ud[-1,0])
        self.angle[-1,-1] = 0.5*(angle_lr[-1,-1]+angle_ud[-1,-1])

    def _calculate_angle_rho(self):
        self.angle_rho = np.arctan2(np.diff(0.5*(self.y_vert[1:,:]+self.y_vert[:-1,:])), \
                                    np.diff(0.5*(self.x_vert[1:,:]+self.x_vert[:-1,:])))

    def calculate_orthogonality(self):
        '''
        Calculate orthogonality error in radians
        '''
        z = self.x_vert + 1j*self.y_vert
        du = np.diff(z, axis=1); du = (du/abs(du))[:-1,:]
        dv = np.diff(z, axis=0); dv = (dv/abs(dv))[:,:-1]
        ang1 = np.arccos(du.real*dv.real + du.imag*dv.imag)
        du = np.diff(z, axis=1); du = (du/abs(du))[1:,:]
        dv = np.diff(z, axis=0); dv = (dv/abs(dv))[:,:-1]
        ang2 = np.arccos(du.real*dv.real + du.imag*dv.imag)
        du = np.diff(z, axis=1); du = (du/abs(du))[:-1,:]
        dv = np.diff(z, axis=0); dv = (dv/abs(dv))[:,1:]
        ang3 = np.arccos(du.real*dv.real + du.imag*dv.imag)
        du = np.diff(z, axis=1); du = (du/abs(du))[1:,:]
        dv = np.diff(z, axis=0); dv = (dv/abs(dv))[:,1:]
        ang4 = np.arccos(du.real*dv.real + du.imag*dv.imag)
        ang = np.mean([abs(ang1), abs(ang2), abs(ang3), abs(ang4)], axis=0)
        ang = (ang-np.pi/2.0)
        return ang

    def mask_polygon(self, polyverts, mask_value=0.0):
        """
        Mask Cartesian points contained within the polygon defined by polyverts
        A cell is masked if the cell center (x_rho, y_rho) is within the
        polygon. Other sub-masks (mask_u, mask_v, and mask_psi) are updated
        automatically.
        mask_value [=0.0] may be specified to alter the value of the mask set
        within the polygon.  E.g., mask_value=1 for water points.
        """

        polyverts = np.asarray(polyverts)
        assert polyverts.ndim == 2, \
            'polyverts must be a 2D array, or a similar sequence'
        assert polyverts.shape[1] == 2, \
            'polyverts must be two columns of points'
        assert polyverts.shape[0] > 2, \
            'polyverts must contain at least 3 points'

        mask = self.mask_rho
        #inside = points_inside_poly(
        #    np.vstack( (self.x_rho.flatten(), self.y_rho.flatten()) ).T,
        #    polyverts)
        path = mpl.path.Path(polyverts)
        inside = path.contains_points(np.vstack( (self.x_rho.flatten(), self.y_rho.flatten()) ).T)
        if np.any(inside):
            self.mask_rho.flat[inside] = mask_value

    def _get_mask_u(self):
        return self.mask_rho[:,1:]*self.mask_rho[:,:-1]

    def _get_mask_v(self):
        return self.mask_rho[1:,:]*self.mask_rho[:-1,:]

    def _get_mask_psi(self):
        return self.mask_rho[1:,1:]*self.mask_rho[:-1,1:]* \
               self.mask_rho[1:,:-1]*self.mask_rho[:-1,:-1]

    def _set_mask_rho(self, mask_rho):
        self.mask_rho = mask_rho

    x = property(lambda self: self.x_vert, None, None, 'Return x_vert')
    y = property(lambda self: self.y_vert, None, None, 'Return x_vert')
    mask = property(lambda self: self.mask_rho, _set_mask_rho, None, 'Return mask_rho')
    mask_u   = property(_get_mask_u, None, None, 'Return mask_u')
    mask_v   = property(_get_mask_v, None, None, 'Return mask_v')
    mask_psi = property(_get_mask_psi, None, None, 'Return mask_psi')


class CGrid_geo(CGrid):
    """
    Curvilinear Arakawa C-grid defined in geographic coordinates
    For a geographic grid, a projection may be specified, or The default
    projection for will be defined by the matplotlib.toolkits.Basemap
    projection:
    proj = Basemap(projection='merc', resolution=None, lat_ts=0.0)
    For a geographic grid, the cell widths are determined by the great
    circle distances. Angles, however, are defined using the projected
    coordinates, so a projection that conserves angles must be used. This
    means typically either Mercator (projection='merc') or Lambert
    Conformal Conic (projection='lcc').
    """
    def _calculate_metrics(self):
        # calculate metrics based on x and y grid
        super(CGrid_geo, self)._calculate_metrics()

        # optionally calculate dx and dy based on great circle distances
        # for more accurate cell sizes.
        if self.use_gcdist:
            geod = pyproj.Geod(ellps=self.ellipse)
            az_forward, az_back, dx = geod.inv(self.lon[:,1:],  self.lat[:,1:], \
                                               self.lon[:,:-1], self.lat[:,:-1])
            self.dx = 0.5*(dx[1:,:]+dx[:-1,:])
            self.pm = 1.0/self.dx
            az_forward, az_back, dy = geod.inv(self.lon[1:,:],  self.lat[1:,:], \
                                               self.lon[:-1,:], self.lat[:-1,:])
            self.dy = 0.5*(dy[:,1:]+dy[:,:-1])
            self.pn = 1.0/self.dy


    def _calculate_derivative_metrics(self):
        if isinstance(self.dy, np.ma.MaskedArray):
            self.dndx = np.ma.zeros(self.dy.shape, dtype='d')
        else:
            self.dndx = np.zeros(self.dy.shape, dtype='d')

        if isinstance(self.dx, np.ma.MaskedArray):
            self.dmde = np.ma.zeros(self.dx.shape, dtype='d')
        else:
            self.dmde = np.zeros(self.dx.shape, dtype='d')

        self.dndx[1:-1,1:-1] = 0.5*(self.dy[1:-1,2:] - self.dy[1:-1,:-2])
        self.dmde[1:-1,1:-1] = 0.5*(self.dx[2:,1:-1] - self.dx[:-2,1:-1])

    def _calculate_angle_rho(self):
        if isinstance(self.lon, np.ma.MaskedArray) or \
           isinstance(self.lat, np.ma.MaskedArray):
            self.angle_rho = np.ma.zeros(self.lon.shape, dtype='d')
        else:
            self.angle_rho = np.zeros(self.lon.shape, dtype='d')

        # calculate metrics based on x and y grid
        super(CGrid_geo, self)._calculate_angle_rho()

        # optionally calculate dx and dy based on great circle distances
        # for more accurate cell sizes.
        if self.use_gcdist:
            geod = pyproj.Geod(ellps=self.ellipse)
            az_forward, az_back, dx = geod.inv(self.lon[:,:-1], self.lat[:,:-1], \
                                               self.lon[:,1:], self.lat[:,1:])

            angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
            self.angle_rho = (90 - angle) * np.pi/180.


    def __init__(self, lon_vert, lat_vert, proj, use_gcdist=True, ellipse='WGS84', \
                    lon_rho=None, lat_rho=None, lon_u=None, lat_u=None, \
                    lon_v=None, lat_v=None, lon_psi=None, lat_psi=None, dx=None, dy=None, \
                    dndx=None, dmde=None, angle_rho=None):

        x, y = proj(lon_vert, lat_vert)
        self.lon_vert = lon_vert
        self.lat_vert = lat_vert
        self.proj = proj

        self.use_gcdist = use_gcdist
        self.ellipse = ellipse

        if lon_rho is None or lat_rho is None or lon_u is None or lat_u is None or \
             lon_v is None or lat_v is None or lon_psi is None or lat_psi is None:

            super(CGrid_geo, self).__init__(x, y)

            self.lon_rho, self.lat_rho = self.proj(self.x_rho, self.y_rho,
                                                   inverse=True)
            self.lon_u, self.lat_u = self.proj(self.x_u, self.y_u, inverse=True)
            self.lon_v, self.lat_v = self.proj(self.x_v, self.y_v, inverse=True)
            self.lon_psi, self.lat_psi = self.proj(self.x_psi, self.y_psi,
                                                   inverse=True)
        else:
            self.lon_rho = lon_rho
            self.lat_rho = lat_rho
            self.lon_u = lon_u
            self.lat_u = lat_u
            self.lon_v = lon_v
            self.lat_v = lat_v
            self.lon_psi = lon_psi
            self.lat_psi = lat_psi
            #calculate cartesian position
            self.x_vert, self.y_vert = proj(lon_vert, lat_vert)
            self.x_rho, self.y_rho = proj(lon_rho, lat_rho)
            self.x_u, self.y_u = proj(lon_u, lat_u)
            self.x_v, self.y_v = proj(lon_v, lat_v)
            self.x_psi, self.y_psi = proj(lon_psi, lat_psi)

        if dx is None or dy is None:
            self._calculate_metrics()
        else:
            self.dx = dx
            self.dy = dy

        self.xl = np.maximum(self.dx[0,:].sum(), self.dx[-1,:].sum())
        self.el = np.maximum(self.dy[:,0].sum(), self.dy[:,-1].sum())

        if dndx is None or dmde is None:
            self._calculate_derivative_metrics()
        else:
            self.dndx = dndx
            self.dmde = dmde

        if angle_rho is None:
            self._calculate_angle_rho()
        else:
            self.angle_rho = angle_rho

        self.f = 2.0 * 7.29e-5 * np.sin(self.lat_rho * np.pi / 180.0)
        self.spherical = 'T'


    def mask_polygon_geo(lonlat_verts, mask_value=0.0):
        lon, lat = list(zip(*lonlat_verts))
        x, y = proj(lon, lat, inverse=True)
        self.mask_polygon(list(zip(x, y)), mask_value)

    lon = property(lambda self: self.lon_vert, None, None, 'Shorthand for lon_vert')
    lat = property(lambda self: self.lat_vert, None, None, 'Shorthand for lat_vert')
    
    
class ROMS_Grid(object):
    """
    grd = ROMS_Grid(hgrid, vgrid)
    ROMS Grid object combining horizontal and vertical grid
    """

    def __init__(self, name, hgrid=CGrid, vgrid=s_coordinate):
        self.name = name
        self.hgrid = hgrid
        self.vgrid = vgrid


class CGrid_TPXO8(object):

    # CGrid object for TPXO8 v1

    def __init__(self, lon_t, lat_t, lon_u, lat_u, lon_v, lat_v, \
                 mask_t, mask_u, mask_v, z_t, z_u, z_v, missing_value, name, xrange, yrange):

        self.name = name

        self.missing_value = -9999.

        self.xrange = xrange
        self.yrange = yrange

        self.z_t = z_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.z_u = z_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.z_v = z_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_t = lon_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_t = lat_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_u = lon_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_u = lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lon_v = lon_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_v = lat_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_t_vert = 0.5 * (lon_t[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lon_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_t_vert = 0.5 * (lat_t[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lat_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.lon_u_vert = 0.5 * (lon_u[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lon_u[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_u_vert = 0.5 * (lat_u[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lat_u[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lon_v_vert = 0.5 * (lon_v[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lon_v[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_v_vert = 0.5 * (lat_v[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
        lat_v[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.mask_t = mask_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_u = mask_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_v = mask_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        ones = np.ones(self.z_t.shape)
        a1 = lat_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] - \
        lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        a2 = lon_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] - \
        lon_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        a3 = 0.5*(lat_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] + \
        lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1])
        a2 = np.where(a2 > 180*ones, a2 - 360*ones, a2)
        a2 = np.where(a2 < -180*ones, a2 + 360*ones, a2)
        a2 = a2 * np.cos(np.pi/180.*a3)
        self.angle = np.arctan2(a1, a2)

def get_nc_CGrid_TPXO8(grdfile, name='TPXO8', \
                       xrange=(4350, 4550), yrange=(6600, 6900), missing_value=-9999):
    """
    grd = get_nc_CGrid_TPXO8(grdfile)
    Load Cgrid object for TPXO8 from netCDF file
    """

    nc = Dataset(grdfile)

    lonh = nc.variables['lon_z'][:]
    lath = nc.variables['lat_z'][:]
    lonu = nc.variables['lon_u'][:]
    latu = nc.variables['lat_u'][:]
    lonv = nc.variables['lon_v'][:]
    latv = nc.variables['lat_v'][:]

    zh = nc.variables['hz'][:]
    zu = nc.variables['hu'][:]
    zv = nc.variables['hv'][:]
    nc.close()

    # land mask
    h_msk = zh!=0
    u_msk = zu!=0
    v_msk = zv!=0

    # longitude from -180 to 180
    lonh[lonh>180] = lonh[lonh>180]-360
    lonu[lonu>180] = lonu[lonu>180]-360
    lonv[lonv>180] = lonv[lonv>180]-360

    lathh, lonhh = np.meshgrid(lath, lonh)
    latuu, lonuu = np.meshgrid(latu, lonu)
    latvv, lonvv = np.meshgrid(latv, lonv)

    # generate tpxo8 grid
    # xrange = [4400, 4600]
    # yrange = [6600, 6900]
    return CGrid_TPXO8(lonhh, lathh, lonuu, latuu, lonvv, latvv, \
                                   h_msk, u_msk, v_msk, \
                                   zh, zu, zv, missing_value, 'TPXO8', xrange, yrange)
                                   
                                                                
def horiz_interp_tide(variable_pgrd, lon_pgrd, lat_pgrd,lon_cgrd, lat_cgrd):
  interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())  
  points = (lat_pgrd.ravel(), lon_pgrd.ravel())	
  ndir=griddata(points, variable_pgrd.ravel(), interp_points, method='linear').reshape(lon_cgrd.shape)
  ff=np.isnan(ndir)     
  if ff.mean() != 0:
    points=(lat_cgrd[~ff], lon_cgrd[~ff])	
    interp_points = (lat_cgrd[ff], lon_cgrd[ff])   
    ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')	
  return ndir