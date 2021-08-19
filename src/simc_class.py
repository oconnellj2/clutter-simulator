"""
    simc_class.py
    =============
    Purpose:
    --------
        Defines objects, allowing new instances of these to be made in
        `simc.py` and `simc_utils.py`.
    Objects:
    ----------------
        Loc: A location in 3D space.
        Vec: Vector in a 3D Cartesian coordinate system.
        Grid: A 2D array of Loc objects, representing a grid.
        Path: A list of loc objects, representing a path.
        Facet: Holds a triangular facet made of thee loc objects.
        Dem: Holds an array of DEM data.
    Contributing:
    -------------
        Univerisity of Arizona | Lunar & Planetary Laboratory
        TAPIR | Terrestrial And Planetary Investigation with Radar
        James O'Connell | Research Assistant
"""
from simc_utils import *


class Loc:
    """
     A location in 3D space, either in cartesian or geodetic coordinates
     The two are switched between often in this program so they are just
     combined. X is the same field as longitude, Y is the same field as
     latitude, Z is the same field as the radius.
    """
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
        self.ol = False
        self.nd = False
    
    def __eq__(self, other):
        return (self.x == other.x and self.y == other.y and self.z == other.z)
    
    def __ne__(self, other):
        return not self.__eq__(other)
   
    def __str__(self):
        return('(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')')

    def __repr__(self):
        return str([self.x, self.y, self.z])
    
    def __add__(self,vec):
        return Loc(self.x + vec.i, self.y + vec.j, self.z + vec.k)
    
    def __radd__(self,vec):
        return self.__add__(vec)
    
    def __sub__(self,vec):
        return Loc(self.x - vec.i, self.y - vec.j, self.z - vec.k)
    
    def __rsub__(self,vec):
        return self.__add__(vec)

    def __getitem__(self, i):
        return self.locList()[i]
  
    def locList(self):
        return [self.x, self.y, self.z]
    
    def copy(self):
        return Loc(self.x, self.y, self.z)
    
    def topix(self,dem):
        """
        This transforms a point to an (x,y) pixel location on a DEM using the
        input geotransform. `IF YOU USE THIS ON ITS OWN WITHOUT THE toground()
        FUNCTION FOR A Pointlist MAKE SURE THAT THE POINT COORDINATES ARE IN THE
        SAME COORDINATE SYSTEM AS THE GEOTRANSFORM.`
        """
        gt = dem.gt
        x = int((gt[0] - self.x)/-gt[1])
        y = int((gt[3] - self.y)/-gt[5])

        ## Bump bad points back onto DEM - THIS IS DANGEROUS!!!!
        #if(x >= dem.data.shape[1]):
        #    x = dem.data.shape[1] - 1

        #if(y >= dem.data.shape[0]):
        #    y = dem.data.shape[0] - 1
        
        #if(x < 0):
        #    x = 0
        
        #if(y < 0):
        #    y = 0

        if(x<0 or y <0 or x >= dem.data.shape[1] or y >= dem.data.shape[0]):
            return Loc(-1,-1,0)
    
        # SPECIFIC FOR MOLA.
        out = Loc(x,y,0)

        return out

    def topix_write(self,gt):
        """
        This transforms a point to an (x,y) pixel location on a DEM using the
        input geotransform. `IF YOU USE THIS ON ITS OWN WITHOUT THE toground()
        FUNCTION FOR A Pointlist MAKE SURE THAT THE POINT COORDINATES ARE IN
        THE SAME COORDINATE SYSTEM AS THE GEOTRANSFORM.`
        """
        
        x = int((self.x - gt[0])/gt[1])
        y = int((self.y - gt[3])/gt[5])
        
        return [x,y]

class Vec:
    """Vector in a 3D Cartesian coordinate system."""
    def __init__(self,l0,l1=None):
        if not l1:
            self.i = l0[0]
            self.j = l0[1]
            self.k = l0[2]
        else:
            self.i = l1.x-l0.x
            self.j = l1.y-l0.y
            self.k = l1.z-l0.z

    def __repr__(self):
        return str([self.i, self.j, self.k])
    
    def __mul__(self, c):
        return Vec((c*self.i, c*self.j, c*self.k))
    
    def __rmul__(self,c):
        return self.__mul__(c)
    
    def __neg__(self):
        return Vec((-self.i, -self.j, -self.k))
    
    def vlist(self):
        return [self.i, self.j, self.k]

    def vlen(self):
        """ Calculates the magnitude of a vector point(value). """
        try:
            return (self.i**2 + self.j**2 + self.k**2)**.5
        except:
            raise SystemExit("\033[91m\033[1m" + str(self.i, self.j, self.k))

    def xprd(self,vec1):
        """ Cross product of two Vec objects."""
        i = self.j*vec1.k-self.k*vec1.j
        j = -(self.i*vec1.k-self.k*vec1.i)
        k = self.i*vec1.j-self.j*vec1.i
        return Vec((i,j,k))
    
    def dprd(self, vec1):
        """ Dot product of two Vec objects. """
        return self.i*vec1.i + self.j*vec1.j + self.k*vec1.k
  
    def unit(self):
        """ Return a unit vector in the direction of self."""
        l = self.vlen()
        return Vec((self.i/l, self.j/l, self.k/l))

class Grid:
    """A 2D array of Loc objects, representing a grid."""
    def __init__(self, side, csys=None, grid=[]):
        self._as_parameter_ = grid
        self.grid = grid
        self.csys = csys
        self.side = side
    
    def __setitem__(self,i,item):
        self.grid[i] = item
    
    def __getitem__(self,i):
        return self.grid[i]
    
    def __len__(self):
        return len(self.grid)

    def __repr__(self):
        return str(self.grid)

    def copy(self):
        """Returns a copy of the Pointlist."""
        return Grid(self.side, self.csys,self.grid[:])
    
    def transform(self,targ):
        """
        Transforms a Grid to another coordinate system, can read Proj4 format
        and WKT.
        """
        # Variables for the transformation.
        grid = self.copy()
        source = osr.SpatialReference()
        target = osr.SpatialReference()

        # Deciding whether the coordinate systems are Proj4 or WKT.
        sc0 = grid.csys[0]
        if sc0 == 'G' or sc0 == 'P':
            source.ImportFromWkt(grid.csys)
        else:
            source.ImportFromProj4(grid.csys)

        tc0 = targ[0]
        if tc0 == 'G' or tc0 == 'P':
            target.ImportFromWkt(targ)
        elif tc0 == '+':
            target.ImportFromProj4(targ)
        else:
            raise SystemExit("\033[91m\033[1mUnrecognized target coordinate system:\n{}".format(targ))

        # The actual transformation.
        transform = osr.CoordinateTransformation(source, target)
        xform = transform.TransformPoint
        for i in range(len(grid)):
            for j in range(len(grid[0])):
                npt = list(xform(grid[i,j].x,grid[i,j].y,grid[i,j].z))
                ol = grid[i,j].ol
                nd = grid[i,j].nd
                grid[i,j] = Loc(npt[0],npt[1],npt[2])
                grid[i,j].ol = ol
                grid[i,j].nd = nd
            

        grid.csys = targ
        return grid
    
    def toground(self, dem, outsys = None):
        """
        Function will get the points on the ground directly below a list of points,
        this is not destructive and returns a new list 
        """
        grd = self.copy() # Copy to store on-ground points.
        origsys = grd.csys

        # Transforming to the DEM coordinate system so the geotransform math works.
        grd = grd.transform(dem.csys)

        # Iterate through the points and get the points below them.
        for i in range(len(grd)):
            for j in range(len(grd[0])):
                zpix = grd[i,j].topix(dem)
                if(zpix.x == -1 and zpix.y == -1):
                    grd[i,j].z = dem.nd
                    grd[i,j].nd = True
                else:
                    v = float(dem.data[zpix.y][zpix.x])
                    grd[i,j].z = v
                    grd[i,j].ol = zpix.ol
                    if(v == dem.nd):
                        grd[i,j].nd = True
                
        # Set coordinate systems for the new lists.
        if outsys:
            grd = grd.transform(outsys)
        else:
            grd = grd.transform(origsys)

        return grd
    
    def facet(self):
        """Generate a list of facets from a Grid."""
        flist = [None]*((len(self.grid)-1)*(len(self.grid[0])-1)*2)  # Init empty facet list.
        cpath = Path(self.csys, [None]*((len(self.grid)-1)*(len(self.grid[0])-1)*2))  # Build empty Path object.
        c = 0  # Iterator for building the facet list.
        ncol = len(self.grid[0])

        # Build list of facets.
        for i in range(len(self.grid)-1):
            for j in range(ncol-1):
                flist[c] = Facet(self.grid[i,j+1], self.grid[i,j], self.grid[i+1,j], j, self.side)
                if(flist[c].l0.nd or flist[c].l1.nd or flist[c].l2.nd):
                    flist[c].nd = True
                if(flist[c].l0.ol or flist[c].l1.ol or flist[c].l2.ol):
                    flist[c].ol = True

                flist[c+1] = Facet(self.grid[i+1,j], self.grid[i+1,j+1], self.grid[i,j+1], j, self.side)
                if(flist[c+1].l0.nd or flist[c+1].l1.nd or flist[c+1].l2.nd):
                    flist[c+1].nd = True
                if(flist[c+1].l0.ol or flist[c+1].l1.ol or flist[c+1].l2.ol):
                    flist[c+1].ol = True

                cpath[c] = flist[c].cen
                cpath[c+1] = flist[c+1].cen
                c += 2
            
        return flist, cpath

class Path:
    """
    A list of loc objects, representing a path. Along with some useful
    metadata.
    """
    def __init__(self, csys = None, pts = []):
        self.pts = pts
        self.csys = csys
        
    def __setitem__(self,i,item):
        self.pts[i] = item

    def __getitem__(self,i):
        return self.pts[i]

    def __len__(self):
        return len(self.pts)

    def __str__(self):
        return str(self.pts)

    def append(self,Loc):
        self.pts.append(Loc)

    def copy(self):
        """Returns a copy of the Pointlist."""
        return Path(self.csys,self.pts[:])

    def avgdiff(self):
        """
        For points in XYZ, returns the average difference between the points in
        the path.
        """
        s = 0
        c = 0
        for i in range(len(self.pts)-1):
            if(self.pts[i] != Loc(0,0,0) and self.pts[i+1] != Loc(0,0,0)):
                xd = (self.pts[i+1].x - self.pts[i].x)**2
                yd = (self.pts[i+1].y - self.pts[i].y)**2
                zd = (self.pts[i+1].z - self.pts[i].z)**2
                s = s + (xd + yd + zd)**.5
                c = c+1
            
        avg = s/c
        return avg

    def save(self, path, bins = None):
        """
        Save CSV products.
        """
        with open(path,'w') as f:
            # Write points and bins into csv file path.
            for i in range(len(self.pts)):
                # Just write points.
                if not bins:
                    f.write("{},{},{}\n".format(self.pts[i].x,self.pts[i].y,self.pts[i].z))
                # Write points and bins.
                elif len(bins) == len(self.pts):
                    f.write("{},{},{},{}\n".format(self.pts[i].x,self.pts[i].y,self.pts[i].z,bins[i]))
                else:
                    raise SystemExit("\033[91m\033[1mERROR: Invalid bins")

    def transform(self,targ):
        """
        Transforms a Pointlist to another coordinate system, can read Proj4 format
        and WKT.
        """
        pts = self.copy()

        source = osr.SpatialReference()
        target = osr.SpatialReference()

        # Deciding whether the coordinate systems are Proj4 or WKT.
        sc0 = pts.csys[0]
        if sc0 == 'G' or sc0 == 'P':
            source.ImportFromWkt(pts.csys)
        else:
            source.ImportFromProj4(pts.csys)

        tc0 = targ[0]
        if tc0 == 'G' or tc0 == 'P':
            target.ImportFromWkt(targ)
        elif tc0 == '+':
            target.ImportFromProj4(targ)
        else:
            raise SystemExit("\033[91m\033[1mUnrecognized target coordinate system:\n{}".format(targ))

        # The actual transformation.
        transform = osr.CoordinateTransformation(source, target)
        xform = transform.TransformPoint
        for i in range(len(pts)):
            npt = list(xform(pts[i].x,pts[i].y,pts[i].z))
            pts[i] = Loc(npt[0],npt[1],npt[2])

        pts.csys = targ
        return pts

    def toground(self, dem, outsys = None):
        """
        Function will get the points on the ground directly below a list of
        points, this is not destructive and returns a new list.
        """
        grd = self.copy() # Copy to store on-ground points.
        origsys = grd.csys

        # Transforming to the DEM coordinate system so the geotransform math works.
        grd = grd.transform(dem.csys)

        # Iterate through the points and get the points below them.
        for i in range(len(grd)):
            zpix = grd[i].topix(dem)
            if(zpix.x == -1 and zpix.y == -1):
                grd[i].z = dem.nd
            else:
                grd[i].z = float(dem.data[zpix.y][zpix.x])

        # Set coordinate systems for the new lists.
        if outsys:
            grd = grd.transform(outsys)
        else:
            grd = grd.transform(origsys)

        return grd

    def getdir(self,i):
        """
        Calculates a unit vector between the point ahead of idx in the list and
        the one behind it has cases to handle the first and last points in the
        list.
        """
        if i == 0:
            d = Vec(self.pts[i], self.pts[i+1])
        if i == len(self)-1 :
            d = Vec(self.pts[i-1], self.pts[i])
        else:
            d = Vec(self.pts[i-1], self.pts[i+1])

        return d.unit()

class Facet:
    """Holds a triangular facet made of thee loc objects."""
    def __init__(self, l0, l1, l2, idx, side):
        self.l0 = l0
        self.l1 = l1
        self.l2 = l2
        self.side = side
        self.area, self.norm = self.farea()
        self.cen = self.center()
        self.idx = idx
        self.ol = False
        self.nd = False

    def __repr__(self):
        return str([self.l0, self.l1, self.l2])

    def __getitem__(self, i):
        return self.facetList()[i]

    def __len__(self):
        return len(self.facetList())

    def facetList(self):
        return [self.l0, self.l1, self.l2]

    def farea(self):
        if(self.side == 'l'):
            v1 = Vec(self.l1,self.l0)
            v2 = Vec(self.l1,self.l2)
            xp = v1.xprd(v2)
        elif(self.side == 'r'):
            v1 = Vec(self.l1,self.l2)
            v2 = Vec(self.l1,self.l0) 
            xp = v1.xprd(v2)
        else:
            raise SystemExit("\033[91m\033[1mInvalid facet side")

        return .5*xp.vlen(), xp

    def center(self):
        x = (self.l0.x + self.l1.x + self.l2.x)/3
        y = (self.l0.y + self.l1.y + self.l2.y)/3
        z = (self.l0.z + self.l1.z + self.l2.z)/3
        return Loc(x,y,z)

class Dem:
    """Holds an array of DEM data and relevant metadata."""
    def __init__(self,dem_path):
        src_dem = gdal.Open(dem_path,GA_ReadOnly)
        self.csys = src_dem.GetProjection()
        self.gt = src_dem.GetGeoTransform()
        self.nd = src_dem.GetRasterBand(1).GetNoDataValue()
        if not self.nd:
            self.nd = -np.finfo('d').max
        self.data = np.array(src_dem.GetRasterBand(1).ReadAsArray())
        src_dem = None

    def __str__(self):
        return str(self.csys)