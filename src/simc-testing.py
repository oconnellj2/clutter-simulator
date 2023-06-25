import numpy as np
import scipy.misc
from osgeo import gdal, osr
from gdalconst import *
import sys
from laspy.file import File
import xml.etree.ElementTree as et

# Data Location
data_path = sys.argv[1]

# XML Info Path
xml_path = sys.argv[2]

# Parse XML and get info
print("Parsing XML info...")
xml = et.parse(xml_path)
root = xml.getroot()
# Coordinate reference system
srs = root.find("header").find("srs").find("wkt").text
# Offset, min, max
xo = float(root.find("header").find("offset").find("x").text)
xmin = float(root.find("header").find("minimum").find("x").text)
xmax = float(root.find("header").find("maximum").find("x").text)
yo = float(root.find("header").find("offset").find("y").text)
ymin = float(root.find("header").find("minimum").find("y").text)
ymax = float(root.find("header").find("maximum").find("y").text)
zo = float(root.find("header").find("offset").find("z").text)
zmin = float(root.find("header").find("minimum").find("z").text)
zmax = float(root.find("header").find("maximum").find("z").text)
# Scale
xs = float(root.find("header").find("scale").find("x").text)
ys = float(root.find("header").find("scale").find("y").text)
zs = float(root.find("header").find("scale").find("z").text)

# Grid Params
dx = 10
dy = 10

# No data value in grid
ndata = -9999

# Import point data
print("Reading point data...")
data_file = File(data_path, mode='r')
x = data_file.X*xs + xo
y = data_file.Y*ys + yo
z = data_file.Z*zs + zo
# np.savetxt('pts.csv',np.vstack([x,y,z]).transpose()[0::10000],fmt="%f",delimiter=",")
n = len(x)

## Do the Gridding ##
print("Mapping points to grid space...")
# Transform x and y from xy space to grid coordinate space
nx = np.floor_divide(x, dx).astype('int32')
ny = np.floor_divide(y, dy).astype('int32')
xlim = nx.min()*dx
ylim = ny.max()*dy
nx = np.subtract(nx, nx.min())
ny = np.subtract(ny, ny.min())
ny = np.subtract(ny.max(), ny)
nz = z.astype('float32')

# Make grid with no data values and a counting grid
gdata = (np.ones((ny.max()+1, nx.max()+1))*ndata).astype('float32')
count = np.zeros((ny.max()+1, nx.max()+1))
# print(nx)
# print(ny)

pct = 0
# Put points into the grid
print("Gridding point data...")
for i in range(n):
	if (int(float(i*100)/n) > pct):
		sys.stdout.write("\r")
		sys.stdout.write(str(int(float(i*100)/n))+"%")
		sys.stdout.flush()
		pct = int(float(i*100)/n)
	if (count[ny[i]][nx[i]] == 0):
		gdata[ny[i]][nx[i]] = nz[i]
		count[ny[i]][nx[i]] = 1
	else:
		count[ny[i]][nx[i]] = count[ny[i]][nx[i]] + 1
		gdata[ny[i]][nx[i]] = gdata[ny[i]][nx[i]] * \
			(1-(1/count[ny[i]][nx[i]])) + nz[i]*(1/count[ny[i]][nx[i]])

sys.stdout.write("\r")
sys.stdout.write("100%\n")

# Make output name
oname = data_path.replace(".las", "." + str(dx) + "m.tif")

# Make and save geotiff w/ fake UTM
print("Saving GeoTIFF...")
drv = gdal.GetDriverByName('GTiff')
ds = drv.Create(oname, int(nx.max()+1), int(ny.max()+1), int(1), GDT_Float32)
ds.SetGeoTransform([xlim, dx, 0, ylim, 0, -dy])
ds.SetProjection(srs)
ds.GetRasterBand(1).WriteArray(gdata)
ds.GetRasterBand(1).SetNoDataValue(ndata)
ds.GetRasterBand(1).SetUnitType("m")
ds = None
