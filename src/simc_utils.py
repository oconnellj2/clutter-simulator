"""
	simc_utils.py
	=============
	Purpose:
	--------
		Organize functions in support of the creation of surface clutter
		simulations to aid in the analysis of airborne and spaceborn sounding
		radar images in simc.py.
	Functions:
	----------
		readConfig:
			configChecks(): Checks values of config file.
		simPrep:
			systems(): Defines XYZ and LLE coordinate systems.
			miscSimPrep(): Opens DEM, reads nav data constructs the shift and nadbin arrays.
			echomapPrep(): Defines empty grid and geotransform.
			dataArrays(): Constructs arrays for sim.
		runSim:
			makeGrid(): Makes grid from transceiver location.
			grids(): Makes left & right grid objects projected to the ground.
			c_sim(): Calculates two way travel time & power in C.
			py_sim(): Python method to calculate two way travel time & power.
			fcenAvg(): Calculates average facet center location along track per trace.
			averages(): Calculates twtt avg and facet center avg per trace.
			checkSides(): Checks left & right side of range bin.
			cluttergram(): Build cluttergrams and echomap_ref.
			printInfo(): Outputs relevant info for each trace.
		saveProducts:
			combined(): Produce combined/adjusted products.
			binary(): Produce binary products.
			sides(): Produce left and right side cluttergram products.
			echomapProd(): Produce echomap / adjusted cluttergram products.
			csvProd(): Produce csv products.
		nav functions:
			getNav_ak(): Get navdat for alaska.
			getNav_fritz(): Get navdat for fritz.
			getNav_geom(): Get navdat for geom.
	Contributing:
	-------------
		Univerisity of Arizona | Lunar & Planetary Laboratory
		TAPIR | Terrestrial And Planetary Investigation with Radar
		James O'Connell | Research Assistant
"""
from ctypes import PyDLL, c_double, c_int
from sys import argv, platform

import numpy as np
from backports import configparser
from gdalconst import GA_ReadOnly
from numpy.ctypeslib import ndpointer
from osgeo import gdal, osr
from PIL import Image, ImageOps
from skimage import transform

from simc_class import *


def configChecks(confDict):
	"""
	Checks values from config file.

	Args:
		confDict(2D dict): Config parameters.
	Raises:
		SystemExit: If config values are invalid for sim to run.
	"""
	at_dist = int(confDict['facet_params']['at_dist'])
	at_step = int(confDict['facet_params']['at_step'])
	ct_dist = int(confDict['facet_params']['ct_dist'])
	ct_step = int(confDict['facet_params']['ct_step'])

	# Checks if commandline args are neccessary.
	if 'nav_path' not in confDict['paths']:
		confDict['paths']['nav_path'] = argv[2]  # Nav path is second arg.
		confDict['paths']['dem_path'] = argv[3]  # DEM path is third arg.

	if at_dist < at_step:
		raise SystemExit(
			"\033[91m\033[1mINVALID PARAM!\nALONG TRACK DIST < ALONG TRACK STEP")
	if ct_dist < ct_step:
		raise SystemExit(
			"\033[91m\033[1mINVALID PARAM!\nCROSS TRACK DIST < CROSS TRACK STEP")
	if at_dist % at_step != 0:
		raise SystemExit(
			"\033[91m\033[1mINVALID PARAM!\nALONG TRACK DIST IS NOT EQUALLY DIVISIBLE BY ALONG TRACK STEP")
	if ct_dist % ct_step != 0:
		raise SystemExit(
			"\033[91m\033[1mINVALID PARAM!\nCROSS TRACK DIST IS NOT EQUALLY DIVISIBLE BY CROSS TRACK STEP")

	# Make output path for unix.
	unix = set(["linux", "Linux", "darwin"])
	if platform in unix:
		confDict['paths']['out_path'] = '../out/'
		navfile = confDict['paths']['nav_path'].split('/')[-1]
	# Make output path for Windows.
	else:
		confDict['paths']['out_path'] = confDict['paths']['out_path'] + '\\'
		navfile = confDict['paths']['nav_path'].split('\\')[-1]

	navname = navfile.split('.')[0]
	confDict['paths']['out_path'] = confDict['paths']['out_path'] + navname + '_'


def systems(confDict, simParams):
	"""
	Defines `XYZ and LLE coordinate systems` for internal calculations into
	simParams.

	Prints information in preperation for the sim to run.

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters to be passed into simRun.
	Notes:
		If you are doing simulations on a body that is not listed here you
		will need to add your own coordinate systems, make sure to add both
		Long/Lat and XYZ.
	"""
	# Long/Lat, xyz, echomap sys.
	if confDict['sim_params']['body'] == 'mars':
		# IAU Mars 2000.
		llesys = '+proj=longlat +a=3396190 +b=3376200 +no_defs'
		# IAU Mars 2000.
		xyzsys = '+proj=geocent +a=3396190 +b=3376200 +no_defs'
		ecosys = '+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=180 +x_0=0 +y_0=0 +a=3396000 +b=3396000 +units=m +no_defs'
	if confDict['sim_params']['body'] == 'moon':
		# IAU Moon 2000.
		llesys = '+proj=longlat +a=1737400 +b=1737400 +no_defs'
		# IAU Moon 2000.
		xyzsys = '+proj=geocent +a=1737400 +b=1737400 +no_defs'
	if confDict['sim_params']['body'] == 'earth':
		# IAU Earth 2000.
		llesys = '+proj=longlat +a=6378140 +b=6356750 +no_defs'
		# IAU Earth 2000.
		xyzsys = '+proj=geocent +a=6378140 +b=6356750 +no_defs'

	print("\n\033[1mPlanet:\033[0m " + confDict['sim_params']['body'])
	print("\033[1mNav:\033[0m " + confDict['paths']['nav_path'].split('/')[-1])
	print("\033[1mDEM:\033[0m " + confDict['paths']['dem_path'].split('/')[-1])
	print("\033[1mGrid:\033[0m " + confDict['facet_params']['at_dist'] +
		  "m along track, " + confDict['facet_params']['ct_dist'] + "m cross track")
	print("\033[1mFacet:\033[0m " + confDict['facet_params']['at_step'] +
		  "m along track, " + confDict['facet_params']['ct_step'] + "m cross track")
	print("\033[1mCoordinate system:\033[0m " + llesys)

	simParams['sys'] = llesys, xyzsys, ecosys


def miscSimPrep(confDict, simParams):
	"""
	Opens DEM, reads in navigation data with specified function, defines
	number of cross track/along track points to make facets and constructs the
	shift and nadbin arrays. 

	Objects are assigned to the simParams dictionary to be used in runSim
	function. 

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters to be passed into simRun.
	"""
	# Open DEM and make DEM object.
	simParams['topo'] = Dem(confDict['paths']['dem_path'])

	# Read in navigation data.
	simParams['navdat'] = eval(confDict['navigation']['navfunc'])(
		confDict['paths']['nav_path'])
	simParams['navdat'].csys = confDict['navigation']['navsys']

	# Define number of along track and cross track points to make facets from.
	simParams['num_ct'] = int(
		confDict['facet_params']['ct_dist'])//int(confDict['facet_params']['ct_step'])
	# Not used.
	simParams['num_at'] = int(
		confDict['facet_params']['at_dist'])//int(confDict['facet_params']['at_step'])

	# This is specific to making FPB sims with areoid referenced products.
	simParams['nad_loc'] = simParams['navdat'].toground(
		simParams['topo'], confDict['navigation']['navsys'])
	# Transform xyzsys.
	simParams['nad_xyz'] = simParams['nad_loc'].transform(simParams['sys'][1])
	simParams['shift'] = [0]*len(simParams['navdat'])
	simParams['nadbin'] = [0]*len(simParams['navdat'])

	c = float(confDict['sim_params']['speedlight'])
	binsize = float(confDict['sim_params']['binsize'])
	datum_sample = int(confDict['datum_params']['datum_sample'])
	for i in range(len(simParams['navdat'])):
		simParams['shift'][i] = (
			simParams['navdat'][i].z)*2/c/binsize - datum_sample
		# Calculate bin of nadir return, if requested.
		if confDict['display_params']['show_nadir'] == 'True':
			simParams['nadbin'][i] = int(
				(simParams['navdat'][i].z-simParams['nad_loc'][i].z)*2/c/binsize - simParams['shift'][i])


def echomapPrep(confDict, simParams):
	"""
	Defines empty grid and geotransform for `echomap_ref`. Creates empty
	arrays for the `echomap_ref` and the `return times`.

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters to be passed into simRun.
	"""
	if confDict['outputs']['echomap_ref'] == 'True':
		# Transform ecosys.
		ecom = simParams['navdat'].transform(simParams['sys'][2])
		# Get grid limits and step size for output echomap_ref.
		maxdist = max(int(confDict['facet_params']['at_dist']), int(
			confDict['facet_params']['ct_dist']))
		exmin = ecom[0].x-maxdist
		exmax = ecom[0].x+maxdist
		eymin = ecom[0].y-maxdist
		eymax = ecom[0].y+maxdist
		for i in range(1, len(ecom)):
			exmin = min(exmin, ecom[i].x-maxdist)
			exmax = max(exmax, ecom[i].x+maxdist)
			eymin = min(eymin, ecom[i].y-maxdist)
			eymax = max(eymax, ecom[i].y+maxdist)

		epsize = max(int(confDict['facet_params']['at_step']), int(
			confDict['facet_params']['ct_step'])) * 12

		# Make empty grid and geotransform for echomap_ref.
		simParams['egt'] = [exmin, epsize, 0, eymax, 0, -epsize]
		simParams['exs'] = int(math.ceil((exmax-exmin)/float(epsize)))
		simParams['eys'] = int(math.ceil((eymax-eymin)/float(epsize)))
		simParams['emap_ref'] = np.zeros((simParams['eys'], simParams['exs']))


def dataArrays(confDict, simParams):
	"""
	Define various arrays to prepare for clutter simulation and `transform nav
	file to XYZ.`

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters to be passed into simRun.
	Returns:
		data(dict): Data structures to hold clutter simulations 
	"""
	data = {}  # Empty Dictionary to hold data for runSim.

	# Empty arrays for the echomap_ref and the return times.
	if confDict['outputs']['echomap'] == 'True':
		data['emap'] = np.zeros(
			(2*simParams['num_ct'], len(simParams['navdat'])))
		data['bins'] = np.ones(
			(2*simParams['num_ct'], len(simParams['navdat'])))*int(confDict['sim_params']['tracelen'])

	# Make empty arrays to hold the right and left side clutter simulations.
	data['rcgram'] = np.zeros(
		[int(confDict['sim_params']['tracelen']), len(simParams['navdat'])])
	data['lcgram'] = np.zeros(
		[int(confDict['sim_params']['tracelen']), len(simParams['navdat'])])
	data['rRed'] = np.zeros(
		[int(confDict['sim_params']['tracelen']), len(simParams['navdat'])])
	data['lRed'] = np.zeros(
		[int(confDict['sim_params']['tracelen']), len(simParams['navdat'])])
	data['tmap'] = [None]*len(simParams['navdat'])

	# Keeping track of the bins data is assigned to.
	# xyzsys.
	data['fret_loc'] = Path(simParams['sys'][1], [
							Loc(0, 0, 0)]*len(simParams['navdat']))
	data['minbin'] = [None]*len(simParams['navdat'])

	# Transform navigation file to XYZ.
	simParams['grdlle'] = simParams['navdat'].toground(
		simParams['topo'], simParams['sys'][0])  # Not used.
	# Transform xyzsys.
	simParams['navxyz'] = simParams['navdat'].transform(simParams['sys'][1])
	simParams['grdxyz'] = simParams['navdat'].toground(
		simParams['topo'], simParams['sys'][1])  # xyzsys.

	data['badf'] = 0

	return data


def simChecks(simParams, data, i):
	"""
	Checks over DEM data, bins, and locations and corrects them.

	Args:
		simParams(dict): Parameters for runSim.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
		i(int): Iterator from runSim.
	"""
	# Check if spacecraft over DEM.
	if simParams['nad_loc'][i].z == simParams['topo'].nd:
		simParams['nad_loc'][i] = Loc(0, 0, 0)
		simParams['nad_xyz'][i] = Loc(0, 0, 0)
		simParams['nadbin'][i] = 0
		data['fret_loc'][i] = Loc(0, 0, 0)
		data['minbin'][i] = 0

	# Check if nav points are identical.
	if simParams['navxyz'][i] == simParams['navxyz'][i-1] and i != 0:
		data['lcgram'][:, i] = data['lcgram'][:, i-1]
		data['rcgram'][:, i] = data['rcgram'][:, i-1]
		data['emap'][:, i] = data['emap'][:, i-1]


def makeGrid(xcvr, unitVec, cvec, at_step, at_dist, ct_step, ct_dist):
	"""
	Make a grid from a transceiver location and a vector orthogonal to
	direction of travel and parallel to the surface.

	Args:
		xcvr(Loc object): Spacecrafts coordinate vector.
		unitVec(Loc object): Approximate direction of travel, unit vector.
		cvec(Vec object): Points to the Right/Left of the aircraft.
		at_step(int): Along track facet dimension.
		at_dist(int): Distance to the front and rear for each nav point.
		ct_step(int): Cross track facet dimension.
		ct_dist(int): Distance to either side.
	Returns:
		Grid(3D numpy array): Grid from a transceiver location and a vector
							orthogonal to direction of travel.
	Notes:
		The xcvr location needs to be in xyz avec/cvec need to be unit vectors.
	"""
	numRows = int(2*(at_dist/at_step)+1)
	nCols = int((ct_dist/ct_step)+1)

	grid = np.empty((numRows, nCols), dtype=object)

	newPostion = xcvr + (at_dist*unitVec)
	# Fill grid.
	for i in range(len(grid)):
		for j in range(len(grid[i])):
			grid[i, j] = newPostion - unitVec*i*at_step + cvec*j*ct_step

	return grid


def grids(confDict, data, simParams, i):
	"""
	Constructs left & right grid objects projected to the ground.

	Args:
		confDict(2D dict): Config parameters.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
		simParams(dict): Parameters for runSim.
		i(int): Iterator from runSim.
	"""
	unitVec = simParams['navxyz'].getdir(
		i)  # Approximate direction of travel, unit vector.
	ground_point = simParams['grdxyz'][i]  # Point on ground at nadir.
	# Vector from nadir to xcvr.
	simParams['grd_vec'] = Vec(ground_point, simParams['navxyz'][i])

	# (Travel Vector X Ground Vector), this will point to the left and right sides of the aircraft.
	cross = unitVec.xprd(simParams['grd_vec'].unit()).unit()
	rvec = cross
	lvec = -cross

	# Make right and left side grids.
	at_step = int(confDict['facet_params']['at_step'])
	ct_step = int(confDict['facet_params']['ct_step'])
	ct_dist = int(confDict['facet_params']['ct_dist'])
	at_dist = int(confDict['facet_params']['at_dist'])
	rgrid = makeGrid(simParams['navxyz'][i], unitVec,
					 rvec, at_step, at_dist, ct_step, ct_dist)
	lgrid = makeGrid(simParams['navxyz'][i], unitVec,
					 lvec, at_step, at_dist, ct_step, ct_dist)

	# Make Grid objects and project them to the ground.
	rgrid = Grid('r', simParams['sys'][1], rgrid)  # xyzsys.
	lgrid = Grid('l', simParams['sys'][1], lgrid)  # xyzsys.
	rgrid = rgrid.toground(simParams['topo'], simParams['sys'][1])  # xyzsys.
	lgrid = lgrid.toground(simParams['topo'], simParams['sys'][1])  # xyzsys.

	# Generate facets from the grids.
	data['rFacet'], data['rcen'] = rgrid.facet()
	data['lFacet'], data['lcen'] = lgrid.facet()


def c_sim(confDict, data, xcvr, side):
	"""
	Calculate two way travel time and power via `Ctypes`.

	Args:
		confDict(2D dict): Config parameters.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
		xcvr(Loc object): Spacecrafts coordinate vector.
		side(str): Indicates which side we're calculating.
	"""
	# Unpack params.
	c = float(confDict['sim_params']['speedlight'])
	flist = data[side + 'Facet']

	# Allows accessing Python API functions.
	slib = PyDLL('../code/simlib.so')

	# Converts parameters to ctypes.
	c_twtt = (c_double * len(flist))(0)
	c_power = (c_double * len(flist))(0)
	c_xcvr = (c_double * 3)(*xcvr.locList())
	c_flist = (c_double * 3 * 3 * len(flist))(*(tuple(tuple(j)
													  for j in i) for i in flist))
	c_flen = c_int(len(flist))
	c_c = c_double(c)

	# Converts return type to array of floats.
	slib.c_pwr.restype = ndpointer(dtype=c_double, shape=(len(flist),))
	slib.c_twtt.restype = ndpointer(dtype=c_double, shape=(len(flist),))

	# Calls C function and returns twtt and power arrays.
	data[side + 'Twtt'] = list(slib.c_twtt(c_flist,
							   c_xcvr, c_twtt, c_c, c_flen))
	data[side + 'Power'] = list(slib.c_pwr(c_flist, c_xcvr, c_power, c_flen))
	data['twtt_min'] = max(data['rTwtt'])


def py_sim(confDict, data, xcvr, side):
	"""
	`NOT IN USE`\n
	Python method to calculate two way travel time and power.

	Args:
		confDict(2D dict): Config parameters.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
		xcvr(Loc object): Spacecrafts coordinate vector.
		side(str): Indicates which side we're calculating.
	"""
	# Unpack params.
	c = float(confDict['sim_params']['speedlight'])
	flist = data[side + 'Facet']

	# Empty arrays to hold the twtt and power output per trace.
	data[side + 'Twtt'] = [None]*len(flist)
	data[side + 'Power'] = [None]*len(flist)

	for i in range(len(flist)):
		vsc = Vec(flist[i].cen, xcvr)
		vsclen = vsc.vlen()  # Vector magnitude.
		data[side + 'Twtt'][i] = (vsclen)*2/c  # Calculate twtt.
		if flist[i].nd:
			data[side + 'Power'][i] = 0
		else:
			ct = vsc.dprd(flist[i].norm)/(vsclen*flist[i].norm.vlen())
			if ct < 0:
				data[side + 'Power'][i] = 0
			else:
				data[side + 'Power'][i] = abs((flist[i].area*ct)**2/vsclen**4)


def fcenAvg(flist):
	"""
	Calculates average facet center location along track per trace.

	Args:
		flist(list): list of facets along track wide.
	Returns:
		x_avg(float): Average center x value in 'flist'.
		y_avg(float): Average center y value in 'flist'.
		z_avg(float): Average center z value in 'flist'.
	"""
	x = y = z = 0
	# Sum x and y values in the facet list.
	for facet in flist:
		facet = facet.cen
		x += facet[0]
		y += facet[1]
		z += facet[2]
	x_avg = x/len(flist)  # x avg
	y_avg = y/len(flist)  # y avg
	z_avg = z/len(flist)  # z avg

	return x_avg, y_avg, z_avg


def averages(confDict, simParams, data, i):
	"""
	Calculates twtt avg and facet center avg per trace and assigns it to `tmap`.

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters of runSim.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
		i(int): Iterator from runSim.
	Raises:
		SystemExit: If left and right facet counts do not equal eachother.
	"""
	# Checks that facet count is equal on both sides of the aircraft.
	if len(data['rTwtt']) != len(data['lTwtt']):
		raise SystemExit("\033[91mFacet count mismatch, exiting")

	if confDict['outputs']['tmap'] == 'True':
		# Holds meta data per trace.
		rdata = [[None]*7 for i in range(2*simParams['num_ct'])]
		# Flip rtwtt and rFacet to produce data from bottom to top.
		twtt = data['rTwtt'][::-1] + data['lTwtt']
		flist = data['rFacet'][::-1] + data['lFacet']

		k = 0
		step = simParams['num_at']*4
		# Assign rdata values for output products.
		for j in range(0, len(twtt), step):
			rdata[k][0] = i  # Trace number.
			rdata[k][1] = simParams['navdat'][i][0]  # Lat.
			rdata[k][2] = simParams['navdat'][i][1]  # Lon.
			rdata[k][3] = sum(twtt[j:j+step]) / step  # Twtt avg along track.
			# Facet center average.
			rdata[k][4], rdata[k][5], rdata[k][6] = fcenAvg(flist[j:j+step])
			k += 1
		data['tmap'][i] = rdata


def checkSides(confDict, simParams, data, box, index, side):
	"""
	Checks left or right side of the range bin(boxl, boxr).

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters used in runSim.
		data(dict): Misc arrays that hold left & right side clutter simulations.
		box(int): Left or right box value.
		index(array): [i, j] iterators.
		side(str): Indicates left or right side data params.
	"""
	i, j = index  # Unpack iterators.
	traceLength = int(confDict['sim_params']['tracelen'])

	if box > -traceLength and box < traceLength:
		data[side + 'cgram'][box][i] = data[side +
											'cgram'][box][i] + data[side + 'Power'][j]
		if data[side + 'Facet'][j].ol:
			data[side + 'Red'][box][i] = 1

		if confDict['outputs']['echomap'] == 'True' or confDict['outputs']['echomap_adjusted'] == 'True':
			# Construct right side emap and data bins.
			if side == 'r':
				data['emap'][simParams['num_ct']-1 + data[side + 'Facet'][j].idx,
							 i] = data['emap'][simParams['num_ct']-1 + data[side + 'Facet'][j].idx, i] + data[side + 'Power'][j]
				data['bins'][simParams['num_ct']-1 + data[side + 'Facet'][j].idx, i] = min(
					data['bins'][simParams['num_ct']-1 + data[side + 'Facet'][j].idx, i], box)
			# Construct left side emap and data bins.
			if side == 'l':
				data['emap'][simParams['num_ct']-1 - data[side + 'Facet'][j].idx,
							 i] = data['emap'][simParams['num_ct']-1 - data[side + 'Facet'][j].idx, i] + data[side + 'Power'][j]
				data['bins'][simParams['num_ct']-1 - data[side + 'Facet'][j].idx, i] = min(
					data['bins'][simParams['num_ct']-1 - data[side + 'Facet'][j].idx, i], box)

		if confDict['outputs']['echomap_ref'] == 'True':
			[xAxis, yAxis] = rcen[j].topix_write(simParams['egt'])
			simParams['emap_ref'][yAxis, xAxis] = simParams['emap_ref'][yAxis,
																		xAxis] + np.log(data[side + 'Power'][j])
	else:
		data['badf'] += 1


def cluttergram(confDict, simParams, data, i):
	"""
	Build cluttergrams and `echomap_ref`. Called from runSim.

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters used in runSim.
		data(dict): Misc arrays to hold left & right side clutter simulations.
		i(int): Iterator from runSim.
	"""
	if confDict['outputs']['echomap_ref'] == 'True':
		# Necessary for gridding.
		rcen = rcen.transform(simParams['sys'][2])  # Transform ecosys.
		lcen = lcen.transform(simParams['sys'][2])  # Transform ecosys.

	data['twtt_min'] = max(data['rTwtt'])
	for j in range(len(data['rTwtt'])):
		boxr = int(data['rTwtt'][j]/float(confDict['sim_params']
				   ['binsize']) - simParams['shift'][i])
		boxl = int(data['lTwtt'][j]/float(confDict['sim_params']
				   ['binsize']) - simParams['shift'][i])

		# Checks right two way travel time is less than the two way travel time minimum.
		# Assigns right data to final products.
		if data['rTwtt'][j] < data['twtt_min']:
			data['twtt_min'] = data['rTwtt'][j]
			data['fret_loc'][i] = data['rFacet'][j].cen
			data['minbin'][i] = boxr

		# Checks left two way travel time is less than the two way travel time minimum.
		# Assigns left data to final products.
		if data['lTwtt'][j] < data['twtt_min']:
			data['twtt_min'] = data['lTwtt'][j]
			data['fret_loc'][i] = data['lFacet'][j].cen
			data['minbin'][i] = boxl

		# Check right side of the range bin.
		checkSides(confDict, simParams, data, boxr, [i, j], 'r')

		# Check left side of the range bin.
		checkSides(confDict, simParams, data, boxl, [i, j], 'l')


def printInfo(simParams, i):
	"""
	Outputs relevant info to the terminal for each trace in the simulation.

	Args:
		simParams(dict): Parameters for the simulation.
		i(int): Iterator from runSim.
	"""
	trace = "\033[92m\033[1m Trace: {} \033[0m\033[1m|".format(i)
	sat = "\033[92m\033[1mSpacecraft elevation = {:.3f}m".format(
		simParams['grd_vec'].vlen())
	percent = 100*(i/len(simParams['navxyz']))
	prog = "\033[0m\033[1m|\033[94m Progress[{:.2f}%]".format(percent)
	print(trace, sat, prog, end='\r')


def combined(confDict, simParams, data, adj):
	"""
	Produce combined/adjusted products.

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters for the simulation.
		data(dict): Misc arrays that hold left/right side clutter sims and
					results.
		adj(str): If producing adjusted products.
	"""
	cscale = data['cgram']*(255.0/data['cgram'].max())
	cstack = np.dstack((cscale, cscale, cscale)).astype('uint8')

	if confDict['display_params']['show_fret'] == 'True':
		for i in range(len(data['minbin'])):
			cstack[data['minbin'][i], i] = data['fre']

	if confDict['display_params']['show_nadir'] == 'True':
		for i in range(len(simParams['nadbin'])):
			if simParams['nadbin'][i] < int(confDict['sim_params']['tracelen']):
				cstack[simParams['nadbin'][i], i] = data['nad']

	cimg = Image.fromarray(cstack)
	cimg = cimg.convert('RGB')
	cimg.save(confDict['paths']['out_path'] + 'combined' + adj + '.png')

	if confDict['outputs']['red'] == 'True':
		one = np.ones(data['cRed'].shape)
		redstack = np.dstack((one, data['cRed'], data['cRed'])).astype('uint8')
		cstack = cstack*redstack
		if confDict['display_params']['show_nadir'] == 'True':
			for i in range(len(simParams['nadbin'])):
				if simParams['nadbin'][i] < int(confDict['sim_params']['tracelen']):
					cstack[simParams['nadbin'][i], i] = data['nad']
		cimg = Image.fromarray(cstack)
		cimg = cimg.convert('RGB')
		cimg.save(confDict['paths']['out_path'] +
				  'combined' + adj + '_red.png')


def binary(confDict, data):
	"""
	Produce binary products.

	Args:
		confDict(2D dict): Config parameters.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
	"""
	cgram = data['rcgram'] + data['lcgram']
	cgram.astype('float32').tofile(
		confDict['paths']['out_path'] + 'combined.img')


def sides(confDict, simParams, data, cgram, side):
	"""
	Produce left and right side cluttergram products.

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters of runSim.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
		cgram(numpy array): Left or Right side cgram to be handled.
		side(str): To indicate which side cgram product is being produced.
	"""
	for i, row in enumerate(cgram):
		for j, pt in enumerate(row):
			cgram[i][j] = 10*np.log10(pt+1)

	cgram = cgram*(255/cgram.max())
	rstack = np.dstack((cgram, cgram, cgram)).astype('uint8')
	if confDict['display_params']['show_fret'] == 'True':
		for i in range(len(data['minbin'])):
			rstack[data['minbin'][i], i] = data['fre']

	if confDict['display_params']['show_nadir'] == 'True':
		for i in range(len(simParams['nadbin'])):
			if simParams['nadbin'][i] < int(confDict['sim_params']['tracelen']):
				rstack[simParams['nadbin'][i], i] = data['nad']

	sideImg = Image.fromarray(rstack)
	sideImg = sideImg.convert('RGB')
	sideImg.save(confDict['paths']['out_path'] + side + '.png')


def echomap(confDict, simParams, data, echoParam, adj):
	"""
	Produce echomap / adjusted cluttergram products.

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters of runSim.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
		echoParam(array): [escale/nemap, yf, yd]
		adj(bool): If producing adjusted products.
	"""
	scale, yDistortion, yDimension = echoParam  # Unpack params for operations.
	# If producing adjusted products.
	if adj:
		adj = '_adjusted'
	else:
		adj = ''

	estack = np.dstack((scale, scale, scale)).astype('uint8')
	for i in range(len(data['bins'])):
		for j in range(len(data['bins'][0])):
			if data['bins'][i, j] == data['minbin'][j]:
				estack[int(i*yDistortion), j] = data['fre']  # Assign echo loc.

	estack[yDimension-1][0:len(simParams['navdat'])-1] = data['nad']
	eimg = Image.fromarray(estack)
	eimg = eimg.convert('RGB')
	eimg.save(confDict['paths']['out_path'] + 'echomap' + adj + '.png')


def csv(confDict, simParams, data):
	"""
	Produce csv products.

	Args:
		confDict(2D dict): Config parameters.
		simParams(dict): Parameters of runSim.
		data(dict): Misc arrays that hold left/right side clutter sims and
				results.
	"""
	# Produce fret.csv product.
	if confDict['outputs']['fret'] == 'True':
		# Transform llesys.
		fret_loc = data['fret_loc'].transform(simParams['sys'][0])
		fret_loc.save(confDict['paths']['out_path'] +
					  'fret.csv', data['minbin'])

	# Produce nadir.csv product.
	if confDict['outputs']['nadir'] == 'True':
		# Transform llesys.
		nad_loc = simParams['nad_loc'].transform(simParams['sys'][0])
		nad_loc.save(confDict['paths']['out_path'] +
					 'nadir.csv', simParams['nadbin'])

	# Produce tmap.csv produt.
	if confDict['outputs']['tmap'] == 'True':
		path = confDict['paths']['out_path'] + 'tmap.csv'
		with open(path, "w") as f:
			for trace in data['tmap']:
				for row in trace:
					f.write("%s\n" % ','.join(str(col) for col in row))


def GetNav_ak(navfile):
	"""Get Nav for Alaska nav file."""
	f = open(navfile, 'r')
	navdat = f.read()
	f.close()
	navdat = navdat.split('\n')
	navdat_ls = Path()
	for i in navdat:
		if len(i) > 0:
			i = i.split(',')
			navdat_ls.append(Loc(float(i[1]), float(i[0]), float(i[2])))

	return navdat_ls


def GetNav_fritz(navfile):
	"""Get nav file for fritz."""
	f = open(navfile, 'r')
	navdat_raw = f.read()
	f.close()
	navdat_raw = navdat_raw.split('\n')
	navdat_raw = filter(None, navdat_raw)
	navdat = Path()
	for i in navdat_raw:
		if len(i) > 0:
			i = i.split(',')
			navdat.append(Loc(float(i[1]), float(i[0]), float(i[2])*1000))

	navdat.csys = '+proj=longlat +a=3396190 +b=3376200 +no_defs'
	# Adjust with areoid ... this relies on specific directory structure.
	aer = Dem('../test/dem/mega_16.tif')
	aer_nadir = navdat.toground(aer)

	for i in range(len(navdat)):
		if aer_nadir[i].z == aer.nd:
			aer_nadir[i].z = aer_nadir[i-1].z

		newr = navdat[i].z - 3396190 - aer_nadir[i].z
		if newr > 400000:
			raise SystemExit("\033[91m\033[1m" +
							 str(i, newr, aer_nadir[i].z, navdat[i].z))

		navdat[i].z = newr

	return navdat


def GetNav_geom(navfile):
	"""Get nav Geom."""
	f = open(navfile, 'r')
	navdat_raw = filter(None, f.read().split('\n'))
	f.close()

	navdat = Path()
	# Build Path object.
	for i in navdat_raw:
		if len(i) > 0:
			i = i.split(',')
			navdat.append(Loc(float(i[3]), float(i[2]), float(i[5])*1000))

	navdat.csys = '+proj=longlat +a=3396190 +b=3376200 +no_defs'
	# Adjust with areoid ... this relies on specific directory structure.
	aer = Dem('../test/dem/mega_16.tif')
	aer_nadir = navdat.toground(aer)

	for i in range(len(navdat)):
		if aer_nadir[i].z == aer.nd:
			aer_nadir[i].z = aer_nadir[i-1].z
		newr = navdat[i].z - 3396190 - aer_nadir[i].z
		navdat[i].z = newr

	return navdat
