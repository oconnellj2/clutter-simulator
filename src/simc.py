"""
    simc.py
    =======
    Purpose:
    --------
        Create surface clutter simulations to aid in the analysis of
        airborne and spaceborn sounding radar images.
    Contributing:
    -------------
        Univerisity of Arizona | Lunar & Planetary Laboratory
        TAPIR | Terrestrial And Planetary Investigation with Radar
        James O'Connell | Research Assistant
"""
from simc_utils import *


def readConfig():
    """
    Reads in config file into dictionary.

    Returns:
        confDict(2D dict): Config parameters.
    """
    config = configparser.RawConfigParser()
    config.read(argv[1])  # Config file path is first arg.
    confDict = {section: dict(config.items(section)) for section in config.sections()}  # Dict of config file.

    # Check users config file params.
    configChecks(confDict)
    
    return confDict

def simPrep(confDict):
    """
    Prepares clutter simulation to run with various helper functions that
    append to `simParams`.

    Args:
        confDict(2D dict): Config parameters.
    Returns:
        simParams(dict): Simulation parameters to be passed into the
                        clutter sim.
        data(dict): Arrays to prepare for clutter simulation and
                    transform nav file to XYZ
    """
    simParams = {}  # Empty dictionary to hold parameters for the sim.
    systems(confDict, simParams)  # Lat/Long, XYZ, echomap systems. [llsys, xyzsys, ecosys].
    miscSimPrep(confDict, simParams)  # DEM, navdat, num_ct, shift, and nadbin.
    echomapPrep(confDict, simParams)  # Echomap_ref
    data = dataArrays(confDict, simParams)  # Left & right side arrays.
    return simParams, data

def runSim(confDict, simParams, data):
    """
    Runs the actual clutter simulator.

    Args:
        confDict(2D dict): Config parameters.
        simParams(dict): Parameters for the simulation.
        data(dict): Arrays to hold clutter simulation data and
                    transform nav file to XYZ
    """
    print('\033[92m\033[1mSimulating {} traces: '.format(len(simParams['navxyz'])))
    # Run sim along navxyz track.
    for i in range(len(simParams['navxyz'])):
        # Checks data, bins, and locations and corrects them.
        simChecks(simParams, data, i)

        # Constructs left and right grid objects projected to the ground.
        grids(confDict, data, simParams, i)

        # Calculate two way travel time and power.
        py_sim(confDict, data, simParams['navxyz'][i], 'r')  # Right side of spacecraft.
        py_sim(confDict, data, simParams['navxyz'][i], 'l')  # Left side of spacecraft.

        # Two way travel time and facet center avgs per trace along track.
        averages(confDict, simParams, data, i)

        # Build cluttergrams and echomap_ref.
        cluttergram(confDict, simParams, data, i)

        # Print relevant info to the terminal for each trace.
        printInfo(simParams, i)

    redc = np.logical_xor(data['rcgram'],data['lcgram'])
    redc = np.logical_or(redc, np.logical_or(data['rRed'],data['lRed']))
    data['cRed'] = np.logical_not(redc)
    print('\033[92mSimulation Complete!                                          ')

def saveProducts(confDict, simParams, data):
    """
    Generates output products to `out_path` directory.
    
    Args:
        confDict(2D dict): Config parameters.
        simParams(dict): Parameters of runSim.
        data(dict): Misc arrays that hold left/right side clutter sims
                    and results.
    """
    # Checks for bad facets.
    if data['badf'] > 0:
        print('\n\033[91m' + str(data['badf']) + ' bad facets')
    else:
        print('\n' + str(data['badf']) + ' bad facets')

    print('\033[92m Saving products...', end="\r")
    data['nad'] = [50, 200, 200]  # Cyan for flight path.
    data['fre'] = [255, 0, 255]  # Pink for echo locations.

    # Handle cgram to prep for combined/adjusted products.
    data['cgram'] = data['rcgram'] + data['lcgram']
    for i,row in enumerate(data['cgram']):
        for j,pt in enumerate(row):
            data['cgram'][i][j] = 10*np.log10(pt+1)

    # Produce combined.png product.
    if confDict['outputs']['combined'] == 'True':
        combined(confDict, simParams, data, '')

    # Produce combined_adjusted.png product.
    if confDict['outputs']['combined_adjusted'] == 'True':
        data['cgram'] = data['cgram']*(255.0/(data['cgram'].max()-.1*data['cgram'].max()))
        data['cgram'] = np.minimum(data['cgram'], 255)
        data['cgram'] = np.maximum(15*(data['cgram'])**.5, data['cgram'])
        data['cgram'] = np.minimum(data['cgram'], 255)
        combined(confDict, simParams, data, '_adjusted')
    
    # Produce combined binary product.
    if confDict['outputs']['binary'] == 'True':
        binary(confDict, data)

    # Produce right.png product.
    if confDict['outputs']['right'] == 'True':
        sides(confDict, simParams, data, data['rcgram'], 'right')
    
    # Produce left.png product.   
    if confDict['outputs']['left'] == 'True':
        sides(confDict, simParams, data, data['lcgram'], 'left')

    # Produce echomap.png product.
    emap = data['emap']
    if confDict['outputs']['echomap'] == 'True':
        escale = emap*(255.0/emap.max())
        echomap(confDict, simParams, data,[escale, 1, simParams['num_ct']], False)

    # Produce echomap_adjusted.png product.
    if confDict['outputs']['echomap_adjusted'] == 'True':
        yDistortion = int(confDict['facet_params']['ct_step']) / simParams['nad_xyz'].avgdiff()
        emap[emap != 0] = emap[emap != 0] - (emap[emap != 0].mean()-emap[emap != 0].std())
        emap = emap*(255.0/(emap[emap != 0].mean()+emap[emap != 0].std()))
        emap = np.minimum(emap, 255)
        emap = np.maximum(0, emap)
        yDimension = int(emap.shape[0]*yDistortion)
        nemap = transform.resize(emap, (yDimension, len(simParams['navdat'])))
        echomap(confDict, simParams, data, [nemap, yDistortion, yDimension//2], True)
    
    # Produce echomap_ref.tif product.
    if confDict['outputs']['echomap_ref'] == 'True':
        driver = gdal.GetDriverByName('GTiff')
        dataSet = driver.Create(confDict['paths']['out_path'] + 'echomap_ref.tif', simParams['exs'], simParams['eys'], 1, GDT_Float32)
        dataSet.SetGeoTransform(simParams['egt'])
        ref = osr.SpatialReference()
        ref.ImportFromProj4(simParams['sys'][2]) # ecosys.
        dataSet.SetProjection(ref.ExportToWkt())
        dataSet.GetRasterBand(1).WriteArray(simParams['emap_ref'].astype('float32'))
        dataSet.GetRasterBand(1).SetNoDataValue(0)
        dataSet.GetRasterBand(1).SetUnitType("m")
        dataSet = None

    # Produce csv products.
    csv(confDict, simParams, data)

def simc():
    """ main() """
    confDict = readConfig()
    simParams, data = simPrep(confDict)
    runSim(confDict, simParams, data)
    saveProducts(confDict, simParams, data)

simc()