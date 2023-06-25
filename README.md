# Multi Clutter Simulator

Builds surface clutter simulations to aid in the analysis of airborne and space-born sounding radar images.

## Description

### What is a clutter simulator?

The purpose of this program is to create surface clutter simulations to aid in the analysis of airborne and space-born sounding radar images. In the context of sounding radar, "surface clutter" refers to power reflected from off-nadir features on the surface of the body being sounded. The power reflected by these off-nadir features can reach the receiver at times similar or identical to when power reflected from subsurface interfaces would return. This creates apparently subsurface features in the radar image that are in truth from surface features, hindering interpretation. With an accurate topographic model and knowledge of spacecraft location over the course the the radar sounding, the power returned from the surface can be roughly simulated to produce an image that can be qualitatively compared to the radar image to help determine if a candidate subsurface reflector is in fact from a surface feature.

### Overview

The MultiSim Clutter Simulator (MultiSim) was written to replace UTSIM, an older clutter simulator written primarily in C to simulate clutter on Mars, specifically for the SHARAD instrument on the MRO. The motivation behind MultiSim was to write a more flexible clutter simulator that could accept many types of input data, and simulate clutter on several different planets.

### Cluttergram datum

Cluttergrams can be produced datumed to the transmitter, or can be
warped/shifted to be datumed to a different elevation.

### *Note: vertical datums*

The navigation data, DEM, and datum information must all be from an identical vertical datum when they hit the program. Coordinate systems are not an issue, these are readily translated by the GDAL library, but vertical datums are often not included with coordinate system information, or are complex shapes (such as a geoid). Because of this, the user must make sure that all of the input data is referenced to the same vertical datum.

## Setup / Dependencies

**Read [here](./environment.yml) for required external libraries:**

Use the package manager [Anaconda](https://docs.anaconda.com/anaconda/install/) to build a sim environment.

```bash
conda env create -f environment.yml
```

## Config File Structure

```conf
[paths]
dem_path = [Path to the DEM to be used in the simulation]
nav_path = [Path to the Nav file to be used in the sim]

[navigation]
navsys = [Proj4 or WKT of nav file coordinate system]
navfunc = [Name of the function to read navigation data]

[sim_params]
speedlight = [Speed of light to use in meters/second]
binsize = [Sampling period of radar data in seconds]
tracelen = [Number of samples in each trace]
body = [Body that is being sounded]

[facet_params]
at_dist = [Distance to the front and rear for each nav point]
ct_dist = [Distance to either side]
at_step = [Along track facet dimension]
ct_step = [Cross track facet dimension]

[datum_params]
datum = [Datum information]
datum_sample = [Sample in radargram to assign datum to]

[display_params]
show_nadir = [Show nadir return in output image]
show_fret = [Show first return in output image]

[outputs]
combined_adjusted = [Output a combined adjusted product]
combined = [Output a combined product]
left = [Outut a left product]
right = [Output a right product]
binary = [Output a binary product]
echomap = [Output a echomap product]
echomap_adjusted = [Output an echomap product]
echomap_ref = [Output an echomap ref product]
red = [Output a red echomap product]
nadir = [Output nadir csv data]
fret = [Output fret csv data]
tmap = [Output range migration csv]
```

**Exceptions:**

* at_dist must be less than at_step.
* ct_dist must be less than ct_step.
* at_dist must be divisible by at_step.
* ct_dist must be divisible by ct_step.

## Running Simulations

Usage:

```bash
python3 ./simc.py <config_file>
```

Example:

```bash
python3 ./simc.py ../config/sharad_fpb.cfg
```

## Contributing

### [**Univerisity of Arizona | Lunar & Planetary Labs**](https://www.lpl.arizona.edu/)

[TAPIR | Terrestrial And Planetary Investigation with Radar](https://tapir.lpl.arizona.edu/)

**James O'Connell | Research Assistant**
