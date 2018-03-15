# GTEDGE3 Codebase

Georgia Tech GTEDGE3 Python-based Codebase for interpreting edge plasmas in tokamak fusion reactors

## Version history

*  1.0 initial commit

## Package dependencies

* Python 2.7
* Matplotlib
* SciPy
* Numpy

## Getting Started

Clone the repository to your local directory. GTEDGE2 comes in two code versions that should be maintained simultaneously: GTEDGE3 and GTEDGE3_cli.

The following input files are needed to run GTEDGE2:

	Inputs\GT3Consts_shotid_timeid.csv				Plasma constant values from DIII-D
	Inputs\GT3NBI_shotid_timeid.csv					Neutral beam data for NBeams from DIII-D
	Inputs\GT3Profs_shotid_timeid.csv				Plasma profiles from DIII-D
	Inputs\GTEDGEsupp_shotid_timeid.csv				Supplementary data from GTEDGE Fortran 77 version. Soon to be deprecated.
	MaxPlasma\shotid_timeid_toplasma				Background plasma data needed for Max Hill code
	MaxPlasma\inputs\shotid_timeid_runid_X.txt			Boundary and plasma data needed for Max Hill code from DIII-D
	MaxPlasma\inputs\pshotid_timeid_runid_ersplrhob_sc.txt		Radial electric field file from DIII-D
## Usage

# GTEDGE3.py
GTEDGE3.py can be run from the command line with 

	python GTEDGE3.py
	
To use this functionality, ensure that the correct shot information is set in the file. The following is an example for DIII-D shot 118890, timeid 1560, runid r90, with NBeams activated and IOL correction applied

	shotid=118890
	timeid=1560
	runid="r90"
	nbRun=True	
	IOL=True
	
# GTEDGE3_cli.py

GTEDGE3 command line version can be run from the command line:

	usage: GTEDGE3_cli.py [-h] [-nbRun [NBRUN]] [-IOL [IOL]] shotid timeid runid

	positional arguments:
	  shotid          DIII-D shot id
	  timeid          DIII-D time id
	  runid           DIII-D run id

	optional arguments:
	  -h, --help      show this help message and exit
	  -nbRun [NBRUN]  Run NBeams (default: True)
	  -IOL [IOL]      Correct for IOL (default: True)

## Usage

# GTEDGE3.py
GTEDGE3.py can be run from command line as follows:

	python GTEDGE.py
With this method, please ensure that

## File structure

	.
	├── Inputs
	│   └── fullplasma     # Storage of data files for shots from DIII-D database. GT2 does not call these directly
	│       ├── 118888     # Example directory of shot 118888 with 2 timeids
	│       │   ├── 1525
	│       │   └── 1570
	├── lib                         # Library of packages
	│   ├── beams			# NBeams directory
	│   │   └── NBeamsMDS
	│   │       ├── 118888		# Example of stored data from DIII-D database for NBeams. GT2 does not call these directly
	│   │       │   ├── 1525
	│   │       │   └── 1570
	│   │       ├── Debug		# Fortran 77 compiler debug folder. Used by GT2.
	│   │       ├── Fastloss	# Unused fastloss calculations
	│   │       │   └── Debug
	│   │       ├── nbeams2gtedge	# Deprecated 
	│   │       ├── OSM		# Deprecated
	│   │       │   └── Debug
	│   │       ├── Release		# Fortran 77 compiler release folder. GT2 does not call these directly
	│   │       └── Test		# Fortran 77 test folder.GT2 does not call these directly
	│   │           ├── Debug
	│   │           └── Release
	│   ├── funcs			# Repository of packages used to perform calculations in GT2
	│   └── graphs			# Repository of codes for making pretty pictures for GT2.
	├── MaxPlasma			# Max Hill's background plasma calculation IOL calculations. Called by GT2. This folder includes shotid_timeid_toplasma files required for code
	│   └── inputs			# Repository of data from DIII-D database. Called by GT2.
	│       └── Originals		# Deprecated
	├── Outputs                     # GTEDGE2 outputs .dat files of the previous run if such a file cannot be found. If found, GTEDGE does not run again on the data. Useful for analysis on results.
	└── RawDIIID			# Raw DIII-D data from servers
	    └── 118890
		   ├── 1515
		   └── 1560
		       └── test

## Authors

* **Jonathan Roveto** - *Conversion from F77 to Python* - [veto1024](https://github.com/veto1024)
* **Maxwell Hill** - *Provided background plasma code and IOL code* 

## Acknowledgments

* Georgia Tech for allowing me to pay them thousands of dollars to do work for them
