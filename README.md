# Planarian Shape Regulation During Growth and Degrowth

Supplementary Data for:

> Mechanistic Regulation of Planarian Shape During Growth and Degrowth<br>
> Jason Ko, Waverly Reginato, and Daniel Lobo<br>
> Submitted, 2023<br>


## Installation

The source code is located in the `src` folder and can be directly executed with MATLAB (The Mathworks, Inc.) from that folder. It has been tested with release R2021a. The simulations use [ROWMAP](https://www.mathematik.uni-halle.de/wissenschaftliches_rechnen/forschung/software/#anchor1303104) to numerically solve the equations. To install it, simply download the file [rowmap.m](http://sim.mathematik.uni-halle.de/software/rowmap/rowmap.m) into the `src` folder.


## Growth and degrowth simulations

To run the simulations, execute any of the following commands in the MATLAB command prompt, replacing any variables in ALL_CAPS with real data:

For growth simulation, using manual parameters, run:
```matlab
manualGrowthSimulation(NEW_SIM_ID)
```

For growth simulation, using the calibrated parameters by the evolutionary algorithm, run:
```matlab
evolvedGrowthSimulation(NEW_SIM_ID)
```

For degrowth simulation, run:
```matlab
degrowthSimulation(NEW_SIM_ID)
```

where `NEW_SIM_ID` is a string identifier for the output directory (e.g. 'sim1').

To run some other parameter set, we recommend the configuration in one of the above files.


#### Note:
As cached simulation data and video data are created, new subdirectories with the results will be automatically created under the `results` directory.


## Evolutionary algorithm
To run the evolutionary algorithm, execute the following commands in a Bash shell:

```shell
JOB_NAME=run1
	
export MATLAB_LOG_DIR=~/$1
mkdir "$MATLAB_LOG_DIR"
 
TMPDIR=$(mktemp -d /tmp/matlab-lobo-XXXX)
trap 'rm -rf "$TMPDIR"' EXIT
 
matlab -nodisplay -r "geneticAlgorithm('$TMPDIR', 36, inf, '$JOB_NAME', 36); exit"
```

where `JOB_NAME` can be any identifying string for the search being run.


These are the most important scripts for running the evolutionary algorithm:

  * geneticAlgorithm.m
  * configureProject.m
  * configureSimulation.m
  * runSimulation.m


After running the evolutionary algorithm, the results can be plot with:

  * plotGA.m


## Image processing pipeline
The following scripts implement the automatic image analysis pipeline:

  * extractShape.m
  * extractShape2.m


