_The discrete part of the dLDS model described in:_
*** NOW AT JMLR ***
**[Noga Mudrik, Yenho Chen, Eva Yezerets, Christopher Rozell, Adam Charles. "Decomposed Linear Dynamical Systems (dLDS) for learning the latent components of neural dynamics". 2024. JMLR.](https://www.jmlr.org/papers/volume25/23-0777/23-0777.pdf)**

Abstract - Mudrik and Chen et al. 2024:

Learning interpretable representations of neural dynamics at a population level is a crucial first step to understanding how neural activity patterns over time relate to perception and behavior. Models of neural dynamics often focus on either low-dimensional projections of neural activity, or on learning dynamical systems that explicitly relate to the neural state over time. We discuss how these two approaches are interrelated by considering dynamical systems as representative of flows on a low-dimensional manifold. Building on this concept, we propose a new decomposed dynamical system model that represents complex nonstationary and nonlinear dynamics of time-series data as a sparse combination of simpler, more interpretable components. The decomposed nature of the dynamics generalizes over previous switched approaches and enables modeling of overlapping and non-stationary drifts in the dynamics. We further present a dictionary learning- driven approach to model fitting, where we leverage recent results in tracking sparse vectors over time. We demonstrate that our model can learn efficient representations and smoothly transition between dynamical modes in both continuous-time and discrete-time examples. We show results on low-dimensional linear and nonlinear attractors to demonstrate that decomposed systems can well approximate nonlinear dynamics. Additionally, we apply our model to _C. elegans_ data, illustrating a diversity of dynamics that is obscured when classified into discrete states.

This repository is subject to the MIT License. 

=================================================================
# A) Installation Instructions and Orientation - MATLAB version

NOTE: Please see the dev branch of this repository for the code used in Yezerets, E., Mudrik, N. & Charles, A.S. Decomposed Linear Dynamical Systems (dLDS) models reveal instantaneous, context-dependent dynamic connectivity in C. elegans. Commun Biol 8, 1218 (2025). https://doi.org/10.1038/s42003-025-08599-3

This code is known to be compatible with MATLAB versions from 2022 through 2023, but may have a few years of backward and forward compatibility, as well (for example, there is a known `cellfun()` call that does not work in version 2020a). Please reach out to our team if you have experience challenges with compatibility with other versions of MATLAB.

In order to use our code, please clone this repository. There is also [a version of this code in Python](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Python-Model) - the MATLAB code was used to work with the _C. elegans_ data.

All of the code, including scripts to run examples and External Packages, are included in the **code** directory. Data used to run the _C. elegans_ scripts is available from [Kato2015 whole brain imaging data](https://osf.io/2395t/). The [SSM](https://github.com/lindermanlab/ssm/blob/master/notebooks/4-Recurrent-SLDS.py) package from the Linderman lab was used to run recurrent switched linear dynamical systems (rSLDS) on the _C. elegans_ data. 

**Key scripts and functions:**

1) [code/script_flexibleMultiNetworkTesting.m](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/code/script_flexibleMultiNetworkTesting.m) - Start here to run a simulated example of two independent systems (no observation matrix - "noobs" option). The parameters here are selected to sample from this simulated data and learn a set of maximally distinct dynamics operators that most efficiently account for the independent states of the two simultaneously active systems. Additional functions called by this script: [generateContMultiSubNetwork.m](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/code/generateContMultiSubNetwork.m) creates the simulated data and [plotMultiNetworkOutput.m](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/code/plotMultiNetworkOutput.m) plots a comparison of the inferred and ground truth coefficients and learned dynamics operators.
2) [code/Dynamics_Learning/bpdndf_dynamics_learning.m]() - Expectation-Maximization algorithm. Observation matrix and dynamics operators get initialized (if not provided to the function) and updated here ([code/Learning_Code/Learning_Lib/dictionary_update()](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/code/Learning_Code/Learning_Lib/dictionary_update.m)). Calls [code/Dynamics_Learning/subfunctions/parallel_bilinear_dynamic_inference](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/code/Dynamics_Learning/subfunctions/parallel_bilinear_dynamic_inference.m) to return inferred coefficients `A_cell` (latent states _x_ in the paper) and `B_cell` (dynamics coefficients _c_). Note that throughout this code, the subscript `b`, as in `lambda_historyb`, is used to refer to the dynamics coefficients. This is a vestige of a previous version of the code.
3) [code/Dynamics_Learning/subfunctions/check_inf_params](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/code/Dynamics_Learning/subfunctions/check_inf_params.m) - See default `inf_opts` parameter settings here.
4) [code/Dynamics_Learning/subfunctions/parallel_bilinear_dynamic_inference](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/code/Dynamics_Learning/subfunctions/parallel_bilinear_dynamic_inference.m) - Returns inferred `A_cell` and `B_cell` as described above, given an observation matrix, dynamics operators, and observed data.
5) [code/script_Celegans.m](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/code/script_Celegans.m) - Example script run on _C. elegans_ data with recommended `inf_opts` parameters. These must be rescaled for each dataset given the average dimensions and distribution of the data. We used the same parameters across multiple worms from one dataset.

=================================================================
# B) Example figures

### Simulation example figure
![image](https://github.com/dLDS-Decomposed-Linear-Dynamics/dLDS-Discrete-Matlab-Model/blob/main/SimulationFigure.png)

### _C. elegans_ example figure from dLDS paper
![image](https://user-images.githubusercontent.com/90283200/171279482-fb59ffa1-8755-475a-a97b-c161afc615a1.png)
