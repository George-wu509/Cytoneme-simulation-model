# Guild to to run create3Dhex_dir to build cytoneme simulation model
cytoneme-mediated morphogen transport model made by George.Wu 

Outline of this document
=================

This document introduces the code structure of create3Dhex_dir.m matlab file which contains everything you will need to build cytoneme simulation model and follow-up analysis. Here is the outline of this guild including:

  * [1: How to simulate the first cytoneme model ](#How-to-simulate-the-first-cytoneme-model)
  * [2: Code structure](#Code-sttructure)
  * [3: Setup parameters](#ch-3-linear-regression)
  * [4: Main functions](#ch-4-support-vector-machines)
  * [5: simulation results](#ch-5-nearest-neighbor-methods)

![image](https://github.com/George-wu509/Cell-3D-segmentation-display-GUI/blob/master/cover/Segshow3D%20cover2.png)

## How-to-simulate-the-first-cytoneme-model
## [1: How to simulate the first cytoneme model ](##How-to-simulate-the-first-cytoneme-model)

To run the first cytoneme model:
step 1: change the parameter values in function p=parameter_set() (line 50)
For example 1:
To test different hypothesis cytoneme model, you can change: p.AN, p.LE, p.PF, p.MA
For example 2:
To change wing disc physical model, you can change: p.nx, p.ny, p.w_dppcent
For example 3:
To choice morphogen, you can change: p.morphogen_init (Dpp =1, Hh=2), p.amount
For example 4:
To display or save figures, you can change: p.img, p.show_every_img, p.save_every_img

step 2: >> create3Dhex_dir(1) (type in matlab command line and ENTER)

step 3: finished!
All results will be stored in cytoneme.mat including p(parameter), re_txt(explanations of re), sta_all(all generation results), sta_save(generation results you displayed), re('number of cytonemes in distance range groups';'#cytoneme events';'# cytonemes';'mean cytoneme distances’).

## [2: Code structure ](##Code-sttructure)

This is what you will see when you collapse the create3Dhex_dir.m file. There are totally 51 functions in this single matlab file.

![image](https://github.com/George-wu509/Cell-3D-segmentation-display-GUI/blob/master/cover/Segshow3D%20cover2.png)

The first function create3Dhex_dir(key, p_pre) is the main function of this cytoneme simulation model. You can choice which script you want to run using the first input: key. p_pre is the per-parameter variable. If you want to run the model using default parameters setup, you don’t need to input p_pre. You can change the default model parameters in the second function parameter_set() before running the model. Script1~7 are script functions.
Functions in ‘’2. main program’’ are important core functions. create_mesh() will create wing disc hexagonally epithelial cell physical model, initial_state() will initiate the morphogen concentration distribution and the model. Update_apical_dpp() will calculate the morphogen transport status.
Example:
To build one cytoneme simulation model and run it, you can type this in matlab command line and run:
> create3Dhex_dir(1);
To batch run several cytoneme simulation models:
> create3Dhex_dir(3)


