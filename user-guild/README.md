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



<a name="How-to-simulate-the-first-cytoneme-model" />
## 1: How to simulate the first cytoneme model 

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

<a name="Code-sttructure" />
## 2: Code structure

This is what you will see when you collapse the create3Dhex_dir.m file. There are totally 51 functions in this single matlab file.

![image](https://github.com/George-wu509/Cell-3D-segmentation-display-GUI/blob/master/cover/Segshow3D%20cover2.png)

The first function create3Dhex_dir(key, p_pre) is the main function of this cytoneme simulation model. You can choice which script you want to run using the first input: key. p_pre is the per-parameter variable. If you want to run the model using default parameters setup, you don’t need to input p_pre. You can change the default model parameters in the second function parameter_set() before running the model. Script1~7 are script functions.
Functions in ‘’2. main program’’ are important core functions. create_mesh() will create wing disc hexagonally epithelial cell physical model, initial_state() will initiate the morphogen concentration distribution and the model. Update_apical_dpp() will calculate the morphogen transport status.
Example:
To build one cytoneme simulation model and run it, you can type this in matlab command line and run:
> create3Dhex_dir(1);
To batch run several cytoneme simulation models:
> create3Dhex_dir(3)

<a name="ch-3-linear-regression" />
## 3: Setup parameters

You can change the default parameter values in p=parameter_set() function. p is a matlab structure variable which contains all model parameters. Basically this is the only place you will need to change in this code. Here are some important parameters to run cytoneme simulation listed below:


<a name="ch-4-support-vector-machines" />
## 4: Main functions

p=create_mesh(p)
This function creates the wing disc hexagonally epithelial cell 2D and 3D physical model with user-defined sizes and properties which you can define in mesh parameters. For example: wing disc size = p.nx(=40) * p.ny(=20) means this simulated wing disc contains 800 hexagonally epithelial cells. Related parameters: p.nx, p.ny, p.r, p.w_dppcent, p.height
sta=initial_state(p)
initial_state function output the initial morphogen concentration on Drosophila wing disc at time = 0. In Dpp morphogen case, Dpp is produced from morphogen-producing cells in Dpp producing center region which controlled by parameter p.w_dppcent. This function will output matlab struct variable sta. sta.dpp is the morphogen concentration distribution on wing disc, and sta.cyto is cytoneme on wing disc.
Related parameters: p.nx, p.ny, p.morphogen_init, p.prod_init, p.w_dppcent
sta=Update_apical_dpp(sta,p)
This function calculate and update the cytoneme growing and morphogen transport in every generation.
Step 1: sub function dpp_grow_reaction(sta_old,p) will calculate the morphogen producing and the morphogen reaction with receptor for every cells and update the variable sta.
Step 2: In every cell in every time generation, function target_position will calculate the probability to grow cytoneme, and select the potential target cell which cytoneme will grow and extend from this cell to that target cell.
Step 3: We assume cytoneme-mediated transport is viewed as a stochastic process. If a random number between 0 to 1 is larger than cytoneme formation probability we get from the previous step, cytoneme will extend from this cell to target cell. If not, there is no cytoneme growth.
step 4: If cytoneme growth is happened, amount=morphogen_amount(d,p) will decide how much morphogen will transport from morphogen-producing cell to morphogen-receiving cell.
      
step 5: calculation process will be applied to every cells on wind disc, and
then update and output new dpp concentration distribution and cytoneme status.
[sta,p]=show_figure(sta,p,ge,sta_save)
This function will display and save result figures. If p.img =0, then the code will not run this function. If p.show_every_img = 1, it will display every figures in the running process from t=0 to t=max generation (p.max_gener). If p.save_every_img = 1, it will save every figures in the running process from t=0 to t=max generation


<a name="ch-5-nearest-neighbor-methods" />
## 5: simulation results


