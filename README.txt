README file - Critical weaknesses in shielding strategies for COVID-19

--------------------------------------------------------------------------------------------

Contact information:
Name: Cameron Smith
Email: cs640@bath.ac.uk

--------------------------------------------------------------------------------------------

This is the README file for the various parts of code that can be found in the ESM.

--------------------------------------------------------------------------------------------

./Workspaces/ (Directory)

Contains all of the workspaces required for plotting the main text and supplementary plots. The numbers and suffixes are required in order to plot the correct workspaces.

--------------------------------------------------------------------------------------------

./herd_immunity_main_execution.m (MATLAB script)

Code for the implementation of the SEIR model described in the paper. Saves the workspace as a .mat file, which can be used for plotting. Requires the c-mex file "quickfind". For more information on mex files please see:

https://uk.mathworks.com/help/matlab/call-mex-files-1.html.

Default parameter values are provided, and given in a comment.

--------------------------------------------------------------------------------------------

./main_text_plots.m (MATLAB script)

Plots the five figures from the main text. The ordering of the workspaces in the directory are important for plotting in the correct order.

--------------------------------------------------------------------------------------------

./supp_plots.m (MATLAB script)

Plots the seven figures from the main text. The ordering of the workspaces in the directory are important for plotting in the correct order.

--------------------------------------------------------------------------------------------

./plot_singular_workspace.m (MATLAB script)

Allows for a single workspace to be plotted. Replace the *** in line 15 with the name of the workspace of the form 'workspace.mat' or, if in another directory, './directory/workspace.mat'.

--------------------------------------------------------------------------------------------

./quickfind.c (C function)

A function to find the index of the next event that has fired using the Gillespie algorithm. 

INPUTs:
-------
a - vector of propensity function values
sum - set to 1
r - uniform random variable between 0 and sum(a)

OUTPUTs:
--------
j - index of the next firing event

--------------------------------------------------------------------------------------------

./quickfind.mex64 (mex file)

Mex file for the quickfind.c function. For more information, see: 

https://uk.mathworks.com/help/matlab/call-mex-files-1.html.
