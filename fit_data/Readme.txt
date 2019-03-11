Readme file for Data fitting Haimerl et al. PNAS (2019)

========================================================

"Data" folder contains the 34 recording sessions analyzed in Fig. 3.
Each ".mat" file has two vector structures Speed of the mouse (cm/s) 
and Velocity of the sequence (neurons/s).

----

1.main_model.m creates the two dimensional diagram in Fig 3B, calculating
the sequence velocity predicted by the model in the two parameter space 
(I_0,\delta). The velocity is calculated integrating the dynam_model.m
equations. Two auxiliary functions are used in the integration of equations.

Basic_model.m : Defines the connectivity matrix and calls the ODE45 routine.
My_events.m   : An event detection function that halts the integration when
		reaching an unstable state where no sequences are generated.

The result of this script is stored in the .mat file PhaseDiagram.mat
with the following data structures 

Vector I_0: 	Vector with the values of I_0 used
Vector delta:	Vector with the values of delta used
Data Matrix:	Matrix where element (i,j) is the resulting velocity of the model
when using the values of the parameters I_0(i) and delta(j).

----

2. Fit_of_Data.m script loads the results from 1.main_model.m stored in 
Phase_diagram.mat. Then, it loads the data from the recording sessions
and performs a minimization procedure to find the best values of 
(I_0, delta) that fits the sequence velocities calculated with the model
and the recorded data. At the end of the script a 34x3 Matrix named "Parameter"
contains the information of the best values of (delta,alpha,beta) defined in
Supplementary material that fitted the data in each recording session and stores 
this matrix in the file "Parameters.mat" for post-processing.

