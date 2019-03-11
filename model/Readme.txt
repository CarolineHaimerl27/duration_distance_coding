Readme file for Computational Model, Haimerl et al. PNAS (2019)

========================================================



call_model.m: 	define parameters and call model simulation and evaluation, 
		graphs: 	facilitation, depression, firing rates and theta input over time (Fig.2.C)
				sequence slope over input power (Fig.2.D)

dynam_model.m : dynamics of firing rate
FRModel.m :	simulate dynamics

Helping functions:
RunSeqSlope.m	compute the slope of firing rate sequences
cont_slop.m:	alternative way to compute the slope of firing rate sequences
SLOPE.m		wrapper for cont_slop function, finds highest activity neuron per
		theta cycle to compute the slope with
ThSlope.m: 	computes the average slope of the theta sequence (formed by the
		subgroup of cells active within a theta cycle)
