#1. Fig2.m
This file contains the code to plot the 1-D manifold and equilibrium points for the slow-fast consumer-resource system given by Eqn. 1.
Setting the Temperature T=287, the code gives the time evolution of the C-R dynamics on the manifold.

##one can change the value of T to obtain trajectories for different temperatures\

#----------------------------------------------------------------------------------------------------------------------------------------------\
#----------------------------------------------------------------------------------------------------------------------------------------------\

#2(a). FIG3_QSSA.m
This files generate the quasi-steady-state attractor for both the cases, setting T<Topt and T>Topt.

#2(b). FIG3_rate.m
This codes is used to obtain the species density along the temperature axis for varying rates of warming, alpha.
# alpha1 is the rate corresponding to the case T<Topt.
# alpha2 is the rate corresponding to the case T>Topt.

#----------------------------------------------------------------------------------------------------------------------------------------------\
#----------------------------------------------------------------------------------------------------------------------------------------------\

#3(a). Fig5_negative.m
This file generates phase portrait of the de-singularized system for T<T_opt 

#3(b). Fig5_positive.m
This file generates phase portrait of the de-singularized system for T>T_opt 

#----------------------------------------------------------------------------------------------------------------------------------------------\
#----------------------------------------------------------------------------------------------------------------------------------------------\

#4(a). Fig6_manifold.m
This codes gives the two-dimensional manifold, the plot for the fold line and the QSSA on the manifold.

#4(b). Fig6_canard_negative.m
This is used to obtain the singular canard upon the two-dimesional surface of the manifold for T<Topt.

#4(c). Fig6_canard_positive.m
This is used to obtain the singular canard upon the two-dimesional surface of the manifold for T>Topt.
#NOTE:  start with code explained in 4(a) and the supeimpose the figures obtained through codes in 4(b) and 4(c) on it.

#----------------------------------------------------------------------------------------------------------------------------------------------\
#----------------------------------------------------------------------------------------------------------------------------------------------\

#5. noise_t_series.m

This .m file generates and stores the time series data of the stochastic system. The time series thus obtained is then used for EWSs analyses (using Rstudio ) followed from reference 1.



#NOTE: We use MATLAB version R2015b or later to run each matlab programme\
#NOTE: R version 3.6.3 or later. Open using any platform such as Rstudio or Anaconda.\

Reference\
--------------
1.Dakos V, Carpenter SR, Brock WA, Ellison AM, Guttal V, Ives AR, K\'e9fi S, Livina V, Seekell DA, van Nes EH, Scheffer M (2012a) Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data. PLoS ONE 7:e41,010\
