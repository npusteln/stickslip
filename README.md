# stickslip
Nonlinear denoising for solid friction dynamics characterization
***************************************************************************
* authors: Nelly Pustelnik and Valerie Vidal                              
* institution: Laboratoire de Physique de l'ENS de Lyon                   
* date: April 30 2019                                                   
* License CeCILL-B                                                        
***************************************************************************
*********************************************************
* RECOMMENDATIONS:                                   	*
* This toolbox is designed to work with Matlab 2017b    *
*********************************************************

------------------------------------------------------------------------------------------------------------------------
DESCRIPTION:

This work investigates the transition between different regimes in 
paper-paper friction. To probe solid friction with the fewest contact 
damage possible, the load applied to the paper samples has to be small. 
This constraint, as well as other experimental limitations, introduces 
noise - in particular, impulsive - on the slider motion and force signals. 
This study proposes variational techniques to process the signal and 
extract the slider velocity. A unique nonlinear filtering technique relying 
on recent nonsmooth optimization algorithms is retained in order to denoise 
both the stick-slip and inertial regimes. Although by definition it leads to 
stepped velocity profiles, it makes possible, in the stick-slip regime, 
to capture the slider velocity asymmetry and, thus, the creep motion before
sliding. A precise estimation of the time spent in the stick regime and in 
motion can thus be inferred, and a criterion is proposed to quantify the 
transition between the stick-slip and inertial regimes. A regime diagram is 
then proposed for the dynamics of this solid frictional system.

------------------------------------------------------------------------------------------------------------------------
SPECIFICATIONS for using DMS toolbox:

The main function for denoising L2L1 is "PD_ChambollePock.m".
The main function to detect t_start,t_stop, tau_m, tau_s is "detect_tstartstop.m"

Two demo files are provided in examples :
- ex1 : denoising + t_start and t_stop detection for a stick-slip regime 
        signal with k = 168 and v = 42.
- ex2 : denoising + t_start and t_stop detection for an inertial regime 
        signal with k = 168 and v = 4300.

The file plot_signal_example_Fig8.m makes it possible to reproduce any panel 
in Figure 8 of the article.

\data_raw contains examples of raw data Fnorm vs. t
\data_Fig8 contains all raw and filtered signals displayed in the panels of Fig.8.

------------------------------------------------------------------------------------------------------------------------
RELATED PUBLICATION:

J. Colas, N. Pustelnik, C. Oliver, P. Abry, J.-C. Geminard, V. Vidal, 
"Nonlinear denoising for solid friction dynamics characterization", 
accepted to Physical Review E, 2019.

------------------------------------------------------------------------------------------------------------------------
