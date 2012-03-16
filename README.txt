********************
*  READ CAREFULLY  *
********************


- The EV_RL.r file contains the algorithm used to estimate both the return levels 
  (for a number of desired return periods) and their confidence intervals pixel per pixel. 
  Also, the algorithm estimates the parameters of the extreme value distributions and their 
  confidence intervals using the Maximum Likelihood Estimator (MLE). 

- The EV_RP.r file contains the algorithm used to estimate the probability of exceedance
  (or return period, if you take the inverse) for some desired return level.
  Also, the algorithm estimates the parameters of the extreme value distributions and their 
  confidence intervals using the Maximum Likelihood Estimator (MLE). 

- The myEVfunctions.r contains all the external modules called by the main scripts (either EV_RL.r 
  or EV_RP.r). Therefore, always place this file in the same folder with the main script.


*****
NOTE:
*****

At this stage, the algorithm is developed on our case study, thus it assumes that the input
dataset is formatted according to certain rules. We look forward to expand this repository 
by adding the code that could be used to read the native satellite image files and format
the spatiotemporal database accordingly. Also, we used the algorithm over monthly minima. Therefore,
you would need to manually adjust the code (looking at our comments throughout the code)
if you are using monthly maxima or, more in general, minima/maxima over a different time
frame such as year, season, etc.



**************************************************************
Contact me for any further information or error messages:
f_tonini@hotmail.com
**************************************************************