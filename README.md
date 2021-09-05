# ztfgkncr
Zwicky Transient Facility (ZTF) Grandma Kilonova Catcer ligthcurve Regression 
## Analyse, treatement, light curbe, and Linear regressions 
of KN-Catcher-ReadyforO4 Data¶ and others transient 
## Summury 

**kn_catcher.py**
    That file contains all the functions to obtain (using the API key of alerte), 
    the updated data of the targets as well as the loading of the data of observations
    provided by Grandma's amateurs (KN-Catcher-ReadyforO4 data). Data such as magnitude, 
    observation time and errors 
    on magnitudes, etc, are extracted and saved in csv format.
    
 **pwl_catcher.py**
    csv file will be used by in this to perform regression using a few methods:
    a Bayesian method with pymultinest (can be installed with conda) and linear regression 
    models, such as scipy.stats, sklearn and Exhaustive Statistics (allowing the determination
    of the maximum likelihood).
 
 **launch.py**
    Runs the first two.   
    In output, a first folder outputs will be created, within which a sub-folder ligthcurve 
    destined to hold the light curves of each event. The creation of a second sub-folder will 
    allow one to save the light curve of each filter as well as the superposition of the 
    different regression models.
