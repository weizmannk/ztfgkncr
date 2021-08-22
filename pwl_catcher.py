import time
import os
import glob
import pandas as pd
from astropy.io import ascii
import numpy as np
import optparse
import corner
import pymultinest
from matplotlib import pyplot as plt
from scipy import stats
from sklearn.linear_model import LinearRegression


# =============================================================================
# returns current working directory of a process.
# =============================================================================
datapath = os.getcwd()
#datapath = os.getcwd() + "/" + "outputs/csv"
# =============================================================================
#  input the csv file
# =============================================================================
csv_file = glob.glob((datapath+'/**/*.csv').split("/")[-1], recursive=True)
save_file_name = [f.split(".")[0] for f in csv_file]

# =============================================================================
# scipy  linear regeression model
# =============================================================================
def scipy_linear_model(t, a, b):
    """[linear regression ]

    :param t: [delay of observation]
    :type t: [array]
    :return: [the regression model for prediction]
    :rtype: [type]
    """
    return a*t + b

# =============================================================================
# sklearn linear regeression model
# =============================================================================
def sklearn_linear_model(t, mag):
    """[Linear Regression model using sklearn]

    :param t: [time of observation in days]
    :type t: [float]
    :param mag: [magnitude of observation]
    :type mag: [float]
    """
    sk_t = np.array(t).reshape((-1, 1))
    result = LinearRegression().fit(sk_t, mag)
    score = result.score(sk_t, mag)
    model = result.predict(sk_t)
    
    return model, score

# =============================================================================
# linear regeression  using Max liklihood bayesain  model
# =============================================================================
def sufficient_statistics(t, mag, magerr):
    """[summary]

    :param t: [observation time delays  ]
    :type t: [float in jd time]
    :param mag: [magnitude of observation]
    :type mag: [float]
    :param magerr: [magnitude error]
    :type magerr: [float]
    :return: [ Matrix ]
    :rtype: [type]
    """
    x = t
    y = mag
    sigma = magerr
    w = 1.0 / sigma**2
    So = np.sum(w)
    Sx = np.sum(w * x)
    Sy = np.sum(w * y)
    Sxy = np.sum(w * x * y)
    Sxx = np.sum(w * x * x)
    Matrix = np.array([[Sxx, Sx], [Sx, So]])
    vector = np.array([Sxy, Sy])
    return Matrix, vector


# =============================================================================
# Read csv file and  create a base directory (folder) to save the plot 
# =============================================================================
def parse_commandline(csv):
    """[Loadind the data of cvs file ]

    :param folder: [output name to save the plots]
    :type folder: [str]
    :param csv: [ intput csv data ]
    :type csv: [type]
    :return: [return a dictionary contain plot directory, cvs file]
    :rtype: [optparse.Values]
    """

    parser = optparse.OptionParser()
    parser.add_option("--doPlots",  action="store_true", default=False)
    parser.add_option("-p", "--plotDir", default= "outputs")
    parser.add_option("-l", "--lightcurve", default=csv)

    opts, args = parser.parse_args()

    return opts

# =============================================================================
# Regeression using the pymultinest (bayesain model)
# =============================================================================
def plaw(t, f0, alpha=1.0, t0=None):
    """
    Power law function:
    f / f0 = (t - t0) ** -alpha
    """
    return f0 * (t - t0) ** -alpha

def pmag(t, f0, alpha=1.0, t0=None):
    """
    Use indices from plaw fit to calculate magnitudes
    -2.5 log(f/f0) = 2.5 * alpha * log(dt)
    m = -2.5 log(f)
    m + 2.5 log(f0) = 2.5 * alpha * log(dt)
    """ 
    return 2.5 * alpha * np.log10(t - t0) - 2.5 * np.log10(f0)

def eflux(emag, flux):
    """
    Error propoagation:
    m - m0 = dm = 2.5 log(1 + df / f0)
    10.0 ** (dm / 2.5) = 1 + df / f0
    df = f0 * (10.0 ** (dm / 2.5) - 1)
    """
    return flux * (10.0 ** (emag / 2.5) - 1)

def myprior(cube, ndim, nparams):

    cube[0] = cube[0]*(tmax - tmin) + tmin
    cube[1] = cube[1]*3.0
    cube[2] = cube[2]*20.0 - 10.0

def myloglike(cube, ndim, nparams):
    
    t0fit = cube[0]
    alpha = cube[1]
    f0 = 10**cube[2]

    mod = pmag(t, f0, alpha=alpha, t0=t0fit)
    
    idx1 = np.where(~np.isnan(y))[0]
    idx2 = np.where(np.isnan(y))[0]

    chisq = -(1/2)*np.sum(((y[idx1] - mod[idx1])**2.) /
                            (dy[idx1]**2.))/(len(t[idx1]) - nparams)

    #gaussprobvals = np.sum(np.log(1-scipy.stats.norm.cdf(dy[i
    # dx2], mod[idx2], 0.1)))

    return chisq

# =============================================================================
# extraction and filters cvs data
# =============================================================================
for i in range(len(csv_file)):
    # Parse command line
    opts = parse_commandline(csv_file[i])
    baseplotDir = opts.plotDir
    if not os.path.isdir(baseplotDir):
        os.makedirs(baseplotDir)
    

    # =============================================================================
    # Read csv file 
    # =============================================================================
    lc = ascii.read(opts.lightcurve, format='csv')
    parameters = ["t0", "alpha", "f0"]
    n_params = len(parameters)

    n_live_points = 1000
    evidence_tolerance = 0.1
    max_iter = -1

    # =============================================================================
    # Using filter to create the data table
    # =============================================================================
    for filter in lc['filters']:
        tab = lc[(lc['filters'] == filter)]
        
        # =========================================================================
        # select data with more than 2 filters i.e len (tab["filter"]) >2 
        # =========================================================================
        if len(tab) > 3:
            tab.rename_column("magerr", "e_mag")
            tab.rename_column("time", "mjd")

            # =====================================================================
            # data selected and reclasse by time increasing0 
            # =====================================================================
            tab.sort('mjd')
            idx = np.where(tab['mag'] > 50)[0]
            tab['mag'][idx] = np.nan

            # ======================================================================
            # time observation  max and min
            # ======================================================================
            tmax = np.min(tab[tab['mag'] < 50]['mjd'])
            tmin = np.min(tab[tab['mag'] < 50]['mjd'])-3

            # ======================================================================
            # print and listed, time, magnitude and instrument 
            # ======================================================================
            print("Data that will be used for the fit:")
            print("MJD, mag, instrument")
            for l in tab:
                print(l['mjd'], l['mag'], l['instrument'])
                print("---")
            
            # ======================================================================
            # create new folder by each filter to put each filter plot in 
            # ======================================================================
            plotDir = os.path.join(baseplotDir, save_file_name[i]+"/"+filter)
            if not os.path.isdir(plotDir):
                os.makedirs(plotDir)

            # ======================================================================
            # redefine the names of  time , magnitude and mag-error  
            # ======================================================================
            t = tab['mjd']
            y = tab['mag']
            dy = tab["e_mag"]
            dy = np.sqrt(dy**2 + 0.1**2)  # np.mean

            # =============================================================================
            # Make a regeression with pymultinest (bayesian method)
            # =============================================================================
            pymultinest.run(myloglike, myprior, n_params, importance_nested_sampling=False, resume=True, verbose=True, sampling_efficiency='parameter',
                            n_live_points=n_live_points, outputfiles_basename='%s/2-' % plotDir, evidence_tolerance=evidence_tolerance, multimodal=False, max_iter=max_iter)

            multifile = os.path.join(plotDir, '2-post_equal_weights.dat')
            data = np.loadtxt(multifile)

            #stat_pred = os.path.abspath(".") + "/" + os.path.join(plotDir,  '2-stats.dat')
            #stat_pred = pd.read_table(stat_pred, sep='\t')
            
            
            # =============================================================================
            # plot corner and save the data in  plotName directory 
            # =============================================================================
            labels = [f"$t_0$", f"$\\alpha$", f"$\\log_{10} f_0$"]
            plotName = "%s/corner.pdf" % (plotDir)
            figure = corner.corner(data[:, :-1], labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={
                                    "fontsize": 24}, label_kwargs={"fontsize": 24}, title_fmt=".2f", smooth=3,  color="coral")
            figure.set_size_inches(14.0, 14.0)
            plt.savefig(plotName)
            plt.close()
            
            # =============================================================================
            # Call for posterio  mag calculate by using bayesian methode with pymultinest
            # =============================================================================
            idx = np.argmax(data[:, -1])
            t_0, alpha, f0 = data[idx, :-1]
            tt = np.linspace(np.min(t), np.max(t), 1000)
            mod = pmag(tt, 10**f0, alpha=alpha, t0=t_0)

            # =============================================================================
            # scipy model
            # =============================================================================
            a, b, r_value, p_value, std = stats.linregress(t, y)
            scipy_model = scipy_linear_model(t, a, b)

            # =============================================================================
            # scikitlearn linear model regeression
            # =============================================================================
            sklearn_model, sk_score = sklearn_linear_model(t, y)

            # =============================================================================
            # sufficient statistics to get the likelihood maximun, bayesian approach 
            # =============================================================================
            M, v = sufficient_statistics(t, y, dy)
            x1, x2 = np.linalg.solve(M, v)
            max_likelihood = x1*t + x2

            # =============================================================================
            # plot the lightcurve and Regression models  
            # =============================================================================
            plotName = "%s/lightcurve.pdf" % (plotDir)
            plt.figure(filter)
            plt.errorbar(t-tmin, y, dy, fmt='o', c='k')
            plt.plot(tt-tmin, mod, 'r-', linewidth=2,
                        label="pymultinest model")
            plt.plot(t - tmin, sklearn_model, c='brown', linewidth=2,
                        label="sklearn linear model pred = "+np.str(np.round(sk_score*100, 2))+"%")
            plt.plot(t - tmin, scipy_model, 'b-', linewidth=2,
                        label="scipy linear model pred = "+np.str(np.round(r_value*100, 2))+"%")
            plt.plot(t - tmin, max_likelihood, 'g-',
                        linewidth=2, label="Max likelihood")

            plt.xlabel(r'Time [days] [t0 = %.5f]' % tmin, fontsize=18)
            plt.ylabel(r'Magnitude', fontsize=18)
            plt.grid()
            plt.legend(loc='best')
            plt.gca().invert_yaxis()
            plt.savefig(plotName)
            plt.close()

        else:
            print("******************************************************")
            print("Attention, the length of Data in this filter:",
                    filter, "is only", len(tab))
            print("Please wait, we go on  for others filter ")

            time.sleep(5)

            print("*******************************************************")
