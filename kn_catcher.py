import matplotlib.pyplot as plt
import numpy as np
from math import isnan
from astropy.table import Table
from astropy.time import Time
import requests
import pandas as pd
import glob
import optparse

def data_upload(api_key, ztf_id):
    """[Get data for ZTF fields and Fink science module outputs.]

    :param api_key: API key for Fink call 
    :type API_key: [str ]
    :param ztf_id:  each ZTF21ab..... folder name in the KN-Catcher-ReadyforO4
    :type ztf_id: [ str ]
    :return: Get only a subset data of the fields provide by each  Grandma KN-Catcher, ZTF_ID
    :rtype: [ json format ]
    """
    r = requests.post(api_key,
                json={
                    'objectId': ztf_id, 
                    'output-format': 'json'
                    }
                    )
    return r
    
def read_json(api_key, ztf_id):
    """[reading json file]
    
    :return: Format output in a DataFrame
    :rtype: [pandas.core.frame.DataFrame]
    """
    return pd.read_json(data_upload(api_key, ztf_id).content)

def data_frame(pdf):
    """[DataFrame Extracting ]

    :param pdf: json file read by pandas in read_json () function 
    :type pdf: [ DataFrame]
    :return:  extracting the time and magnitude data of the observations
    :rtype: [tuple]
    """
    firstdate_ztf = pdf['i:jd'].values[-1]
    times_ztf = pdf['i:jd'] - firstdate_ztf
    mags_ztf = pdf['i:magpsf']
    magerrs_ztf = pdf['i:sigmapsf']
    
    return firstdate_ztf, times_ztf, mags_ztf, magerrs_ztf

def ztf_labeling_filter(pdf):
    """[Labels of ZTF filters]

    :param pdf: json file read by pandas in read_json () function
    :type pdf: [DataFrame]
    :return:the data of g and r filters
    :rtype: [list]
    """
    
    filtdic = {1: 'g', 2: 'r'}
    ztf_filts = []
    for filt in pdf['i:fid']:
        ztf_filts.append(filtdic[filt])

    return ztf_filts

def data_file_name(datapath, ztf_id):
    """[using recursively to find files]

    :param datapath: the filename directory
    :type datapath: [str]
    :return: all target.vot files 
    :rtype: [list]
    """
    return glob.glob(datapath+ztf_id+'/**/*.target.vot', recursive=True)

def data_observation(datapath, ztf_id, pdf):
    """[observation data extraction]
    :param datapath: the filename directory
    :type datapath: [str]
    :param ztf_id: each ZTF21ab..... folder name
    :type ztf_id: [str]
    :param pdf: json file read by pandas
    :type pdf: [DataFrame]
    :return: usernames, filters, date_time,  delays_ztf, magerrs, mags, decs, ras, real_filters, flux)
    :rtype: [tuple]
    """

    flux = []
    ras = []
    decs = []
    mags = []
    magerrs = []
    filters = []
    usernames = []
    delays_ztf = []
    real_filters = []
    date_time = []
    for filename in data_file_name(datapath, ztf_id):

        data = Table.read(filename)
        username = filename.split("_")[2]
        date = filename.split('_')[3]
        date_obs = Time(date.split('T')[0]+' ' +
                        date.split('T')[1].replace('-', ':'))
        real_filters.append(filename.split('_')[4])
        ras.append(float(data['ra']))
        decs.append(float(data['dec']))
        if isnan(float(data['mag_calib'])):
            mags.append(float(data['mag_limit']))
            magerrs.append(0)
        else:
            mags.append(float(data['mag_calib']))
            magerrs.append(float(data['magerr']))
        delays_ztf.append(date_obs.jd - data_frame(pdf)[0])
        filters.append(data['mag_filter_name'][0])
        usernames.append(username)
        date_time.append(
            (date.split('T')[0]+'T' + date.split('T')[1].replace('-', ':')))
        flux.append(float(data['flux']))
    return usernames, filters, date_time,  delays_ztf, magerrs, mags, decs, ras, real_filters, flux

def parse_commandline(): 
    """[Create the output files to save the processed data ]

    :return: [create an output folder]
    :rtype: [optparse.Values]
    """

    parser = optparse.OptionParser()
    parser.add_option("-p", "--plotDir", default= "outputs")
    opts, args = parser.parse_args()

    return opts 

def plot(ztf_id,  df, plotDir, mags_ztf,  figsize=(17, 15), font_plot=24):
    """[Plot ligthcurve]

    :param ztf_id: [figure name]
    :type ztf_id: [str]
    :param df: [selection of data provide by users and ZTF online ]
    :type df: [pandas.core.frame.DataFrame]
    :param plotDir: [plot Directory for saving figure]
    :type plotDir: [str]
    :param mags_ztf: [observation ype(magnitude]
    :type mags_ztf: [pandas.core.series.Series]
    :param figsize: [figure size], defaults to (17, 15)
    :type figsize: tuple, optional
    :param font_plot: [lengend font], defaults to 20
    :type font_plot: int, optional
    :return: [a plot of lightcure  of each ZTF21ab....... foders ]
    :rtype: [ plt.close()]
    """
    
    plt.figure(ztf_id, figsize)

    for filt in np.unique(df.filters.values):
        
        df_meas = df.loc[(df['filters'] == filt) &
                    (df['magerr'] > 0)].sort_values("time", ascending=False)
        df_ul = df.loc[(df['filters'] == filt) &
                (df['magerr'] == 0)].sort_values("time", ascending=False)

        if filt == 'B':
            col = 'b'
            plt.errorbar(df_meas.time, df_meas.mag+0.4,
                            yerr=df_meas.magerr, fmt='o-', label=filt+'+0.4', c=col)
            plt.scatter(df_ul.time, df_ul.mag+0.4, marker='v', c=col)
            col = 'g'
            plt.errorbar(df_meas.time, df_meas.mag+0.2,
                            yerr=df_meas.magerr, fmt='o-', label=filt+'+0.2', c=col)
            plt.scatter(df_ul.time, df_ul.mag+0.2, marker='v', c=col)

        elif filt in ['gmag']:
            col = 'lightgreen'
            plt.errorbar(df_meas.time, df_meas.mag+0.2,
                            yerr=df_meas.magerr, fmt='o-', label=filt+'+0.2', c=col)
            plt.scatter(df_ul.time, df_ul.mag+0.2, marker='v', c=col)

        elif filt == 'R':
            col = 'r'
            plt.errorbar(df_meas.time, df_meas.mag,
                            yerr=df_meas.magerr, fmt='o-', label=filt, c=col)

        elif filt == 'rmag':
            col = 'lightsalmon'
            plt.errorbar(df_meas.time, df_meas.mag,
                            yerr=df_meas.magerr, fmt='o-', label=filt, c=col)
            plt.scatter(df_ul.time, df_ul.mag, marker='v', c=col)

        elif filt == 'I':
            col = 'orange'
            plt.errorbar(df_meas.time, df_meas.mag-0.2,
                            yerr=df_meas.magerr, fmt='o-', label=filt+'-0.2', c=col)
            plt.scatter(df_ul.time, df_ul.mag-0.2, marker='v', c=col)

    plt.errorbar(df[0:len(mags_ztf-1)].loc[df['filters'] == 'g'].time,
                    df[0:len(mags_ztf-1)].loc[df['filters'] == 'g'].mag+0.2,
                    df[0:len(mags_ztf-1)].loc[df['filters'] == 'g'].magerr,
                    fmt='d', label='ZTF g+0.2', c='darkgreen')

    plt.errorbar(df[0:len(mags_ztf-1)].loc[df['filters'] == 'r'].time,
                    df[0:len(mags_ztf-1)].loc[df['filters'] == 'r'].mag,
                    df[0:len(mags_ztf-1)].loc[df['filters'] == 'r'].magerr,
                    fmt='d', label='ZTF r', c='darkred')

    # Figure plot settings
    plt.gca().invert_yaxis()
    plt.legend(loc='best')
    plt.grid(True)
    plt.xlabel(r'Time since first ZTF detection [days]', fontsize=font_plot)
    plt.ylabel(r'absolute magnitude (AB)', fontsize=font_plot)
    plt.xticks(fontsize=font_plot)
    plt.yticks(fontsize=font_plot)
    plt.title('Kilonova-Catcher follow-up of ' + ztf_id, fontsize=font_plot)
    plt.savefig(plotDir+"/"+ztf_id+"_KNC_lc.pdf")

    return plt.close()
