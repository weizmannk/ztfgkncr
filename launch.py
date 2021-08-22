import pandas as pd
import os
import kn_catcher as kn  # the kn_cather.py script 

# =============================================================================
# api_key for Fink call and recharge the local data 
# =============================================================================
api_key = 'http://134.158.75.151:24000/api/v1/objects'

# =============================================================================
# returns current working directory of a process
# =============================================================================
datapath =  os.getcwd()+"/" 

# =============================================================================
# Filters only the ZTF21ab..... folder in /KN-Catcher-ReadyforO4/ 
# with exactly 12 characters
# =============================================================================
ZTF_ID = [f.name for f in os.scandir(datapath) if f.is_dir() and len(f.name)==12]

# =============================================================================
# Parse command line to create folders for  output data and figures
# =============================================================================
opts = kn.parse_commandline()
baseplotDir = opts.plotDir
if not os.path.isdir(baseplotDir):
    os.makedirs(baseplotDir)
    
# =============================================================================
# read json file, 
# =============================================================================
for ztf_id  in ZTF_ID:
    print("Read the file",  ztf_id)
    
    # read update json file data from ZTF 
    pdf = kn.read_json(api_key, ztf_id)
        
    # ==========================================================================
    # This line discards the bad folders that have been loaded  in ZTF_ID
    # ==========================================================================
    if pdf.shape !=(0, 0):
        
        # ======================================================================
        # DataFrame from ZTF pipeline
        # ======================================================================
        firstdate_ztf, times_ztf, mags_ztf, magerrs_ztf = kn.data_frame(pdf)

        # ======================================================================
        # Labels of ZTF filters
        # ======================================================================
        ztf_filts = kn.ztf_labeling_filter(pdf)

        # ======================================================================
        # Data of Observation
        # ======================================================================
        usernames, filters, date_time,  delays_ztf, magerrs, mags, decs, ras, real_filters, flux = kn.data_observation(datapath, ztf_id, pdf)

        # =======================================================================
        # Make a custom pandas dataframe
        # =======================================================================
        #create NA for lack of usernames
        NA =  len(ztf_filts)*["N/A"]
        instrument = len(usernames)*[ztf_id]
        
        ZTF = len(NA)*["ZTF"]
        
        # ========================================================================
        # Selection of data provide by users and ZTF online
        # ========================================================================
        d = {'time': times_ztf.values.tolist() + delays_ztf,
            'mag': mags_ztf.values.tolist() + mags,
            'magerr': magerrs_ztf.values.tolist() + magerrs,
            'filters': ztf_filts + filters,
            'username' : NA + usernames,
            'instrument' : ZTF + instrument 
            }
        
        # ========================================================================
        # Make a custom pandas dataframe 
        # Read data and save it in csv
        # ========================================================================
        df = pd.DataFrame(data=d) 
        df.sort_values("time", ascending=False)

        # ========================================================================
        # Create a csv file and save it under the ztf_id name
        # ========================================================================
        df.to_csv(ztf_id + '.csv')

        """"
        # ========================================================================
        # Create a csv file and save it in cvsDir directory
        # ========================================================================
        #csvDir = os.path.join(baseplotDir, "csv")
        #if not os.path.isdir(csvDir):
        #   os.makedirs(csvDir)  
        #df.to_csv(csvDir +"/"+ ztf_id + '.csv')
        """
        # =============================================================================
        #create a directory to save the lightcurve figures
        # =============================================================================
        plotDir = os.path.join(baseplotDir, "lightcurve")
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        # =============================================================================
        #Plot the lightcurbe 
        # =============================================================================
        kn.plot(ztf_id,  df, plotDir, mags_ztf,  figsize=(20, 12), font_plot=18)

# =============================================================================
# Run the pwl_catcher.py to get the Regressions for each filter 
# =============================================================================
exec(open('pwl_catcher.py').read())
