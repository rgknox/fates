import pandas as pd
import numpy as np
import code  # For development: code.interact(local=dict(globals(), **locals()))
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime

class met_driver:

    def TrainModelRb(self):

        # A = np.array([[1, 2], [3, 4], [5, 6]])
        # b = np.array([3, 7, 11])
        # code.interact(local=dict(globals(), **locals()))
        rb = self.data['r_b']
        rA = np.array([self.data['forc_t'],self.data['forc_t']**2.0,self.data['Rtot'], \
                      self.data['hod'], self.data['u_ref'],self.data['u_ref']**2.0, \
                      self.data['u_ref']*self.data['forc_t']]).transpose()
        rx, rresiduals, rrank, rs = np.linalg.lstsq(rA, rb, rcond=None)
        rb_mod = np.dot(rA,rx)

        tb = self.data['t_veg']
        tA = np.array([self.data['forc_t'],self.data['forc_t']**2.0,self.data['Rtot'], \
                      self.data['hod'], self.data['u_ref'],self.data['u_ref']**2.0, \
                      self.data['u_ref']*self.data['forc_t']]).transpose()
        tx, tresiduals, trank, ts = np.linalg.lstsq(tA, tb, rcond=None)
        tb_mod = np.dot(tA,tx)

        
        fig01,(ax1,ax2) = plt.subplots(1,2,figsize=(9.5,5.5))

        ax1.scatter(rb,rb_mod)
        ax1.set_title('Multiple Regression for\n Boundary Layer Resistance')
        ax1.grid('on')
        ax1.axis('equal')
        ax1.set_xlim([0,80])
        ax1.set_ylim([0,80])
        ax1.set_xlabel('r_b [s/m] (ELM)')
        ax1.set_ylabel('r_b [s/m] (Regression)')

        ax2.scatter(tb-273.14,tb_mod-273.14)
        ax2.set_title('Multiple Regression for\n Leaf Temperature')
        ax2.grid('on')
        ax2.axis('equal')
        ax2.set_xlim([15,45])
        ax2.set_ylim([15,45])
        ax2.set_xlabel('T_v [C] (ELM)')
        ax2.set_ylabel('T_v [C] (Regression)')
        plt.show()

        print("weights (r_b):",rx)
        print("weights (tveg):",tx)
        
        
    def ModelRb(self,tforce,rtot,hod,u_ref):

        # Rb is in s/m
        # Tveg in K
        # hod is 24 hour cycle
        # u_ref = m/s

        rb_x = [ 4.74117599e+00, -1.48710379e-02, -6.09059694e-02, -8.81964792e-01, \
                 -7.46423948e+02,  1.56705310e+00,  2.43424162e+00]

        tveg_x = [ 1.02542670e+00, -9.37678968e-05,  8.46295812e-03,  3.88044372e-02, \
                   -6.39628996e+00,  5.25790309e-02,  2.02054394e-02]

        
        A = np.array([tforce,tforce*tforce,rtot,hod,u_ref,u_ref*u_ref,u_ref*tforce])
        

        rb = np.dot(A,rb_x)

        tveg = np.dot(A,tveg_x)

        
        return rb,tveg

    def ModelSWComponents(sswdn):

        use_simple_model = False
        
        # Model one, simple proportions
        if use_simple_model:
            swvdr = sswdn*0.28
            swndr = sswdn*0.31
            swvdf = sswdn*0.24
            swndf = sswdn*0.17

        else:
            # Model twwo
            # relationship between incoming NIR or VIS radiation and ratio of
            # direct to diffuse radiation calculated based on one year's worth of
            # hourly CAM output from CAM version cam3_5_55

            swndr = sswdn * 0.50
            
            ratio_rvrf =  np.min([0.99,np.max([0.29548 + 0.00504*swndr -1.4957e-05*swndr**2 + 1.4881e-08*swndr**3,0.01])])

            swndr = ratio_rvrf*swndr
            swndf = (1. - ratio_rvrf)*sswdn * 0.50

             swvdr = avstrm%rAttr(sswdn,n) * 0.50_R8
             ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr  &
                  -9.0039e-06_R8*swvdr**2 + 8.1351e-09_R8*swvdr**3,0.01_R8))
             a2x%rAttr(kswvdr,n) = ratio_rvrf*swvdr
             swvdf = avstrm%rAttr(sswdn,n) * 0.50_R8
             a2x%rAttr(kswvdf,n) = (1._R8 - ratio_rvrf)*swvdf
        
        
    def __init__(self,filepath):

        # This routine simply reads a comma delimited text file.
        # It should have no header, and should contain
        # the following 15 columns:
        #
        # yr: year, integer, 4 digits
        # mon: month, integer
        # day: day of month, integer
        # tod: time of day in seconds
        # r_b: boundary layer resistance [s/m]
        # u_ref: wind velocity at reference height (above canopy) [m/s]
        # forc_q: specific humidity above the canopy
        # vpress: vapor pressure in the canopy [Pa]
        # t_veg:  vegetation temperature [K]
        # can_press: air pressure in the canopy [Pa]
        # forc_t: air temperature above canopy [K]
        # thm: mean temperature between air and vegetation [K] (not used)
        # t_can: canopy air temperature [K]
        # Rbeam: downwelling visible beam radiation [W/m2]
        # Rdiff: downwelling visible diffuse radiation [W/m2]
        
        # Not all of these are used. We also convert tod to "hour"od,
        # which is decimal hours

        model_read = False
        
        if(model_read):
            df = pd.read_csv(filepath, delimiter=",", header=None)
            self.data = {}
            self.data['yr'] = np.array(df[0].values)
            self.data['mon'] = np.array(df[1].values)
            self.data['day'] = np.array(df[2].values)
            self.data['hod'] = np.array(df[3].values/3600.0)
            self.data['r_b'] = np.array(df[4].values)
            self.data['vpress'] = np.array(df[7].values)
            self.data['t_veg'] = np.array(df[8].values)
            self.data['can_press'] = np.array(df[9].values)
            self.data['t_can'] = np.array(df[12].values)
            self.data['Rbeam'] = np.array(df[13].values)
            self.data['Rdiff'] = np.array(df[14].values)
            self.data['Rtot'] = self.data['Rbeam'] + self.data['Rdiff']
            self.data['u_ref'] = np.array(df[5].values)
            self.data['forc_t'] = np.array(df[11].values)
            self.ndata = len(self.data['yr'])
            
            #self.ModelRb()
            self.TrainModelRb()

        else:

            # Date_UTC_start,Date_UTC_end,Date_local_start,Date_local_end,SR_W_m2.,SR_flag,Temp_o_C.,T_Flag,RH_%,RH_flag,RA_mm_d,RA_flag,BP_hPa,BP_flag,WS_m_s,WS_flag,WD_deg,WD_flag
            df = pd.read_csv(filepath, delimiter=",")
            self.data = {}
            self.data['can_press'] = np.array(df['BP_hPa'].values*100.)
            npts = len(self.data['can_press'])
            self.data['t_can'] = np.array(df['Temp_o_C.'].values+273.14)
            self.data['Rtot'] = np.array(df['SR_W_m2.'].values)
            self.data['u_ref'] = np.array(df['WS_m_s'].values)
            self.data['yr'] = np.zeros(npts)
            self.data['mon'] = np.zeros(npts)
            self.data['day'] = np.zeros(npts)
            self.data['hod'] = np.zeros(npts)
            self.data['r_b'] = np.zeros(npts)
            self.data['t_veg'] = np.zeros(npts)
            #self.data['Rbeam'] = np.array(df[13].values)
            #self.data['Rdiff'] = np.array(df[14].values)
            #self.data['Rtot'] = self.data['Rbeam'] + self.data['Rdiff']
            #self.data['forc_t'] = np.array(df[11].values)
            for idate in range(npts):
                yr,mo,dy,hod = self.ProcDateStr(df['Date_local_start'].values[idate],df['Date_local_end'].values[idate])
                self.data['yr'][idate] = yr
                self.data['mon'][idate] = mo
                self.data['day'][idate] = dy
                self.data['hod'][idate] = hod
                self.data['r_b'][idate],self.data['t_veg'][idate] = self.ModelRb(self.data['t_can'][idate], \
                                                                                 self.data['Rtot'][idate],  \
                                                                                 self.data['hod'][idate],   \
                                                                                 self.data['u_ref'][idate])
                
            code.interact(local=dict(globals(), **locals()))
        
    def ProcDateStr(self,date_local_start,date_local_end):


        ts1 = pd.Timestamp(date_local_start)
        ts2 = pd.Timestamp(date_local_end)
        ts = ts1 + (ts2-ts1)/2
        yr = ts.year
        mo = ts.month
        dy = ts.day
        hr = ts.hour
        mn = ts.minute

        hod = hr + mn/60.0
        
        return yr, mo, dy, hod
        
        
    def FilterTimes(self,study_period):

        # For Rey-Sanchez, PNM Panama:
        # Create masks for wet season, dry season, morning and afternoon
        #
        # From: https://gml.noaa.gov/grad/solcalc/
        # Solar noon Nov 1 2011:  17:02:57 GMT
        #            Dec 31 2011: 17:22:06 GMT
        #            Feb 1 2013:  17:32:58 GMT
        #            Mar 31 2013: 17:23:33 GMT
        # Sunrise: 11:17
        
        # bfilter is short for binary-filter
    
        if(study_period == 'reysanchez_wetssn_morning'):
            morning = [ (self.data['hod'][iyr] <= 17.25 and self.data['hod'][iyr] > 8) and (self.data['Rtot'][iyr]>0.) for iyr,year in enumerate(self.data['yr'])]
            bfilter = [ (self.data['yr'][iyr] == 2011 and (self.data['mon'][iyr]==11 or self.data['mon'][iyr]==12)) and morning[iyr] for iyr,year in enumerate(self.data['yr'])]

        if(study_period == 'reysanchez_dryssn_morning'):
            morning = [ (self.data['hod'][iyr] <= 17.25 and self.data['hod'][iyr] > 8) and (self.data['Rtot'][iyr]>0.) for iyr,year in enumerate(self.data['yr'])]
            bfilter = [ (self.data['yr'][iyr] == 2013 and (self.data['mon'][iyr]==2 or self.data['mon'][iyr]==3)) and morning[iyr] for iyr,year in enumerate(self.data['yr'])]

        if(study_period == 'reysanchez_wetssn_afternoon'):
            afternoon = [ (self.data['hod'][iyr] > 17.25 or self.data['hod'][iyr] < 8) and (self.data['Rtot'][iyr]>0.) for iyr,year in enumerate(self.data['yr'])]
            bfilter =  [ (self.data['yr'][iyr] == 2011 and (self.data['mon'][iyr]==11 or self.data['mon'][iyr]==12)) and afternoon[iyr] for iyr,year in enumerate(self.data['yr'])]

        if(study_period == 'reysanchez_dryssn_afternoon'):
            afternoon = [ (self.data['hod'][iyr] > 17.25 or self.data['hod'][iyr] < 8) and (self.data['Rtot'][iyr]>0.) for iyr,year in enumerate(self.data['yr'])]
            bfilter = [ (self.data['yr'][iyr] == 2013 and (self.data['mon'][iyr]==2 or self.data['mon'][iyr]==3)) and afternoon[iyr] for iyr,year in enumerate(self.data['yr'])]
        
        if(study_period == 'reysanchez_wetssn'):
            bfilter = [ (self.data['yr'][iyr] == 2011 and (self.data['mon'][iyr]==11 or self.data['mon'][iyr]==12)) for iyr,year in enumerate(self.data['yr'])]

        if(study_period == 'reysanchez_dryssn'):
            bfilter  = [ (self.data['yr'][iyr] == 2013 and (self.data['mon'][iyr]==2 or self.data['mon'][iyr]==3)) for iyr,year in enumerate(self.data['yr'])]

        if(study_period == 'unfiltered'):
            bfilter = [ True for iyr,year in enumerate(self.data['yr'])]

        if(study_period == 'daytime'):
            bfilter = [ (self.data['Rtot'][iyr]>0.) for iyr,year in enumerate(self.data['yr'])]

            
        if(sum(bfilter)<1):
            print('The filtering of met data produced no datapoints')
            exit(2)

            
        # Loop through all met data entries and filter
        for key, val in self.data.items():
            self.data[key] = self.data[key][bfilter]
        
        self.ndata = len(self.data['yr'])
