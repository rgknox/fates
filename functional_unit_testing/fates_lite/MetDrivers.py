import pandas as pd
import numpy as np
import code  # For development: code.interact(local=dict(globals(), **locals()))
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime

eval_cosz = True

def qsat(tempk,rh100):

    # tempk: Air temperature in Kelvin
    # rh100: relative humidity in percent
    # Reference:  Polynomial approximations from:
    #             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
    #             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
    
    tfrz = 273.14
    
    # For water vapor (temperature range 0C-100C)
    a0 =  6.11213476
    a1 =  0.444007856
    a2 =  0.143064234e-01
    a3 =  0.264461437e-03
    a4 =  0.305903558e-05
    a5 =  0.196237241e-07
    a6 =  0.892344772e-10
    a7 = -0.373208410e-12
    a8 =  0.209339997e-15
    # For derivative:water vapor
    b0 =  0.444017302
    b1 =  0.286064092e-01
    b2 =  0.794683137e-03
    b3 =  0.121211669e-04
    b4 =  0.103354611e-06
    b5 =  0.404125005e-09
    b6 = -0.788037859e-12
    b7 = -0.114596802e-13
    b8 =  0.381294516e-16
    # For ice (temperature range -75C-0C)
    c0 =  6.11123516
    c1 =  0.503109514
    c2 =  0.188369801e-01
    c3 =  0.420547422e-03
    c4 =  0.614396778e-05
    c5 =  0.602780717e-07
    c6 =  0.387940929e-09
    c7 =  0.149436277e-11
    c8 =  0.262655803e-14
    # For derivative:ice
    d0 =  0.503277922
    d1 =  0.377289173e-01
    d2 =  0.126801703e-02
    d3 =  0.249468427e-04
    d4 =  0.313703411e-06
    d5 =  0.257180651e-08
    d6 =  0.133268878e-10
    d7 =  0.394116744e-13
    d8 =  0.498070196e-16

    td = np.min([100.0, np.max([-75.0, tempk - tfrz] )])

    if (td >= 0.0):
        vpress_sat = 100.*(a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 + td*(a5 + td*(a6 + td*(a7 + td*a8))))))))
    else:
        vpress_sat = 100.*(c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 + td*(c5 + td*(c6 + td*(c7 + td*c8))))))))

    vpress = rh100*0.01*vpress_sat
    
    return vpress,vpress_sat

    
    

def GetCosZ(ZONE,JULIAN,LAT,LON,LHOUR):

    #SUBROUTINE COSZENITH (LON,LATD,LHOUR,ZONE,JULIAN,CZENITH)
    #
    #     The purpose is to calculate the following:
    #        1)  Day angle (GAMMA)
    #        2)  Solar DEClination
    #        3)  Equation of time
    #        4)  Local apparent time
    #        5)  Hour angle
    #        6)  Cosine of zenith angle
    #
    #     All equations come from "An Introduction to
    #     Solar Radition" By Muhammad Iqbal, 1983.
    #
    # Integers:
    #     ZONE,                ! time zone (1-24) GMT=12
    #     JULIAN               ! julian day
    # Real
    #     CZENITH,             ! cosine of zenith angle (radians)
    #     DECd,                ! solar declination (degrees)
    #     DEC,                 ! solar declination (radians)
    #     et,                  ! equation of time (minutes)
    #     GAMMA,               ! day angle (radians)
    #     LATime,              ! local apparent time
    #     LCORR,               ! longitudical correction
    #     LHOUR,               ! local standard time
    #     LON,                 ! local longitude (deg)
    #     LLAT,                ! local latitude in radians
    #     LATD ,               ! local latitude in degrees
    #     LS,                  ! standard longitude (deg)
    #     OMEGAD,              ! omega in degrees
    #     OMEGA ,              ! omega in radians
    #     PI,                  ! universal constant PI [-]
    #     ZENITH,              ! zenith angle(radians)
    #     ZEND                 ! zenith angle(degress)
    #
    #     Neither ZENITH nor ZEND are necessary for this program.
    #     I originally used them as checks, and left them here in
    #     case anyone else had a use for them.
    #
    #     1)  Day angle GAMMA (radians) page 3

    PI= 3.141592              # universal constant PI
    GAMMA=2.*PI*(JULIAN-1)/365.

    #     2) Solar declination (assumed constant for a 24 hour period)  page 7
    #     in radians
    #
    DEC=(0.006918-0.399912*np.cos(GAMMA)+0.070257*np.sin(GAMMA) \
         -0.006758*np.cos(2.*GAMMA)+0.000907*np.sin(2.*GAMMA) \
         -0.002697*np.cos(3.*GAMMA)+0.00148*np.sin(3.*GAMMA))
    DECd=DEC*(180./PI)
    #
    #     maximum error 0.0006 rad (<3'), leads to error of less than 1/2 degree
    #     in ZENITH angle
    #     3)  Equation of time  page 11

    et=(0.000075+0.001868*np.cos(GAMMA)-0.032077*np.sin(GAMMA) \
        -0.014615*np.cos(2*GAMMA)-0.04089*np.sin(2*GAMMA))*229.18
    #
    #     4) Local apparent time  page 13
    #
    #     LS     standard longitude (nearest 15 degree meridian)
    #     LON     local longitude
    #     LHOUR  local standard time
    #     LATIME local apparent time
    #     LCORR  longitudunal correction (minutes)
    #
    
    LS=((ZONE-1)*15)-180.
    LCORR=4.*(LS-LON)*(-1.)
    LATIME=LHOUR+LCORR/60.+et/60.
    if (LATIME<0.):
        LATIME=LATIME+24.
    if (LATIME>24.):
        LATIME=LATIME-24.

    #     5) Hour angle OMEGA  page 15
    #
    #     hour angle is zero at noon, postive in the morning
    #     It ranges from 180 to -180
    #
    #     OMEGAD is OMEGA in degrees, OMEGA is in radians
    #
    OMEGAD=(LATIME-12.)*(-15.)
    OMEGA=OMEGAD*PI/180.
    #
    #     6)  Zenith angle page 15
    #
    #     CZENITH cosine of zenith angle (radians)
    #     LAT=local latitude in degrees
    #     LLAT=local latitude in radians
    #
    LLAT=LAT*PI/180.
    CZENITH=np.sin(DEC)*np.sin(LLAT)+np.cos(DEC)*np.cos(LLAT)*np.cos(OMEGA)
    CZENITH=np.max([0.,CZENITH])
    ZENITH=np.arcsin(CZENITH)
    ZEND=180.*ZENITH/PI
    
    return CZENITH


class met_driver:

    def EvalCosz(self):

        if(eval_cosz):
            fig02,(ax1) = plt.subplots(1,1,figsize=(5.5,5.5))
            ax1.scatter(self.data['hod'],self.data['cosz'])
            ax1.grid('on')
            ax1.set_xlabel('Hour of Day')
            ax1.set_ylabel('Cos(Z)')
            plt.show()
    
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

    def GetVpressSat(self,tcan,rh,can_press):


        vpress_sat 


        vpress = rh*0.01*vpress_sat
        
        return vpress,vpress_sat
    
    # =============================================================================================
    
    def ModelSWComponents(self,sswdn):

        # This method takes total downwelling short-wave radiation
        # at the land-surface reference point [W/m2], and partitions
        # it into the two bands VIS and NIR, as well as the diffuse and
        # beam components of both, resulting in four different fluxes
        # that should add up to the total provided
        # This method was borrowed from CLM's DATM model:
        # components/cpl7/components/data_comps_mct/datm/src/datm_comp_mod.F90

        # sswdn (in) total downwelling short-wave radiation at reference [w/m2]
        # swvdr (out) direct (beam) downwelling visible radiation at ref [w/m2]
        # swndr (out) direct (beam) downwelling NIR radiation at ref [w/m2]
        # swvdf (out) diffuse downwelling visible radiation at ref [w/m2]
        # swndf (out0 diffuse downwelling NIR radiation at ref [w/m2]

        visfrac = 0.4  # Mean earth visible fraction of sw
        
        
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

            # NIR
            swndr = sswdn * (1.-visfrac)
            ratio_rvrf =  np.min([0.99,np.max([0.29548 + 0.00504*swndr -1.4957e-05*swndr**2 + 1.4881e-08*swndr**3,0.01])])
            swndr = ratio_rvrf*swndr
            swndf = sswdn * 0.50
            swndf = (1. - ratio_rvrf)*swndf

            # VIS
            swvdr = sswdn * visfrac
            ratio_rvrf =   np.min([0.99,np.max([0.17639 + 0.00380*swvdr - 9.0039e-06*swvdr**2 + 8.1351e-09*swvdr**3,0.01])])
            swvdr = ratio_rvrf*swvdr
            swvdf = sswdn * 0.50
            swvdf = (1. - ratio_rvrf)*swvdf

            return swvdr,swndr,swvdf,swndf
        
        
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

        eval_met_forcing = False
        
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

            # CSV HEADERS:
            # Date_UTC_start,Date_UTC_end,Date_local_start,Date_local_end,SR_W_m2.,
            #                SR_flag,Temp_o_C.,T_Flag,RH_%,RH_flag,RA_mm_d,RA_flag,
            #                BP_hPa,BP_flag,WS_m_s,WS_flag,WD_deg,WD_flag
            df = pd.read_csv(filepath, delimiter=",")
            self.data = {}
            self.data['can_press'] = np.array(df['BP_hPa'].values*100.)
            npts = len(self.data['can_press'])
            self.data['t_can']      = np.array(df['Temp_o_C.'].values+273.14)
            self.data['Rtot']       = np.array(df['SR_W_m2.'].values)
            self.data['u_ref']      = np.array(df['WS_m_s'].values)
            self.data['yr']         = np.zeros(npts,dtype=int)
            self.data['mon']        = np.zeros(npts,dtype=int)
            self.data['day']        = np.zeros(npts,dtype=int)
            self.data['hod']        = np.zeros(npts,dtype=int)
            self.data['r_b']        = np.zeros(npts)
            self.data['t_veg']      = np.zeros(npts)
            self.data['visbdn']     = np.zeros(npts)
            self.data['visddn']     = np.zeros(npts)
            self.data['nirbdn']     = np.zeros(npts)
            self.data['nirddn']     = np.zeros(npts)
            self.data['cosz']       = np.zeros(npts)
            self.data['vpress']     = np.zeros(npts)
            self.data['vpress_sat'] = np.zeros(npts)
            for idate in range(npts):
                yr,mo,dy,hod = self.ProcDateStr(df['Date_local_start'].values[idate], \
                                                df['Date_local_end'].values[idate])
                self.data['yr'][idate] = yr
                self.data['mon'][idate] = mo
                self.data['day'][idate] = dy
                self.data['hod'][idate] = hod
                self.data['visbdn'][idate],self.data['nirbdn'][idate], \
                    self.data['visddn'][idate],self.data['nirddn'][idate] = self.ModelSWComponents(self.data['Rtot'][idate])
                self.data['r_b'][idate],self.data['t_veg'][idate] = self.ModelRb(self.data['t_can'][idate], \
                                                                                 self.data['Rtot'][idate],  \
                                                                                 self.data['hod'][idate],   \
                                                                                 self.data['u_ref'][idate])
                self.data['vpress'][idate],self.data['vpress_sat'][idate] = qsat(self.data['t_can'][idate],df['RH_%'].values[idate])

                
            if (eval_met_forcing):

                # Radiation Scattering
                fig01,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(6.5,6.0))

                hods = np.linspace(0,23,24)

                hod_counts = np.zeros(len(hods))
                visbdn24   = np.zeros(len(hods))
                visddn24   = np.zeros(len(hods))
                nirbdn24   = np.zeros(len(hods))
                nirddn24   = np.zeros(len(hods))
                for ih,hod in enumerate(self.data['hod']):
                    ihod = np.argmin(np.abs(hods-hod))
                    hod_counts[ihod] = hod_counts[ihod] + 1.
                    visbdn24[ihod] = visbdn24[ihod] + self.data['visbdn'][ih]
                    visddn24[ihod] = visddn24[ihod] + self.data['visddn'][ih]
                    nirbdn24[ihod] = nirbdn24[ihod] + self.data['nirbdn'][ih]
                    nirddn24[ihod] = nirddn24[ihod] + self.data['nirddn'][ih]
                    
                    
                ax1.scatter(hods,visbdn24/hod_counts)
                ax1.grid('on')
                ax1.set_xlabel('Hour')
                ax1.set_ylabel('VIS Beam [W/m2]')
                ax2.scatter(hods,visddn24/hod_counts)
                ax2.grid('on')
                ax2.set_xlabel('Hour')
                ax2.set_ylabel('VIS Diff [W/m2]')
                ax3.scatter(hods,nirbdn24/hod_counts)
                ax3.grid('on')
                ax3.set_xlabel('Hour')
                ax3.set_ylabel('NIR Beam [W/m2]')
                ax4.scatter(hods,nirddn24/hod_counts)
                ax4.grid('on')
                ax4.set_xlabel('Hour')
                ax4.set_ylabel('NIR Diff [W/m2]')

                plt.show()
        
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
