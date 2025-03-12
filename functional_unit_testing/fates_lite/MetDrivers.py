import pandas as pd
import numpy as np
import code  # For development: code.interact(local=dict(globals(), **locals()))

class met_driver:

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
        
        self.ndata = len(self.data['yr'])

        
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


        if(sum(bfilter)<1):
            print('The filtering of met data produced no datapoints')
            exit(2)

            
        # Loop through all met data entries and filter
        for key, val in self.data.items():
            self.data[key] = self.data[key][bfilter]
        
        self.ndata = len(self.data['yr'])
