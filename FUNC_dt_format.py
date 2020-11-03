#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: Andres F. Zambrano Moreno, 2019
License: GPLv3

This file is part of SSAR_mag_corr.

    SSAR_mag_corr is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SSAR_mag_corr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SSAR_mag_corr.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np

def dt_format(dt_type,exponent,Ms_type_only,M_range_min,M_range_max,within_dt_and_TMD,CAT_shorthand):
    X_NEAREST_Dt=None
    time_threshold=None
    percent=None
    time_threshold_print=None
    dt_lower_end_percent_exponent=None
    Full_CAT_type=None
    
    if(dt_type=='x<Dt<y' or dt_type=='Dt<y'):        
        X_NEAREST_Dt =1 # number of rows to look for an AS after mother event, can be an integer or 'All'(to search whole catalog)
        time_threshold = 10**exponent # time in seconds (specify when choosing thresholding ='time' in line above)
        time_threshold = np.around(time_threshold,decimals=2)
        percent = 0.5 # percentage of upper time to use as lower time treshold when choosing dt_type='x<Dt<y'

        if(within_dt_and_TMD == True and M_range_min == None):
            if (Ms_type_only!=None):
                Full_CAT_type='%s_and_MD_%s_(%s)'%(dt_type,CAT_shorthand,Ms_type_only)
            else:
                Full_CAT_type='%s_and_MD_%s'%(dt_type,CAT_shorthand)
        elif(within_dt_and_TMD == True and M_range_min != None):
            if (Ms_type_only!=None):
                Full_CAT_type='%s_and_MD_M=%s-%s_%s_(%s)'%(dt_type,M_range_min,M_range_max,CAT_shorthand,Ms_type_only)
            else:
                Full_CAT_type='%s_and_MD_M=%s-%s_%s'%(dt_type,M_range_min,M_range_max,CAT_shorthand)
        else:
            Full_CAT_type='%s_%s'%(dt_type,CAT_shorthand)
            
		###shouldn't go to the next two lines. Added as precaution:
        dt_lower_end_percent = 'CHECK FUNC_dt_format.py'
        dt_lower_end_percent_exponent ='CHECK FUNC_dt_format.py' 
        ###
        
        if (dt_type == 'x<Dt<y'):
            time_threshold_print ='%.2f<dt<%.2f'%(dt_lower_end_percent_exponent,exponent) 
            print'\n'
            print'dt_lower_end_percent:\n ',dt_lower_end_percent
            print'-------------------------'
            print'TIME THRESHOLD VALUE', time_threshold_print
            print'----------------------------'
            print'\n'
        elif (dt_type == 'Dt<y'):
            time_threshold_print ='dt<%.2f'%(exponent)
            
    if(dt_type!='x<Dt<y' and dt_type!='Dt<y'):
        time_threshold=None
        time_threshold_print=None
    
    return(X_NEAREST_Dt,Full_CAT_type,time_threshold,time_threshold_print,percent,dt_lower_end_percent_exponent)
