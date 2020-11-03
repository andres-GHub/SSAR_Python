#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Andres F. Zambrano Moreno, 2019
License: GPLv3

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
######################################################################################
# References:[1] J. Davidsen and M. Baiesi 2015, "Self-similar Aftershock Rates" 
# Self-Similar model as a branching process
#-----------------------------------

# By default, python interprets any number that includes a decimal point as 
# a double precision floating point number
#-----------------------------------

#% Always check catalog format first. Reads files without headers and multiple spaces
# Catalog downloads:
#http://service.scedc.caltech.edu/ftp/catalogs/hauksson/Socal_DD/hs_1981_2011_catalog_v01.format
#http://scedc.caltech.edu/research-tools/downloads.html

###################################################
#to use when modules are not intalled in a system:
#import pip
#pip.main(["install","numpy"])
#pip.main(["install","pandas"])
#pip.main(["install","matplotlib"])
#pip.main(["install","scipy"])
#pip.main(["install","math"])
###################################################
import time
import os
import re
import sys
import random
import numpy as np
import matplotlib.pyplot as plt
plt.ioff() #turn off interactive plot

cwd = os.getcwd()

save_to_folder=cwd
load_from_folder='%s/data/SSM/'%cwd

from FUNC_mag_diff_SIMPLE import mag_diff_simple
from FUNC_mag_diff_rand_SIMPLE import mag_diff_rand_simple
from save_load_dP_CDFs import save_dP_arr

from FUNC_dt_format import dt_format
from FUNC_file_name import file_name
from FUNC_load_CAT_type import load_CAT_type_1st
from FUNC_load_CAT_type import load_CAT_type_2nd
from scipy.stats import itemfreq

start_time = time.time()

import datetime
today = datetime.date.today()
today.strftime('%m-%d-%Y')
#################################################
# THIS IS ONLY FORE REFERENCE. WHEN GENERATING CATALOGS THESE ARE THE RELEVANT PARAMETERS FOR THE MODEL
lam       = 200000
T         = 100000.     #PeriodSC_CAT_analysis==True:
Mmin      = 1.0
Mmax      = 7.3
p         = 1.15
c_0       = 210.
tau_0     = 8000. #Value of 10^4 from Ref.[1]
b         = 1.08    #background b value
g         = 0.66
z         = 0.24

b_as = g+z
alpha = z+(p*g)

print'b (background):%s ---- b_as (g+z):%s ---- alpha(z+(p*g)):%s'%(b,b_as,alpha)
#################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINING FUNCTIONS:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             
# Random number generator takes values between 0 and 1:
#   Python uses the Mersenne Twister(a determimistic PRNG) as the core generator.
#    It produces 53-bit precision floats and has a period of 2**19937-1.
#       MT is used often  in Monte-Carlo simulations
def RAND():
    r=random.uniform(0.,1.)
    return r

global CAT_Arr

CAT_Arr=['20181120-115204_FULL_2.0e+05*10^5Back_tau10001_c210_81per']
#Name for the folder:
CAT_shorthand='SSAR-SC'

#Getting name for type of catalog:
split_CAT_Arr=re.split('_|-',CAT_Arr[0])
#print 'len CAT arr name', len(split_CAT_Arr)
if (len (split_CAT_Arr)==6 or len (split_CAT_Arr)==7):
    CAT_type_label=split_CAT_Arr[3]+'_'+split_CAT_Arr[4]
if (len (split_CAT_Arr)==4):
    CAT_type_label=split_CAT_Arr[2]
    
#############################################################################
#############################################################################
# IF USING Southern-California Catalog:
SC_CAT_analysis=False#False
#SC_CAT_analysis=True

#CAT_Arr=['HS17_moded_tm_seconds'] #for unmodified SC CAT (must download from Hauksson et al. reference)
##CAT_shorthand = ['SC_CAT']

#Getting name for type of catalog:
split_CAT_Arr=re.split('_|-',CAT_Arr[0])
print split_CAT_Arr
#print 'len CAT arr name', len(split_CAT_Arr)
if (len (split_CAT_Arr)==6 or len (split_CAT_Arr)==7):
    CAT_type_label=split_CAT_Arr[3]+'_'+split_CAT_Arr[4]
if (len (split_CAT_Arr)==4):
    CAT_type_label=split_CAT_Arr[2]

##########################################################
# Optional: the following is used to find the min. Dts for dm<X and dm>X (need to set histo_minima_dt_dml == True, below)
histo_minima_dt_dml = False
if(histo_minima_dt_dml == True):
    
    min_times_trigger_mag = '6.0'
    CAT_Arr=[]
    CAT_Arr = os.listdir('/home/andres/SSM/SSM/Code/Python/Feb27/data/SSM/1BG/CAT/min_times_%s_%.2EAS/'%(min_times_trigger_mag) )
    print 'LENGTH CAT for 1stAS: ',len(CAT_Arr)
#####################################################

CAT_Arr= list(CAT_Arr)
CATS = len(CAT_Arr)

min_dt_list = []
First_Load = 0

global m_as_1st_list
global percent
global image_counter
percent = 1.0 #initializing to 100%
m_as_1st_list = [0,0,0,0,0,0,0,0]
m_as_1st_list = np.array(m_as_1st_list)
m_as_1st_list = np.reshape(m_as_1st_list, (1,8) )
MAX_for_loop = 2#45000
reps = -1
X_NEAREST=[]
X_NEAREST_Dt=[]

##############################################################
def func_print_cond_type(conditions_list,exponent):
                cond_type_print=['subseq','within_dt','TMD','within_dt_and_TMD']#,'unrelat_subseq','subseq']
                for k in range(0,len(cond_type_print) ):
                    if conditions_list[k]==True:
                        print'\n'
                        print'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                        print'________________________________________________'
                        print 'CONDITIONING TYPE:',cond_type_print[k]
                        print'current time exponent:',exponent
                        print'\n'               
###################################################################      
               
min_range_exponent = 1.0 # value of exponent with base 10 (i.e. 10**(min_range_exponent) ) 
max_range_exponent = 2.8 # value of exponent with base 10 (i.e. 10**(max_range_exponent) )
steps_exponent = 0.2 # step increments in exponent
image_counter=1
for exponent in np.arange(min_range_exponent,max_range_exponent,steps_exponent):
    for divisor_max_mix in np.arange(1,MAX_for_loop,5000):

        for ii in xrange(0,CATS):
        
            catNumb=CAT_Arr[ii]
            CAT_NUMB=CAT_Arr[ii]
            
            TEST_conditioning=False
            
            All_mother_daught_1BG=False #Set to true to only consider only mother-daughter events of first generation
            
            single_tree = []
            FIRST_gen_only=False#True; USE Only with 1BG CATS. Refers to 1st Gen after the mainshock. True will find the next event which is only 1 generation away from the previous event. When using 1BG CAT with 1st gen only events, set this to FALSE
            max_mag_th = 3.6 # maximum magnitude threshold to consider
            Low_Mag_Thr = 1.6 # value of lowest magnitude to keep in catalog (which goes up to max_mag_th=4.0 below)
            N_shuffles = 500 # number of shuffles
            CAT_type = 'loadSSM' #SSM,loadSSM,SC
            Background = 'Full' # 1BG or Full. 1 BG refers to a single background event
            #Shuffle_type (randomizing) is either within list L1 or L2-> jo, for any type of conditioning
            # i.e. in Dt_MD conditioning will only shuffle within the, say, L2 list,
            # or within the full catalog -> lip, randomized mags will come from full catalog
            shuffle_type = 'lip' #'all','jo', 'alternating', 'time_thresh','lip' # time_thresh requires value for dt_threshold #lip uses full cat as input into mag_diff_rand_simple
            lip_consider_dt = False #UNUSED, 
            lip_consider_MD = False #REMOVE
            shuffle_mode = 'Replace' #No_Replace,Replace
            L_Shuffle = 'L2'#'L1' or 'L2' (see FUNC_mag_diff_rand_SIMPLE for details. This shuffles either subsequent events or preceding events) 
            
            print 'Shuffle type: ', shuffle_type
            print 'Number of shuffles: ', N_shuffles
            print 'Shuffle mode: ', shuffle_mode
            
            thresholding_addendum=None#'next_nearest_TrueMD' #'next_nearest_TrueMD' or None
            if (thresholding_addendum=='next_nearest_TrueMD'):
                X_NEAREST='first_subsequent'#'first_subsequent'#'All' or 'first_subsequent' or an INT-> No. of rows to look at after mother
            
            Mag_1BG=0
            if (Background=='1BG'):
                m = re.search('1BG_(.+?)_', CAT_Arr[0])
                if m == None:
                    print'\n'
                    print '***** ERROR: check type of CAT (is it FULL or 1BG?) *****\n'
                    sys.exit()  
                Mag_1BG = np.array( m.group(1) ) #MAGNITUDE of BG event for single mag CAT (reads from CAT_Arr)       
                
            check_trigg_mag = False # REMOVE IN CLEANUP  ,\Option to check waht the magnitude of the triggering event is.\
            Trigg_Mag = 00 # REMOVE IN CLEANUP # \Choose triggering magnitude of mother events which you'd like look at (00 for none)\
            histo_minima_dt_dml = False#SEP13-True #
            histo_m_as_1st = False#SEP13-True
                    
            thresholding = 'subsequent'
			
#			#CONDITIONING -------------------------------:
            subseq            = False # no conditioning
            within_dt         = True  # True or False # condition on events that fall within a certain \delta t value. Used in FUNC_mag_diff_SIMPLE.py
            TMD               = False # True or False #'TMD_cond' in function: FUNC_mag_diff_SIMPLE.py
            within_dt_and_TMD = False # True or False # condition on events that fall within a certain \delta t value and that are mother-daughter events
            MD_1st_gen        = False # SET to true in order to load a CAT in the cwd which has MD events else set to False; this will still look for 1st gen events
            unrelat_subseq    = False # True or False # unrelated subsequent events
            
            if SC_CAT_analysis==False:
                no_background     = False # True or False # Removes all background events starting in func mag_diff_simple
            else:
                no_background     = False
            
            conditions_list = [subseq,within_dt,TMD,within_dt_and_TMD]
            func_print_cond_type(conditions_list,exponent)
            
            consider_all_events           = False # True or False
            mothers_consider_bckgrnd_only = False # True or False
            mothers_consider_AS_only      = True  # True or False # CONSIDER ONLY MOTHER THAT ARE AFTERSHOCKS (ie, NO BACKGROUNDS)
            if(mothers_consider_bckgrnd_only==True):
                Ms_type_only='Ms_bckgrnd_only'
            elif(mothers_consider_AS_only==True):
                Ms_type_only='Ms_AS_only'
            else:
                Ms_type_only=None

            M_range_min = None #1.5 #None or a min mag value
            M_range_max = None #2.0 #None or a max mag value
            if (M_range_max<max_mag_th and M_range_max!=None):
                print 'INSIDE IF BEGINNING: M_range_max<max_mag_th and M_range_max!=None'
                max_mag_th = M_range_max
  
            dt_type ='Dt<y' # OPTIONS: 'x<Dt<y', 'Dt<y' or None
            
            X_NEAREST_Dt,Full_CAT_type,time_threshold,time_threshold_print,percent,dt_lower_end_percent_exponent=\
            dt_format(dt_type,exponent,Ms_type_only,M_range_min,M_range_max,within_dt_and_TMD,CAT_shorthand)
            
            if(TMD==True):
                dt_type=None
                exponent=max_range_exponent#-steps_exponent #so it only loops once
                if(M_range_min!=None and M_range_max!=None):
                    if (Ms_type_only!=None):
                        Full_CAT_type='MD_M=%s-%s_%s_(%s)'%(M_range_min,M_range_max,CAT_shorthand,Ms_type_only) #used when dt_type is None
                    else:
                        Full_CAT_type='MD_M=%s-%s_%s'%(M_range_min,M_range_max,CAT_shorthand)
                else:
                    if MD_1st_gen==True:
                        Full_CAT_type='MD_1st_gen_%s'%(CAT_shorthand)
                    else:
                        Full_CAT_type='MD_%s'%(CAT_shorthand)
                    
            if(subseq==True):
                 Full_CAT_type='uncond_%s'%CAT_shorthand
                 dt_type=None
            
            lower_time_threshold=120 #Xsec#None#percent*time_threshold
            
            inside_file=True#default:False
            File_Label,inside_file_name=file_name(Background,dt_type,X_NEAREST_Dt,percent,Full_CAT_type,time_threshold,reps,catNumb,shuffle_type,shuffle_mode,\
                                                  L_Shuffle,N_shuffles,thresholding,time_threshold_print,Trigg_Mag,single_tree,divisor_max_mix,\
                                                  All_mother_daught_1BG,Mag_1BG)
            
            
            def params_dP (compare=None,shuffle_type = shuffle_type,shuffle_mode=shuffle_mode,group_ID=None,catNumb=CAT_NUMB):
                
                # define ecdf function to be used in Mag_CORR() function:      
                def ecdf(x,norm):    
                    if norm == True:
                        xs = np.sort(x)
                        ys = ( np.arange(1, len(xs)+1) /len(xs) )
                    else:
                        xs = np.sort(x)
                        ys = (np.arange(1, len(xs)+1))
                    return xs, ys
                
                def MAG_CORR(mag_th,CAT_type,Background,shuffle_type,N_shuffles):
                    global catNumb
                    global First_Load
                    global image_counter
                    MAG_TH = mag_th
                    global TimeMagL,TimeMagLT,TimeMagL_RAND
                    #TimeMagL will be the CAT with mags above m_th
                    #TimeMagLT has the following columns(t_as,m_as,tree,gen,leaf,leaf_mom,global_ID,mother_ID)
                    TimeMagL=[]
                    TIME_count = 0
                    INIT = 0
                    First_Load=0
                    
                    for mag_th in np.arange(Low_Mag_Thr,max_mag_th,0.20):
                        mag_th = round(mag_th,2) 
                        print 'mag_th: ',mag_th
                        print '---------------------------------------\n'
                        if INIT==0:
                            #XXXX
                            INIT,First_Load,TimeMagLT,TimeMagL=load_CAT_type_1st(cwd,load_from_folder,INIT,CAT_type,catNumb,Background,\
                                                                  histo_minima_dt_dml,First_Load,mag_th,exponent,min_range_exponent,All_mother_daught_1BG)
                        else: 
                            INIT,First_Load,TimeMagL=load_CAT_type_2nd(cwd,load_from_folder,INIT,CAT_type,catNumb,Background,\
                                                                  histo_minima_dt_dml,First_Load,mag_th,exponent,min_range_exponent,All_mother_daught_1BG,TimeMagLT,TimeMagL)
               
                                                
                        print 'TimeMagLT (total CAT) SIZE -----------------: ', len(TimeMagLT)

                        global TimeMagL_Dir_Off_arr
                        global TimeMagL_Unfiltered
                        TimeMagL_Unfiltered=[]
                        global subseq_arr,within_dt_arr,MD_arr,within_dt_MD_arr
                        global mas_M_arr
                        
                        #converting NaNs (background) value of 'mother event' to -2
                        where_are_NaNs = np.isnan(TimeMagL)  
                        TimeMagL[where_are_NaNs]=-2
                        print 'TML Before any changes: ',np.shape(TimeMagL)
                        TML_for_ratio = TimeMagL
                        
                        # XXXX mag_diff function:                     
                        dml,subseq_arr,within_dt_arr,MD_arr,within_dt_MD_arr=\
                        mag_diff_simple(split_CAT_Arr[1],mag_th,TimeMagL,SC_CAT_analysis,shuffle_type,\
                                        time_threshold,dt_type,MD_1st_gen,\
                                        within_dt_and_TMD,TMD,within_dt,subseq,M_range_min,\
                                        M_range_max,mothers_consider_bckgrnd_only,mothers_consider_AS_only,no_background)
                        
                        within_dt_for_ratio = within_dt_arr
                        
                        if(within_dt_and_TMD == True):
                            mas_M_arr = within_dt_MD_arr

                        if(within_dt == True):
                            mas_M_arr = within_dt_arr

                        if(TMD == True):
                            mas_M_arr = MD_arr

                        if(subseq == True):
                            mas_M_arr=subseq_arr
                        print 'TML AFTER: ', np.shape(TimeMagL)
                        if len(mas_M_arr) == 0:
                            mas_M_arr = np.zeros((1,17))
                        print 'SHAPE mas_M_ARR:',np.shape(mas_M_arr)
                        print'\n'            
                                             
                        def cdfs(data1f):  
                            
                            global x1,y1,data1,cdf1T,cdf2T,cdf1T_norm,cdf2T_norm,L1,L2
                            # calculate, once, mag diff. for unshuffled catalog:  
                            if (histo_minima_dt_dml == True):
                                data1 = dml
                            elif(histo_minima_dt_dml == False):
                                global data2,TimeMagL_RAND
                                global mags_Dir_Off_culled_dml
                                mags_Dir_Off_culled_dml=[]
                                if (INIT == 1):
                                    data1 = dml
                                data2 = []
                      
                                data2,L1,L2,TimeMagL_RAND=\
                                mag_diff_rand_simple(kk,TimeMagL,dml,shuffle_type,shuffle_mode,L_Shuffle,thresholding,SC_CAT_analysis,\
                                                All_mother_daught_1BG,time_threshold,dt_type,lip_consider_dt,lip_consider_MD,within_dt_and_TMD,TMD,within_dt,subseq,\
                                                subseq_arr,within_dt_arr,MD_arr,within_dt_MD_arr,return_only_dml=False)
                       
                            #calculate CDFs:
                                if (INIT == 1):     
                                    cdf1T_x, cdf1T_y = ecdf(data1,norm=False) #plot these to see CDFs
                                    
                                    cdf1T = ecdf(data1,norm=False)
                                    cdf1T_norm = ecdf(data1,norm=True)
                                    cdf1T_norm_x,cdf1T_norm_y =ecdf(data1,norm=True)
                                    
                                cdf2T_x, cdf2T_y = ecdf(data2,norm=False) #plot these to see CDFs
                                
                                cdf2T = ecdf(data2,norm=False)
                                cdf2T_norm = ecdf(data2,norm=True)
                                cdf2T_norm_x,cdf2T_norm_y =ecdf(data2,norm=True)
                                
                                cdf1T = np.transpose(cdf1T)
                                cdf2T = np.transpose(cdf2T)
                                cdf1T_norm = np.transpose(cdf1T_norm)
                                cdf2T_norm = np.transpose(cdf2T_norm)
                                
#                                plt.clf()
#                                plt.plot(cdf1T_x,cdf1T_y/np.max(cdf1T))
#                                plt.plot(cdf2T_x,cdf2T_y/np.max(cdf2T))
#                                plt.axis('tight')
#                                plt.xlim(-2,2)
#                                plt.show()
                                """
                                print 'CDF1\n',np.shape(cdf1T_norm)
                                print'******************************************'
                                print 'CDF2\n',np.shape(cdf2T_norm)
                                """
                                
                        global m0_min
                        global m0_max
                        global m0L
                        global mags_Dir_Off_culled_dml
                        m0_min = -4.0
                        m0_max = 4.0
                        m0 = m0_min
                        
                        def dP(dml,m0,m0_min,m0_max):
                            ProbL=[] 
                            m0L=[]     
                            j=0
                            for m0 in np.arange( m0_min, m0_max+0.10,0.10):
        
                                m0 = round(m0,2)
                                dml = np.array(dml)
                                data1f = np.array(dml[dml<m0])
                                xf,list1f = ecdf(data1f,norm=False)
                                
                                # USING scipystats itemfreq function:
                                freq = itemfreq(xf)
                                # taking the number of occurences (itemfreq has 2 outputs; Column 1 contains sorted, unique values from data
                                #                                , column 2 contains their respective counts.):
                                counts1 = freq[:,1]
                                prob1 = np.sum(counts1/( len(dml) ) )
                                probt = [ prob1 for x in [j] ]
                                dplist = [(probt)]
                                ProbL.extend(dplist)
                                
                                mt = [ m0 for x in [j] ]
                                m0list = [(mt)]
                                m0L.extend(m0list)
                                j+=1   
                                
                            if (INIT==1):
                                global kk
                                kk=1
                                cdfs(data1f)
                                
                            return m0L,ProbL
                
                        if(histo_minima_dt_dml == False):
                            if (INIT == 1):
                                # array of shape n x 1, every row n corresponds to a value of m_0 
                                m0L,ProbL = ( dP(dml,m0,m0_min,m0_max) )
                                               
                        print'\n'
                        if(dt_type=='x<Dt<y'):
                            print 'time exponent (10^t) value for x<Dt<y thresholding: ',exponent
                        if(dt_type=='Dt<y'):
                            print 'time exponent (10^t) value for Dt<y thresholding: ',exponent
                        def shuffles(N):
                            global dml_rand
                            m0 = m0_min
                            
                            for kk in xrange(0,N):

                                print kk+1,
                                #"""
                                def dP_rand(dml_rand,m0,m0_min,m0_max):
                                    ProbL_rand = []                   
                                    j = 0
                                    
                                    for m0 in np.arange(m0_min,m0_max+0.10,0.10):
                                        m0=round(m0,2)
                                        dml_rand = np.array(dml_rand)
                                        data2f =np.array(dml_rand[dml_rand < m0])
                                        x2f,counts2f = ecdf(data2f,norm=False)
                                        
                                        #USING scipystats itemfreq function:
                                        freq = itemfreq(x2f)
                                        #taking the number of occurences (itemfreq has 2 outputs):
                                        counts2 = freq[:,1]
            
                                        prob2 = np.sum(counts2/len(dml_rand) )
                                        probt = [ prob2 for x in [j] ]
                                        dplist = [(probt)]   
                                        ProbL_rand.extend(dplist)                                
                                        j+=1
                                        
                                    return ProbL_rand

                                global Prob2LT
                                if(N==0):
                                    Prob2LT = ( dP_rand(data2,m0,m0_min,m0_max) )
                                else:
                                    data2_Prob2LT,L1,L2,TimeMagL_RAND=\
                                    mag_diff_rand_simple(kk,TimeMagL,dml,shuffle_type,shuffle_mode,L_Shuffle,thresholding,SC_CAT_analysis,\
                                                    All_mother_daught_1BG,time_threshold,dt_type,lip_consider_dt,lip_consider_MD,within_dt_and_TMD,TMD,within_dt,subseq,\
                                                    subseq_arr,within_dt_arr,MD_arr,within_dt_MD_arr,return_only_dml=False)
                   
                                    Prob2LT = ( dP_rand(data2_Prob2LT,m0,m0_min,m0_max) )
                                
                                global Prob2LTF

                                if (kk==0):
                                    Prob2LTF = np.transpose(Prob2LT)

                                else:
                                    Prob2LTF = np.append(Prob2LTF,np.transpose(Prob2LT),axis=0)
                                                                        
                                global std

                                std = np.std(Prob2LTF,axis=0,ddof=0)
        
                        if(histo_minima_dt_dml == False):
                            global std
        
                            shuffles(N_shuffles)
                                    
                            std = np.array(std).reshape((len(std),1))
        
                            mean_Prob2LTF = np.mean(Prob2LTF,axis=0)
                            mean_Prob2LTF = np.array(mean_Prob2LTF).reshape((len(mean_Prob2LTF),1))                  

                            dP_difference = (ProbL - mean_Prob2LTF)

                            dP_all = np.concatenate((m0L,dP_difference), axis = 1)
                            global dP_all_only
                            dP_all_only = dP_all[:,1]
                            dP_all_only = np.reshape(dP_all_only, (len(dP_all_only),1) )
                            dP_all = np.concatenate((dP_all,std),axis = 1)

                            #Saving data:
                            global timestr
                            global file_number
                            global text_file
                            
                            if (TIME_count == 0):
                                timestr = time.strftime("%Y%m%d-%H%M%S")
                                #"FREEZES" time string (timestr) so that all other outputs filenames have same times:
                                timestr = timestr
                            
                            else:
                                timestr = timestr
                            TIME_count+=1
                            
                            file_number = save_dP_arr(save_to_folder,File_Label,TimeMagL,m0,dP_all,mag_th,\
                                                      cdf1T,cdf2T,cdf1T_norm,cdf2T_norm,ProbL,\
                                                      Prob2LTF,mean_Prob2LTF,std,TIME_count,\
                                                      CAT_type,compare,group_ID,inside_file_name)

                            if compare == True:
                                if (inside_file_name != None):
                                    
                                    text_file = open("%s/data/SSM/%s/dP/GROUP_%s/%s/%s_%s/PARAMS_%s.txt"\
                                                 %(save_to_folder,CAT_type,group_ID,inside_file_name,file_number,File_Label,mag_th), "w")
                                else:   
                                    text_file = open("%s/data/SSM/%s/dP/GROUP_%s/%s_%s/PARAMS_%s.txt"\
                                                     %(save_to_folder,CAT_type,group_ID,file_number,File_Label,mag_th), "w")
                            
                            elif compare== False:               
                                if (inside_file_name != None):
                                
                                    text_file = open("%s/data/SSM/%s/dP/%s/%s_%s/PARAMS_%s.txt"\
                                                     %(save_to_folder,CAT_type,inside_file_name,file_number,File_Label,mag_th), "w")
                                else:
                                    
                                    text_file = open("%s/data/SSM/%s/dP/%s_%s/PARAMS_%s.txt"\
                                                     %(save_to_folder,CAT_type,file_number,File_Label,mag_th), "w")
                            
                            text_file.write('Catalog_Number: %s\n' %catNumb)
                            text_file.write("Lower Mag_th: %.1f\n" % MAG_TH)
                            text_file.write("CAT_type: %s\n" % CAT_type)
                            text_file.write("Background type: %s\n" % Background)
                            text_file.write("Number of shuffles: %d\n" %N_shuffles )
                            text_file.write("Shuffling type: %s\n" % shuffle_type)
                            text_file.write('Shuffle Mode: %s\n'%shuffle_mode)
                            text_file.write('Test_conditioning: %s\n'%TEST_conditioning)
                            text_file.write('All_mother_daught_1BG: %s\n'%All_mother_daught_1BG)
                            text_file.write('FIRST_gen_only: %s\n'%FIRST_gen_only)
                            text_file.write('thresholding: %s\n'%thresholding)
                            text_file.write('X_NEAREST: %s\n'%X_NEAREST) # (x nearest event. 1 (i.e. subsequent) by default)
                            text_file.write('\n')
                            
                            text_file.write('dt_type: %s\n'%dt_type) 
                            text_file.write('within_dt_and_TMD: %s\n'%within_dt_and_TMD)
                            text_file.write('within_dt: %s\n'%within_dt)
                            text_file.write('TMD: %s\n'%TMD)
                            text_file.write('unrelat_subseq: %s\n'%unrelat_subseq)  
                            text_file.write('dt_exponent: %s\n'%exponent)
                            if(dt_type=='x<Dt<y'):
                                text_file.write('dt_min (x in x<dt<y): %s\n'%dt_lower_end_percent_exponent)
                                text_file.write('dt_max (y in x<dt<y): %s\n'%10**exponent)
                            
                                            
                            text_file.write('\n')
                            text_file.write('TimeMagL_TMD: %s\n'%len(MD_arr))
                            text_file.write('TimeMagL_within_dt: %s\n'%len(within_dt_arr))
                            text_file.write('TimeMagL_within_dt_TMD: %s\n'%len(within_dt_MD_arr))
                            tmd_full=np.float(len(MD_arr)) / np.float(len(TML_for_ratio))#within_dt_for_ratio))
                            text_file.write('TimeMagL_TMD/full_CAT: %s\n'%(tmd_full) )
                            text_file.write('TimeMagL_within_time/full_CAT: %s\n'%(np.float(len(within_dt_arr)) /len(TML_for_ratio) ) )
                            text_file.write('TimeMagL_within_time_TMD/full_CAT: %s\n'%(np.float(len(within_dt_MD_arr)) /len(TML_for_ratio) ) )
                            
                            text_file.close()
                            
                            if (inside_file_name != None):
                                
                                    text_file = open('%s/data/SSM/%s/dP/%s/%s_%s/RATIOS_%s.txt'\
                                                     %(save_to_folder,CAT_type,inside_file_name,file_number,File_Label,mag_th), "w")
                            else:
                                
                                text_file = open('%s/data/SSM/%s/dP/%s_%s/RATIOS_%s.txt'\
                                                 %(save_to_folder,CAT_type,file_number,File_Label,mag_th), "w")
                                
                            text_file.write('#TMD/full,-- within_Dt/full,-- within_Dt_TMD/full,-- 10**exponent\n')
                            text_file.write('%1.5f\n'%(np.float(len(MD_arr)) / np.float(len(TML_for_ratio))) )
                            text_file.write('%1.5f\n'%(np.float(len(within_dt_arr)) /np.float(len(TML_for_ratio)) ) )
                            text_file.write('%1.5f\n'%(np.float(len(within_dt_MD_arr)) /np.float(len(TML_for_ratio)) ) )
                            text_file.write('%1.5f'%(10.**exponent ) )
                            text_file.close()
                        
                ########################################
                            # SAVE TML, histogram fig and dP fig:
                            path = '%s/data/SSM/%s/dP/%s/%s_%s'%(save_to_folder,CAT_type,inside_file_name,file_number,File_Label)
                            path_save_figs='%s/data/SSM/%s/dP/%s/FIGS'%(save_to_folder,CAT_type,inside_file_name)
                            if not (os.path.exists('%s/histo/m_th_%s'%(path_save_figs,mag_th) and '%s/dP/m_th_%s'%(path_save_figs,mag_th) and '%s/dP/m_th_%s_fixed_frame'%(path_save_figs,mag_th) )):
                                os.makedirs('%s/histo/m_th_%s'%(path_save_figs,mag_th))
                                os.makedirs('%s/dP/m_th_%s'%(path_save_figs,mag_th))
                                os.makedirs('%s/dP/m_th_%s_fixed_frame'%(path_save_figs,mag_th))

                            if(subseq == True and SC_CAT_analysis == False):
                                if (dt_type == None):
                                    np.savetxt('%s/TimeMagL_mth=%.2f_%s.txt'\
                                               %(path,mag_th,CAT_NUMB),TimeMagL,\
                                               fmt='%f %1.2f %i %1.0f %1.0f %1.0f %1.0f %1.0f',delimiter=' ')
                                else:
                                    print 'Error in setting parameters'
                                    break
                            if(subseq == True and SC_CAT_analysis == True):
                                if (dt_type == None):

                                    np.savetxt('%s/TimeMagL_uncond_mth=%.2f_%s.txt'\
                                               %(path,mag_th,CAT_NUMB),TimeMagL,\
                                               fmt='%1f %1.2f',delimiter=' ')
        
                                else:
                                    print 'Error in setting parameters'
                                    break
                                    
                            ################
                            if ((within_dt==True or TMD==True or within_dt_and_TMD==True) and SC_CAT_analysis==False):
                               
                                #daughter_mother_arr has the following column values:
                                ##m_as(0),M(1),dm(2),M_ID(3),m_as_tree(4),M_tree(5),m_as_M_ID(6),m_as_time(7),dt(8):
                                
                                np.savetxt('%s/mas_M_arr_mth=%.2f_%s.txt'%(path,mag_th,CAT_NUMB),mas_M_arr,fmt='%1.2f %1.2f %1.2f %f %i %i %i %i %i %i %i %i %i %i %i %i %i',delimiter=' ')
                            if(within_dt==True and (SC_CAT_analysis==True and (dt_type=='x<Dt<y' or dt_type=='Dt<y') )):

								np.savetxt('%s//mas_M_arr_within_dt=%.1E_mth=%.2f_%s.txt'\
										   %(path,time_threshold,mag_th,CAT_NUMB),mas_M_arr,\
										   fmt='%1.2f %1.2f %1.2f %f %f',delimiter=' ')
                            
                            ################################################
                            
                            plt.clf()
                            fig = plt.figure()
                            time_threshold_print_exponent='$dt<10^{%s}$'%exponent
                            plt.errorbar(dP_all[:,0],dP_all[:,1], yerr = 3*dP_all[:,2],ecolor='m',label='%s$\,(3\sigma)$\n%s'%(mag_th,time_threshold_print_exponent))
                            plt.legend(loc='lower right')
                            plt.grid(True,linestyle='-')
                            plt.savefig('%s/plot_dP_%s.png'%(path,mag_th), bbox_inches='tight')

                            plt.clf()
                            plt.errorbar(dP_all[:,0],dP_all[:,1], yerr = 3*dP_all[:,2],ecolor='m',label='%s$\,(3\sigma)$\n%s'%(mag_th,time_threshold_print_exponent))
                            plt.legend(loc='lower right')
                            plt.grid(True,linestyle='-')

                            plt.savefig('%s/dP/m_th_%s/plot_dP_%04d.png'%(path_save_figs,mag_th,image_counter),dpi=200 )
                            plt.ylim(-0.1,0.1)
                            plt.savefig('%s/dP/m_th_%s_fixed_frame/plot_dP_%04d.png'%(path_save_figs,mag_th,image_counter),dpi=200 )
                            plt.close(fig)
                            ###############################################
                            
                            if ( len(TimeMagL)>0 or len(mas_M_arr>0) ):
                                plt.clf()
                         
                                if(subseq==True):
                                    counts,bins=np.histogram(TimeMagL[:,1],bins=20)
                                    counts_daughters,bins_daughters=np.histogram(mas_M_arr[:,0],bins=20)
                                    bins=bins[1:]

                                    fig = plt.figure()
                                    ax = plt.gca()
                                    ax.scatter(bins,counts,marker='s',color='b',label='$m$')
                                
                                    plt.title('Frequency-magnitude distribution')
                                    ax.set_yscale('log')

                                    y_high = lambda x: 3*np.max(counts) * 10**( -(0.90)*(x-mag_th) )
                                    plt.plot(bins , y_high(bins), '--', color='b',linewidth=2,label = 'b=0.90 (g+z)')
                                    
                                    plt.legend(loc='upper right')

                                    fig.savefig('%s/histogram_L1-L2_%s.png'%(path,mag_th) )#, bbox_inches='tight')
                                    if(dt_type=='Dt<y' or dt_type=='x<Dt<y'):
                                        fig.savefig('%s/histo/m_th_%s/histogram_m_as-M_%04d.png'%(path_save_figs,mag_th,image_counter),dpi=200 )

                                    else:
                                        fig.savefig('%s/histo/m_th_%s/histogram_L1-L2_%04d.png'%(path_save_figs,mag_th,image_counter),dpi=200 )

                                    plt.close(fig)
                                    plt.clf()
                                    
                                #######################
                                if( (All_mother_daught_1BG==False and\
                                     ( thresholding=='subsequent') ) ):
                                    print 'plotting from 1st'
                                    #daughter_mother_arr has the following column values:
                                    #m_as(0),M(1),dm(2),M_ID(3),m_as_tree(4),M_tree(5),m_as_M_ID(6),m_as_time(7)                                    
                                    if(within_dt==True or TMD == True or within_dt_and_TMD==True):
                                        plt.clf()
                                        print 'Histogram for MDs or subseq_dt'
                                        counts_M,bins_M=np.histogram(mas_M_arr[:,1],bins=20)
                                        counts_mas,bins_mas=np.histogram(mas_M_arr[:,0],bins=20)
                                        
                                        bins_M=bins_M[1:]
                                        bins_mas=bins_mas[1:]

                                        fig = plt.figure()
                                        ax = plt.gca()
                                        
                                        ax.scatter(bins_M,counts_M,marker='s',color='b',label='$M$')
                                        ax.scatter(bins_mas,counts_mas,color='g',label='$m_{as}$')
                                        plt.title('Frequency-magnitude distribution')
                                        ax.set_yscale('log')
                                        
                                        y_high = lambda x: 35*np.max(counts_M) * 10**( -(0.90)*(x-mag_th) )
                                        plt.plot(bins_M , y_high(bins_M), '--', color='b',linewidth=2,label = 'b=0.90 (g+z)')
                                        
                                        y_low = lambda x: 1.3*np.max(counts_mas) * 10**( -(0.24)*(x-mag_th) )
                                        plt.semilogy(bins_mas, y_low(bins_mas), '--', color='#FF8C00',linewidth=2,label = 'b=0.24 (z)')
                                        
                                        y_low = lambda x: 1.3*np.max(counts_mas) * 10**( -(0.9)*(x-mag_th) )
                                        plt.semilogy(bins_mas, y_low(bins_mas), '--', color='g',linewidth=2,label = 'b=0.9 (g+z)')
                                        plt.legend(loc='upper right')
                                        
                                        plt.grid(True,linestyle='-')
                                        fig.savefig('%s/histogram_L1-L2_%s.png'%(path,mag_th) )
                                        
                                        fig.savefig('%s/histo/m_th_%s/histogram_m_as-M_%04d.png'%(path_save_figs,mag_th,image_counter),dpi=200 )#, bbox_inches='tight')
                                        plt.close(fig)

                                        plt.clf()
                                    
                                elif(All_mother_daught_1BG==True):
                                    print 'plotting from 3rd'
                                    print 'see end of program for list of things to save'

                                        
                                if (inside_file_name != None):
                                    
                                    print 'SAVING HISTO FILES'
                                  
                                    np.savetxt("%s/data/SSM/%s/dP/%s/%s_%s/L1_mags_%s.txt"\
                                                     %(save_to_folder,CAT_type,inside_file_name,file_number,File_Label,mag_th),L1)
                                    
                                    np.savetxt("%s/data/SSM/%s/dP/%s/%s_%s/L2_mags_%s.txt"\
                                                         %(save_to_folder,CAT_type,inside_file_name,file_number,File_Label,mag_th),L2)
                                                                        
                                else:    
									print 'in else'

                                time_threshold_print_exponent='$dt<10^{%s}$'%exponent
                                
                                plt.errorbar(dP_all[:,0],dP_all[:,1], yerr =3*dP_all[:,2],ecolor='m',label='%s$\,(3\sigma$)\n%s'%(mag_th,time_threshold_print_exponent))
                                plt.legend(loc='lower right')
                               
                MAG_CORR(Low_Mag_Thr,CAT_type,Background,shuffle_type,N_shuffles)
                global image_counter
                image_counter+=1

                print 'END main function\n'
                           
                if (histo_minima_dt_dml == False):
                    
                    if compare == True:
                        path_files = '%s/data/SSM/%s/dP/GROUP_%s/%s'%(save_to_folder,CAT_type,group_ID,file_number)
                   
                    if compare == False:                      
                        path_files = '%s/data/SSM/%s/dP/%s'%(save_to_folder,CAT_type,file_number)
                
                elif(histo_minima_dt_dml == True):
                    print'----- min_dt_list', np.shape(min_dt_list)
                    
                if (reps!=-1):
                    global dP_all_mean
                    
                    if(reps == 0):
                        dP_all_mean = dP_all_only
                    else:
                        dP_all_mean = np.concatenate( (dP_all_mean,dP_all_only), axis = 1)
                    print 'shpae dP-all-mean', np.shape(dP_all_mean)
                    global mean_dP_all
                    mean_dP_all = np.mean(dP_all_mean,axis=1)
                    print'MEAN shape ', np.shape(mean_dP_all)

            params_dP (compare=False,shuffle_type = shuffle_type,shuffle_mode = shuffle_mode,group_ID=0,catNumb=CAT_NUMB)
