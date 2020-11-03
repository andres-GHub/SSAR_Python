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

def mag_diff_rand_simple(kk,TimeMagL,dml,shuffle_type,Shuffle_Mode,L_Shuffle,thresholding,SC_CAT_analysis,\
                         All_mother_daught_1BG,time_threshold,dt_type,lip_consider_dt,lip_consider_MD,within_dt_and_TMD,TMD,within_dt,subseq,\
                         subseq_arr,within_dt_arr,MD_arr,within_dt_MD_arr,return_only_dml):
    
    if kk<1:
        print '########INSIDE MAG_DIFF_RAND_SIMPLE############'
    global dml_rand
    global L_RAND
    global dm_rand
    dm_rand = -99999999.9
    dml_rand=[]
    L_RAND = []
    ########################################################################################        
    global L1
    global L2
        
    if(shuffle_type=='jo' or shuffle_type=='lip'):
        if subseq == True:
            if kk<1:
                print'SUBSEQ'
            arr_to_shuffle=subseq_arr
        elif within_dt == True:
            if kk<1:
                print'within dt'
            arr_to_shuffle = within_dt_arr
        elif TMD == True:
            if kk<1:
                print'TMD'
            arr_to_shuffle = MD_arr  
        elif within_dt_and_TMD == True:
            if kk<1:
                print'within_dt_and_TMD'
            arr_to_shuffle = within_dt_MD_arr               
        else:
            print 'ERROR in thresholding type'
            
        if (kk<1):
            print'---------------------'
            print'SHAPE M-mas ARR',np.shape(arr_to_shuffle)
                
        L1 = arr_to_shuffle[:,1] #MOTHERS or PRECEDING EVENT
        L2 = arr_to_shuffle[:,0] #DAUGHTERS or SUBSEQUENT EVENT
        L1 = np.reshape( L1,(len(L1),1) )
        L2 = np.reshape( L2,(len(L2),1) )
 
        replace_cond=True if Shuffle_Mode=='Replace' else False
        if shuffle_type=='jo':
            if (L_Shuffle =='L1'):
                L_RAND = L1[np.random.choice(L1.shape[0], len(L1), replace=replace_cond),:]
            if (L_Shuffle =='L2'):
                L_RAND = L2[np.random.choice(L2.shape[0], len(L2), replace=replace_cond),:]
        elif shuffle_type =='lip':
            if (L_Shuffle =='L1'):
                indx = np.arange(0,len(L1),dtype=int)
                indx = np.random.choice(indx,size=len(indx),replace=replace_cond)
                L_RAND = TimeMagL[indx]
                L_RAND = L_RAND[:,0]
                L_RAND = np.reshape( L_RAND,(len(L_RAND),1) )
            if (L_Shuffle =='L2'):
                indx = np.arange(0,len(L2),dtype=int)
                indx = np.random.choice(indx,size=len(indx),replace=replace_cond)
                L_RAND = TimeMagL[indx,1]
                L_RAND = np.reshape( L_RAND,(len(L_RAND),1) )
        
        for i in xrange(0,len(L_RAND)):
            #subtracting, L2-L1:
            if (L_Shuffle=='L1'):
                dm_rand = [ (-L_RAND[i,0] + L2[i,0]) for x in [i] ]
                dml_rand.extend(dm_rand)
                    
            if (L_Shuffle=='L2'):
                dm_rand = [ (-L1[i,0] + L_RAND[i,0]) for x in [i] ]
                dml_rand.extend(dm_rand)  
        if(thresholding == 'subsequent_true_daught' or\
           thresholding == 'subsequent' or All_mother_daught_1BG==True):
            if ( len(dml_rand) > len(dml) ):
                dml_rand = dml_rand[0:len(dml)]
                if (kk<1):
                    print 'dml_rand > dml'
            
            elif ( len(dml_rand) < len(dml) ):
                dml = dml[0:len(dml_rand)]
                if (kk<1):
                    print 'dml > dml_rand'
            
            elif ( len(dml_rand) == len(dml) ):
                if (kk<1):
                    print 'Length of lists dml_rand and dml are the same'
        if (kk<1):
            print 'FUNC_mag_diff_rand -- dml_rand --',np.shape(dml_rand)
            print 'FUNC_mag_diff_rand -- TML_RAND --',np.shape(L_RAND)
            
            print'----------------END OF FUNC------------------\n'
            
        if(return_only_dml==True):
            return dml_rand,L_RAND
        else:
            return (dml_rand,L1,L2,L_RAND)
