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

import sys
import numpy as np
import os
cwd = os.getcwd()

def mag_diff_simple(split_CAT_Arr,mag_th,TimeMagL,SC_CAT_analysis,shuffle_type,time_threshold,\
                    dt_type,MD_1st_gen,within_dt_and_TMD,TMD,within_dt,subseq,M_range_min,M_range_max,\
                    mothers_consider_bckgrnd_only,mothers_consider_AS_only,no_background):
    print '\n'
    print 'INSIDE mag_diff_simple'
    subseq_arr=[]
    within_dt_arr=[]
    MD_subseq_repeated_arr=[]
    MD_subseq_arr=[]
    within_dt_MD_subseq_arr=[]
    MD_arr=[]
    within_dt_MD_arr=[]
    
    #TimeMagLT has the following columns:
    #(t_as(0),m_as(1),tree(2),gen(3),leaf(4),leaf_mom(5),global_ID(6),mother_ID(7) )
    if(no_background == True):
        TimeMagL = filter(lambda row:row[7]!=-2, TimeMagL)
    TimeMagL=np.array(TimeMagL)
        
    if time_threshold==None:
        time_threshold=np.max(TimeMagL[:,0])
    #############################################################
    if MD_1st_gen==True:
        i=0
        while(i<len(TimeMagL)-1):
            m_as = TimeMagL[i+1,1]#np.around(TimeMagL[i+1,1],decimals=2)
            M = TimeMagL[i,1]#np.around(TimeMagL[i,1],decimals=2)
            dm = np.float(m_as-M)
            dt = np.float(TimeMagL[i+1,0] - TimeMagL[i,0])
                
            #TimeMagLT has the following columns:
            #(t_as(0),m_as(1),tree(2),gen(3),leaf(4),leaf_mom(5),global_ID(6),mother_ID(7) )
            if(SC_CAT_analysis==True):
                subseq_arr.append( [m_as,M,dm,dt,TimeMagL[i+1,0]] )
            
            elif( SC_CAT_analysis == False ):
                subseq_arr.append([m_as,M,dm,dt,-3,\
                                   TimeMagL[i,2],TimeMagL[i,3],TimeMagL[i,4],\
                                   TimeMagL[i,5],TimeMagL[i,6],TimeMagL[i,7],\
                                   TimeMagL[i+1,2],TimeMagL[i+1,3],TimeMagL[i+1,4],\
                                   TimeMagL[i+1,5],TimeMagL[i+1,6],TimeMagL[i+1,7]])
            i+=1
        print'shape arr subseq MD 1st gen',np.shape(subseq_arr)
        subseq_arr=np.array(subseq_arr)
        
        #Find events within dt threshold for subseq_arr:
        if subseq==False:
            if (within_dt==True or within_dt_and_TMD==True):
                print'TIME THRESHOLD %s(s)'%time_threshold
                mask_within_dt_lower = subseq_arr[:,3]>120
                mask_within_dt_upper = subseq_arr[:,3]<time_threshold
                if dt_type=='Dt<y':
                    mask_within_dt=subseq_arr[:,3]<time_threshold
                elif dt_type=='x<Dt<y':
                    mask_within_dt= (mask_within_dt_lower) & (mask_within_dt_upper)    
                within_dt_arr=subseq_arr[mask_within_dt,:]
                print 'Shape within dt arr',np.shape(within_dt_arr)
                
                
        print'Searching for mother-daughter (1st gen) pairs'

        if split_CAT_Arr == '174212':
            print 'loading MD_arr file for SSAR-SC...'
            MD_arr = np.loadtxt('%s/MD_CAT/MD_CAT_%s.txt'%(cwd,mag_th))
            print 'done.'
        elif split_CAT_Arr == '161822':
            print 'loading MD_arr file for SSAR-LRG1...'
            MD_arr = np.loadtxt('%s/MD_CAT_LRG1/MD_CAT_%s.txt'%(cwd,mag_th))
            print 'done.'
        print 'shape MD_arr',np.shape(MD_arr)
            
#        else:
#            
#            for q in range(0,len(TimeMagL)-1):
#                f=np.where((TimeMagL[:,7] == TimeMagL[q,6]) &\
#                           (TimeMagL[:,3] == TimeMagL[q,3]+1) )[0]
#                if len(f)>=1:
#                    f=f[0]
#                    M=TimeMagL[q,1]
#                    m_as=TimeMagL[f,1]
#                    dm=m_as-M
#                    dt=np.float(TimeMagL[f,0] - TimeMagL[q,0])
#                    
#                    MD_subseq_pair_val=1
#                    MD_arr.append([m_as,M,dm,dt,MD_subseq_pair_val,\
#                                   TimeMagL[f,2],TimeMagL[f,3],TimeMagL[f,4],\
#                                   TimeMagL[f,5],TimeMagL[f,6],TimeMagL[f,7],\
#                                   TimeMagL[q,2],TimeMagL[q,3],TimeMagL[q,4],\
#                                   TimeMagL[q,5],TimeMagL[q,6],TimeMagL[q,7]])
#        
#            MD_arr=np.array(MD_arr)
#            np.savetxt('%s/MD_CAT_LRG1_leo/MD_CAT_%s.txt'%(cwd,mag_th),MD_arr)
#            print 'shape MD_arr',np.shape(MD_arr)
        
        if ( SC_CAT_analysis == False ):
            j=0
#            MD_subseq_arr=MD_arr
            print 'shape MD ARR before',np.shape(MD_arr)
            
            if(mothers_consider_bckgrnd_only == True and (TMD==True or within_dt_and_TMD==True) ):
                print 'considering only mothers that are background events'
                mask_mothers_bckgrnd =MD_arr[:,9]==-2
                MD_arr =MD_arr[mask_mothers_bckgrnd]
            elif(mothers_consider_AS_only == True and (TMD==True or within_dt_and_TMD==True) ):
                print 'considering only mothers that are AS events'
                mask_mothers_AS = MD_arr[:,9]!=-2 #KEEPS ALL COLS EXCEPT THOSE == -2 (mother events)
                MD_arr = MD_arr[mask_mothers_AS]
    
            if (M_range_min!=None and M_range_max!=None):
                print'Keeping only preceding events (M) in range: %s<M<%s'% (M_range_min,M_range_max)
                mask_min_max_M_range = (MD_arr[:,1] > M_range_min) & (MD_arr[:,1] < M_range_max)
                MD_arr = MD_arr[mask_min_max_M_range,:]
        
        if (within_dt_and_TMD==True):
            mask_within_dt_lower = MD_arr[:,3]>120
            mask_within_dt_upper = MD_arr[:,3]<time_threshold
            if dt_type=='Dt<y':
                mask_within_dt = MD_arr[:,3]<time_threshold
            elif dt_type=='x<Dt<y':
                mask_within_dt = (mask_within_dt_lower) & (mask_within_dt_upper)    
            within_dt_MD_arr = MD_arr[mask_within_dt,:]
            print 'Shape within dt MD arr',np.shape(within_dt_MD_arr)
            
#####################################################################################################################
#####################################################################################################################               
    else:
        i=0
        while(i<len(TimeMagL)-1):
            m_as = TimeMagL[i+1,1]#np.around(TimeMagL[i+1,1],decimals=2)
            M = TimeMagL[i,1]#np.around(TimeMagL[i,1],decimals=2)
            dm = np.float(m_as-M)
            dt = np.float(TimeMagL[i+1,0] - TimeMagL[i,0])
            
            if (SC_CAT_analysis == False):
                MD_subseq_pair_val=None
                test_MD = (  (TimeMagL[i+1,7] == TimeMagL[i,6])\
                            and (TimeMagL[i+1,2] == TimeMagL[i,2])\
                            and (TimeMagL[i+1,5] == TimeMagL[i,4])\
                            and (TimeMagL[i+1,3] == TimeMagL[i,3]+1) ) #TESTS FOR 1ST GEN ONLY
                MD_subseq_pair_val=1 if test_MD==True else 0
                
            #TimeMagLT has the following columns:
            #(t_as(0),m_as(1),tree(2),gen(3),leaf(4),leaf_mom(5),global_ID(6),mother_ID(7) )
            if(SC_CAT_analysis==True):
                subseq_arr.append( [m_as,M,dm,dt,TimeMagL[i+1,0]] )
            
            elif( SC_CAT_analysis == False ):
                subseq_arr.append([m_as,M,dm,dt,MD_subseq_pair_val,\
                                   TimeMagL[i,2],TimeMagL[i,3],TimeMagL[i,4],\
                                   TimeMagL[i,5],TimeMagL[i,6],TimeMagL[i,7],\
                                   TimeMagL[i+1,2],TimeMagL[i+1,3],TimeMagL[i+1,4],\
                                   TimeMagL[i+1,5],TimeMagL[i+1,6],TimeMagL[i+1,7]])
    
                #pick out only MD pairs that are subsequent to each other:
                if MD_subseq_pair_val == 1:
                    MD_subseq_repeated_arr.append([m_as,M,dm,dt,MD_subseq_pair_val,\
                                            TimeMagL[i,2],TimeMagL[i,3],TimeMagL[i,4],\
                                            TimeMagL[i,5],TimeMagL[i,6],TimeMagL[i,7],\
                                            TimeMagL[i+1,2],TimeMagL[i+1,3],TimeMagL[i+1,4],\
                                            TimeMagL[i+1,5],TimeMagL[i+1,6],TimeMagL[i+1,7]])
            i+=1
        print'shape arr subseq 2',np.shape(subseq_arr)
    #    print'shape MD repeated arr',np.shape(MD_subseq_repeated_arr) 
        subseq_arr=np.array(subseq_arr)
        MD_subseq_repeated_arr=np.array(MD_subseq_repeated_arr)
        
        #Find events within dt threshold for subseq_arr:
        if subseq==False:
            if (within_dt==True or within_dt_and_TMD==True):
                print'TIME THRESHOLD %s(s)'%time_threshold
                mask_within_dt_lower = subseq_arr[:,3]>120
                mask_within_dt_upper = subseq_arr[:,3]<time_threshold
                if dt_type=='Dt<y':
                    mask_within_dt=subseq_arr[:,3]<time_threshold
                elif dt_type=='x<Dt<y':
                    mask_within_dt= (mask_within_dt_lower) & (mask_within_dt_upper)    
                within_dt_arr=subseq_arr[mask_within_dt,:]
                print 'Shape within dt',np.shape(within_dt_arr)
            
            if ( SC_CAT_analysis == False ):
                j=0
                MD_subseq_arr=MD_subseq_repeated_arr
                print 'shape MD ARR before',np.shape(MD_subseq_arr)
                
                if(mothers_consider_bckgrnd_only == True and (TMD==True or within_dt_and_TMD==True) ):
                    print 'considering only mothers that are background events'
                    mask_mothers_bckgrnd = MD_subseq_arr[:,9]==-2
                    MD_subseq_arr = MD_subseq_arr[mask_mothers_bckgrnd]
                elif(mothers_consider_AS_only == True and (TMD==True or within_dt_and_TMD==True) ):
                    print 'considering only mothers that are AS events'
                    mask_mothers_AS = MD_subseq_arr[:,9]!=-2 #KEEPS ALL COLS EXCEPT THOSE == -2 (mother events)
                    MD_subseq_arr = MD_subseq_arr[mask_mothers_AS]
        
                if (M_range_min!=None and M_range_max!=None):
                    print'Keeping only preceding events (M) in range: %s<M<%s'% (M_range_min,M_range_max)
                    mask_min_max_M_range = (MD_subseq_arr[:,1] > M_range_min) & (MD_subseq_arr[:,1] < M_range_max)
                    MD_subseq_arr=MD_subseq_arr[mask_min_max_M_range,:]
                    
                #Filter out repeated rows: if there is an [mas_id1-M_id1] pair and another row exists with
                # pair [mas_id2-M_id2] (where M_id2=mas_id1), all rows where [M_id2=mas_id1] are removed.
                # it will find the *first* next M_id that matches mas_id, remove that row and then go back to 
                # the beginning of the list to repeat the search. THis will be done for all column values in
                # mas_id:      
                while (j<len(MD_subseq_arr)-1):
                    mas_id = MD_subseq_arr[:,9]
                    M_id = MD_subseq_arr[:,15]
                    indices=np.nonzero(mas_id==M_id[j])
            #        print 'len MD arr',np.shape(MD_subseq_arr)
                    if indices!=[]:
                        MD_subseq_arr=np.delete(MD_subseq_arr,indices,axis=0)
                    else:
                        pass
                    j+=1
                #Find MD within dt threshold:
                print 'shape MD ARR AFTER',np.shape(MD_subseq_arr)
                mask_within_dt_MD = MD_subseq_arr[:,3]<time_threshold
            #    print 'mask within dt md', mask_within_dt_MD[:1000]
                within_dt_MD_subseq_arr=MD_subseq_arr[mask_within_dt_MD,:]
                print'shape within dt MD arr:',np.shape(within_dt_MD_subseq_arr)
                  
    #get dm list:
    if subseq==True:
        dml=subseq_arr[:,2]
    elif within_dt==True:
        dml=within_dt_arr[:,2]
    elif TMD==True:
        if MD_1st_gen==True:
            dml=MD_arr[:,2]
            print'dml len',np.shape(dml)
        else:
            dml=MD_subseq_arr[:,2]
            print'dml len',np.shape(dml)
    elif within_dt_and_TMD==True:
        if MD_1st_gen==True:
            dml=within_dt_MD_arr[:,2]
        else:
            dml=within_dt_MD_subseq_arr[:,2]
    else:
        print '****** ERROR in finding dml in func_mag_diff ***'
        sys.exit()
    
    MD_subseq_arr=MD_arr if MD_1st_gen==True else MD_subseq_arr
    within_dt_MD_subseq_arr=within_dt_MD_arr if MD_1st_gen==True else within_dt_MD_subseq_arr
    
    if (SC_CAT_analysis == False):
        print'returning dml from simp\n'
        return(dml,subseq_arr,within_dt_arr,MD_subseq_arr,within_dt_MD_subseq_arr)
    else:
        print'returning dml from simp\n'
        return(dml,subseq_arr,within_dt_arr,MD_subseq_arr,within_dt_MD_subseq_arr)
