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

def load_CAT_type_1st(cwd,load_from_folder,INIT,CAT_type,catNumb,Background,\
                  histo_minima_dt_dml,First_Load,mag_th,exponent,min_range_exponent,All_mother_daught_1BG):
    if (CAT_type=='SSM'):
        
        print 'MAG_CORR_SSM  MAG_CORR_SSM  MAG_CORR_SSM\n'
        
        if (INIT==0):
            #TimeMagLT has the following columns(t_as,m_as,tree,gen,leaf,leaf_mom,global_ID,mother_ID)
            from Syn_Cat12 import catalog
            CAT_goal='A'
            TimeMagBackLT,TimeMagLT,TimeMagMD,Gen1dfL = catalog(g,z,p,b,b_as,Mmin,Mmax,c_0,tau_0,k_prod,lam,T,CAT_goal)
            TimeMagL = filter(lambda row:row[1]>=mag_th, TimeMagLT)
            INIT=1
        else:
            pass
          
       
    elif (CAT_type=='loadSSM'):
        print '**** LOADING PREVIOUS CAT ****\n'
        
        if (INIT==0):
    
            if (Background == '1BG'):
                
                if(histo_minima_dt_dml == True and exponent == min_range_exponent and All_mother_daught_1BG==False ):
    #                                        TimeMagLT = np.loadtxt('%s/data/SSM/1BG/CAT/min_times_%s_%.1EAS/%s/TimeMagLT.txt'%(cwd,min_times_trigger_mag,catNumb))
                    TimeMagLT = np.loadtxt('%s/1BG/CAT/min_times_%s/%s/TimeMagLT.txt'%(load_from_folder,min_times_trigger_mag,catNumb))
                    
                    TimeMagLT = TimeMagLT[1:]
                    
                elif(exponent == min_range_exponent and All_mother_daught_1BG==False):
                    
    #                                        TimeMagLT = np.loadtxt('%s/data/SSM/1BG/CAT/%s/TimeMagLT.txt'%(cwd,catNumb))
                    TimeMagLT = np.loadtxt('%s/1BG/CAT/%s/TimeMagLT.txt'%(load_from_folder,catNumb))
                    TimeMagLT = TimeMagLT[1:]
                    
                elif(All_mother_daught_1BG==True):
    #                                        TimeMagLT = np.loadtxt('%s/data/SSM/1BG/CAT/%s/LIST_MAGS_AllMD_Back6.0.txt'%(cwd,catNumb))
                    print 'loading All mother-daughter CAT'
    #                                        TimeMagLT = np.loadtxt('%s/data/SSM/1BG/CAT/%s/LIST_MAGS_AllMD_Back7.0.txt'%(cwd,catNumb))
                    TimeMagLT = np.loadtxt('%s/1BG/CAT/%s/LIST_MAGS_AllMD_Back7.0.txt'%(load_from_folder,catNumb))
            if (Background == 'Full' and First_Load == 0):
                
    #                                    TimeMagLT = np.array(np.loadtxt('%s/data/SSM/SSM_Full/CAT/%s/TimeMagLT.txt'%(cwd,catNumb) ) )  
                TimeMagLT = np.array(np.loadtxt('%s/SSM_Full/CAT/%s/TimeMagLT.txt'%(load_from_folder,catNumb) ) ) #,usecols=(0, 10) when wanting only magnitudes from HS
                First_Load  = 1
            #TimeMagL = np.array(filter(lambda row:row[1]>=(mag_th), TimeMagLT) )
    
            INIT=1
        else:
            pass
        
        if(All_mother_daught_1BG==False):
            TimeMagL = np.array(filter(lambda row:row[1]>=(mag_th), TimeMagLT) )
        elif(All_mother_daught_1BG==True):
            #S CURVE: TimeMagL = np.array(filter(lambda row:( row[0]>=(mag_th) and row[1]>=(mag_th) ), TimeMagLT) )
            TimeMagL = np.array(filter(lambda row:1.0<=row[0]<=np.max(TimeMagLT),TimeMagLT))#np.max(TimeMagLT)
    #                                TimeMagL = np.array(filter(lambda row:row[0]>=(mag_th),TimeMagLT))
    
    #                            print'Filtered CAT length: ',len(TimeMagL)
    #                            print'length full CAT: ',len(TimeMagLT)
    else:
        print 'MAG_CORR_CATALOG FOR SC MAG_CORR_CATALOG  MAG_CORR_CATALOG\n'
        
        if (INIT==0):
    
            catNumb = 'hs_1981_2016'
#            catNumb = 'hs_1984_2017_byhand'
            print 'PLOTS FOR MAG CORRELATIOSN OF THE SC CATALOG: \n', catNumb
            TimeMagLT = np.loadtxt('%s/%s.txt'%(cwd,catNumb), usecols=(0, 10) )
            TimeMagL = np.array(filter(lambda row:row[1]>=mag_th, TimeMagLT))
            INIT=1
            
        else:
            pass
        
        TimeMagL = np.array(filter(lambda row:row[0]>=mag_th, TimeMagLT))
        
      
    print 'Mag_Corr THRESHOLD VALUE:',mag_th
    print '\n'
    print 'Length Full CAT: ',len(TimeMagLT)
    print 'Min mag Full CAT:',(np.min(TimeMagLT[:,1]))
    print 'length of CAT after m_th > %s'%mag_th, len(TimeMagL)
    print 'Min mag TimeMagL:',(np.min(TimeMagL[:,1]))

    if INIT==First_Load:
        First_Load+=1
        return(INIT,First_Load,TimeMagLT,TimeMagL)
    
    else:
        print 'ERROR'
        
        
        
        
def load_CAT_type_2nd(cwd,load_from_folder,INIT,CAT_type,catNumb,Background,\
                  histo_minima_dt_dml,First_Load,mag_th,exponent,min_range_exponent,All_mother_daught_1BG,TimeMagLT,TimeMagL):
    if (CAT_type=='SSM'):
        
        print 'MAG_CORR_SSM  MAG_CORR_SSM  MAG_CORR_SSM\n'
        
        if (INIT==0):
            #TimeMagLT has the following columns(t_as,m_as,tree,gen,leaf,leaf_mom,global_ID,mother_ID)
            from Syn_Cat12 import catalog
            CAT_goal='A'
            TimeMagBackLT,TimeMagLT,TimeMagMD,Gen1dfL = catalog(g,z,p,b,b_as,Mmin,Mmax,c_0,tau_0,k_prod,lam,T,CAT_goal)
            TimeMagL = filter(lambda row:row[1]>=mag_th, TimeMagLT)
            INIT=1
        else:
            pass
          
       
    elif (CAT_type=='loadSSM'):
        print '**** LOADING PREVIOUS CAT ****\n'
        
        if (INIT==0):
    
            if (Background == '1BG'):
                
                if(histo_minima_dt_dml == True and exponent == min_range_exponent and All_mother_daught_1BG==False ):
    #                                        TimeMagLT = np.loadtxt('%s/data/SSM/1BG/CAT/min_times_%s_%.1EAS/%s/TimeMagLT.txt'%(cwd,min_times_trigger_mag,catNumb))
                    TimeMagLT = np.loadtxt('%s/1BG/CAT/min_times_%s/%s/TimeMagLT.txt'%(load_from_folder,min_times_trigger_mag,catNumb))
                    
                    TimeMagLT = TimeMagLT[1:]
                    
                elif(exponent == min_range_exponent and All_mother_daught_1BG==False):
                    
    #                                        TimeMagLT = np.loadtxt('%s/data/SSM/1BG/CAT/%s/TimeMagLT.txt'%(cwd,catNumb))
                    TimeMagLT = np.loadtxt('%s/1BG/CAT/%s/TimeMagLT.txt'%(load_from_folder,catNumb))
                    TimeMagLT = TimeMagLT[1:]
                    
                elif(All_mother_daught_1BG==True):
    #                                        TimeMagLT = np.loadtxt('%s/data/SSM/1BG/CAT/%s/LIST_MAGS_AllMD_Back6.0.txt'%(cwd,catNumb))
                    print 'loading All mother-daughter CAT'
    #                                        TimeMagLT = np.loadtxt('%s/data/SSM/1BG/CAT/%s/LIST_MAGS_AllMD_Back7.0.txt'%(cwd,catNumb))
                    TimeMagLT = np.loadtxt('%s/1BG/CAT/%s/LIST_MAGS_AllMD_Back7.0.txt'%(load_from_folder,catNumb))
            if (Background == 'Full' and First_Load == 0):
                
    #                                    TimeMagLT = np.array(np.loadtxt('%s/data/SSM/SSM_Full/CAT/%s/TimeMagLT.txt'%(cwd,catNumb) ) )  
                TimeMagLT = np.array(np.loadtxt('%s/SSM_Full/CAT/%s/TimeMagLT.txt'%(load_from_folder,catNumb) ) ) 
                First_Load  = 1
            #TimeMagL = np.array(filter(lambda row:row[1]>=(mag_th), TimeMagLT) )
    
            INIT=1
        else:
            pass
        
        if(All_mother_daught_1BG==False):
            TimeMagL = np.array(filter(lambda row:row[1]>=(mag_th), TimeMagLT) )
        elif(All_mother_daught_1BG==True):
            #S CURVE: TimeMagL = np.array(filter(lambda row:( row[0]>=(mag_th) and row[1]>=(mag_th) ), TimeMagLT) )
            TimeMagL = np.array(filter(lambda row:1.0<=row[0]<=np.max(TimeMagLT),TimeMagLT))#np.max(TimeMagLT)
    #                                TimeMagL = np.array(filter(lambda row:row[0]>=(mag_th),TimeMagLT))
    
    #                            print'Filtered CAT length: ',len(TimeMagL)
    #                            print'length full CAT: ',len(TimeMagLT)
    else:
        print 'MAG_CORR_CATALOG FOR SC MAG_CORR_CATALOG  MAG_CORR_CATALOG\n'
        
        if (INIT==0):
    
            catNumb = 'hs_1981_2016'
            print 'PLOTS FOR MAG CORRELATIOSN OF THE SC CATALOG: \n', catNumb
            TimeMagLT = np.loadtxt('%s/%s.txt'%(cwd,catNumb), usecols=(0, 10) )
            TimeMagL = np.array(filter(lambda row:row[1]>=mag_th, TimeMagLT))
            INIT=1
            
        else:
            pass
        
        TimeMagL = np.array(filter(lambda row:row[0]>=mag_th, TimeMagLT))
        
      
    print 'Mag_Corr THRESHOLD VALUE:',mag_th
    print '\n'
    print 'Length Full CAT: ',len(TimeMagLT)
    print 'Min mag Full CAT:',(np.min(TimeMagLT[:,1]))
    print 'length of CAT after m_th > %s'%mag_th, len(TimeMagL)
    print 'Min mag TimeMagL:',(np.min(TimeMagL[:,1]))


    return(INIT,First_Load,TimeMagL)
