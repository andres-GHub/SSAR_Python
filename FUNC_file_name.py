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

def file_name(Background,dt_type,X_NEAREST_Dt,percent,Full_CAT_type,time_threshold,reps,catNumb,shuffle_type,shuffle_mode,\
              L_Shuffle,N_shuffles,thresholding,time_threshold_print,Trigg_Mag,single_tree,divisor_max_mix,\
              All_mother_daught_1BG,Mag_1BG):
    if (Background=='1BG'and dt_type=='x<Dt<y'):
        inside_file_name='1BG/%s/M_%s_%sNEAREST_%spercent'%(dt_type,Mag_1BG,X_NEAREST_Dt,percent*100)
    
    elif (Background=='1BG'and dt_type=='Dt<y'):
        inside_file_name='1BG/%s/M_%s_%sNEAREST_%spercent'%(dt_type,Mag_1BG,X_NEAREST_Dt,percent*100)
    
    elif(Background=='1BG'):
        inside_file_name='1BG'
    
    elif(Background=='Full'):
        inside_file_name='Full/%s'%Full_CAT_type    
    ###########################
        
    if (reps != -1):
        File_Label = '%s||%s_%s_%s_%sShfl_%sShfls_Thresh:%s_dt:%s_Trigg_Mag:%.2f_Background:Sing_tree:%s_div_max_mix:%.1f_REP:%d'\
        %(catNumb,Full_CAT_type,shuffle_type,shuffle_mode,L_Shuffle,N_shuffles,thresholding.\
          time_threshold,Trigg_Mag,single_tree,divisor_max_mix,reps)
    
    if (reps == -1 and inside_file_name!=-1):
        File_Label = '%s||_%s_%s_%sShfl_%sShfls_Thresh:%s_dt:%s_%sNEAREST_Trigg_Mag:%.2f_Background:Sing_tree:%s_div_max_mix:%.1f'\
        %(catNumb,shuffle_type,shuffle_mode,L_Shuffle,N_shuffles,thresholding,\
          time_threshold_print,X_NEAREST_Dt,Trigg_Mag,single_tree,divisor_max_mix)
        
    elif(reps==-1 and inside_file_name!=-1 and All_mother_daught_1BG==True):
        File_Label = '%s||All_MD:_%s'%(catNumb,All_mother_daught_1BG)
        
    print 'File_Label\n',File_Label
    
    return (File_Label,inside_file_name)
