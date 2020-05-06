#Kairos Reactor Integral Scaling Program (KRISP), Created by Stephen Louria
#This code generates integral effects test facility design specs based on prototype design \
#    information using the scaling methodology of the scaling topical report for KP-1.
# To execute, you need krisp_input.csv in the same working folder as this code.
#------------------------------------------------------------------------------------------------
import math as ma
import pandas as pd
import numpy as np
import csv

#INPUT CARD

df_pump = pd.read_csv('pump_trip_csv.csv', header = 0) #pump trip output file
df_dhx = pd.read_csv('DRACS_DHX_csv.csv', header = 0) #dhx output file

numColumns_df_pump = len(df_pump.columns)
numRows_df_pump = len(df_pump)

numColumns_df_dhx = len(df_dhx.columns)
numRows_df_pump = len(df_pump)

# df_pump = pd.DataFrame(data = df_pump)
# time_steps_pump = list(df_pump['time']) #making time-steps into a list
# time_steps_pump_len = len(time_steps_pump)

# df_dhx = pd.DataFrame(data = df_dhx)

### future : change to pump if needed
time_steps_dhx = []
with open('DRACS_DHX_csv.csv', newline = '') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        time_steps_dhx += [row['title column']]
### future : change to pump if needed


# time_steps_dhx = list(df_dhx['time']) #making time-steps into a list
# time_steps_dhx = time_steps_dhx[1:]
time_steps_dhx_len = len(time_steps_dhx)
time_steps_dhx_len = int(time_steps_dhx_len)

#OUTPUT CARD

## future : change add more parameters to var_out_param
var_out_param = ['title column', 'time', 'pr_in','pr_out'] 

var_out = pd.DataFrame('-', index=list(range(time_steps_dhx_len-1)), columns=var_out_param) #create dataframe output variables
var_out.columns = var_out_param


df_dhx = df_dhx.set_index('title column') #important for doing the scaling of the prototype dimensions, so each row can be called based on name.
var_out = var_out.set_index('title column')

# NAMING CONVENTION
# p = prototype
# m = model
# r = ratio
# FC = Forced Circulation
# NC = Natural Circulation
# SS = Steady State
# rx = core
# hx = Intermediate Heat exchanger
# dc = Downcomer
# rf = Reflector
# Pandas uses [row#,column#] convention

#SCALING CHOICES, NORMAL OPERATIONS (t=0)
#  1. Force Pr_r = 1, find T_m.
#  2. (T_h_m-T_c_m) is determined.
#  3. Pebble length ratio is chosen
#  4. length and area scaling, Ri and Eu are matched, result in pump head and u_ref_m
#  5. matching F_ratio, integrated loop flow resistance
#  6. matching St#, Re# for loop components (matching heat transfer)
#
#SCALING CHOICES, NATURAL CIRCULATION
#  7. Match Ri Number for natural circulation
#  8a. Match solid structures heat capacity number
#  8b. Component level heat capacity matching
#  9. Match loop geometry number G
#  10. Matching Re, St, Gr for loop components
#  11. Matching Bi# for loop components
#  12. Match solid time constants
#  13. deriving power ratio for core
#------------------------------------------------------------------------------------------------
# SS = 'SS' #krisp og code


for t in list(range(time_steps_dhx_len-1)): 

   	#DHX heater in/ heater out temp
    Heater_in_T = float(df_dhx.iloc[t+1, int('5')])
    Heater_out_T = float(df_dhx.iloc[t+1, int('6')])

    cp_in = np.add(1518, np.multiply(2.82,Heater_in_T))
    miu_in = np.divide(0.13, np.power(Heater_in_T,1.072))
    k_in = np.subtract(0.142, np.multiply(0.00016, Heater_in_T))
    pr_in = np.multiply(cp_in, np.divide(miu_in, k_in))

    cp_out = np.add(1518, np.multiply(2.82,Heater_out_T))
    miu_out = np.divide(0.13, np.power(Heater_out_T,1.072))
    k_out = np.subtract(0.142, np.multiply(0.00016, Heater_out_T))
    pr_out = np.multiply(cp_out, np.divide(miu_out, k_out))

#------------------------------------------------------------------------------------------------
#OUTPUTS
    var_out.iloc[t, 0] = time_steps_dhx[1:][t]
    var_out.iloc[t, 1] = pr_in
    var_out.iloc[t, 2] = pr_out
    ### future : add change more var_out.iloc[t, "no. in int"] = parameter
    # var_out.iloc[53,0] = 'Pr_r' #put the results into the dataframe, and calculate distortions.
    # var_out.iloc[53,t] = Pr_r
    # var_out.iloc[88,t] = (Pr_p - Pr_m)/Pr_p

    # var_out.iloc[54,0] = 'F_r'
    # var_out.iloc[54,t] = F_r
    # var_out.iloc[89,t] = (F_p - F_m)/F_p

    # KRISP change in value by circulation methods
    # if t==1: #forced circulation, normal ops
    #     var_out.iloc[55,0] = 'Eu_FC_r'
    #     var_out.iloc[55,t] = Eu_FC_r
    #     var_out.iloc[90,t] = (Eu_FC_p - Eu_FC_m)/Eu_FC_p

    #     var_out.iloc[56,0] = 'Ri_FC_r'
    #     var_out.iloc[56,t] = Ri_FC_r
    #     var_out.iloc[91,t] = (Ri_FC_p - Ri_FC_m)/Ri_FC_p

    #     var_out.iloc[57,0] = 'St_FC_rx_r'
    #     var_out.iloc[57,t] = St_FC_rx_r
    #     var_out.iloc[92,t] = (St_FC_rx_p - St_FC_rx_m)/St_FC_rx_p

    #     var_out.iloc[58,0] = 'St_FC_dc_r'
    #     var_out.iloc[58,t] = St_FC_dc_r
    #     var_out.iloc[93,t] = (St_FC_dc_p - St_FC_dc_m)/St_FC_dc_p

    #     var_out.iloc[59,0] = 'Re_FC_dc_r'
    #     var_out.iloc[59,t] = Re_FC_dc_r
    #     var_out.iloc[94,t] = (Re_FC_dc_p - Re_FC_dc_m)/Re_FC_dc_p

    #     var_out.iloc[60,0] = 'St_FC_hx_r'
    #     var_out.iloc[60,t] = St_FC_hx_r
    #     var_out.iloc[95,t] = (St_FC_hx_p - St_FC_hx_m)/St_FC_hx_p

    #     var_out.iloc[61,0] = 'Re_FC_hx_r'
    #     var_out.iloc[61,t] = Re_FC_hx_r
    #     var_out.iloc[96,t] = (Re_FC_hx_p - Re_FC_hx_m)/Re_FC_hx_p

    #     var_out.iloc[62,0] = 'Re_FC_rx_r'
    #     var_out.iloc[62,t] = Re_FC_rx_r
    #     var_out.iloc[97,t] = (Re_FC_rx_p - Re_FC_rx_m)/Re_FC_rx_p

    #     var_out.iloc[63,0] = 't_FC_r'
    #     var_out.iloc[63,t] = t_FC_r
    #     var_out.iloc[98,t] = '-'

    #     var_out.iloc[64,0] = 'Bi_FC_dc_r'
    #     var_out.iloc[64,t] = Bi_FC_dc_r
    #     var_out.iloc[99,t] = '-'

    #     var_out.iloc[65,0] = 'Bi_FC_rf_r'
    #     var_out.iloc[65,t] = Bi_FC_rf_r
    #     var_out.iloc[100,t] = '-'

    #     var_out.iloc[66,0] = 'Bi_FC_rx_r'
    #     var_out.iloc[65,t] = Bi_FC_rx_r
    #     var_out.iloc[101,t] = '-'

    #     var_out.iloc[67,0] = 'St_FC_rb_r'
    #     var_out.iloc[67,t] = St_FC_rb_r
    #     var_out.iloc[102,t] = (St_FC_rb_p - St_FC_rb_m)/St_FC_rb_p

    #     var_out.iloc[68,0] = 'Re_FC_rb_r'
    #     var_out.iloc[68,t] = Re_FC_rb_r
    #     var_out.iloc[103,t] = (Re_FC_rb_p - Re_FC_rb_m)/Re_FC_rb_p

    # else: #natural circulation only
    #     var_out.iloc[69,0] = 'Ri_NC_r'
    #     var_out.iloc[69,t] = Ri_NC_r
    #     var_out.iloc[104,t] = (Ri_NC_p - Ri_NC_m)/Ri_NC_p

    #     var_out.iloc[70,0] = 'St_NC_rx_r'
    #     var_out.iloc[70,t] = St_NC_rx_r
    #     var_out.iloc[105,t] = (St_NC_rx_p - St_NC_rx_m)/St_NC_rx_p

    #     var_out.iloc[71,0] = 'St_NC_dc_r'
    #     var_out.iloc[71,t] = St_NC_dc_r
    #     var_out.iloc[106,t] = (St_NC_dc_p -St_NC_dc_m)/St_NC_dc_p

    #     var_out.iloc[72,0] = 'Gr_NC_dc_r'
    #     var_out.iloc[72,t] = Gr_NC_dc_r
    #     var_out.iloc[107,t] = (Gr_NC_dc_p - Gr_NC_dc_m)/Gr_NC_dc_p

    #     var_out.iloc[73,0] = 'Re_NC_dc_r'
    #     var_out.iloc[73,t] = Re_NC_dc_r
    #     var_out.iloc[108,t] = (Re_NC_dc_p - Re_NC_dc_m)/Re_NC_dc_p

    #     var_out.iloc[74,0] = 'Re_NC_rx_r'
    #     var_out.iloc[74,t] = Re_NC_rx_r
    #     var_out.iloc[109,t] = (Re_NC_rx_p - Re_NC_rx_m)/Re_NC_rx_p

    #     var_out.iloc[75,0] = 't_NC_r'
    #     var_out.iloc[75,t] = t_NC_r
    #     var_out.iloc[110,t] = '-'

    #     var_out.iloc[76,0] = 'Bi_NC_dc_r'
    #     var_out.iloc[76,t] = Bi_NC_dc_r
    #     var_out.iloc[111,t] = (Bi_NC_dc_p - Bi_NC_dc_m)/Bi_NC_dc_p

    #     var_out.iloc[77,0] = 'Bi_NC_rf_r'
    #     var_out.iloc[77,t] = Bi_NC_rf_r
    #     var_out.iloc[112,t] = (Bi_NC_rf_p - Bi_NC_rf_m)/Bi_NC_rf_p

    #     var_out.iloc[78,0] = 'Bi_NC_rx_r'
    #     var_out.iloc[78,t] = Bi_NC_rx_r
    #     var_out.iloc[113,t] = (Bi_NC_rx_p - Bi_NC_rx_m)/Bi_NC_rx_p

    #     var_out.iloc[79,0] = 'St_NC_rb_r'
    #     var_out.iloc[79,t] = St_NC_rb_r
    #     var_out.iloc[114,t] = (St_NC_rb_p -St_NC_rb_m)/St_NC_rb_p

    #     var_out.iloc[80,0] = 'Gr_NC_rb_r'
    #     var_out.iloc[110,t] = Gr_NC_rb_r
    #     var_out.iloc[115,t] = (Gr_NC_rb_p - Gr_NC_rb_m)/Gr_NC_rb_p

    #     var_out.iloc[81,0] = 'Re_NC_rb_r'
    #     var_out.iloc[81,t] = Re_NC_rb_r
    #     var_out.iloc[116,t] = (Re_NC_rb_p - Re_NC_rb_m)/Re_NC_rb_p

    #     var_out.iloc[82,0] = 'Fo_NC_rx_r'
    #     var_out.iloc[82,t] = Fo_NC_rx_r
    #     var_out.iloc[117,t] = '-'

    #     var_out.iloc[83,0] = 'Fo_NC_rf_r'
    #     var_out.iloc[83,t] = Fo_NC_rf_r
    #     var_out.iloc[118,t] = '-'

    #     var_out.iloc[84,0] = 'Fo_NC_dc_r'
    #     var_out.iloc[84,t] = Fo_NC_dc_r
    #     var_out.iloc[119,t] = '-'

    #the rest, FC or NC
    # var_out.iloc[85,0] = 'G_r'
    # var_out.iloc[85,t] = G_r
    # var_out.iloc[120,t] = (G_p - G_m)/G_p

    #scale_out.loc[14,'Scaling_Group'] = 'HR_r' #the heat source number is currently not used in the analysis.
    #scale_out.loc[14,t] = HR_r
    #dist_out.loc[14,t] = (HR_p - HR_m)/HR_p

# var_out.iloc[88:125,0] = var_out.iloc[53:90,0] #filling the first column of the distortions table
# var_out.iloc[122,0] = '-' #making a mis-done label due to bug dissapear from the output card

#dim_out.to_csv("model_outputs.csv", sep=',',index=False) #write the outputs to csv files
#scale_out.to_csv("ratios_outputs.csv", sep=',',index=False)
#dist_out.to_csv("distortions_outputs.csv", sep=',',index=False)
var_out.to_csv("pr_test_output.csv", sep=',',index=False) #change index=True to see indexing column in output card.
print('Done.')
#------------------------------------------------------------------------------------------------
#END










































#-
