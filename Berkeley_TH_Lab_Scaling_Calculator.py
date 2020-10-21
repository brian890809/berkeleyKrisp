#Kairos Reactor Integral Scaling Program (KRISP), Created by Stephen Louria
#This code generates integral effects test facility design specs based on prototype design \
#    information using the scaling methodology of the scaling topical report for KP-1.
# To execute, you need krisp_input.csv in the same working folder as this code.
#------------------------------------------------------------------------------------------------
import math as ma
import pandas as pd
import numpy as np
import scipy.constants
import csv
from coolant_params import * #import coolant constuctor

#INPUT CARD

df_pump = pd.read_csv('pump_trip_csv.csv', header = 0) #pump trip output file
df_dhx = pd.read_csv('DRACS_DHX_csv.csv', header = 0) #dhx output file

numColumns_df_pump = len(df_pump.columns)
numRows_df_pump = len(df_pump)

numColumns_df_dhx = len(df_dhx.columns)
numRows_df_dhx = len(df_dhx)

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


# time_steps_dhx = list(df_dhx['time']) #making time-steps into a list
# time_steps_dhx = time_steps_dhx[1:]
time_steps_dhx_len = int(len(time_steps_dhx))
# time_steps_dhx_len = int(time_steps_dhx_len)

#OUTPUT CARD

## future : change add more parameters to var_out_param
var_out_param = ['title column', 'time', 'pr_in','pr_out', 'pr_avg', 'richardson', 're_in', 're_out', 'fr_in', 'fr_out'] 
var_out = pd.DataFrame('-', index=list(range(time_steps_dhx_len-1)), columns=var_out_param) #create dataframe output variables
var_out.columns = var_out_param


df_dhx = df_dhx.set_index('title column') #important for scaling of the prototype dimensions, so each row can be called based on name.
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
gravity = scipy.constants.g
l_ht = 1
#constants:

for t in list(range(time_steps_dhx_len-1)): 
    

   	#DHX heater in/ heater out temp
    Heater_in_T = np.subtract(float(df_dhx.iloc[t+1, int('5')]), 273.15)
    Heater_out_T = np.subtract(float(df_dhx.iloc[t+1, int('6')]), 273.15)
    Heater_avg_T = np.divide(np.add(Heater_in_T, Heater_out_T), 2)

    #DHX velocity in heater
    v_ref = float(df_dhx.iloc[t+1, int('44')])
    pipe_1_v_out = float(df_dhx.iloc[t+1, int('44')])

    #DowTherm A parameters, future : change to flibe([T]) if needed, refer to coolant_params.py

    cp_in = dowtherm_A(Heater_in_T).cp #np.add(1518, np.multiply(2.82,Heater_in_T))
    miu_in = dowtherm_A(Heater_in_T).miu #np.divide(0.13, np.power(Heater_in_T,1.072))
    k_in = dowtherm_A(Heater_in_T).k #np.subtract(0.142, np.multiply(0.00016, Heater_in_T))
    rho_in = dowtherm_A(Heater_in_T).rho #np.subtract(1078, np.multiply(0.85, Heater_in_T))

    cp_out =  dowtherm_A(Heater_out_T).cp #np.add(1518, np.multiply(2.82,Heater_out_T))
    miu_out = dowtherm_A(Heater_out_T).miu #np.divide(0.13, np.power(Heater_out_T,1.072))
    k_out = dowtherm_A(Heater_out_T).k  #np.subtract(0.142, np.multiply(0.00016, Heater_out_T))
    rho_out = dowtherm_A(Heater_out_T).rho #np.subtract(1078, np.multiply(0.85, Heater_out_T))

    rho_mid = np.subtract(1078, np.multiply(0.85, Heater_avg_T))
    delta_rho = np.subtract(rho_out, rho_in)

    #Pr number
    pr_in = np.multiply(cp_in, np.divide(miu_in, k_in))
    pr_out = np.multiply(cp_out, np.divide(miu_out, k_out))
    pr_avg = np.divide(np.add(pr_in, pr_out), 2)

    #Richardson number
    rh = np.negative(np.divide(np.multiply(gravity, np.multiply(delta_rho, l_ht)), np.multiply(rho_mid, np.power(v_ref, 2))))

    #Reynold number UL/miu U velocity scale, L fluid characteristic scale
    re_in = np.divide(np.multiply(rho_in, np.multiply(v_ref, 1.6383)), miu_in) #1.6383 is the length of pipe_1 (rho * v_ref * length) / viscosity 
    re_out = np.divide(np.multiply(rho_out, np.multiply(pipe_1_v_out, 1.6383)), miu_out) #reynold number for output ar pipe one

    #Froude number U^2/gL
    fr_in = np.divide(np.power(v_ref, 2), np.multiply(gravity, 1.6383))
    fr_out = np.divide(np.power(pipe_1_v_out, 2), np.multiply(gravity, 1.6383))
#------------------------------------------------------------------------------------------------
#OUTPUTS
    var_out.iloc[t, 0] = time_steps_dhx[1:][t]
    var_out.iloc[t, 1] = pr_in
    var_out.iloc[t, 2] = pr_out
    var_out.iloc[t, 3] = pr_avg
    var_out.iloc[t, 4] = rh
    var_out.iloc[t, 5] = re_in
    var_out.iloc[t, 6] = re_out
    var_out.iloc[t, 7] = fr_in
    var_out.iloc[t, 8] = fr_out
    ### future : add change more var_out.iloc[t, "no. in int"] = parameter
var_out.to_csv("test_output.csv", sep=',',index=False) #change index=True to see indexing column in output card.
print('Done.')
#------------------------------------------------------------------------------------------------
#END