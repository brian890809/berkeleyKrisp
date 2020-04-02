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
df = pd.read_csv('krisp_input.csv') #dataframe for prototype inputs
numColumns_df = len(df.columns)
numRows_df = len(df)

#OUTPUT CARD
dim_out = pd.DataFrame(0, index=range(numRows_df), columns=range(numColumns_df)) #create dataframe for model outputs
#dim_out.columns=df.columns # give all of the columns the same name
#dim_out.loc[0]=df.loc[0] # give it same row of acronyms
#dim_out['row index']=df['row index'] #first column is the list of variables, same for both

scale_out = pd.DataFrame(0, index=range(34), columns=range(10)) #create dataframe for scaling groups outputted, 16 ratios, 10 time stamps.
scale_out.columns = ['Scaling_Group', 't=SS','t=0','t=1','t=10','t=100','t=1000','t=10000','t=100000','t=600000']
numColumns_scale_out = len(scale_out.columns)
numRows_scale_out = len(scale_out)
for i in range(1,numColumns_scale_out):  #fill the table with '-' (make distinction between unintentional 0 and intentional 0)
    for j in range(1,numRows_scale_out):
        scale_out.iloc[j,i] = '-'

dist_out = pd.DataFrame(0, index=range(34), columns=range(10)) #create dataframe for scaling groups distortions outputted, 16 ratios, 10 time stamps.
dist_out.columns = ['Scaling_Group', 't=SS','t=0','t=1','t=10','t=100','t=1000','t=10000','t=100000','t=600000']
numColumns_dist_out = len(dist_out.columns)
numRows_dist_out = len(dist_out)
for i in range(1,numColumns_dist_out):  #fill the table with '-' (make distinction between unintentional 0 and intentional 0)
    for j in range(1,numRows_dist_out):
        dist_out.iloc[j,i] = '-'

var_out = pd.DataFrame(0, index=range(125), columns=range(31)) #create dataframe output variables
var_out.columns = ['title column', 'time', 'qsource_m','q_RVACS_m','deltaT_m','T_h_m','T_c_m','deltaT_sf_m','u_ref_m', \
'thermal_diff_rv_m','thermal_diff_bl_m','thermal_diff_rx_m','F_m','f_m','ks_rx_m','ks_rv_m','ks_bl_m','ks_rf_m',\
'ehAs_dc_m','ehAs_rx_m','ehAs_hx_m','h_NC_rx_m','h_NC_dc_m','h_NC_rf_m','h_NC_bl_m','ls_rx_m','ls_rv_m','ls_rf_m','ls_bl_m'\
, 'EC_rx_hx_r', 'EC_rx_dc_r']

for i in range(0,31): #columns
    for j in range(0,125): #rows
        if var_out.iloc[j,i] == 0:
            var_out.iloc[j,i] = '-'

#Naming the rows for the output file.
var_out.loc[0,'title column'] = 't=SS'
var_out.loc[1,'title column'] = 't=0'
var_out.loc[2,'title column'] = 't=1'
var_out.loc[3,'title column'] = 't=10'
var_out.loc[4,'title column'] = 't=100'
var_out.loc[5,'title column'] = 't=1000'
var_out.loc[6,'title column'] = 't=10000'
var_out.loc[7,'title column'] = 't=100000'
var_out.loc[8,'title column'] = 't=600000'
var_out.loc[14,'time'] = 'Model Dimensions'
var_out.loc[52,'time'] = 'Scaling Ratios'
var_out.loc[87,'time'] = 'Scaling Distortions'

var_out.loc[15:51,'time'] = df.loc[15:51,'title column'] #naming first column of model dimensions table.

var_out.iloc[14,2] = df.iloc[14,2] #naming top row of model dimensions table.
var_out.iloc[14,3] = df.iloc[14,3]
var_out.iloc[14,4] = df.iloc[14,4]
var_out.iloc[14,5] = df.iloc[14,5]
var_out.iloc[14,6] = df.iloc[14,6]
var_out.iloc[14,7] = df.iloc[14,7]
var_out.iloc[14,8] = df.iloc[14,8]
var_out.iloc[14,9] = df.iloc[14,9]
var_out.iloc[14,10] = df.iloc[14,10]
var_out.iloc[14,11] = df.iloc[14,11]

var_out.loc[0,'time'] = 'SS' #naming first column  of variables table.
var_out.loc[1,'time'] = 't=0'
var_out.loc[2,'time'] = 't=1'
var_out.loc[3,'time'] = 't=10'
var_out.loc[4,'time'] = 't=100'
var_out.loc[5,'time'] = 't=1000'
var_out.loc[6,'time'] = 't=10000'
var_out.loc[7,'time'] = 't=100000'
var_out.loc[8,'time'] = 't=600000'

var_out.iloc[52,2] = 'SS' #naming top row of ratios table.
var_out.iloc[52,3] = 't=0'
var_out.iloc[52,4] = 't=1'
var_out.iloc[52,5] = 't=10'
var_out.iloc[52,6] = 't=100'
var_out.iloc[52,7] = 't=1000'
var_out.iloc[52,8] = 't=10000'
var_out.iloc[52,9] = 't=100000'
var_out.iloc[52,10] = 't=600000'

var_out.iloc[87,2] = 'SS' #naming top row of distortions table.
var_out.iloc[87,3] = 't=0'
var_out.iloc[87,4] = 't=1'
var_out.iloc[87,5] = 't=10'
var_out.iloc[87,6] = 't=100'
var_out.iloc[87,7] = 't=1000'
var_out.iloc[87,8] = 't=10000'
var_out.iloc[87,9] = 't=100000'
var_out.iloc[87,10] = 't=600000'

df = df.set_index('title column') #important for doing the scaling of the prototype dimensions, so each row can be called based on name.
var_out = var_out.set_index('title column')

 #convert all the input dimensions tablulated values to floats instead of strings.
for i in range(1,11): #columns
    for j in range(15,51): #rows
        df.iloc[j,i] = float(df.iloc[j,i])

#convert all the  tablulated values in the input card to floats instead of the default strings.
for i in range(0,34): #columns
    for j in range(1,10): #rows
        df.iloc[j,i] = float(df.iloc[j,i])
for i in range(0,1): #columns
    for j in range(53,numRows_df): #rows
        df.iloc[j,i] = float(df.iloc[j,i])

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
SS = 'SS'
for t in np.geomspace(1,1e8, num=9): #over the transient, loop over time increments (over time stamp 9,0s,1s,10s,100s,1000s,1e4,1e5,6e5)
    if t==1: #need this if statement section because function geomspace cannot start at 0
        t='SS'
    elif t==10:
        t=0
    elif t==1e8:
        t=6e5 #this is the largest time elapsed that we have data for.
    else:
        t=t/100

    if t == 'SS':
        pass
    else:
        t = int(t)  #convert to integer so it doesn't have the '.0' at the end.
    t = str(t)  #need it to be a string for calling using df.loc function of the PANDAS dataframe function set.

    #DETERMINING CHARACTERISTIC RATIOS, AND MODEL-TO-PROTOTYPE RATIOS
    #1 // First step is matching Pr#
    Pr_p = df.loc[t,'8'] / df.loc[t,'9'] # Pr Number at time t
    Pr_r = 1.0 #forcing this to match
    Pr_m = Pr_p*Pr_r

    # for each time t, we have a Pr_m requirement from above. I used table in 'Salt and simulant fluids thermophysical properties.excel' to select a T_mean such that Pr# = 1.
    # see page 38 of scaling report.
    T_m = df.loc[t,'33'] #in K

    t = str(t) #need it to be a string again for using PANDAS dataframe functions like df.loc.

    # MISC PROPERTIES
    deltaP_FC_p = df.loc[t,'13'] #prototypical pump head
    ks_rx_p = df.loc[t,'14'] #from IET_workbook for pebble average
    ks_rx_m = ks_rx_p #for now, I assume we will match this using teflon blended with metal. Can improve power requirements by increasing ks_rx_m
    ks_rx_r = ks_rx_m/ks_rx_p
    rho_rx_p = df.loc[t,'15'] #material properties of fuel pebble
    cp_rx_p = df.loc[t,'16'] #material properties of fuel pebble
    ls_rx_p = df.loc['ls_rx_p','0'] #pebble diameter
    g = df.loc['g','0'] #[m/s2] acceleration of gravity at earth's surface.
    qsource_p = df.loc[t,'30'] #[W] follows decay heat curve.
    deltaT_sf_p = df.loc[t,'10']
    deltaP_p = df.loc[t,'13']
    u_ref_FC_p = df.loc['SS','11'] # [m/s] reference velocity for FC from KP-SAM
    u_ref_NC_p = df.loc[t,'11'] # reference velocity for NC from KP-SAM
    u_ref_FC_hx_p = df.loc[t,'31'] #from KP-SAM data
    u_ref_FC_hx_m = df.loc['u_ref_FC_hx_m','0'] #balance mass flow rate in model
    u_ref_FC_hx_r = u_ref_FC_hx_m/u_ref_FC_hx_p

    por_p = df.loc['por_p','0'] #not needed for now
    ks_rf_p = df.loc['ks_rf_p','0'] #reflector thermal conductivity, ETU-10. Need to add temp dependence when it becomes known.
    ks_rf_m = df.loc['ks_rf_m','0'] #can input this once it is known based on material selection.
    ks_rf_r = ks_rf_m/ks_rf_p
    ls_rf_p = df.loc['Outer_Reflector_Active','6'] #width of reflector next to active core.
    ls_bl_p = df.loc['Barrel_Middle','6'] #width of core barrel.
    h_rf_p = df.loc['h_rf_p','0'] #need convection from flibe to reflector from KP-SAM
    ks_rv_p = df.loc[t,'22'] #vessel SS316 is the solid used
    #ks_rv_m = 1.4
    ks_rv_m = df.loc['ks_rv_m','0'] #W/mK,Pyrex @ 126C. Value from 'Thermal conductivity of pyrex glass: selected values', Carwile 1966.
    ks_bl_m = df.loc['ks_bl_m','0'] #W/mK,Pyrex is assumed.
    #ks_rv_m = 16.2 #W/mK,when SS is assumed.
    rho_rx_m = df.loc['rho_rx_m','0'] #will come from material selection process
    cp_rx_m = df.loc['cp_rx_m','0'] #will come from material selection process

    ks_rv_r = ks_rv_m/ks_rv_p #assume we can match this.

    # HEAT TRANFER OIL PROPERTIES [Downtherm A], dependent on fixed mean temperatures from above.
    rho_m = df.loc[t,'25']
    visc_m = df.loc[t,'26']
    k_m = df.loc[t,'28']
    cp_m = df.loc[t,'27']
    beta_m = df.loc[t,'29']
    kin_visc_m = df.loc[t,'26'] / df.loc[t,'25']
    therm_diff_m = df.loc[t,'28']/(df.loc[t,'25']*df.loc[t,'27'])

    #FLIBE PROPERTIES
    rho_p = df.loc[t,'4']
    visc_p = df.loc[t,'5']
    k_p = df.loc[t,'6']
    cp_p = df.loc[t,'7']
    beta_p = df.loc[t,'12']
    kin_visc_p = df.loc[t,'8']
    therm_diff_p = df.loc[t,'9']
    T_h_p = df.loc[t,'0']
    T_c_p = df.loc[t,'1']
    deltaT_p = df.loc[t,'3']
    T_p = df.loc[t,'2']
    beta_r = beta_m/beta_p
    rho_r = rho_m/rho_p
    cpf_r = cp_m/cp_p
    visc_r = visc_m/visc_p
    k_r = k_m/k_p

    #2  // then we match thermal expansion of fluid (buoyancy) effects by scaling the deltaT to arrive at T_hot and T_cold temps across core.
    deltaT_m = beta_p*(T_h_p-T_c_p)/beta_m
    deltaT_r = deltaT_m/deltaT_p
    T_h_m = T_m + 0.5*deltaT_m
    T_c_m = T_m - 0.5*deltaT_m
    T_h_r = T_h_m/T_h_p
    T_c_r = T_c_m/T_c_p
    deltaT_sf_r = deltaT_r
    deltaT_sf_m = deltaT_sf_p*deltaT_sf_r
    T_r = T_m/T_p #temp ratio based on average temps of flibe and oil
    var_out.loc['t='+t,'deltaT_m'] = deltaT_m
    var_out.loc['t='+t,'T_h_m'] = T_h_m
    var_out.loc['t='+t,'T_c_m'] = T_c_m
    var_out.loc['t='+t,'deltaT_sf_m'] = deltaT_sf_m

    #3  // pebble length scaling can be done independently from facility-wide length scaling. But this is not used for ITF-0, which has no physical pebbles
    l_peb_p = df.loc['l_peb_p','0'] # Pebble diameter.
    l_peb_r = df.loc['l_peb_r','0']  # currently the same as the height ratio, although this can change to change Re# in pebble bed.
    l_peb_m = l_peb_p*l_peb_r

    if t=='SS': # this whole section is for normal operation (steady-state) only.

        #4a System-level scaling of Ri#, Fr#, Eu#, reference velocity. Here we establish the facility-wide length and area scaling to match these important phenomena.
        length_p = df.loc['length_p','0'] #length of vessel
        length_r = df.loc['length_r','0'] #this is the length scale for all components, free to choose.
        #area_r = 1/11.787 #this is driven by downcomer flow area.
        #area_r = 1/8.89 #this is if DC OD=8.858in
        area_r = df.loc['area_r','0'] #this is 1/35.9, from DC OD=4.72in
        core_scaling_factor = df.loc['core_scaling_factor','0'] #0.4*0.9. The 0.4 comes from the porosity, accounts for the lack of pebbles in the model. 0.9 accounts for smaller volumes of upper and lower core sections, given model is one long pipe.
        A_ref_r = ((length_r)**2)*(area_r) # the ((0.47)**2) comes from 1/2 length for geometric scaling, and area_r for further direct area scaling.
        #print('Check 2 Passed')
        A_ref_p = df.loc['Core_Active_Section','2']
        A_ref_m = A_ref_p*A_ref_r
        n_r = df.loc['n_r','0'] #from scaling report, number of parallel flow channels. may have to update this.

        Ri_FC_p = g*beta_p*(deltaT_p)*length_p/(u_ref_FC_p**2)
        Ri_FC_r = 1 #assume we can match
        Ri_FC_m = Ri_FC_p*Ri_FC_r
        u_ref_FC_r = (length_r)**.5
        u_ref_FC_m = u_ref_FC_p*u_ref_FC_r
        Fr_FC_r = 1 #assume we can match if Ri is matched

        t_FC_r = length_r/u_ref_FC_r #use ref vel of core to determine time scaling for forced circulation

        Eu_FC_p = deltaP_p/(rho_p*(u_ref_FC_p**2)) #used reference vel. of core for Eu#
        Eu_FC_r = 1.0 #matched if pump head ratio (Ho_r) matches length ratio (length_r).
        Eu_FC_m = Eu_FC_p*Eu_FC_r
        Ho_r = length_r
        deltaP_FC_r = rho_r*length_r
        Ho_p = deltaP_FC_p/(rho_p*g) #this is based on the KP-SAM preliminary data for pressure drop head.
        Ho_m = Ho_p*Ho_r
        var_out.loc['t='+t,'u_ref_m'] = u_ref_FC_m

        #6a Core (rx) heat transfer and fluid dynamics scaling

        Dh_rx_p = (df.loc['Core_Active_Section','4'])  # hydraulic diameter in active core region
        Dh_rx_r = (A_ref_r**0.5) # Dh = A/P and P = pi*D so Dh_r = sqrt(A_r)
        Dh_rx_m = Dh_rx_p*Dh_rx_r

        ehAs_rx_r = rho_r*cpf_r*u_ref_FC_r*A_ref_r # eq. 102 of scaling report
        ehAs_rx_p = df.loc['e_rx_p','0']*df.loc['h_rx_p','0']*(2*np.pi*(df.loc['Core_Active_Section','2']))*(df.loc['Core_Active_Section','3']+df.loc['Core_Lower_Section','3']+df.loc['Core_Upper_Section','3']) #e=1,h = 500, assumes entire core has radius of active core (over-estimation of As)
        ehAs_rx_m = ehAs_rx_p*ehAs_rx_r

        St_FC_rx_p = (ehAs_rx_p*deltaT_sf_p/(rho_p*u_ref_FC_p*A_ref_p*cp_p*deltaT_p))
        St_FC_rx_r = 1.0
        St_FC_rx_m = St_FC_rx_p*St_FC_rx_r

        Re_FC_rx_p = rho_p*u_ref_FC_p*l_peb_p/visc_p # Re# for core is a function of pebble diameter, which can be different from system-wide Dh scaling.
        Re_FC_rx_m = rho_m*u_ref_FC_m*l_peb_m/visc_m # although, it will be a function of l_peb, and not l_peb directly. This is an approx. to allow for...
        #... differences between system-wide Dh scaling and pebble diameter scaling
        Re_FC_rx_r = Re_FC_rx_m/Re_FC_rx_p

        #6b Intermediate Heat Exchanger (hx) heat transfer and fluid dynamics scaling
        Dh_hx_r = A_ref_r**0.5 # hydraulic diameter of hx tube
        #Dh_hx_r = .46 #select as independent variable to match Re_FC_hx_r
        Dh_hx_p = (df.loc['Intermediate_HX','4'])  # diameter of HX tube, Dh
        Dh_hx_m = Dh_hx_p*Dh_hx_r

        ehAs_hx_r = rho_r*cpf_r*u_ref_FC_r*n_r*(Dh_hx_r**2)
        ehAs_hx_p = df.loc['e_hx_p','0']*df.loc['h_hx_p','0']*2*np.pi*((df.loc['Intermediate_HX','1']/np.pi)**.5)*(df.loc['Intermediate_HX','3']) # e = 1, h = 500, surface area from circumfrence and length
        ehAs_hx_m = ehAs_hx_p*ehAs_hx_r

        St_FC_hx_p = (ehAs_hx_p*deltaT_sf_p/(rho_p*u_ref_FC_hx_p*(np.pi/4)*(Dh_hx_p**2)*cp_p*deltaT_p))
        St_FC_hx_r = 1.0 #assumed in scaling report by constraining (hAs_rx_r).
        St_FC_hx_m = St_FC_hx_p*St_FC_hx_r

        Re_FC_hx_p = rho_p*u_ref_FC_p*Dh_hx_p/visc_p
        Re_FC_hx_r = rho_r*u_ref_FC_r*Dh_hx_r/visc_r
        Re_FC_hx_m = Re_FC_hx_p*Re_FC_hx_r

        #6c Downcomer (dc) heat transfer and fluid dynamics scaling
        Dh_dc_r = length_r #even though Dh is a function of flow area, in the DC, the channel width length scale drives the flow behavior (Re#).
        Dh_dc_p = (df.loc['Downcomer_Upper_Section','4'])
        Dh_dc_m = Dh_dc_p*Dh_dc_r

        ehAs_dc_r = rho_r*cpf_r*u_ref_FC_r*(Dh_dc_r**2)
        ehAs_dc_p = df.loc['e_dc_p','0']*df.loc['h_dc_p','0']*(2*np.pi*1.876)*(df.loc['Downcomer_Middle_Section','3']+df.loc['Downcomer_Upper_Section','3']+df.loc['Downcomer_Lower_Section','3']+df.loc['Downcomer_Upper_Upper','3'])
        # e=1, h=500,downcomer radius = 1.876 from 'KP-SAM inputs checking'.
        ehAs_dc_m = ehAs_dc_p*ehAs_dc_r

        St_FC_dc_p = (ehAs_dc_p*deltaT_sf_p/(rho_p*u_ref_FC_p*(np.pi/4)*(Dh_dc_p**2)*cp_p*deltaT_p))
        St_FC_dc_r = 1.0 #assumed in scaling report by constraining (hAs_rx_r).
        St_FC_dc_m = St_FC_dc_p*St_FC_dc_r

        Re_FC_dc_p = rho_p*u_ref_FC_p*Dh_dc_p/visc_p
        Re_FC_dc_r = rho_r*u_ref_FC_r*Dh_dc_r/visc_r
        Re_FC_dc_m = Re_FC_dc_p*Re_FC_dc_r

        #6d Reflector Region bypass channel (rb = reflector bypass) heat transfer and fluid dynamics scaling
        Dh_rb_r = A_ref_r**0.5
        Dh_rb_p = (df.loc['Reflector_Bypass','4'])
        Dh_rb_m = Dh_rb_p*Dh_rb_r

        ehAs_rb_r = rho_r*cpf_r*u_ref_FC_r*(Dh_rb_r**2)
        ehAs_rb_p = df.loc['e_rb_p','0']*df.loc['h_rb_p','0']*np.pi*df.loc['Reflector_Bypass','4']*(df.loc['Reflector_Bypass','3']) #e=1,h=500
        ehAs_rb_m = ehAs_rb_p*ehAs_rb_r

        St_FC_rb_p = (ehAs_rb_p*deltaT_sf_p/(rho_p*u_ref_FC_p*(np.pi/4)*(Dh_rb_p**2)*cp_p*deltaT_p))
        St_FC_rb_r = 1.0 #assumed in scaling report by constraining (hAs_rx_r).
        St_FC_rb_m = St_FC_rb_p*St_FC_rb_r

        Re_FC_rb_p = rho_p*u_ref_FC_p*Dh_rb_p/visc_p
        Re_FC_rb_r = rho_r*u_ref_FC_r*Dh_rb_r/visc_r
        Re_FC_rb_m = Re_FC_rb_p*Re_FC_rb_r

    #NATURAL CIRCULATION phenomena, relevant to transient part of tests

    #Right now these are all constant for all time because h_p is constant. If we have h vary in time, than these will change. But this would be very difficult \
    # ... to  change in the experiment.
    var_out.loc['t='+t,'ehAs_dc_m'] = ehAs_dc_m
    var_out.loc['t='+t,'ehAs_hx_m'] = ehAs_hx_m
    var_out.loc['t='+t,'ehAs_rx_m'] = ehAs_rx_m

    #5 Matching F# in the loop.
    # during forced circulation, form loss >> friction loss, so only form factor K is considered to match F.

    sumK = 0.0
    for i in range(15,50):
        if df.iloc[i,1] == 0: #if the component has no area (aka a solid)
            pass
        elif  df.iloc[i,3] == 0: #if the component has no length, same idea...
            pass
        else:
            sumK += 0.5*(df.iloc[i,9]*(A_ref_p/df.iloc[i,1])**2)  # #9 are the Ks, index 1 are the areas. See Table 2 of topical report.

    # There is no F_NC_P because it is assumed given the normal ops matching, aka the matching above for SS can be applied to the transient well (as per the scaling report).

    sumf = 0.0
    for i in range(15,50):
        if df.iloc[i,1] == 0: #if the component has no area (aka a solid)
            pass
        elif  df.iloc[i,3] == 0: #if the component has no length, same idea...
            pass
        else:
            sumf += 0.5*(df.iloc[i,10]*(df.iloc[i,3]/(2.0*df.iloc[i,2]))*(A_ref_p/df.iloc[i,1])**2)  # Row 1 is Areas, #2 is radius, #3 is length.
    f_p = sumf # this is for matching total major loss for the loop
    f_r = 1.0
    f_m = f_p*f_r

    F_p = sumK + sumf # these are split up because currently the topical report only accounts for minor loss since it assumes K>>f during SS.
    F_r = 1.0 #assume it can be matched because we can just adjust the needle valves to match F_p.
    F_m = F_p*F_r

    var_out.loc['t='+t,'F_m'] = F_m
    var_out.loc['t='+t,'f_m'] = f_m

    #7a Matching system-wide nat. circ. behavior using Ri# and time and reference velocity scaling
    Ri_NC_r = 1 # we can assume this if we match the decay curve with our power input into facility during nat. circulation phase (at the reduced scale).
    Ri_NC_p = (g*beta_p*deltaT_p*length_p/((u_ref_NC_p)**2)) #this length scale should be [centerline of DC - centerline of core], which is about 0.
    # also, the temperature difference across the core is assumed to be similar to the normal ops case (deltaT_p)
    # although, this number is never needed, since Ri# matching is assumed already.
    Ri_NC_m = Ri_NC_p*Ri_NC_r
    u_ref_NC_r = u_ref_FC_r #eq 106 of scaling report
    u_ref_NC_m = u_ref_NC_p*u_ref_NC_r
    t_NC_r = t_FC_r #eq 106 of scaling report
    var_out.loc['t='+t,'u_ref_m'] = u_ref_NC_m

    #8a System-wide heat capacity matching
    # this section is dependent on model material selection, so the model materials are assumed. (pyrex for vessel and steel for barrel)

    #Material property data for prototype (SS3016H for vessel/barrel, and ETU-10 for reflector)
    cp_rf_p = df.loc[t,'18']
    cp_rv_p = df.loc[t,'23']
    cp_bl_p = df.loc[t,'23']
    rho_rf_p = df.loc[t,'19']
    rho_rv_p = df.loc[t,'21']
    rho_bl_p = df.loc[t,'21']

    #Additional needed material property data for model (pyrex for vessel/barrel, and ETU-10 for reflector (rough approx))
    cp_rf_m = df.loc['cp_rf_m','0'] # ETU-10 (unimportant for now since it isn't used)
    rho_rf_m = df.loc['rho_rf_m','0'] # ETU-10 (unimportant for now since it isn't used)
    cp_rv_m = df.loc['cp_rv_m','0'] #pyrex
    cp_bl_m = df.loc['cp_bl_m','0'] #stainless steel
    rho_rv_m = df.loc['rho_rv_m','0'] #pyrex
    rho_bl_m = df.loc['rho_bl_m','0'] #stainless steel

    #Calculating the volumes of all of the KP-SAM components based on their dimensions in the KP-SAM input card
    vol_pb_s_p = (df.loc['Core_Active_Section','1'])*(df.loc['Core_Active_Section','3'])*(1-por_p) + (df.loc['Core_Upper_Section','1'])*(df.loc['Core_Upper_Section','3'])*(1-por_p) + (df.loc['Core_Lower_Section','1'])*(df.loc['Core_Lower_Section','3'])*(1-por_p)
    vol_pb_f_p = (df.loc['Core_Active_Section','1'])*(df.loc['Core_Active_Section','3'])*(por_p) + (df.loc['Core_Upper_Section','1'])*(df.loc['Core_Upper_Section','3'])*(por_p) + (df.loc['Core_Lower_Section','1'])*(df.loc['Core_Lower_Section','3'])*(por_p)
    vol_rf_s_p = 2*np.pi*(df.loc['Outer_Reflector_Lower','2'])*(df.loc['Outer_Reflector_Lower','6'])*(df.loc['Outer_Reflector_Lower','3']) + 2*np.pi*(df.loc['Outer_Reflector_Active','2'])*(df.loc['Outer_Reflector_Active','6'])*(df.loc['Outer_Reflector_Active','3']) + \
               2*np.pi*(df.loc['Outer_Reflector_Upper','2'])*(df.loc['Outer_Reflector_Upper','6'])*(df.loc['Outer_Reflector_Upper','3']) + 2*np.pi*(df.loc['Outer_Reflector_Upper_Upper','2'])*(df.loc['Outer_Reflector_Upper_Upper','6'])*(df.loc['Outer_Reflector_Upper_Upper','3'])
    vol_rf_f_p = (df.loc['Pipe_from_Lower_Plenum_to_Active_Core','1'])*(df.loc['Pipe_from_Lower_Plenum_to_Active_Core','3'])+(df.loc['Pipe_from_Active_Core_to_Upper_Plenum','1'])*(df.loc['Pipe_from_Active_Core_to_Upper_Plenum','3'])\
               +(df.loc['Reflector_Bypass','1'])*(df.loc['Reflector_Bypass','3'])
    vol_plenums_s_p = 0.8967 #lower plenum solid accounts for lower vessel head mass, from KP-SAM input card (as V_wall)
    vol_plenums_f_p = (df.loc['Lower_Plenum','5']) + (df.loc['Upper_Plenum','5'])
    vol_bl_s_p = 2*np.pi*(df.loc['Barrel_Lower','2'])*(df.loc['Barrel_Lower','6'])*(df.loc['Barrel_Lower','3'])+ 2*np.pi*(df.loc['Barrel_Middle','2'])*(df.loc['Barrel_Middle','6'])*(df.loc['Barrel_Middle','3']) + 2*np.pi*(df.loc['Barrel_Upper','2'])*(df.loc['Barrel_Upper','6'])*(df.loc['Barrel_Upper','3'])\
               + 2*np.pi*(df.loc['Barrel_Upper_Upper','2'])*(df.loc['Barrel_Upper_Upper','6'])*(df.loc['Barrel_Upper_Upper','3'])
    vol_bl_f_p = 0.0 #there is no fluid associated with barrel component
    vol_vs_s_p = 2*np.pi*(df.loc['Reactor_Vessel_Lower','2'])*(df.loc['Reactor_Vessel_Lower','6'])*(df.loc['Reactor_Vessel_Lower','3'])+ 2*np.pi*(df.loc['Reactor_Vessel_Middle','2'])*(df.loc['Reactor_Vessel_Middle','6'])*(df.loc['Reactor_Vessel_Middle','3'])\
              + 2*np.pi*(df.loc['Reactor_Vessel_Upper','2'])*(df.loc['Reactor_Vessel_Upper','6'])*(df.loc['Reactor_Vessel_Upper','3'])+ 2*np.pi*(df.loc['Reactor_Vessel_Upper_Upper','2'])*(df.loc['Reactor_Vessel_Upper_Upper','6'])*(df.loc['Reactor_Vessel_Upper_Upper','3'])
    vol_vs_f_p = 0.0 #there is no fluid associated with vessel component
    vol_dc_f_p = (df.loc['Downcomer_Lower_Section','1'])*(df.loc['Downcomer_Lower_Section','3']) + (df.loc['Downcomer_Middle_Section','1'])*(df.loc['Downcomer_Middle_Section','3']) + (df.loc['Downcomer_Upper_Section','1'])*(df.loc['Downcomer_Upper_Section','3']) \
               + (df.loc['Downcomer_Upper_Upper','1'])*(df.loc['Downcomer_Upper_Upper','3'])
    vol_dc_s_p = 0.0 #there is no solid associated with downcomer channel. See vessel and barrel components.
    vol_IHXs_s_p = (((df.loc['Intermediate_HX','1'])/np.pi)**.5)*(df.loc['Intermediate_HX','3']) #rough approx, may want to be more accurate based on # tubes.
    vol_IHXs_f_p = (df.loc['Intermediate_HX','1'])*(df.loc['Intermediate_HX','3'])
    vol_ppool_s_p = 1.104 # from KP-SAM input card (as V_wall), this volume of steel accounts for upper vessel head mass
    vol_ppool_f_p = (df.loc['Pump_Pool','5'])
    vol_legs_s_p = 0.0
    vol_legs_f_p = (df.loc['Pipe_from_Intermediate_HX_to_Downcomer','1'])*(df.loc['Pipe_from_Intermediate_HX_to_Downcomer','3']) + (df.loc['Pipe_From_Pump_to_Intermediate_HX','1'])*(df.loc['Pipe_From_Pump_to_Intermediate_HX','3'])

    # prototype total fluid volume
    tot_vol_f_p = vol_pb_f_p + vol_rf_f_p + vol_plenums_f_p + vol_dc_f_p + vol_IHXs_f_p + vol_legs_f_p
    #tot_vol_f_p = 2.93E+01 #value from 'inputs-checking_KP-SAM'

    #Prototype energy absorption capacity of components (energy per degree celcius change) based on volumes calculated above
    energy_core_solid_p = (vol_pb_s_p*rho_rx_p*cp_rx_p)+(vol_rf_s_p*rho_rf_p*cp_rf_p)+(vol_plenums_s_p*rho_rv_p*cp_rv_p) # pebbles + reflector + lower plenum +upper plenum
    energy_core_fluid_p = (vol_pb_f_p*rho_p*cp_p)+(vol_rf_f_p*rho_p*cp_p)+(vol_plenums_f_p*rho_p*cp_p) # pebble bed + ref bypass + lower plenum +upper plenum
    energy_DC_solid_p = (vol_vs_s_p*rho_rv_p*cp_rv_p) + (vol_bl_s_p*rho_bl_p*cp_bl_p) # vessel + barrel
    energy_DC_fluid_p = (vol_dc_f_p*rho_p*cp_p) #DC channel fluid
    energy_IHXs_solid_p = vol_IHXs_s_p*rho_rv_p*cp_rv_p #assume SS316H
    energy_IHXs_fluid_p = (vol_IHXs_f_p*rho_p*cp_p) #fluid in IHXs
    energy_ppool_solid_p = vol_ppool_s_p*rho_rv_p*cp_rv_p #assume SS316H, this accounts for upper vessel head mass
    energy_ppool_fluid_p = vol_ppool_f_p*rho_p*cp_p # assume that fluid is in half pump pool volume only
    energy_legs_solid_p = vol_legs_s_p*rho_rv_p*cp_rv_p #assume SS316H #KP-SAM has adiabatic boundary, no solid walls for hot/cold legs
    energy_legs_fluid_p = (vol_legs_f_p*rho_p*cp_p) #fluid in hot/cold legs

    # model total fluid volume, data from 'KP-SAM Relative Volumes' spreadsheet.
    tot_vol_f_m = 8.82E-02 + 5.87E-02 + 2.46E-02 + 5.98E-02 + 2*(2.76E-03)

    #Model energy absorption capacity of components (energy per degree celcius change) based on volumes and heat capacities
    energy_core_solid_m = (0*rho_rx_m*cp_rx_m)+(0*rho_rf_m*cp_rf_m) # + (wall of core model) # pebbles + reflector + walls
    energy_core_fluid_m = (8.82E-02*rho_m*cp_m) #core fluid
    energy_DC_solid_m = (3.60E-01*rho_rv_m*cp_rv_m) + (7.44E-02*rho_bl_m*cp_bl_m) # vessel + barrel
    energy_DC_fluid_m = (2.46E-02*rho_m*cp_m) #DC channel fluid
    energy_IHXs_solid_m = (2.35*np.pi*0.002546*0.001)*rho_bl_m*cp_bl_m # all IHX tubing estimate, assume length = 2.35m, Dh = 0.002546, pipe thickness = 1mm. Steel matieral (same as barrel)
    energy_IHXs_fluid_m = 2*(5.98E-02*rho_m*cp_m) # two IHXs worth of fluid
    energy_pumppool_solid_m = np.pi*((0.4497*4/(np.pi*0.564))**0.5)*0.002*0.564*rho_bl_m*cp_bl_m # solid volume estimate derived from model output card fluid volume and length, material assumed same as barrel, assume 2mm thick tank
    energy_pumppool_fluid_m = (5.87E-02*rho_m*cp_m) #fluid in pump pool, half full
    energy_plenums_solid_m = np.pi*((0.1275*4/(np.pi*0.24534))**0.5)*0.24534*.002*rho_bl_m*cp_bl_m + np.pi*((0.00438*4/np.pi/0.0658)**0.5)*0.0658*0.002*rho_bl_m*cp_bl_m #wall of plenum tanks, derived from model output card fluid volume and length, assume 2mm thick tank
    energy_legs_solid_m = 0.000152466*rho_rv_m*cp_rv_m + 0.000152466*rho_rv_m*cp_rv_m # hot/cold leg piping solids. Assume hot and cold leg have same length. 1NPS,Sch5. Solid volume = (pi*1.8048/4)*((.0334^2)-((.0334-.001651)^2)) = 0.000152466
    energy_legs_fluid_m = (0.00142883)*rho_m*cp_m + (0.00142883)*rho_m*cp_m #fluid in hot & cold legs. Volume = pi*(Di^2)*L/4 = 0.00142883m3

    sumEnergy = 0
    sumEnergy += energy_core_solid_p + energy_core_fluid_p #core including plenums
    sumEnergy += energy_ppool_solid_p + energy_ppool_fluid_p #pump pool
    sumEnergy += energy_DC_solid_p + energy_DC_fluid_p #DC
    sumEnergy += energy_IHXs_solid_p + energy_IHXs_fluid_p
    sumEnergy += energy_legs_solid_p + energy_legs_fluid_p #Hot/Cold legs
    HS_p = (sumEnergy)/(tot_vol_f_p*rho_p*cp_p)

    sumEnergy = 0
    sumEnergy += energy_core_solid_m + energy_core_fluid_m #Core including plenum tanks
    sumEnergy += energy_pumppool_solid_m + energy_pumppool_fluid_m #pump pool
    sumEnergy += energy_DC_solid_m + energy_DC_fluid_m #DC
    sumEnergy += energy_IHXs_solid_m + energy_IHXs_fluid_m # two IHXs worth of fluid
    sumEnergy += energy_legs_solid_m + energy_legs_fluid_m #Hot and Cold legs, two of each
    HS_m = (sumEnergy)/(tot_vol_f_m*rho_m*cp_m)

    HS_r = HS_m/HS_p #solid structure heat capacity number. We want to match this as per scaling report Eq 39. This is system-wide only.

    # 8b. Component-level heat capacity matching.
    #We try to match relative energy capacity between solid and fluid in both the downcomer and core in this section

    #for prototype

    #ratio of solid/fluid energy capacity in core
    core_sf_ratio_p = energy_core_solid_p / energy_core_fluid_p

    #ratio of solid/fluid energy capacity in DC
    dc_sf_ratio_p = energy_DC_solid_p / energy_DC_fluid_p

    #for model

    #ratio of solid/fluid energy capacity in core
    core_sf_ratio_m = energy_core_solid_m / energy_core_fluid_m

    #ratio of solid/fluid energy capacity in DC
    dc_sf_ratio_m = energy_DC_solid_m / energy_DC_fluid_m

    # model/prototype ratios
    # model-to-prototype ratio of solid/fluid energy capacity ratio for Core
    EC_rx_r = core_sf_ratio_m/core_sf_ratio_p

    # model-to-prototype ratio of solid/fluid energy capacity ratio for Downcomer
    EC_dc_r = dc_sf_ratio_m/dc_sf_ratio_p

    # We also want to try to match relative energy capacity between core/downcomer and core/IHXs.

    #prototype

    #ratio of core/downcomer energy capacity
    EC_rx_dc_p = (energy_core_solid_p + energy_core_fluid_p)/(energy_DC_solid_p + energy_DC_fluid_p)
    #ratio of core/IHXs energy capacity
    EC_rx_hx_p = (energy_core_solid_p + energy_core_fluid_p)/(energy_IHXs_solid_p + energy_IHXs_fluid_p)

    #model

    #ratio of core/downcomer energy capacity
    EC_rx_dc_m = (energy_core_solid_m + energy_core_fluid_m)/(energy_DC_solid_m + energy_DC_fluid_m)
    #ratio of core/IHXs energy capacity
    EC_rx_hx_m = (energy_core_solid_m + energy_core_fluid_m)/(energy_IHXs_solid_m + energy_IHXs_fluid_m)

    #important model/prototype relative solid energy capacity ratios to monitor. These are 1 ideally.
    EC_rx_dc_r = EC_rx_dc_m / EC_rx_dc_p
    EC_rx_hx_r = EC_rx_hx_m / EC_rx_hx_p

    var_out.loc['t='+t,'EC_rx_dc_r'] = EC_rx_dc_r
    var_out.loc['t='+t,'EC_rx_hx_r'] = EC_rx_hx_r

    #9a System-level scaling of loop geometry number, for preserving relative volumes and residence times (in-vessel 'loop')

    #loop is as follows: lower core, active core, upper core (defueling chute), pipe from active core to upper plenum,... //
    #pipe from upper plenum to diode, pipe from diode to downcomer, downcomer, pipe from lower plenum to active core.

    sumLengths = 0.0 #this section is matching loop geometry number. loop is for in-vessel NC, not the whole PHTS.
    sumAreas = 0.0
    sumLA_p = 0.0
    referenceRatio = length_r/A_ref_r

    # row 3 is Lengths. Row 1 is Areas.
    sumLA_p += df.loc['Core_Lower_Section','3']/df.loc['Core_Lower_Section','1'] #core_lower
    sumLA_p += df.loc['Pipe_from_Lower_Plenum_to_Active_Core','3']/df.loc['Pipe_from_Lower_Plenum_to_Active_Core','1'] #pipe from lower plenum to active core
    sumLA_p += df.loc['Core_Active_Section','3']/(df.loc['Core_Active_Section','1']*core_scaling_factor) #active core, must reduce area to account for the pebbles and reflector bypass considerations of the model active core area measurement
    sumLA_p += df.loc['Core_Upper_Section','3']/df.loc['Core_Upper_Section','1'] #upper core
    sumLA_p += df.loc['Pipe_from_Active_Core_to_Upper_Plenum','3']/df.loc['Pipe_from_Active_Core_to_Upper_Plenum','1'] #pipe from active core to upper plenum
    sumLA_p += df.loc['Pipe_from_Pump_Pool_to_Fluidic_Doide','3']/df.loc['Pipe_from_Pump_Pool_to_Fluidic_Doide','1'] #pipe from upper plenum to diode
    sumLA_p += df.loc['Pipe_from_Fluidic_Diode_to_Downcomer','3']/df.loc['Pipe_from_Fluidic_Diode_to_Downcomer','1'] #pipe from diode to downcomer
    sumLA_p += df.loc['Downcomer_Upper_Upper','3']/df.loc['Downcomer_Upper_Upper','1'] #downcomer upper upper
    sumLA_p += df.loc['Downcomer_Upper_Section','3']/df.loc['Downcomer_Upper_Section','1'] #downcomer upper
    sumLA_p += df.loc['Downcomer_Middle_Section','3']/df.loc['Downcomer_Middle_Section','1'] #downcomer middle
    sumLA_p += df.loc['Downcomer_Lower_Section','3']/df.loc['Downcomer_Lower_Section','1'] #downcomer lower

    sumLA_m = referenceRatio*sumLA_p # this must be true to match eq 108. Design requirement

    #Calculating loop geometry number (in-vessel 'loop')
    G_r = 1.0 #given above forcing of sumLA_m. G = Loop geomtry ratio
    sumG = 0.0
    sumG += (df.loc['Core_Lower_Section','1']/length_p)*(A_ref_p/df.loc['Core_Lower_Section','3'])  # row 3 is Lengths. Row 1 is Areas.
    sumG += (df.loc['Pipe_from_Lower_Plenum_to_Active_Core','1']/length_p)*(A_ref_p/df.loc['Pipe_from_Lower_Plenum_to_Active_Core','3'])
    sumG += (df.loc['Core_Active_Section','1']/length_p)*(A_ref_p/df.loc['Core_Active_Section','3'])
    sumG += (df.loc['Core_Upper_Section','1']/length_p)*(A_ref_p/df.loc['Core_Upper_Section','3'])
    sumG += (df.loc['Pipe_from_Active_Core_to_Upper_Plenum','1']/length_p)*(A_ref_p/df.loc['Pipe_from_Active_Core_to_Upper_Plenum','3'])
    sumG += (df.loc['Pipe_from_Pump_Pool_to_Fluidic_Doide','1']/length_p)*(A_ref_p/df.loc['Pipe_from_Pump_Pool_to_Fluidic_Doide','3'])
    sumG += (df.loc['Pipe_from_Fluidic_Diode_to_Downcomer','1']/length_p)*(A_ref_p/df.loc['Pipe_from_Fluidic_Diode_to_Downcomer','3'])
    sumG += (df.loc['Downcomer_Upper_Upper','1']/length_p)*(A_ref_p/df.loc['Downcomer_Upper_Upper','3'])
    sumG += (df.loc['Downcomer_Upper_Section','1']/length_p)*(A_ref_p/df.loc['Downcomer_Upper_Section','3'])
    sumG += (df.loc['Downcomer_Middle_Section','1']/length_p)*(A_ref_p/df.loc['Downcomer_Middle_Section','3'])
    sumG += (df.loc['Downcomer_Lower_Section','1']/length_p)*(A_ref_p/df.loc['Downcomer_Lower_Section','3'])
    G_p = sumG
    G_m = G_p*G_r

    #var_out.loc['t='+t,'G_m'] = G_m
    #var_out.loc['t='+t,'sumLA_m'] = sumLA_m # this is the sum of the lengths/areas for loop components.

    #10a Core (rx) component heat transfer and fluid dynamics scaling
    St_NC_rx_r = 1.0
    St_NC_rx_p = (ehAs_rx_p*deltaT_sf_p/(rho_p*u_ref_NC_p*A_ref_p*cp_p*deltaT_p))
    St_NC_rx_m = St_NC_rx_p*St_NC_rx_r

    Re_NC_rx_p = rho_p*u_ref_NC_p*l_peb_p/visc_p
    Re_NC_rx_m = rho_m*u_ref_NC_m*l_peb_m/visc_m
    Re_NC_rx_r = Re_NC_rx_m/Re_NC_rx_p

    Gr_NC_rx_r = beta_r*deltaT_r*(Dh_rx_r**3)/((visc_r/rho_r)**2)
    Gr_NC_rx_p = (g*beta_p*deltaT_p*(Dh_rx_p**3)/((visc_p/rho_p)**2))
    Gr_NC_rx_m = Gr_NC_rx_p*Gr_NC_rx_r

    #10c Downcomer (dc) component heat transfer and fluid dynamics scaling
    St_NC_dc_r = 1.0
    St_NC_dc_p = (ehAs_dc_p*deltaT_sf_p/(rho_p*u_ref_NC_p*(np.pi/4)*(Dh_dc_p**2)*cp_p*deltaT_p))
    St_NC_dc_m = St_NC_dc_p*St_NC_dc_r

    Re_NC_dc_p = rho_p*u_ref_NC_p*Dh_dc_p/visc_p  #for downcomer width of channel is from length scaling.
    Re_NC_dc_m = rho_m*u_ref_NC_m*Dh_dc_m/visc_m
    Re_NC_dc_r = Re_NC_dc_m/Re_NC_dc_p

    Gr_NC_dc_r = beta_r*deltaT_r*(Dh_dc_r**3)/((visc_r/rho_r)**2)  # we assume deltaT_r during NC is the same as during FC (suspect)
    Gr_NC_dc_p = (g*beta_p*deltaT_p*(Dh_dc_p**3)/((visc_p/rho_p)**2)) #assuming deltaT during NC is the same as during FC (suspect)
    Gr_NC_dc_m = Gr_NC_dc_p*Gr_NC_dc_r

    #10d Reflector Bypass (rb) component heat transfer and fluid dynamics scaling
    St_NC_rb_r = 1.0
    St_NC_rb_p = (ehAs_rb_p*deltaT_sf_p/(rho_p*u_ref_NC_p*(np.pi/4)*(Dh_rb_p**2)*cp_p*deltaT_p))
    St_NC_rb_m = St_NC_rb_p*St_NC_rb_r

    Re_NC_rb_p = rho_p*u_ref_NC_p*Dh_rb_p/visc_p
    Re_NC_rb_m = rho_m*u_ref_NC_m*Dh_rb_m/visc_m
    Re_NC_rb_r = Re_NC_rb_m/Re_NC_rb_p

    Gr_NC_rb_r = beta_r*deltaT_r*(Dh_rb_r**3)/((visc_r/rho_r)**2)  # we assume deltaT_r during NC is the same as during FC (suspect)
    Gr_NC_rb_p = (g*beta_p*deltaT_p*(Dh_rb_p**3)/((visc_p/rho_p)**2)) #assuming deltaT during NC is the same as during FC (suspect)
    Gr_NC_rb_m = Gr_NC_rb_p*Gr_NC_rb_r

    #11a Core (rx) part 2- component solid heat transfer
    h_NC_rx_p = df.loc['h_NC_rx_p','0'] #need heat transfer coefficient from pebble to flibe from KP-SAM.
    h_NC_rx_m = df.loc['h_NC_rx_m','0']#assume it is the same in facility.
    h_NC_rx_r = h_NC_rx_m/h_NC_rx_p #assume we can match this for now. Can change this.

    Bi_NC_rx_r = 1.0 #we assume this, drives ls_rx_m and/or material selection
    Bi_NC_rx_p = (h_NC_rx_p*ls_rx_p/ks_rx_p)
    Bi_NC_rx_m = Bi_NC_rx_p*Bi_NC_rx_r #need to ensure we match this.

    ls_rx_m = (h_NC_rx_p*ls_rx_p/ks_rx_p)*(ks_rx_m/h_NC_rx_m) #pebble length scale (diameter) to match Bi#
    ls_rx_r = ls_rx_m/ls_rx_p

    var_out.loc['t='+t,'h_NC_rx_m'] = h_NC_rx_m
    var_out.loc['t='+t,'ls_rx_m'] = ls_rx_m
    var_out.loc['t='+t,'ks_rx_m'] = ks_rx_m

    #11b Reflector (rf) - component solid heat transfer
    #choice of model material will dictate ls.
    h_NC_rf_p = df.loc['h_NC_rf_p','0'] #need, heat transfer coefficient from flibe to reflector from KP-SAM.
    h_NC_rf_m = df.loc['h_NC_rf_m','0'] #can input this once it is determined
    h_NC_rf_r = h_NC_rf_m/h_NC_rf_p

    Bi_NC_rf_r = 1.0 #assume we want to match this.
    Bi_NC_rf_p = (h_rf_p*ls_rf_p/ks_rf_p)
    Bi_NC_rf_m = Bi_NC_rf_p*Bi_NC_rf_r

    ls_rf_r = ks_rf_r/h_NC_rf_r  # eq 112
    ls_rf_m = ls_rf_p*ls_rf_r #thickness of reflector design output

    var_out.loc['t='+t,'ls_rf_m'] = ls_rf_m
    var_out.loc['t='+t,'h_NC_rf_m'] = h_NC_rf_m
    var_out.loc['t='+t,'ks_rf_m'] = ks_rf_m

    #11c: Barrel (bl) - component solid heat transfer

    h_NC_bl_p = df.loc['h_NC_bl_p','0'] #need, heat transfer coefficient from flibe to barrel from KP-SAM.
    h_NC_bl_m = df.loc['h_NC_bl_m','0'] #can input this once it is determined
    ks_bl_p = ks_rv_p #using SS316 for both
    ks_bl_r = ks_rv_r
    h_NC_bl_r = h_NC_bl_m/h_NC_bl_p

    Bi_NC_bl_r = 1.0 #assume we want to match this.
    Bi_NC_bl_p = (h_NC_bl_p*ls_bl_p/ks_bl_p)
    Bi_NC_bl_m = Bi_NC_bl_p*Bi_NC_bl_r

    ls_bl_r = ks_bl_r/h_NC_bl_r  # eq 112
    ls_bl_m = ls_bl_p*ls_bl_r #thickness of reflector design output

    var_out.loc['t='+t,'ls_bl_m'] = ls_bl_m
    var_out.loc['t='+t,'h_NC_bl_m'] = h_NC_bl_m
    var_out.loc['t='+t,'ks_bl_m'] = ks_bl_m

    #11d Downcomer vessel wall (dc) component solid heat transfer

    h_NC_dc_p = df.loc['h_NC_dc_p','0'] # heat transfer coefficient from reflector to downcomer salt in KP-SAM
    h_NC_dc_m = df.loc['h_NC_dc_m','0'] #can input this once it is determined. Assume matched for now.
    h_NC_dc_r = h_NC_dc_m/h_NC_dc_p

    #these are included in downcomer section, but are for the vessel wall.
    ls_rv_p = df.loc['Reactor_Vessel_Upper','6'] #solid thickness (vessel)
    ls_rv_r = 2*ks_rv_r/h_NC_dc_r  # eq 112, matching Bi# for vessel heat coupling to downcomer fluid
    ls_rv_m = ls_rv_p*ls_rv_r #output length scale (thickness of vessel)
    Bi_NC_dc_p = (h_NC_dc_p*ls_rv_p/ks_rv_p)
    Bi_NC_dc_m = (h_NC_dc_m*ls_rv_m/ks_rv_m)
    Bi_NC_dc_r = Bi_NC_dc_m/Bi_NC_dc_p #should be 1.0 since we want to match this.

    var_out.loc['t='+t,'h_NC_dc_m'] = h_NC_dc_m
    var_out.loc['t='+t,'ks_rv_m'] = ks_rv_m
    var_out.loc['t='+t,'ls_rv_m'] = ls_rv_m

    #print(Bi_NC_dc_r)
    #print(" Bi# for vessel wall")
    #print(ls_rv_m)
    #print("thickness of vessel")

    #12 Solid time constant matching (See Eq 114), important to preserve component transient heat conduction in solids relative to system-wide time scale

    therm_diff_rx_r = ((ls_rx_r)**2)*u_ref_NC_r/length_r
    therm_diff_rv_r = ((ls_rv_r)**2)*u_ref_NC_r/length_r
    therm_diff_rf_r = ((ls_rf_r)**2)*u_ref_NC_r/length_r
    therm_diff_bl_r = ((ls_bl_r)**2)*u_ref_NC_r/length_r

    therm_diff_rx_p = ((ls_rx_p)**2)*u_ref_NC_p/length_p
    therm_diff_rv_p = ((ls_rv_p)**2)*u_ref_NC_p/length_p
    therm_diff_rf_p = ((ls_rf_p)**2)*u_ref_NC_p/length_p
    therm_diff_bl_p = ((ls_bl_p)**2)*u_ref_NC_p/length_p

    therm_diff_rx_m = therm_diff_rx_p*therm_diff_rx_r
    therm_diff_rv_m = therm_diff_rv_p*therm_diff_rv_r
    therm_diff_rf_m = therm_diff_rf_p*therm_diff_rf_r
    therm_diff_bl_m = therm_diff_bl_p*therm_diff_bl_r

    var_out.loc['t='+t,'thermal_diff_rx_m'] = therm_diff_rx_m # model core thermal diffusivity requirement (not used)
    var_out.loc['t='+t,'thermal_diff_rv_m'] = therm_diff_rv_m # model vessel thermal diffusivity requirement
    var_out.loc['t='+t,'thermal_diff_bl_m'] = therm_diff_bl_m # model barrel thermal diffusivity requirement

    # if the real and outputted thermal diffusivity above do not match, may want to adjust material selection or solid thickness to improve this.
    thermal_diff_bl_m_real = ks_bl_m/(rho_bl_m*cp_bl_m) #actual model thermal diffusivity based on df. selection
    thermal_diff_rv_m_real = ks_rv_m/(rho_rv_m*cp_rv_m) #actual model thermal diffusivity based on df. selection
    thermal_diff_rx_m_real = ks_rx_m/(rho_rx_m*cp_rx_m) #actual model thermal diffusivity based on df. selection (not used)

    #13 // core power scaling section
    #qsource_r = ks_rx_r*deltaT_r/(ls_rx_r**2) #eq 118 , no longer used

    rhocp_eff_p = rho_p*cp_p*por_p + rho_rx_p*cp_rx_p*(1-por_p)
    rhocp_eff_m = rho_m*cp_m*por_p + rho_rx_m*cp_rx_m*(1-por_p)
    rhocp_eff_r = rhocp_eff_m/rhocp_eff_p
    rho_rx_r = rho_rx_m/rho_rx_p
    cp_rx_r = cp_rx_m/cp_rx_p

    qsource_r = rho_r*cpf_r*deltaT_r*u_ref_FC_r*A_ref_r #heating power scaling ratio (Eq. 94)
    qsource_m = qsource_p*qsource_r
    var_out.loc['t='+t,'qsource_m'] = qsource_m
    Qp_FC_r = deltaP_FC_r*u_ref_FC_r*A_ref_r #pumping power ratio

    #RVACS heat removal scaling section.
    # Heat removal due to the RVACS during the transient will be simulated using forced circulation convection from the outer vessel wall.

    q_RVACS_p = df.loc[t,'32'] # [W]  heat removal by RVACS in KP-SAM during transient.
    emiss_r = 1.0 #emissivity. KP-SAM uses 0.8. Assume we can match it since we will be using a convection hx.
    A_h_r = ((A_ref_r)**(1/2))*length_r #area of hot surface, pi*diameter*length
    view_r = 1.0 #view factor ratio. Assume we match this. (0.7 for KP-SAM) Assume we can match it since we will be using a convection hx.
    T_rv_p = df.loc[t,'20']
    T_tcp_p = df.loc['T_tcp_p','0'] #assume that RVACS panels are at constant temp of boiling coolant at 100C after short transient.
    #since we are including rad HT between two surfaces, it is 0 initially since the temps are about the same. After water//
    #rushes into panels, they cool to about 100C approx.
    #deltaTfourth_p = (T_rv_p**4) - (T_tcp_p**4) #tcp = Thermosyphon panel. Data from KP-SAM.
    #deltaTfourth_m = ((T_rv_p*T_r)**4) - ((T_tcp_p*T_r)**4) #for model T_tcp_m is the heat sink temp.
    #deltaTfourth_r = deltaTfourth_m/deltaTfourth_p
    q_RVACS_m = q_RVACS_p*qsource_m/qsource_p #scaled RVACS heat removal requirement.

    var_out.loc['t='+t,'q_RVACS_m'] = q_RVACS_m

    #MISC RATIOS that may or may not be desired in the future
    por_r = 1.0 #matching porosity.

    # prototype heat transfer coefficients for each component. None of these are known, but estimates can be found.
    h_FC_rx_p = df.loc['h_FC_rx_p','0']
    h_FC_rf_p = df.loc['h_FC_rf_p','0']
    h_FC_dc_p = df.loc['h_FC_dc_p','0']
    h_FC_rx_r = 1.0
    h_FC_rf_r = 1.0
    h_FC_dc_r = 1.0
    h_FC_rx_m = h_FC_rx_p*h_FC_rx_r
    h_FC_rf_m = h_FC_rf_p*h_FC_rf_r
    h_FC_dc_m = h_FC_dc_p*h_FC_dc_r

    #these variables are not outputted, but may be useful at some point.
    Fo_NC_rx_r = therm_diff_rx_r*t_NC_r/(ls_rx_r**2)
    Fo_NC_rf_r = therm_diff_rf_r*t_NC_r/(ls_rf_r**2)
    Fo_NC_dc_r = therm_diff_rv_r*t_NC_r/(ls_rv_r**2)

    Bi_Fo_NC_rx_r = Bi_NC_rx_r*Fo_NC_rx_r
    Bi_Fo_NC_rf_r = Bi_NC_rf_r*Fo_NC_rf_r
    Bi_Fo_NC_dc_r = Bi_NC_dc_r*Fo_NC_dc_r

    Fo_FC_rx_r = therm_diff_rx_r*t_FC_r/(ls_rx_r**2)
    Fo_FC_rf_r = therm_diff_rf_r*t_FC_r/(ls_rf_r**2)
    Fo_FC_dc_r = therm_diff_rv_r*t_FC_r/(ls_rv_r**2)

    Bi_FC_rx_r = (h_FC_rx_r*ls_rx_r/ks_rx_r)
    Bi_FC_rf_r = (h_FC_rx_r*ls_rx_r/ks_rx_r)
    Bi_FC_dc_r = (h_FC_rx_r*ls_rx_r/ks_rx_r)

    Bi_Fo_FC_rx_r = Bi_FC_rx_r*Fo_FC_rx_r
    Bi_Fo_FC_rf_r = Bi_FC_rf_r*Fo_FC_rf_r
    Bi_Fo_FC_dc_r = Bi_FC_dc_r*Fo_FC_dc_r

#------------------------------------------------------------------------------------------------
#MODEL DIMENSIONS AND CHARACTERISTIC RATIO CALCULATION. (using the ratios determined in the above calculations to generate model dimensions)
    #this section fills the model_outputs table with scaled dimensions. The inputs come from the KP-SAM input card geometry.

    for i in range(15,51): #Scaling area
        (var_out.iloc[i,1]) = (df.iloc[i,1])*A_ref_r #scaling all component areas with the system-wide area scaling
        if i == 3: #if component is the active core, reducing the area further to account for the lack of pebbles in ITF-0 (the 0.40), the non-prototypical geometry of the ITF-0 core, and the reflector bypass and inflow/outflow fluid volumes (the 0.9).
            (var_out.iloc[i,1]) = (df.iloc[i,1])*A_ref_r*core_scaling_factor
    for i in range(15,51):#Scaling radius for each component (derived from area scaling)
        if i == 3: #if component is the active core
            (var_out.iloc[i,2]) = (df.iloc[i,2])*((A_ref_r*core_scaling_factor)**.5)
        else:
            (var_out.iloc[i,2]) = (df.iloc[i,2])*(A_ref_r**0.5) #based on the area scaling since radius is a function of area
    for i in range(15,51):#Scaling length (height) for all components
        (var_out.iloc[i,3]) = (df.iloc[i,3])*length_r
    for i in range(15,51):#Scaling hydraulic diameter
        (var_out.iloc[i,4]) = (df.iloc[i,4])*Dh_rx_r #this is separate from the area scaling, comes from matching flow behavior (Re#, St#) in core.
    for i in range(15,51):#Scaling volume
        (var_out.iloc[i,5]) = (df.iloc[i,5])*(length_r*A_ref_r) #volume is a function of length (height) and area scaling.
    for i in range(15,51):#Scaling width
        if i == 27: #if component is the barrel lower, scale width as per Bi# scaling, not system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*ls_bl_r
        elif i == 28: #if component is the barrel middle, scale width as per Bi# scaling, not system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*ls_bl_r
        elif i == 29: #if component is the barrel upper, scale width as per Bi# scaling, not system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*ls_bl_r
        elif i == 34: #if component is the barrel upper upper, scale width as per Bi# scaling, not system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*ls_bl_r
        elif i == 12: #if component is the vessel lower (downcomer wall), scale width as per Bi# scaling, not system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*ls_rv_r
        elif i == 13: #if component is the vessel middle (downcomer wall), scale width as per Bi# scaling, not system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*ls_rv_r
        elif i == 14: #if component is the vessel upper (downcomer wall), scale width as per Bi# scaling, not system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*ls_rv_r
        elif i == 36: #if component is the vessel upper upper (downcomer wall), scale width as per Bi# scaling, not system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*ls_rv_r
        elif i == 9: #if component is lower downcomer channel, scale with length scale (47%).
            (var_out.iloc[i,6]) = (df.iloc[i,6])*length_r
        elif i == 10: #if component is middle downcomer channel, scale with length scale
            (var_out.iloc[i,6]) = (df.iloc[i,6])*length_r
        elif i == 11: #if component is upper downcomer channel, scale with length scale
            (var_out.iloc[i,6]) = (df.iloc[i,6])*length_r
        elif i == 35: #if component is upper upper downcomer channel, scale with length scale
            (var_out.iloc[i,6]) = (df.iloc[i,6])*length_r
        else: #for the remaining non-special cases, scale with the system-wide Dh scaling
            (var_out.iloc[i,6]) = (df.iloc[i,6])*Dh_rx_r
    for i in range(15,51):#Scaling porosity
        (var_out.iloc[i,7]) = (df.iloc[i,7])*por_r #we match porosity
    for i in range(15,51):#Scaling temperature
        (var_out.iloc[i,8]) = (df.iloc[i,8])*T_r
    for i in range(15,51):#Scaling form loss coefficient
        (var_out.iloc[i,9]) = (df.iloc[i,9])*F_r #we match this
    for i in range(15,51):#Scaling friction loss coefficent
        (var_out.iloc[i,10]) = (df.iloc[i,10])*F_r #we match this

    #need to check if the requirement from section #9a is still met
    sumLA2_m = 0.0
    sumLA2_m += var_out.iloc[15,3]/var_out.iloc[15,1] #core_lower
    sumLA2_m += var_out.iloc[16,3]/var_out.iloc[16,1] #pipe from lower plenum to active core
    sumLA2_m += var_out.iloc[17,3]/(var_out.iloc[17,1]*(core_scaling_factor)) #active core
    sumLA2_m += var_out.iloc[19,3]/var_out.iloc[19,1] #upper core
    sumLA2_m += var_out.iloc[18,3]/var_out.iloc[18,1] #pipe from active core to upper plenum
    sumLA2_m += var_out.iloc[39,3]/var_out.iloc[39,1] #pipe from upper plenum to diode
    sumLA2_m += var_out.iloc[40,3]/var_out.iloc[40,1] #pipe from diode to downcomer
    sumLA2_m += var_out.iloc[49,3]/var_out.iloc[49,1] #downcomer UU
    sumLA2_m += var_out.iloc[23,3]/var_out.iloc[23,1] #downcomer U
    sumLA2_m += var_out.iloc[24,3]/var_out.iloc[24,1] #downcomer M
    sumLA2_m += var_out.iloc[25,3]/var_out.iloc[25,1] #downcomer L

    sumLA_m = int(sumLA_m) #make them integers to avoid false flag due to tiny fractional differences from python memory
    sumLA2_m = int(sumLA2_m)

    if sumLA2_m == sumLA_m:
        pass
    else:
        pass
        print("Sum of loop length/area doesn't match requirement.")

    #In this section, I recalculate the ratios that were initally set to 1 to ensure consistency.
    Pr_r = visc_r*cpf_r/k_r
    Ri_FC_r = beta_r*(deltaT_r)*length_r/(u_ref_FC_r**2)
    Fr_FC_r =u_ref_FC_r/((length_r))**.5
    Eu_FC_r = deltaP_FC_r/(rho_r*(u_ref_FC_r**2))
    F_r =1
    f_r =1
    St_FC_rx_r =(ehAs_rx_r*deltaT_sf_r/(rho_r*u_ref_FC_r*A_ref_r*cpf_r*deltaT_r))
    St_FC_hx_r =(ehAs_hx_r*deltaT_sf_r/(rho_r*u_ref_FC_r*(np.pi/4)*(Dh_hx_r**2)*cpf_r*deltaT_r)) #use rx u_ref
    St_FC_dc_r =(ehAs_dc_r*deltaT_sf_r/(rho_r*u_ref_FC_r*(Dh_dc_r**2)*cpf_r*deltaT_r))
    Ri_NC_r= (beta_r*deltaT_r*length_r/((u_ref_NC_r)**2))
    G_r =G_m/G_p
    St_NC_rx_r =(ehAs_rx_r*deltaT_sf_r/(rho_r*u_ref_NC_r*A_ref_r*cpf_r*deltaT_r))
    St_NC_dc_r =(ehAs_dc_r*deltaT_sf_r/(rho_r*u_ref_NC_r*(Dh_dc_r**2)*cpf_r*deltaT_r))
    Bi_NC_rx_r =(h_NC_rx_r*ls_rx_r/ks_rx_r)
    Bi_NC_rf_r =(h_NC_rf_r*ls_rf_r/ks_rf_r)
    Bi_NC_dc_r =(h_NC_dc_r*ls_rv_r/ks_rv_r)

#------------------------------------------------------------------------------------------------
#OUTPUTS (this section generates the ratios and distortions onto an output card)
    if t == 'SS':
        pass
    else:
        t = int(t)

    if t == 'SS':
        t = 1
    elif t == 0:
        t = 2
    elif t == 1:
        t = 3
    elif t == 10:
        t = 4
    elif t == 100:
        t = 5
    elif t == 1000:
        t = 6
    elif t == 10000:
        t = 7
    elif t == 100000:
        t = 8
    elif t == 600000:
        t = 9
    t = int(t)

    var_out.iloc[53,0] = 'Pr_r' #put the results into the dataframe, and calculate distortions.
    var_out.iloc[53,t] = Pr_r
    var_out.iloc[88,t] = (Pr_p - Pr_m)/Pr_p

    var_out.iloc[54,0] = 'F_r'
    var_out.iloc[54,t] = F_r
    var_out.iloc[89,t] = (F_p - F_m)/F_p

    if t==1: #forced circulation, normal ops
        var_out.iloc[55,0] = 'Eu_FC_r'
        var_out.iloc[55,t] = Eu_FC_r
        var_out.iloc[90,t] = (Eu_FC_p - Eu_FC_m)/Eu_FC_p

        var_out.iloc[56,0] = 'Ri_FC_r'
        var_out.iloc[56,t] = Ri_FC_r
        var_out.iloc[91,t] = (Ri_FC_p - Ri_FC_m)/Ri_FC_p

        var_out.iloc[57,0] = 'St_FC_rx_r'
        var_out.iloc[57,t] = St_FC_rx_r
        var_out.iloc[92,t] = (St_FC_rx_p - St_FC_rx_m)/St_FC_rx_p

        var_out.iloc[58,0] = 'St_FC_dc_r'
        var_out.iloc[58,t] = St_FC_dc_r
        var_out.iloc[93,t] = (St_FC_dc_p - St_FC_dc_m)/St_FC_dc_p

        var_out.iloc[59,0] = 'Re_FC_dc_r'
        var_out.iloc[59,t] = Re_FC_dc_r
        var_out.iloc[94,t] = (Re_FC_dc_p - Re_FC_dc_m)/Re_FC_dc_p

        var_out.iloc[60,0] = 'St_FC_hx_r'
        var_out.iloc[60,t] = St_FC_hx_r
        var_out.iloc[95,t] = (St_FC_hx_p - St_FC_hx_m)/St_FC_hx_p

        var_out.iloc[61,0] = 'Re_FC_hx_r'
        var_out.iloc[61,t] = Re_FC_hx_r
        var_out.iloc[96,t] = (Re_FC_hx_p - Re_FC_hx_m)/Re_FC_hx_p

        var_out.iloc[62,0] = 'Re_FC_rx_r'
        var_out.iloc[62,t] = Re_FC_rx_r
        var_out.iloc[97,t] = (Re_FC_rx_p - Re_FC_rx_m)/Re_FC_rx_p

        var_out.iloc[63,0] = 't_FC_r'
        var_out.iloc[63,t] = t_FC_r
        var_out.iloc[98,t] = '-'

        var_out.iloc[64,0] = 'Bi_FC_dc_r'
        var_out.iloc[64,t] = Bi_FC_dc_r
        var_out.iloc[99,t] = '-'

        var_out.iloc[65,0] = 'Bi_FC_rf_r'
        var_out.iloc[65,t] = Bi_FC_rf_r
        var_out.iloc[100,t] = '-'

        var_out.iloc[66,0] = 'Bi_FC_rx_r'
        var_out.iloc[65,t] = Bi_FC_rx_r
        var_out.iloc[101,t] = '-'

        var_out.iloc[67,0] = 'St_FC_rb_r'
        var_out.iloc[67,t] = St_FC_rb_r
        var_out.iloc[102,t] = (St_FC_rb_p - St_FC_rb_m)/St_FC_rb_p

        var_out.iloc[68,0] = 'Re_FC_rb_r'
        var_out.iloc[68,t] = Re_FC_rb_r
        var_out.iloc[103,t] = (Re_FC_rb_p - Re_FC_rb_m)/Re_FC_rb_p

    else: #natural circulation only
        var_out.iloc[69,0] = 'Ri_NC_r'
        var_out.iloc[69,t] = Ri_NC_r
        var_out.iloc[104,t] = (Ri_NC_p - Ri_NC_m)/Ri_NC_p

        var_out.iloc[70,0] = 'St_NC_rx_r'
        var_out.iloc[70,t] = St_NC_rx_r
        var_out.iloc[105,t] = (St_NC_rx_p - St_NC_rx_m)/St_NC_rx_p

        var_out.iloc[71,0] = 'St_NC_dc_r'
        var_out.iloc[71,t] = St_NC_dc_r
        var_out.iloc[106,t] = (St_NC_dc_p -St_NC_dc_m)/St_NC_dc_p

        var_out.iloc[72,0] = 'Gr_NC_dc_r'
        var_out.iloc[72,t] = Gr_NC_dc_r
        var_out.iloc[107,t] = (Gr_NC_dc_p - Gr_NC_dc_m)/Gr_NC_dc_p

        var_out.iloc[73,0] = 'Re_NC_dc_r'
        var_out.iloc[73,t] = Re_NC_dc_r
        var_out.iloc[108,t] = (Re_NC_dc_p - Re_NC_dc_m)/Re_NC_dc_p

        var_out.iloc[74,0] = 'Re_NC_rx_r'
        var_out.iloc[74,t] = Re_NC_rx_r
        var_out.iloc[109,t] = (Re_NC_rx_p - Re_NC_rx_m)/Re_NC_rx_p

        var_out.iloc[75,0] = 't_NC_r'
        var_out.iloc[75,t] = t_NC_r
        var_out.iloc[110,t] = '-'

        var_out.iloc[76,0] = 'Bi_NC_dc_r'
        var_out.iloc[76,t] = Bi_NC_dc_r
        var_out.iloc[111,t] = (Bi_NC_dc_p - Bi_NC_dc_m)/Bi_NC_dc_p

        var_out.iloc[77,0] = 'Bi_NC_rf_r'
        var_out.iloc[77,t] = Bi_NC_rf_r
        var_out.iloc[112,t] = (Bi_NC_rf_p - Bi_NC_rf_m)/Bi_NC_rf_p

        var_out.iloc[78,0] = 'Bi_NC_rx_r'
        var_out.iloc[78,t] = Bi_NC_rx_r
        var_out.iloc[113,t] = (Bi_NC_rx_p - Bi_NC_rx_m)/Bi_NC_rx_p

        var_out.iloc[79,0] = 'St_NC_rb_r'
        var_out.iloc[79,t] = St_NC_rb_r
        var_out.iloc[114,t] = (St_NC_rb_p -St_NC_rb_m)/St_NC_rb_p

        var_out.iloc[80,0] = 'Gr_NC_rb_r'
        var_out.iloc[110,t] = Gr_NC_rb_r
        var_out.iloc[115,t] = (Gr_NC_rb_p - Gr_NC_rb_m)/Gr_NC_rb_p

        var_out.iloc[81,0] = 'Re_NC_rb_r'
        var_out.iloc[81,t] = Re_NC_rb_r
        var_out.iloc[116,t] = (Re_NC_rb_p - Re_NC_rb_m)/Re_NC_rb_p

        var_out.iloc[82,0] = 'Fo_NC_rx_r'
        var_out.iloc[82,t] = Fo_NC_rx_r
        var_out.iloc[117,t] = '-'

        var_out.iloc[83,0] = 'Fo_NC_rf_r'
        var_out.iloc[83,t] = Fo_NC_rf_r
        var_out.iloc[118,t] = '-'

        var_out.iloc[84,0] = 'Fo_NC_dc_r'
        var_out.iloc[84,t] = Fo_NC_dc_r
        var_out.iloc[119,t] = '-'

    #the rest, FC or NC
    var_out.iloc[85,0] = 'G_r'
    var_out.iloc[85,t] = G_r
    var_out.iloc[120,t] = (G_p - G_m)/G_p

    #scale_out.loc[14,'Scaling_Group'] = 'HR_r' #the heat source number is currently not used in the analysis.
    #scale_out.loc[14,t] = HR_r
    #dist_out.loc[14,t] = (HR_p - HR_m)/HR_p

    if t == 1:
        t = 'SS'
    elif t == 2:
        t = 0
    elif t == 3:
        t = 1
    elif t == 4:
        t = 10
    elif t == 5:
        t = 100
    elif t == 6:
        t = 1000
    elif t == 7:
        t = 10000
    elif t == 8:
        t = 100000
    elif t == 9:
        t = 600000

var_out.iloc[88:125,0] = var_out.iloc[53:90,0] #filling the first column of the distortions table
var_out.iloc[122,0] = '-' #making a mis-done label due to bug dissapear from the output card

print('Done.')

#dim_out.to_csv("model_outputs.csv", sep=',',index=False) #write the outputs to csv files
#scale_out.to_csv("ratios_outputs.csv", sep=',',index=False)
#dist_out.to_csv("distortions_outputs.csv", sep=',',index=False)
var_out.to_csv("krisp_output.csv", sep=',',index=False) #change index=True to see indexing column in output card.
#------------------------------------------------------------------------------------------------
#END










































#-
