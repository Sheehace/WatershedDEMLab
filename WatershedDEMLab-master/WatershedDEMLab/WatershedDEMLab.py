# CatchmentLab
# Author: Chris Sheehan

#%% Imports and paths

# Print status
print('Importing libraries, capturing paths, importing input parameters, and creating output directories...')

# Libraries
import inspect
import os
import os.path
from os import path
import sys
# sys.path.append("C:\Program Files\QGIS 3.28.10\apps\Python39\Lib")
import numpy as np
from numpy import nan, isnan
import pandas as pd
from landlab import RasterModelGrid, imshow_grid
from landlab.components import StreamPowerEroder, LinearDiffuser, FlowAccumulator, ChannelProfiler, PriorityFloodFlowRouter
from landlab.io.esri_ascii import read_esri_ascii, write_esri_ascii
from matplotlib import pyplot as plt
from osgeo import osr, gdal

# Capture paths
script_name = os.path.basename(__file__)
script_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
model_directory = script_directory.replace('\\model', '')
input_directory = model_directory + '\\input'

# Import inputs
Params = pd.read_csv(input_directory + '\Input.csv')

# Create output directories
if path.exists(str(model_directory)+'/Output') == False:
    os.mkdir(str(model_directory)+'/Output')
if path.exists(str(model_directory)+'/Output//topographic__elevation__png') == False:
    os.mkdir(str(model_directory)+'/Output//topographic__elevation__png')
if path.exists(str(model_directory)+'/Output//topographic__elevation__pdf') == False:
    os.mkdir(str(model_directory)+'/Output//topographic__elevation__pdf')
if path.exists(str(model_directory)+'/Output//topographic__elevation__asc') == False:
    os.mkdir(str(model_directory)+'/Output//topographic__elevation__asc')
if path.exists(str(model_directory)+'/Output//topographic__elevation__tif') == False:
    os.mkdir(str(model_directory)+'/Output//topographic__elevation__tif')
if path.exists(str(model_directory)+'/Output//dzdt__png') == False:
    os.mkdir(str(model_directory)+'/Output//dzdt__png')
if path.exists(str(model_directory)+'/Output//dzdt__pdf') == False:
    os.mkdir(str(model_directory)+'/Output//dzdt__pdf')
if path.exists(str(model_directory)+'/Output//dzdt__asc') == False:
    os.mkdir(str(model_directory)+'/Output//dzdt__asc')
if path.exists(str(model_directory)+'/Output//dzdt__tif') == False:
    os.mkdir(str(model_directory)+'/Output//dzdt__tif')
if path.exists(str(model_directory)+'/Output//channel_map__png') == False:
    os.mkdir(str(model_directory)+'/Output//channel_map__png')
if path.exists(str(model_directory)+'/Output//channel_map__pdf') == False:
    os.mkdir(str(model_directory)+'/Output//channel_map__pdf')
    
# Print space
print(' ')

#%% Retrieve spin-up parameters

# Print status
print('Retrieving model parameters...')

# Assign inputs
# DEM_path = [os.path.dirname(Params.DEM_path[0]) + '\\' + os.path.basename(Params.DEM_path[0])]
DEM_path = os.path.abspath(Params.DEM_path[0])
no_data_value = Params.no_data_value[0]
m_sp = Params.m_sp[0]
n_sp = Params.n_sp[0]
k_sp = Params.k_sp[0]
k_hs = Params.k_hs[0]
u = Params.u[0]
dt = Params.dt[0]     
tmax = Params.tmax[0]
export_interval = Params.export_interval[0]

# Print space
print(' ')

#%% Create grid

# Print status
print('Creating model grid...')

# Import DEM
(mg, zr) = read_esri_ascii(DEM_path, name='topographic__elevation')

# Handle model DEM dimensions and non-value nodes
no_data_nodes = np.where(mg.at_node['topographic__elevation'] == no_data_value)
no_data_nodes = no_data_nodes[0]
mg.at_node['topographic__elevation'][no_data_nodes] = np.nan
mg.set_nodata_nodes_to_closed(mg.at_node['topographic__elevation'], np.nan)
number_of_rows = mg.number_of_cell_rows + 2
number_of_columns = mg.number_of_cell_columns + 2
nodes = np.arange(0, np.size(mg.at_node['topographic__elevation']))

# Find pour point
pffr0 = PriorityFloodFlowRouter(mg, depression_handler = 'fill')
pffr0.run_one_step()
pour_point_node = np.where(mg.at_node['drainage_area'] == np.max(mg.at_node['drainage_area']))
pour_point_node = pour_point_node[0]
# imshow_grid(mg, 'drainage_area', grid_units=('m', 'm'), var_name="Elevation (m)", cmap='terrain', allow_colorbar=True)
# imshow_grid(mg, 'topographic__elevation', grid_units=('m', 'm'), var_name="Elevation (m)", cmap='terrain', allow_colorbar=True)

# Create node keys
mg.add_zeros('node', 'previous__topographic__elevation')
mg.add_zeros('node', 'dzdt')
mg.add_zeros('node', 'erosionrate')

#%% Block 2: The Grid Boundaries

# ENTER VARIABLES ############################################################
##############################################################################

# EDGES
East = 4
North = 4
West = 4
South = 4

# CORNERS
Northeast = 4
Northwest = 4
Southwest = 1
Southeast = 4

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------
if (East < 1) or (North < 1) or (West < 1) or (South < 1) or (Northeast < 1) or (Northwest < 1) or (Southwest < 1) or (Southeast < 1):
    raise Warning("Boundary values must be an integer between 1 and 4")
if (East > 4) or (North > 4) or (West > 4) or (South > 4) or (Northeast > 4) or (Northwest > 4) or (Southwest > 4) or (Southeast > 4):
    raise Warning("Boundary values must be an integer between 1 and 4")
if (type(East) != int) or (type(North) != int) or (type(West) != int) or (type(South) != int) or (type(Northeast) != int) or (type(Northwest) != int) or (type(Southwest) != int) or (type(Southeast) != int):
    raise Warning("Boundary values must be an integer between 1 and 4")
mg.set_status_at_node_on_edges(right=East, top=North, left=West, bottom=South)
mg.at_node.keys()
['topographic__elevation']
if Northeast == 1:
    mg.status_at_node[-1] = mg.BC_NODE_IS_FIXED_VALUE
if Northeast == 2:
    mg.status_at_node[-1] = mg.BC_NODE_IS_FIXED_GRADIENT
if Northeast == 3:
    mg.status_at_node[-1] = mg.BC_NODE_IS_LOOPED 
if Northeast == 4:
    mg.status_at_node[-1] = mg.BC_NODE_IS_CLOSED
if Northwest == 1:
    mg.status_at_node[(number_of_rows * number_of_columns) - number_of_columns] = mg.BC_NODE_IS_FIXED_VALUE
if Northwest == 2:
    mg.status_at_node[(number_of_rows * number_of_columns) - number_of_columns] = mg.BC_NODE_IS_FIXED_GRADIENT
if Northwest == 3:
    mg.status_at_node[(number_of_rows * number_of_columns) - number_of_columns] = mg.BC_NODE_IS_LOOPED 
if Northwest == 4:
    mg.status_at_node[(number_of_rows * number_of_columns) - number_of_columns] = mg.BC_NODE_IS_CLOSED
if Southwest == 1:
    mg.status_at_node[0] = mg.BC_NODE_IS_FIXED_VALUE
if Southwest == 2:
    mg.status_at_node[0] = mg.BC_NODE_IS_FIXED_GRADIENT
if Southwest == 3:
    mg.status_at_node[0] = mg.BC_NODE_IS_LOOPED 
if Southwest == 4:
    mg.status_at_node[0] = mg.BC_NODE_IS_CLOSED
if Southeast == 1:
    mg.status_at_node[number_of_columns - 1] = mg.BC_NODE_IS_FIXED_VALUE
if Southeast == 2:
    mg.status_at_node[number_of_columns - 1] = mg.BC_NODE_IS_FIXED_GRADIENT
if Southeast == 3:
    mg.status_at_node[number_of_columns - 1] = mg.BC_NODE_IS_LOOPED 
if Southeast == 4:
    mg.status_at_node[number_of_columns - 1] = mg.BC_NODE_IS_CLOSED
  




#%% BLOCK 3: ENTER LITHOLOGIC PARAMETERS 

# ENTER VARIABLES ############################################################
##############################################################################

m_sp = 0.5
n_sp = 1

Ksp_1 = 0.001
Khs_1 = 0.01

Ksp_2 = 0.001
Khs_2 = 0.01

Ksp_3 = 0.001
Khs_3 = 0.01

Ksp_4 = 0.001
Khs_4 = 0.01

Ksp_5 = 0.0001
Khs_5 = 0.01

##############################################################################
# ENTER VARIABLES ############################################################

# SELECT MODE ****************************************************************
#*****************************************************************************

# EXACTLY ONE OF THESE FOUR MODES MUST BE TRUE (NO MORE, NO LESS)

# SPATIALLY_UNIFORM MODE
Spatially_Uniform = True

# SPATIALLY_ZONED MODE
Spatially_Zoned = False
K_Edge_2 = 100
K_Edge_3 = 200
K_Edge_4 = 300
K_Edge_5 = 15000

# TILTED_ROCKS MODE
Tilted_Rocks = False
tilted_layer_thickness = 100
tilted_sediment_lithology = 1
dip = 5
dip_direction = 'E'

# FOLDED_ROCKS MODE
Folded_Rocks = False
folded_layer_thickness = 100
folded_sediment_lithology = 1
curvature_x = 0.003
curvature_y = 0.001
Fold_Type = 'Anticline'

#*****************************************************************************
# SELECT MODE ****************************************************************

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------
if (type(Spatially_Uniform) != bool) or (type(Spatially_Zoned) != bool) or (type(Tilted_Rocks) != bool) or (type(Folded_Rocks) != bool):
    raise Warning("The following variables must be either True or False: Spatially_Uniform, Spatially_Zoned, Tilted_Rocks, Folded_Rocks")
if (np.count_nonzero([Spatially_Uniform, Spatially_Zoned, Tilted_Rocks, Folded_Rocks]) > 1):
    raise Warning("Only one of the following variables can be True: Spatially_Uniform, Spatially_Zoned, Tilted_Rocks, Folded_Rocks. The other three must be False")
if (np.count_nonzero([Spatially_Uniform, Spatially_Zoned, Tilted_Rocks, Folded_Rocks]) < 1):
    raise Warning("At least one of the following variables must be True: Spatially_Uniform, Spatially_Zoned, Tilted_Rocks, Folded_Rocks. The other three must be False")
if Spatially_Uniform == True:
    print("Erodibility mode = Spatially_Uniform")
    Ksp = np.ones(mg.number_of_nodes) * Ksp_1 
    Khs = np.ones(mg.number_of_nodes) * Khs_1
if Spatially_Zoned == True: 
    print("Erodibility mode = Spatially_Zoned")
    Ksp = np.ones(mg.number_of_nodes) * Ksp_1 
    Khs = np.ones(mg.number_of_nodes) * Khs_1
    if K_Edge_2 > 0:
        Ksp[np.where(mg.node_x>K_Edge_2)] = Ksp_2 
        Khs[np.where(mg.node_x>K_Edge_2)] = Khs_2
    if K_Edge_3 > 0:
        Ksp[np.where(mg.node_x>K_Edge_3)] = Ksp_3 
        Khs[np.where(mg.node_x>K_Edge_3)] = Khs_3
    if K_Edge_4 > 0:
        Ksp[np.where(mg.node_x>K_Edge_4)] = Ksp_4 
        Khs[np.where(mg.node_x>K_Edge_4)] = Khs_4
    if K_Edge_5 > 0:
        Ksp[np.where(mg.node_x>K_Edge_5)] = Ksp_5 
        Khs[np.where(mg.node_x>K_Edge_5)] = Khs_5   
if Tilted_Rocks == True:
    print("Erodibility mode = Tilted_Rocks")
    attrs = {'K_sp': {0: Ksp_1, 1: Ksp_2, 2: Ksp_3, 3: Ksp_4, 4: Ksp_5}, 'K_hs': {0: Khs_1, 1: Khs_2, 2: Khs_3, 3: Khs_4, 4: Khs_5}}
    z0s = np.arange(0, tilted_layer_thickness * 5, tilted_layer_thickness)
    z0s[-1] += dxy*number_of_rows
    ids = np.arange(0,5,1)
    dip_rad = dip * 0.0174533                                                                       
    multiplier = math.tan(dip_rad)
    if dip_direction == 'E':
        function = lambda x, y: (multiplier * x)
        lith = LithoLayers(mg, z0s, ids, function=function, attrs=attrs, x0=0, y0=0,)
    if dip_direction == 'W':
        function = lambda x, y: (multiplier * -x) 
        lith = LithoLayers(mg, z0s, ids, function=function, attrs=attrs, x0=((number_of_columns*dxy)-1), y0=0,)
    if dip_direction == 'N':
        function = lambda x, y: (multiplier * y) 
        lith = LithoLayers(mg, z0s, ids, function=function, attrs=attrs, x0=0, y0=0,)
    if dip_direction == 'S':
        function = lambda x, y: (multiplier * -y) 
        lith = LithoLayers(mg, z0s, ids, function=function, attrs=attrs, x0=0, y0=((number_of_rows*dxy)-1))
    if dip_direction == 'NE':
        function = lambda x, y: (multiplier * x) + (multiplier * y)
        lith = LithoLayers(mg, z0s, ids, function=function, attrs=attrs, x0=0, y0=0)
    if dip_direction == 'SW':
        function = lambda x, y: (multiplier * -x) + (multiplier * -y)
        lith = LithoLayers(mg, z0s, ids, function=function, attrs=attrs, x0=((number_of_columns*dxy)-1), y0=((number_of_rows*dxy)-1))
    if dip_direction == 'NW':
        function = lambda x, y: (multiplier * -x) + (multiplier * y)
        lith = LithoLayers(mg, z0s, ids, function=function, attrs=attrs, x0=((number_of_columns*dxy)-1), y0=0)
    if dip_direction == 'SE':
        function = lambda x, y: (multiplier * x) + (multiplier * -y)
        lith = LithoLayers(mg, z0s, ids, function=function, attrs=attrs, x0=0, y0=((number_of_rows*dxy)-1))
    lith.rock_id = tilted_sediment_lithology
if Folded_Rocks == True:
    print("Erodibility mode = Folded_Rocks")
    curvature_x_rad = (curvature_x * 0.0174533)                                                                     
    multiplier_x = (math.tan(curvature_x_rad)) / 2
    curvature_y_rad = (curvature_y * 0.0174533)                                                                     
    multiplier_y = (math.tan(curvature_y_rad)) / 2
    attrs = {'K_sp': {0: Ksp_1, 1: Ksp_2, 2: Ksp_3, 3: Ksp_4, 4: Ksp_5}, 'K_hs': {0: Khs_1, 1: Khs_2, 2: Khs_3, 3: Khs_4, 4: Khs_5}}
    z0s = np.arange(0, folded_layer_thickness * 5, folded_layer_thickness)
    z0s[-1] += dxy*number_of_rows
    ids = np.arange(0,5,1)
    if Fold_Type == 'Syncline':
        anticline_func = lambda x, y: ( (multiplier_x * (x**2)) + (multiplier_y * (y**2)) ) *-1
    if Fold_Type == 'Anticline':
        anticline_func = lambda x, y: (multiplier_x * (x**2)) + (multiplier_y * (y**2)) 
    lith = LithoLayers(mg, z0s, ids, x0=(number_of_columns*dxy)/2, y0=(number_of_rows*dxy)/2, function=anticline_func, attrs=attrs)
    lith.rock_id = folded_sediment_lithology   
Lith_dict = {"Ksp_1": Ksp_1}
if Spatially_Zoned == True:
    Lith_dict.update({"Minimum_Block_Ksp": np.min(Ksp)})
if Tilted_Rocks == True:
    Lith_dict.update({"Minimum_Tilted_Ksp": min(lith._attrs['K_sp'].values())})
if Folded_Rocks == True:
    Lith_dict.update({"Minimum_Tilted_Ksp": min(lith._attrs['K_sp'].values())})
Min_Ksp = min(Lith_dict.values())    
    
    
    


#%% BLOCK 4: TECTONIC PARAMETERS

# ENTER VARIABLES ############################################################
##############################################################################

U_1 = 0.001

##############################################################################
# ENTER VARIABLES ############################################################

# OPTIONAL FUNCTIONS A *******************************************************
#*****************************************************************************

# YOU CAN USE BOTH OF THESE FUNCTIONS AT THE SAME TIME OR NEITHER AT ALL.
# YOU CAN ALSO USE EITHER (OR BOTH) OF THEM WITH ANY OF THE FUNCTIONS IN
# OPTIONAL FUNCTIONS B.

Periodic_Uplift = False                                                                                                                              
U_pulse = 0.002                                       
Period = 2000                        
U_pulse_duration = 300              

Normal_Fault = False                                                  
Dip_Angle = 60                  # Degrees                              
Recurrence = 1000                # Years 
Slip = 2                        # meters
x1 = 10000
y1 = 0
x2 = 10000
y2 = 1000
                                                                                
#*****************************************************************************
# OPTIONAL FUNCTIONS A *******************************************************

# OPTIONAL FUNCTIONS B *******************************************************
#*****************************************************************************

# YOU CAN ONLY USE ONE OF THESE THREE FUNCTIONS AT A TIME (OR NONE AT ALL).

Uplifted_Blocks = False
U_2 = 0.002
U_Edge_2 = 4000
U_3 = 0.003
U_Edge_3 = 8000
U_4 = 0.002
U_Edge_4 = 12000
U_5 = 0.001
U_Edge_5 = 16000

Tilting = False
Tilting_High_Uplift_Rate = 0.005
Tilting_Direction = 'S'

Gradient = False
Gradient_High_Uplift_Rate = 0.005
Gradient_Axis = 'N-S'

#*****************************************************************************
# OPTIONAL FUNCTIONS B *******************************************************

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------
if (type(Periodic_Uplift) != bool) or (type(Normal_Fault) != bool) or (type(Uplifted_Blocks) != bool) or (type(Tilting) != bool) or (type(Gradient) != bool):
    raise Warning("The following variables must be either True or False: Periodic_Uplift, Normal_Fault, Uplifted_Blocks, Tilting, Gradient")
if (np.count_nonzero([Uplifted_Blocks, Tilting, Gradient]) > 1):
    raise Warning("Only one of the following variables can be True: Uplifted_Blocks, Tilting, Gradient. Otherwise they must all be False")
U = np.ones(mg.number_of_nodes) * U_1
if Uplifted_Blocks == True:
    print("Uplift function = Uplifted_Blocks")
    if U_Edge_2 > 0:
        U[np.where(mg.node_x>U_Edge_2)] = U_2 
    if U_Edge_3 > 0:
        U[np.where(mg.node_x>U_Edge_3)] = U_3 
    if U_Edge_4 > 0:
        U[np.where(mg.node_x>U_Edge_4)] = U_4
    if U_Edge_5 > 0:
        U[np.where(mg.node_x>U_Edge_5)] = U_5 
if Tilting == True:    
    print("Uplift function = Tilting")
    if Tilting_Direction == 'S':
        U = U_1 + ( (Tilting_High_Uplift_Rate - U_1) / (number_of_rows - 1) ) * ( mg.node_y / dxy )
    if Tilting_Direction == 'N':
        U = U_1 + ( (Tilting_High_Uplift_Rate - U_1) / (number_of_rows - 1) ) * ( mg.node_y[::-1] / dxy )
    if Tilting_Direction == 'E':
        U = U_1 + ( (Tilting_High_Uplift_Rate - U_1) / (number_of_columns - 1) ) * ( mg.node_x / dxy )
    if Tilting_Direction == 'W':
        U = U_1 + ( (Tilting_High_Uplift_Rate - U_1) / (number_of_columns - 1) ) * ( mg.node_x[::-1] / dxy ) 
if Gradient == True:
    print("Uplift function = Gradient")
    if Gradient_Axis == 'N-S':
        node_array = np.arange(0, int((number_of_rows / 2)*dxy), dxy)
        gradient = ( (Gradient_High_Uplift_Rate - U_1) / ((number_of_rows / 2)-1) )
        q = 0
        for i in node_array:
            U[np.where(mg.node_y == i)] = U_1 + gradient * ( node_array[q] / dxy )
            q += 1
        node_array_2 = np.arange( int(((number_of_rows / 2)*dxy)), (mg.node_y[-1]+dxy), dxy )
        q = int((number_of_rows / 2)-1)
        for i in node_array_2:
            U[np.where(mg.node_y == i)] = U_1 + ( (Gradient_High_Uplift_Rate - U_1) / ((number_of_rows / 2)-1) ) * ( node_array[q] / dxy )
            q-= 1      
    if Gradient_Axis == 'E-W':
        node_array = np.arange(0, int((number_of_columns / 2)*dxy), dxy)
        gradient = ( (Gradient_High_Uplift_Rate - U_1) / ((number_of_columns / 2)-1) )
        q = 0
        for i in node_array:
            U[np.where(mg.node_x == i)] = U_1 + gradient * ( node_array[q] / dxy )
            q += 1
        node_array_2 = np.arange( int(((number_of_columns / 2)*dxy)), (mg.node_x[-1]+dxy), dxy )
        q = int((number_of_columns / 2)-1)
        for i in node_array_2:
            U[np.where(mg.node_x == i)] = U_1 + ( (Gradient_High_Uplift_Rate - U_1) / ((number_of_columns / 2)-1) ) * ( node_array[q] / dxy )
            q-= 1
U_prepulse[mg.core_nodes] = U[mg.core_nodes]
U_postpulse[mg.core_nodes] = U[mg.core_nodes] + U_pulse
if Normal_Fault == True:
    print("Uplift function = Normal_Fault")
    fault_dict = {'faulted_surface': 'topographic__elevation',
                  'fault_dip_angle': Dip_Angle,
                  'fault_trace': {'y1': y1,
                                  'x1': x1,
                                  'y2': y2,
                                  'x2': x2},
                  'include_boundaries': True}
    nf = NormalFault(mg, **fault_dict)
    nf.faulted_nodes.reshape(mg.shape)
U_dict = {"U_1": U_1}
if np.count_nonzero([Periodic_Uplift, Normal_Fault, Uplifted_Blocks, Tilting, Gradient]) == 0:
    print("No uplift functions active")
if Periodic_Uplift == True: 
    print("Uplift function = Periodic_Uplift")
    Time_averaged_U = ( (U_pulse * U_pulse_duration) + (U_1 * (Period - U_pulse_duration)) ) / Period
    U_dict.update({"Time_averaged_Uplift_Rate": Time_averaged_U})
if Normal_Fault == True:
    Faulted_U = (Slip / Recurrence) * math.sin(Dip_Angle * (math.pi / 180))
    U_dict.update({"Faulted_Uplift_Rate": Faulted_U})
if Uplifted_Blocks == True:    
     U_dict.update({"Maximum_Block_Uplift_Rate" : np.max(U)})
if Tilting == True:
    U_dict.update({"Tilting_High_Uplift_Rate": Tilting_High_Uplift_Rate})
if Gradient == True:
    U_dict.update({"Gradient_High_Uplift_Rate": Gradient_High_Uplift_Rate})
Max_U = max(U_dict.values())





#%% BLOCK 5: RESET THE MODEL TIME AND THE PLOTTING TIMERS TO 0

total_time = 0 
Plot_Ticker = 0
Export_DEM_Ticker = 0
U_Ticker = 0
Fault_Ticker = 0               





#%% BLOCK 6: TIME PARAMETERS

# ENTER VARIABLES ############################################################
##############################################################################

dt = 100       
tmax = 1E5  

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------
t = np.arange(0, tmax, dt) 





#%% BLOCK 7: PLOTTING OPTIONS

# ENTER VARIABLES ############################################################
##############################################################################

Directory = 'C:/Users/csheehan1/Desktop/'
Export_format = 'png'

Plot_interval = 1E3
Export_DEM_Interval = 5E3
number_of_watersheds = 1
min_drainage_area = 10000
main_channel_only = False

DEM_Image = True
Slope_Area = False
Channel_Profile = False
Channel_Map = True
Ksn_Profile = False
Ksn_Map = False
Chi_Profile = False
Chi_Map = False
Erosion_Rate_Profile = False
Erosion_Rate_Map = False
DZ_DT_Profile = False
DZ_DT_Map = True
Ksp_Profile = False
Ksp_Map = False
Uplift_Profile = False
Uplift_Map = False
Timeseries = True

Export_DEM = True

Terrain_3D = False

# CHOOSE FROM ONE OF THESE THREE COLOR SCALING OPTIONS. IF "Color_Scaling_Prescribed" 
# IS SELECTED, ENTER A VALUE FOR "Max_color_scale"
Color_Scaling_Automatic = False
Color_Scaling_Updated = False
Color_Scaling_Prescribed = True
Max_color_scale = 360 

# CHOOSE FROM ONE OF THESE TWO ELEVATION EXAGGERATION OPTIONS. IF "Elevation_Exaggeration_Prescribed" 
# IS SELECTED, ENTER A CALUE FOR "Exaggeration_factor"
Elevation_Exaggeration_Automatic = False
Elevation_Exaggeration_Prescribed = True
Exaggeration_factor = 10

##############################################################################
# ENTER VARIABLES ############################################################

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------
Plot_Ticker = 0
Export_DEM_Ticker = 0





#%% BLOCK 8: TIME LOOP

#OPERATORS-->DO_NOT_EDIT_ANYTHING_BELOW_THIS_LINE-----------------------------
if path.exists(str(Directory)+'/TerrainSandbox') == False:
    os.mkdir(str(Directory)+'/TerrainSandbox')
ceiling = (Max_U / Min_Ksp) * 10
frr = FlowAccumulator(mg, flow_director='D8')
pffr = PriorityFloodFlowRouter(mg, depression_handler = 'fill')
if Spatially_Zoned == True or Spatially_Uniform == True:
    spr = StreamPowerEroder(mg, K_sp=Ksp, m_sp=m_sp, n_sp=n_sp)
    dfn = LinearDiffuser(mg, linear_diffusivity=Khs)
if Tilted_Rocks == True or Folded_Rocks == True:
    spr = StreamPowerEroder(mg, K_sp='K_sp', m_sp=m_sp, n_sp=n_sp)
    dfn = LinearDiffuser(mg, linear_diffusivity='K_hs')
previous_zr[mg.core_nodes] = zr[mg.core_nodes]
for ti in t:
    previous_zr_plot = np.ones(np.size(mg.nodes)) * zr
    if Periodic_Uplift == True:
        if U_Ticker < Period:
            U_Ticker = U_Ticker + dt
        if U_Ticker >= Period:
            U_Ticker = 0
        if U_Ticker < Period-U_pulse_duration:
            U[mg.core_nodes] = U_prepulse[mg.core_nodes]
        if U_Ticker > Period-U_pulse_duration:
            U[mg.core_nodes] = U_postpulse[mg.core_nodes]                
    if Normal_Fault == True:
        Fault_Ticker += dt
        if Fault_Ticker < Recurrence:
            #U_Plot[mg.nodes] = U[mg.nodes]
            []
        else:
            nf.run_one_earthquake(dz=Slip)
            zr.reshape(mg.shape)
            previous_zr_plot.reshape(mg.shape)
            #U_Plot[mg.nodes] = U[mg.nodes]
            #U_Plot[np.where(nf.faulted_nodes == True)] += Slip
            Fault_Ticker = 0
    zr[mg.core_nodes] += U[mg.core_nodes]*dt     
    U_Plot = ((np.ones(np.size(mg.nodes)) * zr) - previous_zr_plot) / dt     
    dfn.run_one_step(dt)                                
    # frr.run_one_step()  
    pffr.run_one_step()                  
    spr.run_one_step(dt)
    if Tilted_Rocks == True or Folded_Rocks == True:
        dz_ad = np.zeros(mg.size('node'))
        dz_ad[mg.core_nodes] = U[mg.core_nodes]*dt
        lith.dz_advection=dz_ad
        lith.run_one_step()
    dz_dt[mg.core_nodes] = (zr[mg.core_nodes] - previous_zr[mg.core_nodes]) / dt
    erosion_rate[mg.core_nodes] = (previous_zr[mg.core_nodes] - (zr[mg.core_nodes] - (U[mg.core_nodes]*dt))) / dt
    previous_zr[mg.core_nodes] = zr[mg.core_nodes]   
    
    
    if total_time == 0:
        time = total_time
        mean_elev = np.mean(zr)
        max_elev = np.max(zr)
        uplift_rate = U_1
        spe = Ksp_1
    else:
        time = np.append(time, total_time)
        mean_elev = np.append(mean_elev, np.mean(zr))
        max_elev = np.append(max_elev, np.max(zr))
        uplift_rate = np.append(uplift_rate, U_1)    
        spe = np.append(spe, Ksp_1) 
    
    
    total_time += dt     
    print(total_time)
    Plot_Ticker += dt
    Export_DEM_Ticker += dt
    if Plot_Ticker == Plot_interval:
        print('Exporting figures... Please be patient!')
        if main_channel_only == True:
            prf = ChannelProfiler(mg, number_of_watersheds=number_of_watersheds, main_channel_only=True, minimum_channel_threshold=min_drainage_area)
        if main_channel_only == False:
            prf = ChannelProfiler(mg, number_of_watersheds=number_of_watersheds, main_channel_only=False, minimum_channel_threshold=min_drainage_area)
        if total_time == Plot_Ticker:
            sf = SteepnessFinder(mg, reference_concavity=m_sp/n_sp, min_drainage_area=min_drainage_area) # NEED TO FIX! Currently breaks if you don't export images during first run of Cell 8 but then export them during later runs of Cell 8
            cf_check = 'cf' in locals()
            if cf_check == False:
                cf = ChiFinder(mg, reference_concavity=m_sp/n_sp, min_drainage_area=min_drainage_area)
        if DEM_Image == True:
            if path.exists(str(Directory)+'/TerrainSandbox/DEM_Image') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/DEM_Image')
            plt.ioff()
            fig = plt.figure(1)         
            imshow_grid(mg, 'topographic__elevation', grid_units=('m', 'm'), var_name="Elevation (m)", cmap='terrain', allow_colorbar=True)
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/DEM_Image/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if Slope_Area == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Slope_Area') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Slope_Area')
            plt.ioff()
            fig = plt.figure(2)
            for i, outlet_id in enumerate(prf.data_structure):
                for j, segment_id in enumerate(prf.data_structure[outlet_id]):
                    if j == 0:
                        label = "channel {i}".format(i=i + 1)
                    else:
                        label = '_nolegend_'
                    segment = prf.data_structure[outlet_id][segment_id]
                    profile_ids = segment["ids"]
                    color = segment["color"]
                    plt.loglog(mg.at_node["drainage_area"][profile_ids][2 : -1], mg.at_node["topographic__steepest_slope"][profile_ids][2 : -1], '.', color=color, label=label)
            
            plt.ylim([1E-4, 1E-1])
            plt.xlim([1E5, 1E10])
            
            #plt.legend(loc="lower left")
            plt.xlabel("Drainage Area (m^2)")
            plt.ylabel("Channel Slope [m/m]")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Slope_Area/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if Channel_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Channel_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Channel_Profile')
            plt.ioff()
            fig = plt.figure(3)
            prf.plot_profiles(field='topographic__elevation', xlabel='Distance Upstream (m)', ylabel='Elevation (m)')
            
            plt.xlim([0, 70000])
            plt.ylim([0, 360])
            
            title_text = '$Year$='+str(total_time)
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Channel_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=300)
            plt.close(fig)
        if Channel_Map == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Channel_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Channel_Map')
            plt.ioff()
            fig = plt.figure(4)
            prf.plot_profiles_in_map_view()
            title_text = '$Year$='+str(total_time)
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Channel_Map/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=300)
            plt.close(fig)
        if Ksn_Profile == True:
            prf.run_one_step()
            sf.calculate_steepnesses()
            if path.exists(str(Directory)+'/TerrainSandbox/Ksn_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Ksn_Profile')
            plt.ioff()
            fig = plt.figure(5)
            for i, outlet_id in enumerate(prf.data_structure):
                for j, segment_id in enumerate(prf.data_structure[outlet_id]):
                    if j == 0:
                        label = "channel {i}".format(i=i + 1)
                    else:
                        label = '_nolegend_'
                    segment = prf.data_structure[outlet_id][segment_id]
                    profile_ids = segment["ids"]
                    distance_upstream = segment["distances"]
                    color = segment["color"]
                    plt.plot(distance_upstream[2 : -1], mg.at_node["channel__steepness_index"][profile_ids][2 : -1], '.', color=color, label=label)
            
            plt.xlim([0, 70000])
            plt.ylim([0, 250])
            
            plt.xlabel("Distance Upstream (m)")
            plt.ylabel("Steepness Index")
            plt.legend(loc="upper left")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Ksn_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=300)
            plt.close(fig)
        if Ksn_Map == True:
            prf.run_one_step()
            sf.calculate_steepnesses()
            if path.exists(str(Directory)+'/TerrainSandbox/Ksn_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Ksn_Map')
            plt.ioff()
            fig = plt.figure(6)
            #imshow_grid(mg, "channel__steepness_index", grid_units=("m", "m"), var_name="Steepness Index", cmap="jet")
            Ksn = np.ones(np.size(mg.nodes)) * mg.at_node["channel__steepness_index"]
            Ksn[np.where(mg.node_x >= ((dxy * number_of_columns) - (dxy * 2)))] = 0
            Ksn[np.where(mg.node_y >= ((dxy * number_of_rows) - (dxy * 2)))] = 0
            Ksn[np.where(mg.node_x <= dxy * 2)] = 0
            Ksn[np.where(mg.node_y <= dxy * 2)] = 0
            imshow_grid(mg, Ksn, grid_units=("m", "m"), var_name="Steepness Index", cmap="jet")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Ksn_Map/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=300)
            plt.close(fig)
        if Chi_Profile == True:
            prf.run_one_step()
            cf.calculate_chi()
            if path.exists(str(Directory)+'/TerrainSandbox/Chi_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Chi_Profile')
            plt.ioff()
            fig = plt.figure(7)
            for i, outlet_id in enumerate(prf.data_structure):
                for j, segment_id in enumerate(prf.data_structure[outlet_id]):
                    if j == 0:
                        label = "channel {i}".format(i=i + 1)
                    else:
                        label = '_nolegend_'
                    segment = prf.data_structure[outlet_id][segment_id]
                    profile_ids = segment["ids"]
                    color = segment["color"]
                    plt.plot(mg.at_node["channel__chi_index"][profile_ids], mg.at_node["topographic__elevation"][profile_ids], color=color, label=label)
            plt.xlabel("Chi Index (m)")
            plt.ylabel("Elevation (m)")
            plt.legend(loc="lower right")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Chi_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=300)
            plt.close(fig)
        if Chi_Map == True:
            prf.run_one_step()
            cf.calculate_chi()
            if path.exists(str(Directory)+'/TerrainSandbox/Chi_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Chi_Map')
            plt.ioff()
            fig = plt.figure(8)
            imshow_grid(mg, "channel__chi_index", grid_units=("m", "m"), var_name="Chi Index", cmap="jet")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Chi_Map/'+str(total_time)+'.'+Export_format,  format=Export_format,  dpi=300)
            plt.close(fig)
        if Erosion_Rate_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Erosion_Rate_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Erosion_Rate_Profile')
            plt.ioff()
            fig = plt.figure(3)
            prf.plot_profiles(field='erosion_rate', xlabel='Distance Upstream (m)', ylabel='Erosion Rate (m/yr)')
            title_text = '$Year$='+str(total_time)
            plt.title(title_text)
            plt.grid(linestyle='--')
            axes = plt.gca()
            axes.set_xlim([dxy,None])
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Erosion_Rate_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if Erosion_Rate_Map == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Erosion_Rate_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Erosion_Rate_Map')
            plt.ioff()
            fig = plt.figure(8)
            imshow_grid(mg, "erosion_rate", grid_units=("m", "m"), var_name="Erosion Rate (m/yr)", cmap="jet")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Erosion_Rate_Map/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if DZ_DT_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/DZDT_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/DZDT_Profile')
            plt.ioff()
            fig = plt.figure(3)
            prf.plot_profiles(field='elevational_change', xlabel='Distance Upstream (m)', ylabel='Rate of Elevational Change (m/yr)')
            title_text = '$Year$='+str(total_time)
            plt.title(title_text)
            plt.grid(linestyle='--')
            axes = plt.gca()
            axes.set_xlim([dxy,None])
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/DZDT_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if DZ_DT_Map == True:
            #prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/DZDT_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/DZDT_Map')
            plt.ioff()
            fig = plt.figure(8)
            imshow_grid(mg, "elevational_change", grid_units=("m", "m"), var_name="Rate of Elevational Change (m/yr)", cmap="seismic_r", symmetric_cbar = True, vmin = -0.001, vmax = 0.001) # vmin = np.max(dz_dt[mg.core_nodes])*-1, vmax = np.max(dz_dt[mg.core_nodes])
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/DZDT_Map/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if Ksp_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Ksp_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Ksp_Profile')
            plt.ioff()
            fig = plt.figure(117)
            if Spatially_Uniform == True or Spatially_Zoned == True:
                prf.plot_profiles(field=Ksp, xlabel='Distance Upstream (m)', ylabel='Ksp')
            if Tilted_Rocks == True or Folded_Rocks == True:
                prf.plot_profiles(field='K_sp', xlabel='Distance Upstream (m)', ylabel='Ksp')
            title_text = '$Year$='+str(total_time) 
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Ksp_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if Ksp_Map == True:
            #prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Ksp_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Ksp_Map')
            plt.ioff()
            fig = plt.figure(8)
            if Spatially_Uniform == True or Spatially_Zoned == True:
                imshow_grid(mg, Ksp, grid_units=("m", "m"), var_name="Ksp", cmap="jet")
            if Tilted_Rocks == True or Folded_Rocks == True:
                imshow_grid(mg, 'K_sp', grid_units=("m", "m"), var_name="Ksp", cmap="PuOr")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Ksp_Map/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)  
        if Uplift_Profile == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Uplift_Profile') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Uplift_Profile')
            plt.ioff()
            fig = plt.figure(118)
            prf.plot_profiles(field=U_Plot, xlabel='Distance Upstream (m)', ylabel='Uplift Rate (m/yr)')
            title_text = '$Year$='+str(total_time) 
            plt.title(title_text)
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Uplift_Profile/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if Uplift_Map == True:
            prf.run_one_step()
            if path.exists(str(Directory)+'/TerrainSandbox/Uplift_Map') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Uplift_Map')
            plt.ioff()
            fig = plt.figure(300)
            imshow_grid(mg, U_Plot, grid_units=("m", "m"), var_name="Uplift Rate (m/yr)", cmap="jet")
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Uplift_Map/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
        if Terrain_3D == True:
            if path.exists(str(Directory)+'/TerrainSandbox/Terrain_3D') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Terrain_3D')
            plt.ioff()
            fig = plt.figure(9)
            ax = plt.axes(projection ='3d') 
            Z = mg.at_node['topographic__elevation'].reshape(mg.shape)
            if Color_Scaling_Updated == True:
                color = cm.terrain((Z-Z.min())/(Z.max()-Z.min()))                  # USE THIS LINE IF YOU WANT THE COLOR SCALING TO UPDATE EACH TIMESTEP ACCORDING TO THE MAX ELEVATION 
            if Color_Scaling_Prescribed == True:
                color = cm.terrain((Z-Z.min())/(Max_color_scale - Z.min()))                     # USE THIS LINE IF YOU WANT THE COLOR SCALING TO STAY CONSTANT THROUGH THE RUN. ENTER DESIRED MAX COLOR ELEVATION WHERE Z.max() SHOULD BE IN THE LINE ABOVE
            if Color_Scaling_Automatic == True:
                color = cm.terrain((Z-Z.min())/(ceiling-Z.min()))                   # USE THIS LINE IF YOU WANT TO AUTOSCALE THE RUN USING "ceiling". NOTE THAT IF THE Ksp OR U PARAMETERS ARE CHANGED MID-RUN, THE SCALING WILL ALSO CHANGE.
            surf = ax.plot_surface(mg.node_x.reshape(mg.shape), mg.node_y.reshape(mg.shape), Z, rstride=1, cstride=1, facecolors=color, edgecolor='black', linewidth=1, antialiased=True)
            ax.view_init(elev=30, azim=-120)                                    # USE THIS LINE TO SET THE FIGURE EYELEVEL AZIMUTH
            ax.set_xlabel('meters')
            ax.set_ylabel('meters')
            ax.set_zlabel('Elevation')
            title_text = '$Year$='+str(total_time)  
            plt.title(title_text)
            plt.tight_layout()
            if Elevation_Exaggeration_Prescribed == True:
                ax.set_zlim(0, dxy * number_of_rows / Exaggeration_factor)
            if Elevation_Exaggeration_Automatic == True:
                ax.set_zlim(0, ceiling * 10)
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Terrain_3D/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
            
            
        if Timeseries == True:
            if path.exists(str(Directory)+'/TerrainSandbox/Timeseries') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Timeseries')
            #
            if path.exists(str(Directory)+'/TerrainSandbox/Timeseries/elev') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Timeseries/elev')
            plt.ioff()
            fig = plt.figure(118)
            plt.plot(time, mean_elev)
            plt.plot(time, max_elev)
            #title_text = '$Year$='+str(total_time) 
            #plt.title(title_text)
            
            plt.xlim([0, 1E6])
            plt.ylim([0, 360])
            
            plt.xlabel('Years')
            plt.ylabel('Elevation (m)')
            plt.grid(linestyle='--')
            plt.legend(['Mean elevation', 'Max elevation'])
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Timeseries/elev/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
            #
            if path.exists(str(Directory)+'/TerrainSandbox/Timeseries/uplift_rate') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Timeseries/uplift_rate')
            plt.ioff()
            fig = plt.figure(118)
            plt.plot(time, uplift_rate)
            #title_text = '$Year$='+str(total_time) 
            #plt.title(title_text)
            
            #plt.xlim([0, 1E6])
            #plt.ylim([0, 360])
            
            plt.xlabel('Years')
            plt.ylabel('Uplift rate (m / yr)')
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Timeseries/uplift_rate/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
            #
            if path.exists(str(Directory)+'/TerrainSandbox/Timeseries/spe') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Timeseries/spe')
            plt.ioff()
            fig = plt.figure(118)
            plt.plot(time, spe)
            #title_text = '$Year$='+str(total_time) 
            #plt.title(title_text)
            
            #plt.xlim([0, 1E6])
            #plt.ylim([0, 360])
            
            #plt.xlabel('Years')
            #plt.ylabel('Uplift rate (m / yr)')
            plt.grid(linestyle='--')
            plt.tight_layout()
            fig.savefig(str(Directory)+'/TerrainSandbox/Timeseries/spe/'+str(total_time)+'.'+Export_format,  format=Export_format, dpi=300)
            plt.close(fig)
            #
            
        Plot_Ticker = 0
    if Export_DEM == True:
        if total_time == Export_DEM_Ticker:
            if path.exists(str(Directory)+'/TerrainSandbox') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox')
        if Export_DEM_Ticker == Export_DEM_Interval:
            if path.exists(str(Directory)+'/TerrainSandbox/Export_DEM') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Export_DEM')
            write_esri_ascii(str(Directory)+'/TerrainSandbox/Export_DEM/'+str(total_time)+'.asc', mg, names='topographic__elevation')
            
        if Export_DEM_Ticker == Export_DEM_Interval:
            if path.exists(str(Directory)+'/TerrainSandbox/Export_DZDT') == False:
                os.mkdir(str(Directory)+'/TerrainSandbox/Export_DZDT')
            write_esri_ascii(str(Directory)+'/TerrainSandbox/Export_DZDT/'+str(total_time)+'.asc', mg, names='elevational_change')
            
            
            Export_DEM_Ticker = 0
print('')
print('Complete!')




