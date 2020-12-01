# -*- coding: utf-8 -*-
"""
Created on: xxxx
@author   : Anon

NAME
    Farm_Evalautor.py
    
PYTHON VERSION   
    3.7.3 
    
DESCRIPTION
    Calculates Annual Energy Production (AEP) of a Wind Farm
    ============================================================    
    
    Farm_Evalautor.py is a python file that calculates AEP (GWh)
    of a certain arrangement of wind turbines in a wind farm, under 
    given annual wind conditions. 
    
    The code in this script for wake-effect modeling is based on
    standard Jensen (PARK) model. 
    I. Katic, J. Hojstrup and N. Jensen, "A simple model for cluster 
    efficiency," in European Wind Energy Association Conference and 
    Exhibition, 1986.
    
    As its inputs, the code takes three data files containing info 
    about:
    - Turbine Locations
    - Turbine Power Curve
    - Annual Wind Conditions
    
PACKAGE LIST
    You may need to install the packages Tqdm and Shapely in your
    python distribution. These are not pre-installed. 
    ============================================================
    Packages Used:
    Numpy
    Pandas
    Numba
    Tqdm
    Shapely
    math (built-in)
"""

# Module List
import numpy  as np
import pandas as pd    
from   numba  import njit                      # For some speed gains                    
from   math   import radians as DegToRad       # Degrees to radians Conversion

from shapely.geometry import Point             # Used in constraint checking
from shapely.geometry.polygon import Polygon

from tqdm import tqdm                          # For Progressbar
# Use from tqdm import tqdm_notebook as tqdm for running in Jupyter Notebook  
# Use from tqdm import tqdm for running in IDE Environ 



def getTurbLoc(turb_loc_file_name):
    """ 
    -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
    
    Returns x,y turbine coordinates
    
    :Called from
        main function
    
    :param
        turb_loc_file_name - Turbine Loc csv file location
        
    :return
        2D array
    """
    
    df = pd.read_csv(turb_loc_file_name, sep=',')
    turb_coords = df.to_numpy(dtype = np.float32)
    return(turb_coords)


def loadPowerCurve(power_curve_file_name):
    """
    -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
    
    Returns a 2D numpy array with information about
    turbine thrust coeffecient and power curve of the 
    turbine for given wind speed
    
    :called_from
        main function
    
    :param
        power_curve_file_name - power curve csv file location
        
    :return
        Returns a 2D numpy array with cols Wind Speed (m/s), 
        Thrust Coeffecient (non dimensional), Power (MW)
    """
    powerCurve = pd.read_csv(power_curve_file_name, sep=',')
    powerCurve = powerCurve.to_numpy(dtype = np.float32)
    return(powerCurve)
    

def binWindResourceData(wind_data_file_name):
    r"""
    -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
    
    Loads the wind data. Returns a 2D array with shape (36,15). 
    Each cell in  array is a wind direction and speed 'instance'. 
    Values in a cell correspond to probability of instance
    occurence.  
    
    :Called from
        main function
        
    :param
        wind_data_file_name - Wind Resource csv file  
        
    :return
        2D array with estimated probabilities of wind instance occurence. 
        Along: Row-direction (drct), Column-Speed (s)
        
                      |0<=s<2|2<=s<4| ...  |26<=s<28|28<=s<30|
        |_____________|______|______|______|________|________|
        | drct = 360  |  --  |  --  |  --  |   --   |   --   |
        | drct = 10   |  --  |  --  |  --  |   --   |   --   |
        | drct = 20   |  --  |  --  |  --  |   --   |   --   |
        |   ....      |  --  |  --  |  --  |   --   |   --   |
        | drct = 340  |  --  |  --  |  --  |   --   |   --   |
        | drct = 350  |  --  |  --  |  --  |   --   |   --   |        
    """
    
    # Load wind data. Then, extracts the 'drct', 'sped' columns
    df = pd.read_csv(wind_data_file_name)
    wind_resource = df[['drct', 'sped']].to_numpy(dtype = np.float32)
    
    # direction 'slices' in degrees
    slices_drct   = np.roll(np.arange(10, 361, 10, dtype=np.float32), 1)
    ## slices_drct   = [360, 10.0, 20.0.......340, 350]
    n_slices_drct = slices_drct.shape[0]
    
    # speed 'slices'
    slices_sped   = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 
                        18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0]
    n_slices_sped = len(slices_sped)-1

    
    # placeholder for binned wind
    binned_wind = np.zeros((n_slices_drct, n_slices_sped), 
                           dtype = np.float32)
    
    # 'trap' data points inside the bins. 
    for i in range(n_slices_drct):
        for j in range(n_slices_sped):     
            
            # because we already have drct in the multiples of 10
            foo = wind_resource[(wind_resource[:,0] == slices_drct[i])] 

            foo = foo[(foo[:,1] >= slices_sped[j]) 
                          & (foo[:,1] <  slices_sped[j+1])]
            
            binned_wind[i,j] = foo.shape[0] 
    
    wind_inst_freq   = binned_wind/np.sum(binned_wind)
    
    return(wind_inst_freq)


def rotatedFrame(turb_coords, wind_drct):
    """
    -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
    
    Rotates euclidean coordinates to downwind-crosswind coordinates.
    Rotate the axes such that the wind flow direction aligns with the
    positive x-axis. 
    
    :called from
        partAEP
        
    :param
        turb_coords   - 2D array. turbine euclidean x,y coordinates
        wind_drct     - wind direction in degrees. 
                        [{'N':0},{'E':90},{'S':180},{'W':270}]
        
    :return    
        Rotated turbine locations. 
    """
    
    # so that the wind flow direction aligns with the +ve x-axis.
    wind_drct = wind_drct-90
    
    # Convert inflow wind direction from degrees to radians
    wind_drct = DegToRad(wind_drct)

    # Contants for coordinate transformation 
    cos_dir = np.cos(wind_drct)
    sin_dir = np.sin(wind_drct)
    
    # Coordinate Transformation. Rotate coordinates to downwind, crosswind coordinates 
    rotate_coords      =    np.zeros((turb_coords.shape[0],2), dtype=np.float32)
    rotate_coords[:,0] =    (turb_coords[:,0] * cos_dir) - (turb_coords[:,1] * sin_dir)
    rotate_coords[:,1] =    (turb_coords[:,0] * sin_dir) + (turb_coords[:,1] * cos_dir)
    
    return(rotate_coords)


@njit            # jensenParkWake is numba compiled function. Gives speed gains. 
def jensenParkWake(n_turbs, turb_diam, rotate_coords, power_curve, wind_sped):
    """
    -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
    
    For each turbine, returns the total velocity deficit due to the wakes from 
    all contributing upstream turbines for given particular wind instance i.e. 
    wind_sped and drct (indirectly represented by rotated rotate_coords). 
    Uses Jensen's (PARK) model. 
    
    :called from
        partAEP
        
    :param
        n_turbs      - Total number of turbines
        turb_diam    - Diameter of the turbine (m)
        rotate_coords- Rotated turbine coordinates
        power_curve  - For estimating thrust coeff (C_t) 
        wind_sped    - Wind Speed (m/s)
        
    :return
        Velocity deficit for each turbine. returns array
    """
    # turbine radius
    turb_rad = turb_diam/2

    # we use power_curve data as look up to estimate the thrust coeff.
    # of the turbine for the corresponding closest matching wind speed
    idx_foo  = np.argmin(np.abs(power_curve[:,0] - wind_sped))
    C_t      = power_curve[idx_foo,1]
    
    # Wake decay constant kw for the offshore case
    kw = 0.05 
    
    # velocity deficit suffered by each turbine for this particular wind instance
    # impact_on_ibyj - placeholder to calc vel deficit from all turbs on i                 
    impact_on_ibyj = np.zeros((n_turbs,n_turbs), dtype=np.float32) 
        
    # i - target turbine
    for i in range(n_turbs):   
        
        # looping over all other turbs to check their effect
        for j in range(n_turbs):             
            
            # Calculate the x-dist and the y-offset 
            # (wrt downwind/crosswind coordinates)
            x = rotate_coords[i,0] - rotate_coords[j,0]   
            y = rotate_coords[i,1] - rotate_coords[j,1]   
            
            # Naturally, no wake effect of turbine on itself
            if i!=j: 
                
                # either j not an upstream turbine or wake not happening 
                # on i because its outside of the wake region of j                    
                if x<=0 or np.abs(y) > (turb_rad + kw*x):  
                    impact_on_ibyj[i,j] = 0.0
                
                # otherwise, at target i, wake is happening due to j
                else:                               
                    impact_on_ibyj[i,j] = (1-np.sqrt(1-C_t))*((turb_rad/(turb_rad + kw*x))**2)
            
    # Calculate Total vel deficit from all upstream turbs, using sqrt of sum of sqrs
    sped_deficit = np.sqrt(np.sum(impact_on_ibyj**2, axis = 1))

    return(sped_deficit)


def partAEP(n_turbs, turb_diam, turb_coords, power_curve, wind_drct, wind_sped):
    """
    -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
    
    Returns total power produced by all the turbine for this particular 
    'wind instance' represented by wind_drct, wind_sped.
    
    :called from
        totalAEP
        
    :param
        n_turbs      - Total number of turbines
        turb_diam    - Diameter of the turbine (m)
        turb_coords  - 2D array turbine x,y coordinates
        power_curve  - For estimating thrust coeff, C_t
        wind_drct    - wind direction in degrees
        wind_sped    - Wind Speed (m/s)
        
    :return
        total power produced for single wind instance  
    """
    # For given wind_drct rotate coordinates to downwind/crosswind 
    rotate_coords = rotatedFrame(turb_coords, wind_drct)
    
    # Use the jensen park wake model to calc speed deficits by wake
    sped_deficit = jensenParkWake(n_turbs, turb_diam, rotate_coords, power_curve, wind_sped)
    
    # Placeholder for storing power output of turbines
    turb_pwr = np.zeros(n_turbs, dtype=np.float32)

    # calculate the individual turbine power for effective wind speed
    for i in range(n_turbs):
        
        # Effective windspeed due to the happening wake
        wind_sped_eff = wind_sped*(1.0 - sped_deficit[i])
        
        # we use power_curve data as look up to estimate the power produced
        # by the turbine for the corresponding closest matching wind speed
        idx_foo = np.argmin(np.abs(power_curve[:,0] - wind_sped_eff))
        pwr     = power_curve[idx_foo,2]
        
        turb_pwr[i] = pwr
           
    # Sum the power from all turbines for this wind instance
    power = np.sum(turb_pwr)
    
    return power


def totalAEP(turb_diam, turb_coords, power_curve, wind_inst_freq):
    """
    -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
    
    Calculates the wind farm AEP
    
    :called from
        main 
        
    :param
        turb_diam       - Diameter of the turbine (m)
        turb_coords     - 2D array turbine euclidean coordinates
        power_curve     - For estimating thrust coeff (C_t)
        wind_inst_freq  - 2D array with estimated probabilities of 
                          wind instance occurence.(return from func 
                          binWindResourceData)
                          
    :return
        wind farm AEP in Gigawatt Hours, GWh (float)
    """
    # number of turbines
    n_turbs        =   turb_coords.shape[0]
    assert n_turbs ==  50, "Error! Number of turbines is not 50."
    
    # direction 'slices' in degrees
    slices_drct   = np.roll(np.arange(10, 361, 10, dtype=np.float32), 1)
    ## slices_drct   = [360, 10.0, 20.0.......340, 350]
    n_slices_drct = slices_drct.shape[0]
    
    # speed 'slices'
    slices_sped   = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 
                        18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0]
    n_slices_sped = len(slices_sped)-1

    # Power produced by the wind farm from each wind instance
    farm_pwr = np.zeros((wind_inst_freq.shape), dtype=np.float32)
    
    # Looping over every wind instance and calc Power
    # set disable=True for disabling progressbar
    for i in tqdm(range(n_slices_drct), disable=False):
        for j in range(n_slices_sped): 
            
            # take the mid value as effective speed
            wind_drct =  slices_drct[i]
            wind_sped = (slices_sped[j] + slices_sped[j+1])/2
            
            pwr  = partAEP(n_turbs, turb_diam, turb_coords, 
                           power_curve, wind_drct, wind_sped)
            
            farm_pwr[i,j] = pwr
    
    # multiply the respective values with the wind instance probabilities 
    farm_pwr = wind_inst_freq * farm_pwr
    
    # now sum all values
    farm_pwr = np.sum(farm_pwr)
    
    # multiply farm_pwr by hours in a year
    year_hours = 365.*24.
    AEP = year_hours*farm_pwr
    
    # Convert MWh to GWh
    AEP = AEP/1e3
    return(AEP)


def checkConstraints(turb_coords, turb_diam):
    """
    -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
    
    Checks if the turbine configuration satisfies the two
    constraints:(i) perimeter constraint,(ii) proximity constraint 
    Prints which constraints are violated if any. Note that this 
    function does not quantifies the amount by which the constraints 
    are violated if any. 
    
    :called from
        main 
        
    :param
        turb_coords - 2d np array containing turbine x,y coordinates
        turb_diam   - Diameter of the turbine (m)
    
    :return
        None. Prints messages.   
    """
    bound_clrnc      = 50
    prox_constr_viol = False
    peri_constr_viol = False
    
    # create a shapely polygon object of the wind farm
    farm_peri = [(0, 0), (0, 4000), (4000, 4000), (4000, 0)]
    farm_poly = Polygon(farm_peri)
    
    # checks if for every turbine perimeter constraint is satisfied. 
    # breaks out if False anywhere
    for turb in turb_coords:
        turb = Point(turb)
        inside_farm   = farm_poly.contains(turb)
        correct_clrnc = farm_poly.boundary.distance(turb) >= bound_clrnc
        if (inside_farm == False or correct_clrnc == False):
            peri_constr_viol = True
            break
    
    # checks if for every turbines proximity constraint is satisfied. 
    # breaks out if False anywhere
    for i,turb1 in enumerate(turb_coords):
        for turb2 in np.delete(turb_coords, i, axis=0):
            if  np.linalg.norm(turb1 - turb2) < 4*turb_diam:
                prox_constr_viol = True
                break
    
    # print messages
    if  peri_constr_viol  == True  and prox_constr_viol == True:
          print('Somewhere both perimeter constraint and proximity constraint are violated\n')
    elif peri_constr_viol == True  and prox_constr_viol == False:
          print('Somewhere perimeter constraint is violated\n')
    elif peri_constr_viol == False and prox_constr_viol == True:
          print('Somewhere proximity constraint is violated\n')
    else: print('Both perimeter and proximity constraints are satisfied !!\n')
        
    return()
    

if __name__ == "__main__":
    
    # Turbine Specifications.
    # -**-SHOULD NOT BE MODIFIED-**-
    turb_specs    =  {   
                         'Name': 'Anon Name',
                         'Vendor': 'Anon Vendor',
                         'Type': 'Anon Type',
                         'Dia (m)': 100,
                         'Rotor Area (m2)': 7853,
                         'Hub Height (m)': 100,
                         'Cut-in Wind Speed (m/s)': 3.5,
                         'Cut-out Wind Speed (m/s)': 25,
                         'Rated Wind Speed (m/s)': 15,
                         'Rated Power (MW)': 3
                     }
    turb_diam      =  turb_specs['Dia (m)']
    
    
    # Turbine x,y coordinates
    turb_coords    =  getTurbLoc(r'..\Shell_Hackathon Dataset\turbine_loc_test.csv')
    
    # Load the power curve
    power_curve    =  loadPowerCurve('..\Shell_Hackathon Dataset\power_curve.csv')
    
    # Pass wind data csv file location to function binWindResourceData.
    # Retrieve probabilities of wind instance occurence.
    wind_inst_freq =  binWindResourceData(r'..\Shell_Hackathon Dataset\Wind Data\wind_data_2007.csv')
    
    # check if there is any constraint is violated before we do anything. Comment 
    # out the function call to checkConstraints below if you desire. Note that 
    # this is just a check and the function does not quantifies the amount by 
    # which the constraints are violated if any. 
    checkConstraints(turb_coords, turb_diam)
     
    print('Calculating AEP......')
    AEP = totalAEP(turb_diam, turb_coords, power_curve, wind_inst_freq) 
    print('Total power produced by the wind farm is: ', "%.12f"%(AEP), 'GWh')
