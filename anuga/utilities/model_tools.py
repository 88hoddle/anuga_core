"""
 FILE :  MODEL_BUILD_TOOLS_04.py
 DATE:  28/05/2013
 
 Last Change:  Can now define directories by using underscore _ 
 so mannings = 0_05



This script is meant to work with a set of standardised directories so that an ANUGA model script
can be kept extremely minimalist, this module provides functions that will semi-automate the population
of polygons and attributes to:

1. define MESH Refinement
2. define Variation in ROughness in the model
3. define the additional of buildings of various heights
4. define rainfall spatial variation via polygons

To ensure these functions operate it is important to use the
STANDARD_MODEL_CREATE_SCRIPT that will create a set of standardised directories to store model data

The concept is to populate this standard directory structure with the required 'csv' files 
describing the required polygons.

MESH REFINE
EG: c:\Model\Data\Mesh_Refine\10000 and c:\Model\Data\Mesh-Refine\1000
Contain the polygons that describe the mesh refine for triangles no larger 
than 10,000m2 and 1,000m2

BUILDINGS
Similarly buildings within numbered building directories will signify the building height
EG: c:\Model\Data\Buildings\8  and c:\Model\Data\Buildings\12
will contain polygons to describe buildings that are 8 and 12 metres high

ROUGHNESS
The directory structure will be based on identifying for example the roughness directory
under which numbered directories will signify the rougness value.
EG: c:\Model\Data\Roughness\100 and c:\Model\Data\Roughness\10
will contain polygons top describe surface roughness in the model of 0.100 and 0.010

RAINFALL
EG: c:\Model\Data\Rainfall_Polys\Gauge_name_01 and c:\Model\Data\Rainfall_Polys\Gauge_name_02
Will list the polygons that apply rainfall described by the raingauge files in the directory name
options: Either the directory Name is the Raingauge file name, or the Raingauge file (with extension tms) 
is located in the directory and associated with the polygons in the same directory.... or
the raingauge is kept in another directory under FORCEFUNC and there is a pointer to it ???

TO BE RESOLVED 

METHOD OF CALLING:
Identify the directory in which sub directories will appear that have names as numeric values

USAGE FOR REFINE (example):
mesh_refine_directory =join('POLYS','Cycleway')
mesh_refine_defined = get_REFINE_polygon_value_list(mesh_refine_directory) # Returns List of Polygons and Value from Directory Name
print 'Mesh Refine Definitions in Directory....'
print mesh_refine_defined
# Need to Reverse the Lists so that mesh Refine is applied from Coarse to Fine. As default directory read is from smallest to largest number....
mesh_refine_defined.reverse()
print 'Mesh Refine REVERSED HERE....'
#mesh_refine_reversed = mesh_refine_defined.reverse()
print mesh_refine_defined


buildings_directory =join('02POLYS','04_BLDGS')
buildings_defined = get_BUILDING_polygon_value_list(buildings_directory)

such that all polygons under a directory named 45, under directory:- 04_BLDGS  will be assigned a height of 4.5m

Need to add the following lines to Scripts:


from anuga.utilities.model_tools import get_polygon_list_from_files
from anuga.utilities.model_tools import get_polygon_dictionary
from anuga.utilities.model_tools import get_REFINE_polygon_value_list
from anuga.utilities.model_tools import get_ROUGHNESS_polygon_value_list
from anuga.utilities.model_tools import get_BUILDING_polygon_value_list
"""

import os
import glob
import numpy
from anuga.geometry.polygon import read_polygon
from anuga import Boyd_box_operator
from anuga import Boyd_pipe_operator
from anuga import Weir_orifice_trapezoid_operator
from anuga import Inlet_operator



# ---------------------------------------------------------------------------------------------------------

def get_polygon_from_single_file(Rfile):
    fid = open(Rfile)
    lines = fid.readlines()
    fid.close()
    polylist = []
    polygon = []
    polycount = 0
    for line in lines:
        fields = line.split(',')
        #print line
        if line in ('\n', '\r\n'): # Must be a blank Line....
            # Found a line without INcorrect data, assume this signifies the start of a new polygon
            polycount+=1
            #print 'Polygon '+str(polycount)
            polylist.append(polygon)
            #print polygon
            polygon =[]
        else:
            polygon.append([float(fields[0]), float(fields[1])])
    polylist.append(polygon)
    return polylist
# ---------------------------------------------------------------------------------------------------------



def get_polygons_from_Mid_Mif(Rfile):
    """Create List of Polygons from a Directory with a File containing Multiple-Polygons
       
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """
    #print 'in get_polygon_from_Mid_Mif -line 126'
    fid = open(Rfile)
    lines = fid.readlines()
    fid.close()
    # Got the Multi Poly file now process
    polylist = []
    polygon = []
    check_pts_list=[]
    Poly_count=0
    Poly_line_count=0
    Points_in_Poly = 100
    total_lines_in_file= len(lines)
    #print "total number of lines in the Polygons FILE is: ",total_lines_in_file

# ==================================================== FOR LOOP ===========================        
    for i, line in enumerate(lines): 
        if line.strip().startswith('Region'):
            Poly_line_count=0
            check_pts_list=[]
            if Poly_count==0:
                pass
            else:
                polylist.append(polygon)
                #outfid.close()
                #print polygon
                polygon=[]
            Poly_count+=1
            # Create A poly File for each polygon defined in the Multi-Poly file
            #print 'Polygon #',Poly_count
            #poly_write_file="Poly_"+str(Poly_count)+".csv"
            #outfid = open(poly_write_file, 'w')
            #raw_input('Check Output... -line 155')  
            # Instead write to a List
        elif line.strip().startswith('    Pen'):
            pass
        elif line.strip().startswith('    Brush'):
            pass
        else:
            Poly_line_count+=1
            if Poly_line_count > 1 and Poly_line_count <= (Points_in_Poly+1) and Poly_count<>0:
                #print line, #Points_in_Poly,#Poly_line_count
                fields = line.split(' ')
                if line in check_pts_list and Poly_line_count <> Points_in_Poly+1:   # Get rid of any doubled up points NOTE this gets rid of last line !!!
                    #print Poly_line_count, Points_in_Poly+1
                    pass
                else:
                    #outfid.write("%.3f,%.3f\n" % (float(fields[0]),float(fields[1])))
                    polygon.append([float(fields[0]),float(fields[1])])
                    check_pts_list.append(line)
            elif Poly_line_count==1 and Poly_count<>0:
                # read number of points in poly
                #print 'line=',line
                Points_in_Poly=int(line)
                #print 'Number Points in Poly =',Points_in_Poly
                #print polygon
                #raw_input('Check Poly...')
    #print 'End For Loop....-line178'
    polylist.append(polygon)
    #print polylist
    #outfid.close()          
    return polylist          

# ---------------------------------------------------------------------------------------------------------
def get_polygon_list_from_files(dir):
    """Read all polygons found in specified dir and return them in a list
       Called by:
       get_polygon_dictionary
       Purpose:
       To fill a list with all of the polygons read under a specified directory
       CALLS:
       anuga.utilities.polygon.read_polygon
    """
    
    #print 'Reading polygon files from ' + dir
    #print 'This will check the file for Multiple Polygons or read mutiple files with a single polygon per file...' # Need to read files with multiple polys also....
    polylist = []
    for filename in os.listdir(dir):
        Rfile = dir +'/'+filename
        #print Rfile
        
        if Rfile[-4:] == '.svn':  # wHAT DOES THIS DO ??
            continue
        if Rfile[-4:] == '.csv':
            #print 'CSV File'
            polylist.extend(get_polygon_from_single_file(Rfile))
            #polylist.append(polys)
        if Rfile[-4:] == '.mif':
            #print 'MIF File ...'
            #polys = get_polygons_from_Mid_Mif(Rfile)
            polylist.extend(get_polygons_from_Mid_Mif(Rfile))
        #print filename
        #print Rfile
        #raw_input('Hold check file...- line 211')
    #print polylist
    #raw_input('hold at polylist.. -line 213')
    return polylist

# ---------------------------------------------------------------------------------------------------------
def get_polygon_dictionary(dir):
    """Create dictionary of polygons with directory names 
       indicating associated attribute values 
       Called by:
       get_polygon_value_list
       Purpose:
       To Fill a Dictionary with sets of poygons and attribute, from a list of polygons 
       and using the directory name as the attribute
       For Example used to read Mesh Size Directory 1500, using all polygons in the directory
       to create mesh refinement to 1500m2
       CALLS:
       get_polygon_list_from_files
    """
    
    try:
        attribute_values = os.listdir(dir)  # Create the Attribute from the Directory Name
    except:
        msg = 'Directory %s was not found' % dir
        raise Exception(msg)
    D = {}   # Create Empty Dictionary
    for a in attribute_values:
        # How to read a file with multiple polygons ??
        D[a] = get_polygon_list_from_files(os.path.join(dir, a)) # Fill Item [a] in the Dictionary with FIle name and attribute
    return D
# ---------------------------------------------------------------------------------------------------------

# ---- GENERIC POLYGON VALUE LIST Generator
def get_polygon_value_list(dir):
    """Create list of multiple Polygons attributed with a value
       Where the values are obtained from sub directory names based on number and decimal at underscore
       So: Passing Directory ROUGHNESS containing, subs, 0_015, and 0_06 for example
       
       Called by:
       User ANUGA Model SCRIPT
       Purpose:
       CALLS:
       get_polygon_dictionary
    These lists can either be used as interior regions in mesh refinement or as input to Polygon_function
    """
    
    #print 'Read directories of polygons and attributing DIR NAME to polygon'    
    #print 'Naming convention uses the underscore as decimal point eg:0_015, 1000_0'
    D = get_polygon_dictionary(dir)
    polygon_value_list = []
    for key in D:
        try:
            numb_bits = key.split('_')
            attribute = float(numb_bits[0]+'.'+numb_bits[1])
            #print 'Polygon Attribute = ' + str(attribute)
        except:
            print 'Non numerical attributes not yet implemented. I got %s' % key
            return []
        for polygon in D[key]:
            # Create polygon-value pair and append to list for this dir
            pair = [polygon, attribute]
            polygon_value_list.append(pair)
    #print polygon_value_list
    return polygon_value_list
# ---------------------------------------------------------------------------------------------------------


def read_polygon_dir(weight_dict, directory, filepattern='*.csv'):
    """
    In a directory directory looks at all files matching filepattern
    and returns a list of tuples consisting of polygon and a weight 
    """
    pattern = os.path.join(directory, filepattern)
    files = glob.glob(pattern)

    # check that the dictionary contains *all* the files
    
    errors = []
    for f in files:
        try:
            _ = weight_dict[f]
        except KeyError:
            errors.append(f)
            
    if errors:
        msg = ''
        for f in errors:
            msg = msg + ', ' + f
        raise KeyError, 'Files not defined in dictionary: %s' % msg[2:]

    # now get the result list
    result = []
    for f in files:
        result.append((read_polygon(f), weight_dict[f]))
    return result



#define a function with without an attribute
def read_hole_dir_multi_files_with_single_poly(directory, filepattern='*.csv'):
    """
    Looks in a directory, and reads all .csv files as polygons
    and returns a list of polygon 
    """
    pattern = os.path.join(directory, filepattern)
    files = glob.glob(pattern)

    # now get the result list
    result = []
    for f in files:
        result.append(read_polygon(f))
    return result
    
    
    
    
    
# Define a function to read Single File with Multi-polygons
def read_multi_poly_file(multi_P_file):
    """
    Reads a file with multiple polygons, formatted as 
    x,y
    x,y
    x,y
    
    x,y
    x,y
    x,y ...
    
    I.e each poly is defined by x,y position of vertices. New polygon starts after
    a space.
    
    Returns a list of polygons 
    """
    delimiter = ','
    fid = open(multi_P_file)
    lines = fid.readlines()
    fid.close()
    polygon = []
    polygons = []
    for line in lines:
        fields = line.split(delimiter)
        try:
            polygon.append([float(fields[0]), float(fields[1])])
        except:
            # Found a line without correct data, assume this signifies the start of a new polygon
            polygons.append(polygon)
            polygon = []
        
    # Pickup the last polygon
    polygons.append(polygon)

    #print len(polygons)

    #print polygons    
    return polygons



#define a function with without an attribute
def read_hole_dir_single_file_with_multi_poly(directory, filepattern='*.csv'):
    """
    Looks in a directory, and reads 1 .csv file
    containing muliple polygons
    and returns a list of polygon 
    """
    pattern = os.path.join(directory, filepattern)
    files = glob.glob(pattern)

    # now get the result list
    result = []
    for f in files: # For the 1 file
        result.append(read_multi_poly_file(multi_P_file)) # Get the multiple Polygons
    return result


# Define a function to read Single File with Multi-polygons and attribute a value
def read_multi_poly_file_value(multi_P_file,attribute):
    """
    Reads a file with multiple polygons, formatted as 
    x,y
    x,y
    x,y
    
    x,y
    x,y
    x,y ...
    
    I.e each poly is defined by x,y position of vertices. New polygon starts after
    a space.
    
    Returns a list of tuples (polygon, attribute) 
    """
    delimiter = ','
    fid = open(multi_P_file)
    lines = fid.readlines()
    fid.close()
    polygon = []
    polygon_value_list = []
    for line in lines:
        fields = line.split(delimiter)
        try:
            polygon.append([float(fields[0]), float(fields[1])])
        except:
            # Found a line without correct data, assume this signifies the start of a new polygon
            pair = [polygon, attribute] # create polygon , value pair....
            polygon_value_list.append(pair) # add it to the list....
            polygon = []
    # Pickup the last polygon
    pair = [polygon, attribute]
    #print '================================='
    polygon_value_list.append(pair)
    #print len(polygon_value_list)
    #print polygon_value_list
    return polygon_value_list 



# Define a function to read Culvert and Bridge data from Files in Directory
def Create_culvert_bridge_Operator(domain,culvert_bridge_file):
    """This script reads in culvert and bridge data files
    and populates Operator parameters.    
    
    """
    #print culvert_bridge_file
    globals={}
    locals={}    
    
    execfile(culvert_bridge_file, globals, locals)
    #print locals
    #if 'height' and 'z1' and 'z2' in locals:
    if 'z1' and 'z2' in locals:
        culvert = Weir_orifice_trapezoid_operator(domain, **locals)
    elif 'diameter' in locals:
        culvert = Boyd_pipe_operator(domain, **locals)
    elif 'height' in locals:
        culvert = Boyd_box_operator(domain, **locals)
    else:
        raise Exception, 'Cant create culvert'



#-----------------------------------------------------------------------------------------
#          FUNCTION FOR BLOCKAGE
#-----------------------------------------------------------------------------------------
def get_WCC_2016_Blockage_factor(Structure,Event,Scenario, long_result=False, verbose=True):
    """
    If the Structure has a single dimension it is assumed a Diameter (hence Pipe)
    Else it is a Box with w & h
    The Event is grouped 1,2,5= small, 10,20 = med, 50,100,pmp = large
    The Scenario can be a Design or Risk Management Outlook
    Based on these there are two arrays containng the Blockage Factor to be applied
    
    
    --------------------------------------------------------------------
    2017 - 02 - 03
    
    Wollongong Blockage Factor Calculator
    
    Author: Rudy Van Drie
    
    --------------------------------------------------------------------
    Class P1. 
    Pipes 1.2 m internal diameter or smaller.  
    
    Class P2. 
    Pipes  greater  than  1.2  m  internal  diameter.
    
    Class B1
    Box culverts or bridges with a diagonal opening less than 1.5 m,  
    and a width or height less than 0.9 m. 
    
    Class B2. 
    Box  culverts  or  bridges  with a diagonal opening of more than or equal to 1.5 m, less than 3 m 
    and minimum dimension of 0.9 m for both width and height. 
    >= 0.9m w and h
    Class 3. 
    Box culverts or bridges with a diagonal opening of more than or equal 
    to  3  m,  less  than  6  m,  
    and  a  minimum  dimension  of  1.2  m  for  both  width and height.   
    
    Class 4. 
    Box culverts or bridges with a diagonal opening greater than or equal 
    to 6 m, and a minimum dimension of 2.5 m for both width and height.   
    
        CLASSP1   Diam =< 1.2m
        CLASSP2   Diam > 1.2m
    
        CLASSB1:    diag < 1.5m and W or H < 0.9m 
                    
        CLASSB2:    diag >= 1.5m AND diag < 3.0m AND both W and H >= 0.9m
                                    
        CLASSB3:    diag >= 3.0m AND diag < 6.0m AND both W and H >= 1.2m 
    
        CLASSB4:    diag >= 6.0m AND W and H >= 2.5m 
    
    
    DESIGN BLOCKAGE FACTORS
                      CLP1,CLP2
    event,            CLB1,CLB2,CLB3,CLB4
    1,2,5,small,      0.35,0.25,0.15,0.00 
    10,20,medium,     0.50,0.40,0.30,0.05 
    50,100,pmp,large, 0.70,0.50,0.40,0.10 
    
    RISK MANAGEMENT BLOCKAGE FACTORS
                      CLP1,CLP2
    event,            CLB1,CLB2,CLB3,CLB4
    1,2,5,small,      0.60,0.50,0.35,0.05 
    10,20,medium,     0.75,0.65,0.50,0.10 
    50,100,pmp,large, 0.95,0.75,0.60,0.15    
    
    """
    
    # REQUIRED DATA FOR small,medium,large and class 1,2,3,4 for two Scenarios
    BF_DES = [[0.35,0.25,0.15,0.00],[0.50,0.40,0.30,0.05],[0.70,0.50,0.40,0.10]]
    BF_RMN = [[0.60,0.50,0.35,0.05],[0.75,0.65,0.50,0.10],[0.95,0.75,0.60,0.15]]
    
    
    if len(Structure) > 1:# ====== FOR BOX =================
        h = float(Structure[0])
        w = float(Structure[1])
        diag = (w**2+h**2)**0.5
                
        if diag >= 6.00 and w >= 2.5 and h >= 2.5:
            BF_clss = 'CLASS B4'                                    
            cclass = 3
        elif diag >= 3.0 and w >= 1.2 and h >= 1.2:
            BF_clss = 'CLASS B3'            
            cclass = 2
        elif diag >= 1.5 and w >= 0.9 and h >= 0.9:
            BF_clss = 'CLASS B2'            
            cclass = 1
        elif diag < 1.5 or w < 0.9 or h < 0.9:
            BF_clss = 'CLASS B1'
            cclass = 0
    else:   # ====== FOR PIPE ================
        d = float(Structure[0])
        if d < 1.2:
            diag = d
            BF_clss =  'CLASS P1'
            cclass = 0
        else:
            diag = d
            BF_clss =  'CLASS P2'
            cclass = 1
            
    if Event in [1,2,5]: 
        Ev_Row = 0
        Ev_mag = 'Small'
    elif Event in [10,20]: 
        Ev_Row = 1 
        Ev_mag = 'Medium'
    elif Event in [50,100,9999]: 
        Ev_Row = 2
        Ev_mag = 'Large'
    
    if Scenario == 'D':
        Scenario = 'DESIGN'
        BF = BF_DES[Ev_Row][cclass]
    elif Scenario == 'R':
        Scenario = 'RISKMAN'
        BF = BF_RMN[Ev_Row][cclass]

    if verbose:
        print '       Importing Culverts'
        print '   Culvert Size ([H,W] or [d]): ', Structure
        print '                    Event Size: ', Ev_mag
        print '             Blockage Scenario: ', Scenario
        print '               Blockage Factor: ', BF
        print ''

    if long_result:
        return(Scenario, Ev_mag,BF_clss,diag,BF)
    else:
        return(BF)



def get_WCC_2002_Blockage_factor(Structure, verbose=True):
    """
    If the Structure has a single dimension it is assumed a Diameter (hence Pipe)
    Else it is a Box with w & h

    --------------------------------------------------------------------
    2017 - 06 - 22
    
    Wollongong Blockage Factor Calculator for 2002 Blockage Policy
    
    Author: Petar Milevski
    
    --------------------------------------------------------------------
    For all design storm events    
     
    if diag >= 6.0m, blockage factor = 0.25, otherwise blockage factor is 100%

    """
    
    if len(Structure) > 1:# ====== FOR BOX =================
        h = float(Structure[0])
        w = float(Structure[1])
        diag = (w**2+h**2)**0.5
                
    else:   # ====== FOR PIPE ================
        d = float(Structure[0])
        diag = d
     
    if diag >= 6.0:
        BF = 0.25
    else:
        BF = 1.0

    if verbose:
        print '       Importing Culverts'
        print '   Culvert Size ([H,W] or [d]): ', Structure
        print '                      Diagonal: ', diag
        print '               Blockage Factor: ', BF
        print ''

    return(BF)

#!/usr/bin/python
"""
Thing TO DO:
Think how to implement a Parallel Run, may need to create a run file  for each run & call it.... via MPI.... etc


Various Additional Tools To Allow a GENERIC Call to be made to ANUGA



"""
#------------------------------------------------------------------------------------------------------------
def replace_in_file(filename, key, new_value):
    """
    This Script with replace an item in a file based on a Key and a Value and a Separator
    Key = Value where default separator is '='
    Example:
    Blockage = 0.55
    
    """
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    for i, line in enumerate(lines):
        if line.split('=')[0].strip(' \n') == key: #<<< The Key is Blockage....
            lines[i] = key + ' = ' + str(new_value) + '\n'
    f = open(filename, "w")
    f.write("".join(lines))
    f.close()
    return()
#------------------------------------------------------------------------------------------------------------
def set_ANUGA_CULVERT_Blockage(Culvert_Bridge_directory,Event,Blockage_Type):
    """
    CALLS: def replace_in_file(filename, key, new_value):
    
    This Function allows the setting of CULVERT Blockage based on Various Blockage METHODS including
    - All Unblocked
    - All Fully Blocked
    - All Set to the same Global Value
    - WCC2002 Method
    - WCC 2016 Design Method
    - WCC 2016 Risk Method
    - Read File of Pre-Specified Blockage Scenarios (Files should have list of Culvert Filenames and Blockage)
    
    You need to Pass the DIR of where the Culvert Files are located, and Event Flag & the Blockage Type
    
    """
    import os
    import anuga
    print '-------------------------------------------------------------------------------'
    print 'In set ANUGA CULVERT Blockage...'
    print 'Valid Options include: unblocked, fullyblocked, WCC2016_D, WCC2016_R, WCC2002,'
    print '                        a fileRef or a Global Value eg: 0.65'
    print 'The fileRef must have a complete list of Culverts and Blocakge Values you wish'
    print '    to apply.     Eg: Crown_Str_1500pipe, 0.5'
    if Culvert_Bridge_directory is not None:
        for culvert_bridge_file in os.listdir(Culvert_Bridge_directory): # Go through the Culvert files in the Directory and ADD them
            culvert_name = os.path.basename(culvert_bridge_file[:-3])
            try:
                float(Blockage_Type)
                label='"'+culvert_name+'_'+str(Event)+'y'+'_Glob_Blk_'+str(Blockage_Type)+'"'
            except: 
                label='"'+culvert_name+'_'+str(Event)+'y'+'_'+Blockage_Type+'"'
                if Blockage_Type not in ['unblocked','fullyblocked','WCC2016_D','WCC2016_R','WCC2002']:
                    # Must be reading from a file....
                    print 'Reading from File: '+Blockage_Type
                    label='"'+culvert_name+'_'+str(Event)+'y'+'_Blkd_from_File"'
                    binfid = open(Blockage_Type, 'r')
                    lines = binfid.readlines()
                    print lines
                    # Find the Culvert....
                    print culvert_name
                    rightline = [element for element in lines if culvert_name in element] # ... Find the line in file with reference to Culvert
                    print rightline
                    BF = float(rightline[0].split(',')[1])
            Scenario = None
            width = None
            height = None
            Need_BF = True
            file_includes_blockage = False
            file_includes_label = False
            infid = open(anuga.join(Culvert_Bridge_directory,culvert_bridge_file), 'r')
            lines = infid.readlines()
            infid.close()
            try:
                BF = float(Blockage_Type) # <<<<<<<----------------  SET GLOBAL Blockage for ALL Files
            except: 
                if 'WCC2016_R' in Blockage_Type : Scenario = 'R'
                if 'WCC2016_D' in Blockage_Type : Scenario = 'D'
            for line in lines:
                if line.startswith('diameter'):
                    diameter = float(line.split('=')[1])
                    Structure = str(diameter)          
                    insert_index = lines.index(line)
                    if Scenario is not None:
                        BF=anuga.get_WCC_2016_Blockage_factor([diameter],Event,Scenario) 
                if line.startswith('height'):                    
                    height = float(line.split('=')[1].split('#')[0].strip())
                    hgt_index = lines.index(line)
                if line.startswith('width'):
                    width = float(line.split('=')[1].split('#')[0].strip())
                    wdth_index = lines.index(line)
                if width is not None and height is not None and Need_BF:
                    if Scenario is not None:
                        BF=anuga.get_WCC_2016_Blockage_factor([height,width],Event,Scenario)  
                    insert_index = max(hgt_index,wdth_index)
                    Structure = [height,width]
                    Need_BF = False
            if Blockage_Type == 'unblocked': BF = 0.0
            if Blockage_Type == 'fullyblocked': BF = 1.0
            if Blockage_Type == 'WCC2002': 
                print Structure
                print Structure[1]
                BF = anuga.get_WCC_2002_Blockage_factor(Structure)

            #------- GOT NEW BF --------------------------------
            # Blockage can be placed any where in the file....
            for line in lines: # ----- NOW REPLACE LINES if they EXIST
                if line.startswith('blockage'):
                    replace_in_file(anuga.join(Culvert_Bridge_directory,culvert_bridge_file),'blockage',BF) # <<<<----- Replace the Blockage Factor
                    file_includes_blockage = True
                if line.startswith('label'):
                    replace_in_file(anuga.join(Culvert_Bridge_directory,culvert_bridge_file),'label',label) # <<<<----- Replace the Blockage Factor
                    print label
                    file_includes_label = True
            #---- Sometimes File MAy not have a label or blockage flag.... Need to ADD them...
            if not file_includes_blockage:
                # Add to file
                #print 'GOT HERE.......'
                #raw_input()
                lines.insert(insert_index,'blockage = '+str(BF))
                # rewrite the File
            if not file_includes_label:
                # Add to file
                lines.insert(0,'label = '+label)
                # Rewrite the file
            #--- End If
        #---- END for ------------
        print 'Blockage set to :'+str(BF)
    # ---- END if--------------
    #raw_input('check....')
    return()
#------------------------------------------------------------------------------------------------------------
def Get_Model_Rect_Ext(Model_Ext_Poly):
    """
    Call Using:
        W,E,S,N = Get_Model_Rect_Ext(Model_Ext_Poly)
        To return the rectangular extents of any polygon to use as a model Rect Extent
    """
    print 'Get Model Extents...'
    # GET MODEL EXTENTS
    #    Read a Polygon File to extract the Maximum X & Y to set the bounding polygon , E,W,N,S
    bounding_polygon_file = Model_Ext_Poly
    infid = open(bounding_polygon_file)
    i=0
    W=999999999.9
    E=-99999999.9
    S=999999999.9
    N=-99999999.9
    
    for line in infid.readlines():
        fields = line.split(',')
        i=i+1
        print i,fields[0],fields[1]

        W=min(W,float(fields[0]))
        E=max(E,float(fields[0]))
        S=min(S,float(fields[1]))
        N=max(N,float(fields[1]))  
    return(W,E,S,N)
#---------------------------------------------------------------------
    
def mesh_refine_defined(mesh_refine_directory):
    """
    Creates List of Polgons and Values of Mesh Refinement Target Size
    """
    import anuga
    # ADD INTERIOR REFINE REGIONS
    print 'REFINE DIRECTORY TO BE USED:'
    print mesh_refine_directory
    mesh_refine_defined = anuga.get_polygon_value_list(mesh_refine_directory)
    print mesh_refine_defined
    print 'MESH REFINE Target Sizes...'
    # Mesh_Target_Size = zip(*mesh_refine_defined)[1]  # puts it into a tuple not list...
    Mesh_Target_Size =[]
    for row in mesh_refine_defined:
        Mesh_Target_Size.append(row[1])
    print Mesh_Target_Size
    Mesh_Target_Size = Mesh_Target_Size.sort()
    print Mesh_Target_Size
    
    mesh_refine_defined.reverse()   
    #mesh_refine_reversed = mesh_refine_defined.reverse()
    
    return(mesh_refine_defined)
#---------------------------------------------------------------------

def build_Depth_Vary_N_Func(filename):
    """
    FACTORY RD Material File for ANUGA Depth Varying Operator

    #1
    Roadways Impervious
    h = [-0.1, 0.010, 0.1,   2.0,   99.0]
    n = [0.05, 0.015, 0.012, 0.011, 0.01]

    """
    from scipy import interpolate
    frict_functions = {}
    f_func = []
    LandCatNums = []
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    print 'Depth Varying Roughness applied from: '+lines[1]
    ProcessCat = False
    for line in lines:
        if line.startswith('#'):
            LandCat = int(line.split('#')[1].strip())
            print LandCat
            LandCatNums.append(LandCat)
            ProcessCat = True
        elif ProcessCat:
            #print line
            if line.startswith('h'):
                h = line.split('=')[1].replace(" ","").translate(None, "[]").split(',')
                h = [ float(x) for x in h ]
                print h
            if line.startswith('n'):
                n = line.split('=')[1].replace(" ","").translate(None, "[]").split(',')
                n = [ float(x) for x in n ]
                print n
                ProcessCat = False
                f_func.append(interpolate.interp1d(h, n))
        #---- End If
    for item in LandCatNums:
        frict_functions.update({item:f_func[item-1]}) # Create Dictionary of Friction FUnctions 
        # THis Should Create :    frict_functions = {3:f_fun3, 4:f_fun4, 8:f_fun8, 13:f_fun13}
    return(frict_functions)
# --------------------------------------------------------------------------------------------------
def plot_fric_funcs(filename):
    """
    File Descriptor on First Line in File
    
    #1    Land Use Category Number
    Roadways Impervious    Label or Descriptor of Land Use type
    h = [-0.1, 0.010, 0.1,   2.0,   99.0]   Depths
    n = [0.05, 0.015, 0.012, 0.011, 0.01]   Roughness values
    {Blank Line}
    #2  Next land Use etc...
    """
    from matplotlib import pyplot as pyplot
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    print 'Depth Varying Roughness applied from: '+lines[1]
    ProcessCat = False
    Get_title = False
    LandCat = []
    Roughness_Titles = []
    SUPTitle = lines[1].strip()
    for line in lines:
        if line.startswith('#'):
            LandCat =line.split('#')[1].strip()
            Get_title = True
        elif Get_title:
            Roughness_Title = line.strip()
            Get_title = False
            ProcessCat = True
        elif ProcessCat:
            #print line
            if line.startswith('h'):
                h = line.split('=')[1].replace(" ","").translate(None, "[]").split(',')
                h = [ float(x) for x in h ]
                print h
            if line.startswith('n'):
                n = line.split('=')[1].replace(" ","").translate(None, "[]").split(',')
                n = [ float(x) for x in n ]
                print n
                ProcessCat = False
                #pyplot.ion()
                #pyplot.clf() # Clear any Previous Plots...
                pyplot.suptitle(SUPTitle)
                pyplot.title('Adopted Depth Varying Roughness Functions')
                pyplot.ylabel('Depth of Flow (m)')
                pyplot.xlabel("Manning's Roughness")
                pyplot.plot(n,h, linewidth = 1.5,label = LandCat+'-'+Roughness_Title)  # Plot the Time Series
                pyplot.legend(loc='upper right')
                pyplot.yscale('log')
                pyplot.xscale('log')
                #pyplot.savefig(outfile_name+'_Qplot'+Last_Plot_format_used )
    pyplot.show()
    #raw_input('Plot done  Check location !!') 
    return() 
# --------------------------------------------------------------------------------------------------    
def get_LAND_USE_from_Mid_Mif(Rfile):
    """
    This function reads a MID/MIF and extracts the Closed polygons and associated Land Use Catergory
    USED to ASSIGN Depth Varying ROughness
    
    The MID File contains the Value of the ATTRIBUTE, in Columns
    The MIF File contains the Description of the Columns in the MID file & Contains the Polygons
    
    Users need to add an ATTRIBUTE called LANDUSECAT as a integer for this to work.....!!!
    
    Columns 5
      NAME Char(20)
      LAYER Char(17)
      PERIMETER Char(9)
      ENCLOSED_AREA Char(9)
      LANDUSECAT Integer
    """
    print 'In..... get_LAND_USE_from_Mid_Mif'
    #Open the MIF and Find the LANDUSECAT Attribute
    try:
        fmif = open(Rfile+'.mif')
    except:
        fmif = open(Rfile+'.MIF')
        
    linesmif = fmif.readlines()
    fmif.close()
    try:
        fmid = open(Rfile+'.mid')
    except:
        fmid = open(Rfile+'.MID')
    linesmid = fmid.readlines()
    fmid.close()
    MID_Attribute_COL = int(linesmif[4].split(' ')[1])
    Category_Poly_List = []  #Group Polygons with Attribute
    Category_List = []
    PolyCount = 0
    # GET THE LAND USE CATEGORIES AND ASSOCIATE WITH A SET OF POLYGONS
    for line in linesmid:
        if int(line.split(',')[MID_Attribute_COL-1]) in Category_List:
            App_indx = [i for i, lst in enumerate(Category_Poly_List) if int(line.split(',')[MID_Attribute_COL-1]) in lst] # get Index of Category in the Poly List
            print Category_Poly_List[App_indx[0]]
            Category_Poly_List[App_indx[0]].append(PolyCount) # Now Add the additional Poly to Existing Category
        else: 
            Category_List.append(int(line.split(',')[MID_Attribute_COL-1]))
            Category_Poly_List.append([int(line.split(',')[MID_Attribute_COL-1]),PolyCount])
        PolyCount +=1
    print Category_List
    print Category_Poly_List # Got the Multi Poly file now process to put Polys in a List
    polylist = []
    polygon = []
    check_pts_list=[]
    Poly_count=0
    Poly_line_count=0
    Points_in_Poly = 100 # Why is this set to 100 ???
    total_lines_in_file= len(linesmif)
    #print "total number of lines in the Polygons FILE is: ",total_lines_in_file

# ==================================================== FOR LOOP ===========================        
    for i, line in enumerate(linesmif): # Read MIF File and PUT Polys in a List
        if line.strip().startswith('Region'): # Must be Closed POLYS
            Poly_line_count=0
            check_pts_list=[]
            if Poly_count==0:
                pass
            else:
                polylist.append(polygon)
                polygon=[]
            Poly_count+=1
        elif line.strip().startswith('    Pen'):
            pass
        elif line.strip().startswith('    Brush'):
            pass
        else:
            Poly_line_count+=1
            if Poly_line_count > 1 and Poly_line_count <= (Points_in_Poly+1) and Poly_count<>0:
                fields = line.split(' ')
                if line in check_pts_list and Poly_line_count <> Points_in_Poly+1:   # Get rid of any doubled up points NOTE this gets rid of last line !!!
                    pass
                else:
                    polygon.append([float(fields[0]),float(fields[1])])
                    check_pts_list.append(line)
            elif Poly_line_count==1 and Poly_count<>0:
                Points_in_Poly=int(line)
    polylist.append(polygon)
    print len(polylist)
    return (Category_Poly_List,polylist)
#-------------------------------------------------------------------------------
def cleanup_files_after_RUN(outname,Results_DIR,Debug):
    """
    End of Run Doing any Processing between runs here....
    Can Move Files, etc...
    Need anuga outname, and a Results DIR Name
    """
    import os
    import glob
    import shutil
    #-----------------------------------------------------------------------------------------------------------------------------------

    print 'About to copy files to RESULTS DIRS... from SWW outfile'
    Scenario_DIR = outname
    if Results_DIR == None: Results_DIR = "10_RESULTS"
    print os.path.join(Results_DIR,Scenario_DIR )
    if not os.path.exists(os.path.join(Results_DIR)):
        os.makedirs(os.path.join(Results_DIR))
        if Debug: print 'RESULTS DIR Created...'
    if not os.path.exists(os.path.join(Results_DIR,Scenario_DIR )):
        os.makedirs(os.path.join(Results_DIR,Scenario_DIR ))
        if Debug: print 'SCENARIO & STORM RESULTS DIR Created...'
    if os.path.exists(os.path.join(Results_DIR,Scenario_DIR )):
        # Copy Files as required
        print 'DIRS Available... now copy....'
        for file in glob.glob('*.log'):
            print(file)
            if file == 'anuga.log':
                pass
            else:
                shutil.copy(file, os.path.join(Results_DIR,Scenario_DIR))
                os.remove(file)
        for file in glob.glob('*.SWW'):
            print(file)
            shutil.copy(file, os.path.join(Results_DIR,Scenario_DIR))
            os.remove(file)
        if Debug: print 'ALL FILES COPIES TO RESULTS DIR '
        #os.remove() #will remove a file.
        #os.rmdir()  #will remove an empty directory.
        #shutil.rmtree() #will delete a directory and all its contents.
            # Copy Files here
    return()
#--------------------------------------------------------------------------------------------------------

#++++++++++++++++++++++++++++++++++++++++++++ ANUGA RUN SCRIPT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def RUN_ANUGA_GEN(Model_Ext_Poly,
                  basename,
                  mesh_refine_directory,
                  Roughness_Polygons_directory,
                  interior_holes,
                  Depth_VarN_File,
                  Depth_VarN_Map,
                  Culvert_Bridge_directory,
                  SingleQ4Event,
                  Rainfall_Amount_DIR,
                  Rainfall_Apply_DIR,
                  outname,
                  Model_set_flow_algorithm,
                  Evolve_Model,
                  Debug,
                  PlotStuff):
    """
    
    
    """
    import shutil
    import anuga, numpy, time, os, glob
    from anuga.operators.rate_operators import Polygonal_rate_operator
    from anuga import Rate_operator
    from anuga import file_function, Polygon_function, read_polygon, create_mesh_from_regions, Domain, Inlet_operator
    import anuga.utilities.spatialInputUtil as su
    
    """
    import  anuga.parallel
    from anuga import distribute, myid, numprocs, finalize, barrier
    from anuga.parallel.parallel_operator_factory import Inlet_operator, Boyd_box_operator, Boyd_pipe_operator
    """
    
    
    
    from scipy import interpolate
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    #  IF USING INTERNAL POLIES FOR BUILDINGS USE THIS,
    # VERY IMPORTANT NOTE: just remember if using holes in domain and rainfall, you need to factor up the rainfall by the area of holes as you will loose rainfall from the domain from the building holes
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    
    #If you have individual polygon files, use this method
    #BuildingDictionary = 'Model/Buildings'
    
    #if you have a single file with all polygons with all buildings or holes in it separated by a single space use this
    #BuildingDictionary = 'Model/Soundbarrier/Soundbarrier.csv'
    
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # FILENAMES, MODEL DOMAIN and VARIABLES
    # See ARR_Design_Rain for event and losses options
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    """
    import sys
    Scenario = sys.argv[1]
    Blockage = sys.argv[2]
    """
    print 'In RUN_ANUGA_GEN....'
    
    
    print 'About to Create...'+outname
    meshname = anuga.join(os.path.dirname(basename),'Domain.msh')
    
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # ENTER DOMAIN COORDINATES
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    
    W,E,S,N = Get_Model_Rect_Ext(Model_Ext_Poly)
    maximum_triangle_area = 20.0			# mesh size outside the catchment boundary
    minimum_storable_height = 0.02			# minimum depth which is stored
    tide = 0.0				                # outlet boundary water level
    base_friction = 0.055 		            # this sets the roughness of the roads
    alpha = 0.99				            # smoothing of mesh (0 is not smooth 1 is smoothest)
    yieldstep = 60.0                        # timestep at which model run is saved
    finaltime = 3600.0                      # time at end of model run
    verbose = True
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # SETUP DOMAIN ONLY ON PROCESSOR 0
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    
    #interior_holes = anuga.get_polygon_list_from_files(dir)
    #interior_holes = anuga.read_multi_poly_file(holes)
    print 'Holes added to Mesh....'
    
    if anuga.parallel.myid == 0:
        #------------------------------------------------------------------------------
        # CREATING MESH
        #------------------------------------------------------------------------------
        
        bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]
        #bounding_polygon = anuga.read_polygon('Model/Bdy/Catchment.csv')
        print 'Getting Mesh Refinement....'
        interior_regions = mesh_refine_defined(mesh_refine_directory)
    
        
        #If you have individual polygon files in a directory, use this method
        #interior_holes = anuga.read_hole_dir_multi_files_with_single_poly(BuildingDictionary, filepattern='*.csv')
        
        #if you have a single file with multiple polygons
        #interior_holes = anuga.read_multi_poly_file(BuildingDictionary)
    
        
        create_mesh_from_regions(bounding_polygon,
            boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
            maximum_triangle_area=maximum_triangle_area,
            interior_regions=interior_regions,
            interior_holes=interior_holes,
            #breaklines=riverWalls.values(),
            filename=meshname,
            use_cache=False, 
            verbose=True)
     
        domain = Domain(meshname, use_cache=False, verbose=True) # Create the Domain
        domain.set_name(outname)
        #domain.set_datadir(outname)    
        print domain.statistics()
        
        #------------------------------------------------------------------------------
        # APPLY MANNING'S ROUGHNESSES
        #------------------------------------------------------------------------------  
        print 'Applying mannings roughness'
        
        friction_list = None
        if Roughness_Polygons_directory is not None:
            friction_list = anuga.get_polygon_value_list(Roughness_Polygons_directory)
        if not friction_list:
            domain.set_quantity('friction', base_friction)# SETS a Default BASE Friction.....
        else:
            domain.set_quantity('friction', anuga.Polygon_function(friction_list, default=base_friction, geo_reference=domain.geo_reference))
            # Sets Fixed Roughness over a series of Polygons
        
        #if Depth_Var_File is not None:
        # Depth Varying Set by Polygon
        # Single Polygon
        """
        depthVArNPoly = bounding_polygon
        h = [-0.001, 0.0500, 0.20, 1.000, 3.00, 99.00]
        n = [0.250, 0.160, 0.090,  0.060, 0.050, 0.035]
        f_fun1 = interpolate.interp1d(h, n)
        op = anuga.operators.set_friction_operators.Polygonal_depth_friction_operator(domain,friction=f_fun1,polygon=depthVArNPoly,verbose=True)
        """
        #----------------------------------------------------------------------
        #
        # Depth Varying Set by a series of Polygons for a series of Land Uses:
        #
        #----------------------------------------------------------------------
        """
        This approach uses the Function:   build_Depth_Vary_N_Func("Depth_Vary_Materials_List.txt")
        to create Multiple Friction Functions
        which are then applied to a series of Polygons using the Function:  get_LAND_USE_from_Mid_Mif('simpe_EG')
        The Polygons are in a MID/MIF and Must contain the ATTRIBUTE 'LANDUSECAT' as an integer 
        this links the Polygons to the related friction function
        """
        if Depth_VarN_File == None:
            pass
        else:
            print 'Get Friction Function Data for Roughness...' 
            print Depth_VarN_File
            frict_functions = build_Depth_Vary_N_Func(Depth_VarN_File) 
            if PlotStuff :plot_fric_funcs(Depth_VarN_File) 
            print 'Keys:'
            print frict_functions.keys()
            # Builds a Dictionary... of Fiction Functions, frict_functions = {3:f_fun3, 4:f_fun4, 8:f_fun8, 13:f_fun13}
            Category_Poly_List,polylist = get_LAND_USE_from_Mid_Mif(Depth_VarN_Map )
            print len(polylist),Category_Poly_List
            print 'Create Friction Operators...'
            
            for i in Category_Poly_List:
                All_indices = []
                print 'Category: '+str(i[0]) # This is the category
                print 'Friction Func: '
                print '---'
                for j in range(len(i)):
                    if j == 0:pass # No PolyGons ??
                    elif j == 1: # Only 1 Polygon
                        print 'Poly: '+str(i[j]) # These are the Polygons associated with that Category
                        print polylist[i[j]]
                        
                        region = anuga.Region(domain, polygon=polylist[i[j]])
                        #print region.indices
                        All_indices.append(region.indices)          
                    elif j > 1:
                        print 'Poly: '+str(i[j]) # These are the Polygons associated with that Category
                        print polylist[i[j]]
                        region = anuga.Region(domain, polygon=polylist[i[j]])
                        #print region.indices
                        All_indices.append(region.indices)
                # CReate Depth Varying Roughness Operator using the collected Indices
                op = anuga.Depth_friction_operator(domain,friction=frict_functions[i[0]], indices=All_indices[0]) 
                # Using Polygons directly only allows a SINGLE poly to be passed.......
                #op = anuga.operators.set_friction_operators.Polygonal_depth_friction_operator(domain,frict_functions[i[0]],polygon=polylist[i[j]],verbose=True)  
                print '+++++'
            #------------------------------------------------------------------------------
            #          Depth Varying ROughness Operators Completed
            #------------------------------------------------------------------------------        
        
        
        
        #------------------------------------------------------------------------------   
        # Set a Initial Water Level over the Domain
        #------------------------------------------------------------------------------
        
        domain.set_quantity('stage', tide)
        domain.set_quantity('elevation', filename=basename+'.csv', use_cache=True, verbose=True, alpha=alpha) #change csv to pts if using ascii
        
    else:
    	domain = None
    
    if anuga.parallel.myid == 0 and verbose: print 'DISTRIBUTING DOMAIN'
    domain = anuga.parallel.distribute(domain)
    domain.set_minimum_storable_height(minimum_storable_height)
    
    #domain.riverwallData.create_riverwalls(riverWalls) 
    if anuga.parallel.myid == 0 and verbose: print 'CREATING INLETS'  
    
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # ENTER CULVERT DATA, if you have more than one culvert, copy and paste block below and edit data
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # OR IF YOU HAVE MANY CULVERTS YOU CAN READ CULVERT INFO FROM A DIRECTORY, CULVERT INFO MUST BE IN A .py FILE
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    print 'SETUP Culvert Scenario...'
    Event = 20
    Scenario = 'D'
   
    if Culvert_Bridge_directory is not None:
        for culvert_bridge_file in os.listdir(Culvert_Bridge_directory): # Go through the Culvert files in the Directory and ADD them
            print Culvert_Bridge_directory+'/'+culvert_bridge_file
            # Need method to extract data and set Blockage............ and set the Log File name.... ???
            anuga.Create_culvert_bridge_Operator(domain,anuga.join(Culvert_Bridge_directory,culvert_bridge_file)) # Adds all the files as Culverts
    
    
    """
    # STEVE's Method.... 
    if Culvert_Bridge_directory is not None:
        print Culvert_Bridge_directory
        culverts = {}
        for culvert_bridge_file in os.listdir(Culvert_Bridge_directory):
            print Culvert_Bridge_directory+'/'+culvert_bridge_file
            print '---------'
            
            culverts[culvert_bridge_file] = anuga.Create_culvert_bridge_Operator(domain,anuga.join(Culvert_Bridge_directory,culvert_bridge_file)) # Adds all the files as Culverts
            print '========='
            print culverts[culvert_bridge_file]
            culverts[culvert_bridge_file].set_culvert_blockage(get_WCC_2016_Blockage_factor([height,width],Event,Scenario)) # <<<--------- SET BLOCKAGE Factor HERE !!!!
    print culverts
    raw_input('Check Create OP..')
    """
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # APPLY RAINFALL
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    if Rainfall_Amount_DIR is not None:
        print 'Set Up Rainfall.......'
        for filename in os.listdir(Rainfall_Amount_DIR):
            if filename[-3:]== 'tms':
                Rainfile = anuga.join(Rainfall_Amount_DIR,filename)
                RainApply_poly = anuga.join(Rainfall_Apply_DIR,filename[:-3]+'csv')
                print RainApply_poly 
                print Rainfile
                polygon = anuga.read_polygon(RainApply_poly)
                rainfall = anuga.file_function(Rainfile, quantities='rate')
                op1 = Polygonal_rate_operator(domain, rate=rainfall, factor=1.0e-3, polygon=polygon, default_rate = 0.0)
    else: print 'NOT APPLYING RAINFALL FOR THIS RUN !!!.....'
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # APPLY FLOW ACROSS A LINE
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    if SingleQ4Event is not None:
        for Data in range(len(SingleQ4Event)):
            print SingleQ4Event[Data][0],SingleQ4Event[Data][1]
            Qline = SingleQ4Event[Data][1]
            Q = SingleQ4Event[Data][0]
            anuga.Inlet_operator(domain, Qline, Q)
    
    # Flows Could be Hydrographs ???           
                   
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # SETUP BOUNDARY CONDITIONS
    #----------------------------------------------------------------------------------------------------------------------------------------------------
        
    print 'Available boundary tags', domain.get_boundary_tags()
        
    Br = anuga.Reflective_boundary(domain)
    Bd = anuga.Dirichlet_boundary([tide,0,0])
    
    domain.set_boundary({'interior': Br, 'exterior': Bd, 'west': Br, 'south': Bd, 'north': Br, 'east': Bd})
    domain.set_quantities_to_be_stored({'elevation':1,'stage':2,'xmomentum':2,'ymomentum':2,'friction':2})  
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    # EVOLVE SYSTEM THROUGH TIME
    #----------------------------------------------------------------------------------------------------------------------------------------------------
    domain.set_flow_algorithm(Model_set_flow_algorithm)
    domain.print_algorithm_parameters()
    #domain.set_flow_algorithm('2_0')
    if anuga.parallel.myid == 0 and verbose: print 'EVOLVE'
        
    t2 = time.time()
    if Evolve_Model:    
        for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime): #change run time and timstep saved here
            domain.print_operator_timestepping_statistics() 
            if anuga.parallel.myid == 0:
                domain.write_time()
    
        if anuga.parallel.myid == 0:
            print 'Number of processors %g ' %anuga.parallel.numprocs
            print 'That took %.2f seconds' %(time.time()-t0)
            print 'Communication time %.2f seconds'%domain.communication_time
            print 'Reduction Communication time %.2f seconds'%domain.communication_reduce_time
            print 'Broadcast time %.2f seconds'%domain.communication_broadcast_time
        
        domain.sww_merge(delete_old=True)
    
        anuga.parallel.finalize()
    print '+++++++++ ANUGA SCENARIO RUN COMPLETE +++++++++++++++'
#------===========  END OF ANUGA RUN ======================----------------------------    
