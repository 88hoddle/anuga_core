import os.path
import sys

from anuga.utilities.system_tools import get_pathname_from_package
from anuga.geometry.polygon_function import Polygon_function
        
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.abstract_2d_finite_volumes.quantity import Quantity

import anuga

from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
from anuga.structures.inlet_operator import Inlet_operator
                            
#from anuga.culvert_flows.culvert_routines import boyd_generalised_culvert_model
     
from math import pi, pow, sqrt

import numpy as num





"""test_culvert_link

This test exercises the culvert and checks values outside rating curve
are dealt with       
"""

path = get_pathname_from_package('anuga.culvert_flows')    

length = 60.
width = 15.

dx = dy = 0.5          # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx),
                                               int(width/dy),
                                               len1=length, 
                                               len2=width)
domain = anuga.Domain(points, vertices, boundary)   
domain.set_starttime(10)
domain.set_name()                 # Output name
domain.set_default_order(2)



#----------------------------------------------------------------------
# Setup initial conditions
#----------------------------------------------------------------------

def topography(x, y):
    """Set up a weir
    
    A culvert will connect either side
    """
    # General Slope of Topography
    z = -x/1000
    #z = 0
    
    N = len(x)
    for i in range(N):

       # Sloping Embankment Across Channel
        if 5.0 < x[i] < 10.1:
            z[i] +=  0.5*(x[i] - 5.0)     # Sloping Segment  U/S Face
        if 10.0 < x[i] < 12.1:
           z[i] += 5                    # Flat Crest of Embankment
        if 12.0 < x[i] < 14.5:
            z[i] += 2.5 - 1.0 * (x[i] - 12.0) # Sloping D/S Face
                   
    return z


domain.set_quantity('elevation', topography) 
domain.set_quantity('friction', 0.01)         # Constant friction 
domain.set_quantity('stage',
                    expression='elevation')   # Dry initial condition

filename = os.path.join(path, 'example_rating_curve.csv')

end_point0 = num.array([5.0, 7.5])
end_point1 = num.array([43.0, 7.5])

Boyd_pipe_operator(domain,
                   losses=1.5,
                   end_points=[end_point0, end_point1],
                   diameter=2.5,
                   apron=0.5,
                   use_momentum_jet=True, 
                   use_velocity_head=False,
                   manning=0.013,
                   verbose=False)

#-----------------------------------------------------------------------
# Setup boundary conditions
#-----------------------------------------------------------------------

# Inflow based on Flow Depth and Approaching Momentum
Bi = anuga.Dirichlet_boundary([2.0, 0.0, 0.0])
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall

domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})


##-----------------------------------------------------------------------
## Evolve system through time
##-----------------------------------------------------------------------

for t in domain.evolve(yieldstep=1.0, finaltime=100.0):
    domain.write_time()
