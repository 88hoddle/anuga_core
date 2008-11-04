#!/usr/bin/env python

import sys

from types import ListType, TupleType

import exceptions

class NoTrianglesError(exceptions.Exception): pass
import anuga.mesh_engine.mesh_engine_c_layer as triang
#import anuga.mesh_engine.list_dic as triang

from Numeric import array, Float, Int32, reshape

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.anuga_exceptions import ANUGAError
    
def generate_mesh(points=None,
                  segments=None,holes=None,regions=None,
                  pointatts=None,segatts=None,
                  mode=None, dummy_test=None):
    """
    pointatts can be a list of lists.

    generatedtriangleattributelist is used to represent tagged regions.
    #FIXME (DSG-DSG): add comments
    """
    #FIXME (DSG-DSG): Catch parameters that are lists,
    #instead of lists of lists
    # check shape[1] is 2 etc

    if points is None:
        points = []

    if segments is None:
        segments = []

    if holes is None:
        holes = []
        
    if regions is None:
        regions = []

    if dummy_test is None:
        dummy_test  = []
        
    try:
        points =  ensure_numeric(points, Float)
    except ValueError:
        msg = 'ERROR: Inconsistent points array.'
        raise ANUGAError, msg
    if points.shape[1] <>2:
        msg = 'ERROR: Bad shape points array.'
        raise ANUGAError, msg

    #print "pointatts",pointatts 
    # This is after points is numeric
    if pointatts is None or pointatts == []:
        pointatts = [[] for x in range(points.shape[0])]
        
    try:
        # If Int is used, instead of Int32, it fails in Linux
        segments = ensure_numeric(segments, Int32)
        
    except ValueError:
        msg = 'ERROR: Inconsistent segments array.'
        raise ANUGAError, msg
    
    # This is after segments is numeric
    if segatts is None or segatts == []:
        segatts = [0 for x in range(segments.shape[0])]
        
    try:
        holes = ensure_numeric(holes, Float)
    except ValueError:
        msg = 'ERROR: Inconsistent holess array.'
        raise ANUGAError, msg

   
    regions = add_area_tag(regions)
    try:
        regions = ensure_numeric(regions, Float)
    except  (ValueError, TypeError):
        msg = 'ERROR: Inconsistent regions array.'
        raise ANUGAError, msg
        
    if not regions.shape[0] == 0 and regions.shape[1] <= 2:
        msg = 'ERROR: Bad shape points array.'
        raise ANUGAError, msg
    
    try:
        pointatts = ensure_numeric(pointatts, Float)
    except (ValueError, TypeError):
        msg = 'ERROR: Inconsistent point attributes array.'
        raise ANUGAError, msg

    if pointatts.shape[0] <> points.shape[0]:
        msg = """ERROR: Point attributes array not the same shape as
        point array."""
        raise ANUGAError, msg
    if len(pointatts.shape) == 1:
        pointatts = reshape(pointatts,(pointatts.shape[0],1))
    
    try:
        segatts = ensure_numeric(segatts, Int32)
    except ValueError:
        msg = 'ERROR: Inconsistent point attributes array.'
        raise ANUGAError, msg
    if segatts.shape[0] <> segments.shape[0]:
        msg = """ERROR: Segment attributes array not the same shape as
        segment array."""
        raise ANUGAError, msg
    
    #print "mode", mode
    if mode.find('n'):
        #pass
        mode = 'j' + mode
        # j- Jettisons vertices that are not part of the final
        #    triangulation from the output .node file (including duplicate
        #    input vertices and vertices ``eaten'' by holes).  - output a
        #    list of neighboring triangles
        # EG handles lone verts!
            
    #print "points",points 
    #print "segments", segments
    #print "segments.shape", segments.shape
    #print "holes", holes
    #print "regions", regions
    #print "pointatts", pointatts
    #print "segatts", segatts
    #print "mode", mode
    #print "yeah"
    # .copy()
    trianglelist, pointlist, pointmarkerlist, pointattributelist, triangleattributelist, segmentlist, segmentmarkerlist, neighborlist = triang.genMesh(points,segments,holes,regions,
                          pointatts,segatts, mode)
    mesh_dict = {}
    # the values as arrays
    mesh_dict['generatedtrianglelist'] = trianglelist
    mesh_dict['generatedpointlist'] = pointlist
    #print "mesh engine mesh_dict['generatedpointlist']", mesh_dict['generatedpointlist']
    # WARNING - generatedpointmarkerlist IS UNTESTED
    mesh_dict['generatedpointmarkerlist'] = pointmarkerlist
    mesh_dict['generatedpointattributelist'] = pointattributelist
    mesh_dict['generatedsegmentlist'] = segmentlist 
    mesh_dict['generatedsegmentmarkerlist'] =  segmentmarkerlist
    mesh_dict['generatedtriangleneighborlist'] = neighborlist
    mesh_dict['qaz'] = 1 #debugging

    #mesh_dict['triangleattributelist'] = triangleattributelist
    #print "mesh eng generatedtrianglelist", trianglelist
    #print "mesh eng mesh_dict['triangleattributelist'] ",mesh_dict['triangleattributelist']
    #print "mesh eng mesh_dict['generatedtriangleattributelist'] ", mesh_dict['generatedtriangleattributelist']   
    if True: 
        mesh_dict['generatedtriangleattributelist'] = triangleattributelist

        if mesh_dict['generatedtriangleattributelist'].shape[1] == 0:
            mesh_dict['generatedtriangleattributelist'] = None
            
        if mesh_dict['generatedpointattributelist'].shape[1] == 0:
            mesh_dict['generatedpointattributelist'] = None
            
        if mesh_dict['generatedtriangleneighborlist'].shape[1] == 0:
            mesh_dict['generatedtriangleneighborlist'] = None
            
        if trianglelist.shape[0] == 0:
            # There are no triangles.
            # this is used by urs_ungridded2sww
            raise NoTrianglesError
    #print "mesh eng mesh_dict['generatedtriangleattributelist'] ", mesh_dict['generatedtriangleattributelist']  
    a = mesh_dict['generatedtriangleattributelist']
    #print 'mesh_dict', mesh_dict
    # the structure of generatedtriangleattributelist is an list of
    # list of integers.  It is transformed into a list of list of
    # strings later on.  This is then inputted into an triangle
    # object.  The triangle object outputs a list of strings.  Note
    # the subtle change!  How should I handle this?  For my code, when
    # the list of list of integers is transformed, transform it into a
    # list of strings, not a list of list of strings.
    
    return mesh_dict

def add_area_tag(regions):
    """
    So, what is the format?
    A list with
    [x,y,region_tag,area] OR [x,y,region_tag]
    if it's [x,y,region_tag], add a 4th element, value of 0.0.
    """
    if isinstance(regions, ListType):
        for i, region in enumerate(regions):
            if len(region) == 3:
                if isinstance(region, TupleType):
                    #FIXME: How do you convert a tuple to a list?
                    # I can do it a stupid way..
                    tuple = region[:]
                    regions[i] = []
                    for j in tuple:
                        regions[i].append(j)
                    regions[i].append(0.0)
                else:
                    regions[i].append(0.0)
                    
                # let ensure numeric catch this..    
                #len(region) <= 2:
                #msg = 'ERROR: Inconsistent regions array.'
                #raise msg
            #elif
    return regions

if __name__ == "__main__":
    pass 