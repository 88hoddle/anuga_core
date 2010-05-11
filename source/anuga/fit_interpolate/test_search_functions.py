#!/usr/bin/env python


import unittest
from mesh_quadtree import MeshQuadtree, compute_interpolation_values

from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.geometry.quad import Cell
from anuga.geometry.aabb import AABB
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.geometry.polygon import is_inside_polygon, is_inside_triangle    

import numpy as num

class Test_search_functions(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def NOtest_that_C_extension_compiles(self):
        FN = 'search_functions_ext.c'
        try:
            import search_functions_ext
        except:
            from compile import compile

            try:
                compile(FN)
            except:
                raise 'Could not compile %s' %FN
            else:
                import search_functions_ext



    def test_off_and_boundary(self):
        """test_off: Test a point off the mesh
        """

        points, vertices, boundary = rectangular(1, 1, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        root = MeshQuadtree(mesh)
        root.set_last_triangle()

        found, s0, s1, s2, k = root.search_fast([-0.2, 10.7])
        assert found is False

        found, s0, s1, s2, k = root.search_fast([0, 0])
        assert found is True
        
        
    def test_small(self):
        """test_small: Two triangles
        """

        points, vertices, boundary = rectangular(1, 1, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        root = MeshQuadtree(mesh)
        root.set_last_triangle()

        x = [0.2, 0.7]
        found, s0, s1, s2, k = root.search_fast(x)
        assert k == 1 # Triangle one
        assert found is True        
        
    def test_bigger(self):
        """test_bigger
        
        test larger mesh
        """

        points, vertices, boundary = rectangular(4, 4, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        root = MeshQuadtree(mesh)
        root.set_last_triangle()

        for x in [[0.6, 0.3], [0.1, 0.2], [0.7,0.7],
                  [0.1,0.9], [0.4,0.6], [0.9,0.1],
                  [10, 3]]:
            
            found, s0, s1, s2, k = root.search_fast(ensure_numeric(x))                                   
                                                           
            if k >= 0:
                V = mesh.get_vertex_coordinates(k) # nodes for triangle k
                assert is_inside_polygon(x, V)
                assert found is True
                #print k, x
            else:
                assert found is False                

        

    def test_large(self):
        """test_larger mesh and different quad trees
        """

        points, vertices, boundary = rectangular(10, 12, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        #Test that points are arranged in a counter clock wise order
        mesh.check_integrity()

        

        root = MeshQuadtree(mesh)
        root.set_last_triangle()
        #print m, root.show()

        for x in [[0.6, 0.3], [0.1, 0.2], [0.7,0.7],
                  [0.1,0.9], [0.4,0.6], [0.9,0.1],
                  [10, 3]]:
            
            found, s0, s1, s2, k = root.search_fast(x)

            if k >= 0:
                V = mesh.get_vertex_coordinates(k) # nodes for triangle k
                assert is_inside_triangle(x, V, closed=True)
                assert is_inside_polygon(x, V)
                assert found is True
            else:
                assert found is False                

        
            if k == 0: return    

    def test_underlying_function(self):
        """test_larger mesh and different quad trees
        """

        points, vertices, boundary = rectangular(2, 2, 1, 1)
        mesh = Mesh(points, vertices, boundary)

        root = MeshQuadtree(mesh)
        root.set_last_triangle()

        # One point
        x = ensure_numeric([0.5, 0.5])

        triangles = root._trilist_from_data(root.search(x))
    
        found, sigma0, sigma1, sigma2, k = \
               root._search_triangles_of_vertices(triangles, x)

        if k >= 0:
            V = mesh.get_vertex_coordinates(k) # nodes for triangle k
            assert is_inside_polygon(x, V)
            assert found is True
        else:
            assert found is False                

        

        # More points    
        for x in [[0.6, 0.3], [0.1, 0.2], [0.7,0.7],
                  [0.1,0.9], [0.4,0.6], [0.9,0.1],
                  [10, 3]]:
                
            triangles = root._trilist_from_data(root.search(x))

            #print x, candidate_vertices
            found, sigma0, sigma1, sigma2, k = \
                   root._search_triangles_of_vertices(triangles,
                                                 ensure_numeric(x))
            if k >= 0:
                V = mesh.get_vertex_coordinates(k) # nodes for triangle k
                assert is_inside_polygon(x, V)
                assert found is True
            else:
                assert found is False

                

    def expanding_search(self):
        """test_larger mesh and different quad trees
        """
        
        p0 = [2,1]
        p1 = [4,1]
        p2 = [4.,4]
        p3 = [2,4]
        p4 = [5,4]

        p5 = [-1,-1]
        p6 = [1,-1]
        p7 = [1,1]
        p8 = [-1,1]

        points = [p0,p1,p2, p3,p4,p5,p6,p7,p8]
        #
        vertices = [[0,1,2],[0,2,3],[1,4,2],[5,6,7], [5,7,8]]
        mesh = Mesh(points, vertices)

        # Don't do this, want to control the max and mins
        #root = build_quadtree(mesh, max_points_per_cell=4)
    

        root = Cell(-3, 9, -3, 9,
                    max_points_per_cell = 4)
        #Insert indices of all vertices
        root.insert( range(mesh.number_of_nodes) )

        #Build quad tree and return
        root.split()
        
        # One point
        #x = [3.5, 1.5]
        x = [2.5, 1.5]
        element_found, sigma0, sigma1, sigma2, k = root.search_fast(x)
        # One point
        x = [3.00005, 2.999994]
        element_found, sigma0, sigma1, sigma2, k = root.search_fast(x)
        assert element_found is True
        assert k == 1
        
        
    def test_compute_interpolation_values(self):
        """test_compute_interpolation_values
        
        Test that interpolation values are correc
        
        This test used to check element_found as output from 
        find_triangle_compute_interpolation(triangle, n0, n1, n2, x)
        and that failed before 18th March 2009.
        
        Now this function no longer returns this flag, so the test
        is merely checknig the sigmas.
        
        """
        
        triangle = num.array([[306951.77151059, 6194462.14596986],
                              [306952.58403545, 6194459.65001246],
                              [306953.55109034, 6194462.0041216]])


        n0 = [0.92499377, -0.37998227]
        n1 = [0.07945684,  0.99683831]
        n2 = [-0.95088404, -0.30954732]
        
        x = [306953.344, 6194461.5]
        
        # Test that point is indeed inside triangle
        assert is_inside_polygon(x, triangle, 
                                 closed=True, verbose=False)
        assert is_inside_triangle(x, triangle, 
                                  closed=True, verbose=False)                                 
        
        sigma0, sigma1, sigma2 = \
            compute_interpolation_values(triangle, n0, n1, n2, x)
            
        msg = 'Point which is clearly inside triangle was not found'
        assert abs(sigma0 + sigma1 + sigma2 - 1) < 1.0e-6, msg

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_search_functions, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
    
