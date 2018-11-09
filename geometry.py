"""
this file deals with creating appropriate gmsh file(s)
"""

def create_geometry(fileDir, H, W, L, dHf, dHs):
    
    # open the .geo file for writing
    fid = open(fileDir + 'geometry.geo', 'wt')

    # define constants
    print("L = %f;" % (L), file=fid)
    print("H = %f;" % (H), file=fid)
    print("W = %f;" % (W), file=fid)
    print("dHf = %f;" % (dHf), file=fid)
    print("dHs = %f;" % (dHs), file=fid)

    string = """
Mesh.Algorithm = 8;
Point(1) = {0.0, W, 0.0, dHf};
Point(2) = {0.0, W - H, 0.0, dHf};
Point(3) = {  L, W - H, 0.0, dHf};
Point(4) = {  L, W, 0.0, dHf};
Point(5) = {0.0, 0.0, 0.0, dHs};
Point(6) = {  L, 0.0, 0.0, dHs};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 2};
Line(6) = {3, 6};
Line(7) = {6, 5};
Line Loop(8) = {1, 2, 3, 4};
Line Loop(9) = {5, 2, 6, 7};

Plane Surface(1) = {8};
Plane Surface(2) = {9};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
    """

    print("%s" % (string), file=fid)

    fid.close()

def create_geometry_quad(fileDir, H, W, L, dHf, dHs):
    
    # open the .geo file for writing
    fid = open(fileDir + 'geometry_quad.geo', 'wt')
    
    # define constants
    print("L = %f;" % (L), file=fid)
    print("H = %f;" % (H), file=fid)
    print("W = %f;" % (W), file=fid)
    print("dHf = %f;" % (dHf), file=fid)
    print("dHs = %f;" % (dHs), file=fid)
    
    string = """
Mesh.Algorithm = 8;
Point(1) = {0.0, W, 0.0, dHf};
Point(2) = {0.0, W - H, 0.0, dHf};
Point(3) = {  L, W - H, 0.0, dHf};
Point(4) = {  L, W, 0.0, dHf};
Point(5) = {0.0, 0.0, 0.0, dHs};
Point(6) = {  L, 0.0, 0.0, dHs};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 2};
Line(6) = {3, 6};
Line(7) = {6, 5};
Line Loop(8) = {1, 2, 3, 4};
Line Loop(9) = {5, 2, 6, 7};

Plane Surface(1) = {8};
Plane Surface(2) = {9};
Recombine Surface{1};
Recombine Surface{2};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
"""
    
    print("%s" % (string), file=fid)
    
    fid.close()



