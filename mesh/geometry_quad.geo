L = 10.000000;
H = 0.050000;
W = 5.000000;
dHf = 0.025000;
dHs = 0.375000;

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

