a = 10;
b = 8;

// base corner points
Point(1) = {0,0,0};
Point(2) = {a,0,0};
Point(3) = {a,b,0};
Point(4) = {0,b,0};

// base lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// base surface
Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};

// extra work for creating structured mesh
Transfinite Line {1,3} = 30;
Transfinite Line {2,4} = 30;
Transfinite Surface {1} = {1,2,3,4};
Recombine Surface {1};


// extrude base surface in z direction for create 
// dummy 3d geometry for OpenFoam
Extrude {0, 0, 1} {
   Surface{1}; Layers{1}; Recombine;
}


// label boundary patches to get named paches for OpenFoam
Physical Surface("left") = {17};
Physical Surface("right") = {25};
Physical Surface("top") = {13};
Physical Surface("bottom") = {21};
Physical Surface("frontAndBack") = {1,26};


// Don't forget to add a physical volume at the end
Surface Loop(1) = {26, 13, 1, 17, 21, 25};
Physical Volume("volume") = {1};
