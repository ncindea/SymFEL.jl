// Gmsh project created on Mon Jul  6 11:21:02 2020
//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Surface("square", 1) = {1};
//+
Transfinite Surface {1};
Recombine Surface {1};

//+
Physical Curve("top", 3) = {3};
//+
Physical Curve("right", 2) = {2};
//+
Physical Curve("bottom", 1) = {1};
//+
Physical Curve("left", 4) = {4};
//+
Characteristic Length {1, 2, 3, 4} = 0.05;
