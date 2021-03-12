// Gmsh project created on Mon Jul  6 11:21:02 2020
//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
//+
Transfinite Surface {1};
Recombine Surface {1};

Extrude {0,0,1} {
  Surface{1}; Layers{20}; Recombine;
}

Physical Surface("bottom", 1) = {1};
Physical Surface("left", 2) = {2};
Physical Surface("face", 3) = {3};
Physical Surface("right", 4) = {4};
Physical Surface("back", 5) = {5};
Physical Surface("top", 6) = {6};

Physical Volume("cube, 1") = {1};


Characteristic Length {1, 2, 3, 4} = 0.05;
