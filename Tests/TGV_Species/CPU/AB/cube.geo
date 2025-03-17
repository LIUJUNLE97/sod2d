pi = 3.14159265;
n = 12;
h = 2*pi/n;

//+ Defiine the geometry
Point(1) = {0.0,0.0,0.0,h};
Point(2) = {2*pi,0.0,0.0,h};
Point(3) = {2*pi,2*pi,0.0,h};
Point(4) = {0.0,2*pi,0.0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Surface(1) = {1};

Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, 2*pi} {
  Surface{1}; Layers{n}; Recombine;
}

//+ Set every surface to periodic, using surf(1) as master
Physical Surface("Periodic") = {21,13,17,25,26,1};

//+ Set the fluid region
Physical Volume("fluid") = {1};

//+ Set the output mesh file version
Mesh.MshFileVersion = 2.2;

//+ Options controlling mesh generation
Mesh.ElementOrder = 4; //+ Set desired element order beetween 3 and 8 (1,2 not supported in SOD)
Mesh 3;                //+ Volumetric mesh

//+ Generate the periodicity between surface pairs
Periodic Surface {17} = {25} Translate {2*pi, 0, 0};
Periodic Surface {21} = {13} Translate {0, 2*pi, 0};
Periodic Surface {26} = {1}  Translate {0, 0, 2*pi};
