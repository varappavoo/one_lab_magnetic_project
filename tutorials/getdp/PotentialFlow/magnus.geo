Include "magnus_common.pro";

A = (BoxSize-1)/2.;
B = A;
R = 0.5;
lc  = A/10;


If( Flag_Object == 0 ) // Cylinder

lcr = 0.1;
Point(1) = { 2*R, 0.0, 0.0, lcr};
Point(2) = {   R,  -R, 0.0, lcr};
Point(3) = { 0.0, 0.0, 0.0, lcr};
Point(4) = {   R,   R, 0.0 , lcr};
Point(5) = {   R, 0.0, 0.0 , lcr};
Circle(1) = {1,5,2};
Circle(2) = {2,5,3};
Circle(3) = {3,5,4};
Circle(4) = {4,5,1};

// Points to be connected with the outer boundary
PtA = 4;
PtB = 2; 

Else // naca airfoil

lca = 0.03;
Include "nacaAirFoil.geo";
PtA = 121;
PtB = 81; 

EndIf

Point(306)  = { 0, B, 0, lc};
Point(307)  = { 0,-B, 0, lc};

Line(5) = { PtA, 306 };
Line(6) = { PtB, 307 };

Point(308) = {-A,-B, 0, lc};
Point(309) = {-A, B, 0, lc};
Point(310) = { A+1, B, 0, lc};
Point(311) = { A+1,-B, 0, lc};

Line( 7) = { 306, 309 };
Line( 8) = { 309, 308 };
Line( 9) = { 308, 307 };
Line(10) = { 307, 311 };
Line(11) = { 311, 310 };
Line(12) = { 310, 306 };

Line Loop(21) = { 7, 8, 9, -6, 2, 3, 5 }; // aire a gauche
Line Loop(22) = { 10, 11, 12, -5, 4, 1, 6 };

Plane Surface(24) = { 21 };
Plane Surface(25) = { 22 };

Physical Surface("Fluid", 2) = { 24, 25 };
Physical Line("UpStream", 10) = { 8 };
Physical Line("DownStream", 11) = { 11 };
Physical Line("Airfoil", 12) = { 1 ... 4 };
Physical Line("Wake", 13) = { 5 };



//Line Loop(20) = {1,2,3,4};
//Plane Surface(23) = {20};
//Physical Surface("Cylinder", 1) = {23};
//Physical Line("Outer", 10) = { 7 ... 12 };

