
// Onelab parameters

Flag_Regularization = 
  DefineNumber[0, Name "Options/Regularize field", Choices {0,1} ];
Flag_QFromFile = 
  DefineNumber[0, Name "Options/Heat source from file", Choices {0,1} ];
Flag_ConstraintWin2 = 
  DefineNumber[0, Name "Options/Constraint in Window2",
	       Choices {0="Fixed temperature", 1="Fixed heat flux"} ];
MeshRefinement = 
  DefineNumber[1, Name "Options/0Mesh refinement", 
	       Help "Choose 1 for a coarse mesh and 0.1 for a fine mesh."];

// Geometrical dimensions

mm=1e-3; // mm to m conversion factor

dx_Brick=100*mm; dy_Brick= 50*mm;
e_layer = 1*mm;
dx_Win1 = 20*mm; dy_Win1 = 20*mm;
If( Flag_Regularization )
  dx_Win1 += e_layer;
  dy_Win1 += e_layer;
EndIf

dx_Win2 = 20*mm; dy_Win2 = 20*mm;
xc_Win1 = 25*mm; yc_Win1 = dy_Brick/2;
xc_Win2 = 75*mm; yc_Win2 = dy_Brick/2;

// Element sizes

s=1;
lc_Brick = dy_Brick/20 *s;
lc_Win1 = dx_Win1/10 *s;
lc_Win2 = dx_Win2/10 *s;

