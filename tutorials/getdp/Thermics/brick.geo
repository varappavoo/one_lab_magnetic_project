/* -------------------------------------------------------------------
   File "brick.geo"

   This file is the geometrical description used by GMSH to produce
   the file "brick.msh".
   ------------------------------------------------------------------- */

/* Gmsh options can be set directly in the .geo file.
   Setting "Solver.AutoMesh" to 2 ensures that GMSH systematically 
   regenerates the mesh at each model execution, and does not reuse
   the mesh on disk if it exists (which is the default option).
   This option is needed in this model, because the interactive
   model parameters "Flag_Regularization" has an effect 
   on the definition of the regions, which makes remeshing mandatory. */

Solver.AutoMesh = 2; 

Include "brick_common.pro";

/* The (very) simple geometry of this thermal model only contains
   rectangles. It is the opportunity to illustrate
   the function definition capabilities of Gmsh. */

Function Def_Rectangle
p1=newp; Point(newp)={xc_, yc_, 0, lc_};
p2=newp; Point(newp)={xc_+dx_, yc_, 0, lc_};
p3=newp; Point(newp)={xc_+dx_, yc_+dy_, 0, lc_};
p4=newp; Point(newp)={xc_, yc_+dy_, 0, lc_};

l1=newl; Line(newl)={p1,p2};
l2=newl; Line(newl)={p2,p3};
l3=newl; Line(newl)={p3,p4};
l4=newl; Line(newl)={p4,p1};

lines_[] = {l1, l2, l3, l4};

ll1 = newll; Line Loop(ll1) = {lines_[]};
s_ = news; Plane Surface(news) = {ll1, ll_Holes_[]};
Return  // end of Function Def_Rectangle

/* Note that the code above is not parsed before it is called. 
   Good practice for readability (but by no means mandatory)
   is to distinguish function variables from other variables.  
   This is done here with a trailing "_" in the variable name.  */


// Window 1 (with optional layer of thickness "e_layer")
dx_ = dx_Win1; 
dy_ = dy_Win1;
lc_ = lc_Win1;
xc_ = xc_Win1-dx_/2; yc_=yc_Win1-dy_/2;
ll_Holes_[]={};
Call Def_Rectangle;
s_Win1 = s_;
l_Win1[] = lines_[];

If( !Flag_Regularization )
  ll_HolesForBrick[] += ll1;
Else
  ll_HolesForLayer[] += ll1;
  dx_ = dx_Win1+2*e_layer; 
  dy_ = dy_Win1+2*e_layer;
  lc_ = lc_Win1;
  xc_ = xc_Win1-dx_/2; yc_=yc_Win1-dy_/2;

  ll_Holes_[] = {ll_HolesForLayer[]};
  Call Def_Rectangle;
  s_Win1_Layer = s_;
  l_Win1_Layer[] = lines_[];
  ll_HolesForBrick[] += ll1;
EndIf

// Window 2
dx_ = dx_Win2; 
dy_ = dy_Win2;
lc_ = lc_Win2;
xc_ = xc_Win2-dx_/2; yc_=yc_Win2-dy_/2;
ll_Holes_[]={};
Call Def_Rectangle;
s_Win2 = s_;
l_Win2[] = lines_[];
ll_HolesForBrick[] += ll1;

// Brick
dx_ = dx_Brick; dy_ = dy_Brick;
lc_ = lc_Brick;
xc_ = 0; yc_ = 0;
ll_Holes_[] = {ll_HolesForBrick[]};
Call Def_Rectangle;
s_Brick = s_;
l_Brick[] = lines_[];

//Printf("l_Brick", l_Brick[]);


// Physical regions

Physical Surface("Brick", 100) = {s_Brick};
Physical Surface("Window1", 111) = {s_Win1};
Physical Surface("Window2", 112) = {s_Win2};
If( Flag_Regularization )
  Physical Surface("LayerWindow1", 115) = {s_Win1_Layer};
EndIf

Physical Line("Surface1", 201) = { l_Brick[{3}] };
Physical Line("Surface2", 202) = { l_Brick[{1}] };
Physical Line("Surface3", 203) = { l_Brick[{0}], l_Brick[{2}] };

Physical Line("SurfWindow1", 211) = {l_Win1[]};
Physical Line("SurfWindow2", 212) = {l_Win2[]};

