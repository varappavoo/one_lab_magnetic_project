/* -------------------------------------------------------------------
   Tutorial 3 : linear elastic model of a wrench

   Features:
   - "Grad u" GetDP specific formulation for linear elasticity
   - first and second order elements 
   - triangular and quadrangular elements

   To compute the solution interactively from the Gmsh GUI:
       File > Open > wrench.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

/* Linear elasticity with GetDP:

   GetDP has a peculiar way to deal with linear elasticity.
   Instead of vector field "u = Vector[ ux, uy, uz ]",
   the displacement field is regarded as two (2D case) or 3 (3D case) scalar fields.
   Unlike conventional formulations, GetDP's formulation is written 
   in terms of the gradient "Grad u" of the displacement field,
   which is a non-symmetric tensor, and the needed symmetrization 
   (to define the strain tensor and relate it to the stress tensor) 
   is done through the constitutive relationship (Hooke law). 
   The reason for this unusual formulation is to be able to use also for elastic problems
   the powerful geometrical and homological kernel of GetDP,
   which relies on the operators Grad, Curl and Div.

   The "Grad u" formulation entails a small increase of assembly work
   but makes in counterpart lots of geometrical features
   implemented in GetDP (change of coordinates, function spaces, etc...)
   applicable to elastic problems out-of-the-box, 
   since the scalar fields { ux, uy, uz } have exactly the same geometrical properties
   as, e.g. a scalar electrip potential or a temperature field. 
*/


Include "wrench2D_common.pro";

Young = 1e9 * 
  DefineNumber[ 200, Name "Material/Young modulus [GPa]"];
Poisson = 
  DefineNumber[ 0.3, Name "Material/Poisson coefficient []"];
AppliedForce = 
  DefineNumber[ 100, Name "Material/Applied force [N]"];

// Approximation of the maximum deflection by an analytical model: 
// Deflection = PL^3/(3EI) with I = Width^3*Thickness/12
Deflection = 
  DefineNumber[4*AppliedForce*((LLength-0.018)/Width)^3/(Young*Thickness)*1e3, 
	       Name "Solution/Deflection (analytical) [mm]", ReadOnly 1];

Group {
  // Physical regions: Give explicit labels to the regions defined in the .geo file 
  Wrench = Region[ 1 ];
  Grip = Region[ 2 ];
  Force = Region[ 3 ];

  // Abstract regions:
  Vol_Elast_Mec = Region[ { Wrench } ];
  Vol_Force_Mec = Region[ { Wrench } ];
  Sur_Clamp_Mec = Region[ { Grip } ];
  Sur_Force_Mec = Region[ { Force } ];
/* 
  Signification of the abstract regions:
  Vol_Elast_Mec       Elastic domain
  Vol_Force_Mec       Region with imposed volumic force
  Sur_Force_Mec       Surface with imposed surface traction
  Sur_Clamp_Mec       Surface with imposed zero displacements (all components)
*/
}

Function {
  /* Material coefficients. 
     No need to define them regionwise here ( E[{Wrench}] = ... ; )
     as there is only one region in this model. */
  E[] = Young;
  nu[] = Poisson;
  /* Components of the volumic force applied to the region "Vol_Force_Mec"
     Gravity could be defined here ( force_y[] = 7000*9.81; ) ; */
  force_x[] = 0;
  force_y[] = 0;
  /* Components of the surface traction force applied to the region "Sur_Force_Mec" */
  pressure_x[] = 0;
  pressure_y[] = -AppliedForce/(SurfaceArea[]*Thickness); // downward vertical force
}


/* Hooke law

   The material law 

   sigma_ij = C_ijkl epsilon_ij 

   is represented in 2D by 4 2x2 tensors C_ij[], i,j=1,2 
   depending on the Lamé coefficients of the isotropic linear material,

   lambda = E[]*nu[]/(1.+nu[])/(1.-2.*nu[]);
   mu = E[]/2./(1.+nu[]);

   as follows

   EPC:  a[] = E/(1-nu^2)        b[] = mu     c[] = E nu/(1-nu^2)
   EPD:  a[] = lambda + 2 mu     b[] = mu     c[] = lambda
    3D:  a[] = lambda + 2 mu     b[] = mu     c[] = lambda

   respectively for the 2D plane strain (EPD), 2D plane stress (EPS) and 3D cases.
*/

Function {
  Flag_EPC = 1; 
  If(Flag_EPC) // Plane stress
    a[] = E[]/(1.-nu[]^2);
    c[] = E[]*nu[]/(1.-nu[]^2);
  Else // Plane strain or 3D
    a[] = E[]*(1.-nu[])/(1.+nu[])/(1.-2.*nu[]);
    c[] = E[]*nu[]/(1.+nu[])/(1.-2.*nu[]);
  EndIf
  b[] = E[]/2./(1.+nu[]);

  C_xx[] = Tensor[ a[],0  ,0  ,    0  ,b[],0  ,    0  ,0  ,b[] ];
  C_xy[] = Tensor[ 0  ,c[],0  ,    b[],0  ,0  ,    0  ,0  ,0   ];

  C_yx[] = Tensor[ 0  ,b[],0  ,    c[],0  ,0  ,    0  ,0  ,0   ];
  C_yy[] = Tensor[ b[],0  ,0  ,    0  ,a[],0  ,    0  ,0  ,b[] ];
}

/* Clamping boundary condition */
Constraint {
  { Name Displacement_x;
    Case {
      { Region Sur_Clamp_Mec ; Type Assign ; Value 0; }
    }
  }
  { Name Displacement_y;
    Case {
       { Region Sur_Clamp_Mec ; Type Assign ; Value 0; }
    }
  }
}

/* As explained above, the displacement field is discretized 
   as two scalar fields "ux" and "uy", which are the spatial components
   of the vector field "u" in a fixed Cartesian coordinate system. 

   Boundary conditions like

   ux = ... ;
   uy = ... ;

   translate naturally into Dirichlet constraints 
   on the scalar "ux" and "uy" FunctionSpaces. 
   Conditions like, however,

   u . n = ux Cos [th] + uy Sin [th] = ... ;

   are less naturally accounted for within the "Grad u" formulation.
   Fortunately, they are rather uncommon. 

   Finite element shape (triangles or quadrangles) makes no difference
   in the definition of the FunctionSpaces. The appropriate shape functions to be used
   are determined by GetDP at a much lower level on basis of the information
   contained in the *.msh file. 

   Second order elements, on the other hand, are implemented in the hierarchical fashion
   by adding to the first order node-based shape functions 
   a set of second order edge-based functions
   to complete a basis for 2d order polynomials on the reference element. 
 */


// Domain of definition of the "ux" and "uy" FunctionSpaces
Group {
  Dom_H_u_Mec = Region[ { Vol_Elast_Mec, 
			  Sur_Force_Mec, 
			  Sur_Clamp_Mec} ];
}

Flag_Degree = 
  DefineNumber[ 0, Name "Geometry/Use degree 2 (hierarch.)", Choices{0,1}, Visible 1];
FE_Degree = ( Flag_Degree == 0 ) ? 1 : 2; // Convert flag value into polynomial degree

FunctionSpace {
  { Name H_ux_Mec ; Type Form0 ;
    BasisFunction {
      { Name sxn ; NameOfCoef uxn ; Function BF_Node ;
        Support Dom_H_u_Mec ; Entity NodesOf[ All ] ; }
     If ( FE_Degree == 2 )
        { Name sxn2 ; NameOfCoef uxn2 ; Function BF_Node_2E ;
          Support Dom_H_u_Mec; Entity EdgesOf[ All ] ; }
     EndIf
    }
    Constraint {
      { NameOfCoef uxn ;
        EntityType NodesOf ; NameOfConstraint Displacement_x ; }
      If ( FE_Degree == 2 )
         { NameOfCoef uxn2 ;
	   EntityType EdgesOf ; NameOfConstraint Displacement_x ; }
     EndIf
    }
  }
  { Name H_uy_Mec ; Type Form0 ;
    BasisFunction {
      { Name syn ; NameOfCoef uyn ; Function BF_Node ;
        Support Dom_H_u_Mec ; Entity NodesOf[ All ] ; }
     If ( FE_Degree == 2 )
        { Name syn2 ; NameOfCoef uyn2 ; Function BF_Node_2E ;
          Support Dom_H_u_Mec; Entity EdgesOf[ All ] ; }
     EndIf
    }
    Constraint {
      { NameOfCoef uyn ;
        EntityType NodesOf ; NameOfConstraint Displacement_y ; }
      If ( FE_Degree == 2 )
      { NameOfCoef uyn2 ;
        EntityType EdgesOf ; NameOfConstraint Displacement_y ; }
      EndIf
    }
  }
}


Jacobian {
  { Name Vol;
    Case {
      { Region All; Jacobian Vol; }
    }
  }
  { Name Sur;
    Case {
      { Region All; Jacobian Sur; }
    }
  }
}

/* Adapt the number of Gauss points to the polynomial degree of the finite elements
is as simple as this: */
If (FE_Degree == 1)
Integration {
  { Name Gauss_v;
    Case {
      {Type Gauss;
        Case {
	  { GeoElement Line       ; NumberOfPoints  3; }
          { GeoElement Triangle   ; NumberOfPoints  3; }
          { GeoElement Quadrangle ; NumberOfPoints  4; }
        }
      }
    }
  }
}
ElseIf (FE_Degree == 2)
Integration {
  { Name Gauss_v;
    Case {
      {Type Gauss;
        Case {
	  { GeoElement Line       ; NumberOfPoints  5; }
          { GeoElement Triangle   ; NumberOfPoints  7; }
          { GeoElement Quadrangle ; NumberOfPoints  7; }
        }
      }
    }
  }
}
EndIf

Formulation {
  { Name Elast_u ; Type FemEquation ;
    Quantity {
      { Name ux  ; Type Local ; NameOfSpace H_ux_Mec ; }
      { Name uy  ; Type Local ; NameOfSpace H_uy_Mec ; }
    }
    Equation {
      Galerkin { [ -C_xx[] * Dof{d ux}, {d ux} ] ;
        In Vol_Elast_Mec ; Jacobian Vol ; Integration Gauss_v ; }
      Galerkin { [ -C_xy[] * Dof{d uy}, {d ux} ] ;
        In Vol_Elast_Mec ; Jacobian Vol ; Integration Gauss_v ; }
      Galerkin { [ -C_yx[] * Dof{d ux}, {d uy} ] ;
        In Vol_Elast_Mec ; Jacobian Vol ; Integration Gauss_v ; }
      Galerkin { [ -C_yy[] * Dof{d uy}, {d uy} ] ;
        In Vol_Elast_Mec ; Jacobian Vol ; Integration Gauss_v ; }

      Galerkin { [ force_x[] , {ux} ];
        In Vol_Force_Mec ; Jacobian Vol ; Integration Gauss_v ; }
      Galerkin { [ force_y[] , {uy} ];
        In Vol_Force_Mec ; Jacobian Vol ; Integration Gauss_v ; }
      
      Galerkin { [ pressure_x[] , {ux} ];
        In Sur_Force_Mec ; Jacobian Sur ; Integration Gauss_v ; }
      Galerkin { [ pressure_y[] , {uy} ];
        In Sur_Force_Mec ; Jacobian Sur ; Integration Gauss_v ; }
    }
  }
}

Resolution {
  { Name Elast_u ;
    System {
      { Name Sys_Mec ; NameOfFormulation Elast_u; }
    }
    Operation {
      InitSolution [Sys_Mec];
      Generate[Sys_Mec];
      Solve[Sys_Mec];
      SaveSolution[Sys_Mec] ;
    }
  }
}

PostProcessing {
  { Name Elast_u ; NameOfFormulation Elast_u ;
    PostQuantity {
      { Name u ; Value { Term { [ Vector[ {ux}, {uy}, 0 ]]; 
	      In Vol_Elast_Mec ; Jacobian Vol ; } } }
      { Name uy ; Value { Term { [ 1e3*{uy} ]; 
	      In Vol_Elast_Mec ; Jacobian Vol ; } } }
      { Name sig_xx ; Value { Term { 
	    [ CompX[  C_xx[]*{d ux} + C_xy[]*{d uy} ] ]; 
	    In Vol_Elast_Mec ; Jacobian Vol ; } } }
      { Name sig_xy ; Value { Term { 
	    [ CompY[  C_xx[]*{d ux} + C_xy[]*{d uy} ] ]; 
	    In Vol_Elast_Mec ; Jacobian Vol ; } } }
      { Name sig_yy ; Value { Term { 
	    [  CompY [ C_yx[]*{d ux} + C_yy[]*{d uy} ] ]; 
	    In Vol_Elast_Mec ; Jacobian Vol ; } } }
    }
  }
}



PostOperation {
  { Name pos; NameOfPostProcessing Elast_u;
    Operation {

      If(FE_Degree == 1) 
	Print[ sig_xx, OnElementsOf Wrench, File "sigxx.pos" ];
        Print[ u, OnElementsOf Wrench, File "u.pos" ];
      Else
	Print[ sig_xx, OnElementsOf Wrench, File "sigxx2.pos" ];
        Print[ u, OnElementsOf Wrench, File "u2.pos" ];
      EndIf
      Echo[ StrCat["l=PostProcessing.NbViews-1; ", 
      		   "View[l].VectorType = 5; ",
      		   "View[l].ExternalView = l; ",
      		   "View[l].DisplacementFactor = 200; ",
      		   "View[l-1].IntervalsType = 3; "
		   ],
             File "tmp.geo", LastTimeStepOnly] ;
      //Print[ sig_yy, OnElementsOf Wrench, File "sigyy.pos" ];
      //Print[ sig_xy, OnElementsOf Wrench, File "sigxy.pos" ];
      Print[ uy, OnPoint{probe_x, probe_y, 0}, 
	     File > "deflection.pos", Format TimeTable,
	     SendToServer "Solution/Deflection (computed) [mm]", Color "AliceBlue" ];
    }
  }
}


// Tell Gmsh which GetDP commands to execute when running the model.
DefineConstant[
  R_ = {"Elast_u", Name "GetDP/1ResolutionChoices", Visible 0},
  P_ = {"pos", Name "GetDP/2PostOperationChoices", Visible 0},
  C_ = {"-solve -pos -v2", Name "GetDP/9ComputeCommand", Visible 0}
];

