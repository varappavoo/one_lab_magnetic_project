/* -------------------------------------------------------------------
   Tutorial 2 : magnetostatic field of an electromagnet

   Features:
   - Infinite ring geometrical transformation
   - Parameters shared by Gmsh and GetDp, and Onelab parameters
   - FunctionSpaces for the 2D vector potential formulation

   To compute the solution in a terminal:
       getdp electromagnet -solve MagSta_a
       getdp electromagnet -pos Map_a

   To compute the solution interactively from the Gmsh GUI:
       File > Open > electromagnet.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

/* Electromagnetic fields expand to infinity.
   The corresponding boundary condition can be imposed rigorously
   by means of a gometrical transformation that maps a ring (or shell) of finite elements
   to the complementary of its interior.
   As this is a mere geometric transformation,
   it is enough in the model description to attribute a special jacobian
   to the ring region ("AirInf"). See Jacobian{} section below.
   With this information, GetDP is able to deal with the correct transformation
   of all quantities involved in the model.

   The special jacobian "VolSphShell" has parameters.
   There are 2 parameters in this case, "Val_Rint" and "Val_Rext",
   which represent the inner and outer radii of the transformed ring region
   and whose value must match those used
   in the geometrical description of the model (.geo file).
   This is a typical case where Gmsh and GetDP must consistently share parameter values.
   To ensure consistency in all cases, common parameters are defined
   is a specific file "electromagnet_common.pro",
   which is included in both the .geo and .pro file of the model.

   Besides sharing parameters between Gmsh and GetDP,
   it is also useful to share some parameters (not all) with the user of the model,
   i.e., to make them editable in the GUI before running the model.
   Such variables are called Onelab variables (because the sharing mechanism
   between the model and the GUI uses the Onelab interface).
   Onelab parameters are defined with a "DefineNumber" statement,
   which can be invoked in the .geo, .pro, or _common.pro files.
 */


Group {
  // Physical regions:
  Air    = Region[ 101 ];   Core   = Region[ 102 ];
  Ind    = Region[ 103 ];   AirInf = Region[ 111 ];

  Surface_ht0 = Region[ 1100 ];
  Surface_bn0 = Region[ 1101 ];
  Surface_Inf = Region[ 1102 ];


  /* Abstract regions :
     The purpose of abstract regions is to allow a generic definition of
     the FunctionSpace, Formulation and PostProcessing fields
     with no reference to model-specific Physical regions.
     We will show in a later tutorial how abstract formulations can then be isolated
     in geometry independent template files, thanks to an appropriate declaration mechanism
     (using DefineConstant[], DefineGroup[] and DefineFunction[]).

     The abstract regions in this model have the following interpretation:
     - Vol_Nu_Mag  = region where the term [ nu[] * Dof{d a} , {d a} ] is assembled
     - Vol_Js_Mag  = region where the term [ - Dof{js} , {a} ] is assembled
     - Vol_Inf_Mag = region where the infinite ring geometric transformation is applied
     - Sur_Dir_Mag = Homogeneous Dirichlet part of the model's boundary;
     - Sur_Neu_Mag = Homogeneous Neumann part of the model's boundary;
  */
  Vol_Nu_Mag  = Region[ {Air, AirInf, Core, Ind} ];
  Vol_Js_Mag  = Region[ Ind ];
  Vol_Inf_Mag = Region[ {AirInf} ];
  Sur_Dir_Mag = Region[ {Surface_bn0, Surface_Inf} ];
  Sur_Neu_Mag = Region[ {Surface_ht0} ];
}

Function {
  mu0 = 4.e-7 * Pi;
  murCore = DefineNumber[100, Name "Model parameters/Mur core",
			 Help "Magnetic relative permeability of Core"];

  nu [ Region[{Air, Ind, AirInf}] ]  = 1. / mu0;
  nu [ Core ]  = 1. / (murCore * mu0);

  NbTurns = 1000 ;
  Current = DefineNumber[0.01, Name "Model parameters/Current",
			 Help "Current injected in coil [A]"];
  Js_fct[ Ind ] = -NbTurns*Current/SurfaceArea[];
  /* The minus sign is to have the current in -e_z direction,
     so that the magnetic induction field is in +e_y direction */
}

/* In the 2D approximation, the magnetic vector potential A and the current density Js
   are not scalars, but vectors with a z-component only:
   A  = Vector [ 0, 0, az(x,y,t) ]
   Js = Vector [ 0, 0, jsz(x,y,t) ]
   In order to compute derivatives and apply geometric transformations adequately,
   GetDP needs this information.
   Regarding discretization, now, A is node-based, whereas Js is region-wise constant.
   Considering all this, one ends then up with the FunctionSpaces
   "Hcurl_a_Mag_2D" and "Hregion_j_Mag_2D" as they are defined below.

   The function space "Hregion_j_Mag_2D" provides one basis function,
   and hence one degree of freedom, per physical region in the abstract region "Vol_Js_Mag".
   The constraint "SourceCurrentDensityZ" fixes all these dofs,
   so the FunctionSpace "Hregion_j_Mag_2D" is fully fixed and has no FE unknowns.
   One could thus have replaced it by a simple function
   and the Galerkin term would have been

   Galerkin { [ Vector[ 0,0,-Js_fct[] ] , {a} ]; In Vol_Js_Mag;
              Jacobian Vol; Integration Int; }

   instead of

   Galerkin { [ - Dof{js} , {a} ]; In Vol_Js_Mag;
               Jacobian Vol; Integration Int; }

   Thechosen implementation below is however more effeicient
   as it avoids evaluating repeatedly the function Js_fct[] during assembly.
 */


Constraint {
  { Name Dirichlet_a_Mag;
    Case {
      { Region Sur_Dir_Mag ; Value 0.; }
    }
  }
  { Name SourceCurrentDensityZ;
    Case {
      { Region Vol_Js_Mag ; Value Js_fct[]; }
    }
  }
}

Group {
  Dom_Hcurl_a_Mag_2D = Region[ {Vol_Nu_Mag, Sur_Neu_Mag} ];
}
FunctionSpace {
  { Name Hcurl_a_Mag_2D; Type Form1P; // Magnetic vector potential A
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
        Support Dom_Hcurl_a_Mag_2D ; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef ae; EntityType NodesOf;
        NameOfConstraint Dirichlet_a_Mag; }
    }
  }

  { Name Hregion_j_Mag_2D; Type Vector; // Electric current density Js
    BasisFunction {
      { Name sr; NameOfCoef jsr; Function BF_RegionZ;
        Support Vol_Js_Mag; Entity Vol_Js_Mag; }
    }
    Constraint {
      { NameOfCoef jsr; EntityType Region;
        NameOfConstraint SourceCurrentDensityZ; }
    }
  }

}

Include "electromagnet_common.pro";
Val_Rint = rInt; Val_Rext = rExt;
Jacobian {
  { Name Vol ;
    Case { { Region Vol_Inf_Mag ;
             Jacobian VolSphShell {Val_Rint, Val_Rext} ; }
           { Region All ; Jacobian Vol ; }
    }
  }
}

Integration {
  { Name Int ;
    Case { {Type Gauss ;
	Case { { GeoElement Triangle    ; NumberOfPoints  4 ; }
	       { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
	}
      }
    }
  }
}

Formulation {
  { Name Magnetostatics_a_2D; Type FemEquation;
    Quantity {
      { Name a ; Type Local; NameOfSpace Hcurl_a_Mag_2D; }
      { Name js; Type Local; NameOfSpace Hregion_j_Mag_2D; }
    }
    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ]; In Vol_Nu_Mag;
                 Jacobian Vol; Integration Int; }
      Galerkin { [ -Dof{js} , {a} ]; In Vol_Js_Mag;
                 Jacobian Vol; Integration Int; }
    }
  }
}

Resolution {
  { Name MagSta_a;
    System {
      { Name Sys_Mag; NameOfFormulation Magnetostatics_a_2D; }
    }
    Operation {
      Generate[Sys_Mag]; Solve[Sys_Mag]; SaveSolution[Sys_Mag];
    }
  }
}

PostProcessing {
  { Name MagSta_a_2D; NameOfFormulation Magnetostatics_a_2D;
    Quantity {
      { Name a;
        Value {
          Local { [ {a} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name az;
        Value {
          Local { [ CompZ[{a}] ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name b;
        Value {
          Local { [ {d a} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name h;
        Value {
          Local { [ nu[] * {d a} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name js;
        Value {
          Local { [ {js} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
    }
  }
}

e = 1.e-5;
h = 0.02;
p1 = {e,h,0};
p2 = {0.25-e,h,0}; // horizontal cut through model, just above x-axis.

PostOperation {

  { Name Map_a; NameOfPostProcessing MagSta_a_2D;
    Operation {
      Echo[ Str["l=PostProcessing.NbViews-1;",
		"View[l].IntervalsType = 1;",
		"View[l].NbIso = 40;"],
	    File "tmp.geo", LastTimeStepOnly] ;
      Print[ a, OnElementsOf Dom_Hcurl_a_Mag_2D, File "a.pos" ];
      Print[ js, OnElementsOf Dom_Hcurl_a_Mag_2D, File "js.pos" ];
      Print[ az, OnElementsOf Dom_Hcurl_a_Mag_2D, File "az.pos" ];
      Print[ b, OnElementsOf Dom_Hcurl_a_Mag_2D, File "b.pos" ];
      Print[ b, OnLine{{List[p1]}{List[p2]}} {50}, File "by.pos" ];
    }
  }
}
