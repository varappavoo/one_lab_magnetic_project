/* -------------------------------------------------------------------
   Tutorial 5 : thermal problem with contact resistances

   Features:
   - Contact resistances: scalar FunctionSpace with a surface discontinuity
   - Region with a uniform temperature (infinite thermal conductivity)
   - Computation of the heat flux through surfaces
   - Import a source field from a file

   To compute the solution in a terminal: 
     getdp brick.pro -solve Thermal_T -pos Map_T

   To compute the solution interactively from the Gmsh GUI:
       File > Open > brick.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

/* This model is a rectangular brick with two windows, 
   where various kinds of thermal constraints can be set.
   Dirichlet, Neumann and convection boundary conditions are imposed
   on different parts of the surface of the brick. 
   The model is rather academic but it
   demonstrates some useful high-level GetDP features. 

   Governing eqations are 

   div ( -lambda[] grad T ) = Q    in Vol_Lambda_The 
   -lambda[] grad T . n = qn = 0   on Sur_Neumann_The

   Contact thermal resistance:
   First, it is shown how to implement contact thermal resistances
   with surface elements. The surface elements are associated with a 
   thickness and a thermal conductivity (typically much lower than that 
   of surrounding regions). The implementation takes advantage
   of the powerful FunctionSpace definition in GetDP. 

   With the flag "Flag_Regularization", the contact surface
   can be, for the sake of comparison, replaced
   by a thin volume conducting region. 
   
   Thermal "electrode":
   The floating potential idea (introduced in tutorial 4)
   is reconsidered here in a thermal context to represent a region 
   with a very large thermal conductivity where, consequently, 
   the temperature field is uniform (exactly like the electric potential
   is uniform on an electrode).
   The dual quantity of this uniform temperature "T_electrode" [K]
   (which is the "associated global quantity" in GetDP language) 
   is the heat flux "Q_electrode" [W] injected in the electrode
   by the agent that maintains the temperature equal
   to the prescribed value.    
   
   The value of Q_electrode is a by-product of the system resolution 
   provided the term 

   GlobalTerm { [-Dof{Q_electrode} , {T_electrode} ] ; In Tfloating_The ; }

   is present in the "Resolution" section. This term triggers the writing
   in the linear system of a supplementary equation associated 
   with the global basis function BF{T_electrode}.
   All integrations are automatically done by getDP,  
   and the value of Q_electrode is obtained in postprocessing 
   with the PostOperation

   Print[ Q_electrode, OnRegion Tfloating_The, ... ]

   Heat flux through surfaces:
   The purpose of a thermal simulation usually goes beyond
   the mere calculation of a temperature distribution.
   One is in general also interested in evaluating 
   the heat flux q(S) through some specific surface S:

   q(S) = ( -lambda[] grad T . n )_S
  
   This quantity cannot be computed from the temperature
   distribution available on the surface S only. 
   As heat flux is related with the gradient of temperature 
   in the direction normal to the surface, its computation relies on 
   the temperature distribution in a neighborhood of the surface. 
   This means that volume elements in contact with the considered surface
   need be involved in the computation.
   To achieve this with getDP, a good method proceeds
   by the definition a smooth auxiliary function g(S),
   with g(S)=1 on S, and g(S)=0 outside a finite neighborhood of S.
   Typically, g(S) is the sum of the shape functions of the nodes on S.
   Let w(S) be the support of g(S), 
   and let dw(S) denote the boundary of w(S). 
   We then have, just adding and substracting dw(S)
   to the surface of integration S

   q(S) = ( -lambda[] grad T . n g(S) )_{ dw(S) - ( dw(S)-S ) }.
  
   dw(S) being a boundary, Stokes theorem can be invoked and,
   after an integration by part one ends up with

   q(S) = ( -lambda[] grad T . grad g(S) )_w(S)
        - ( -lambda[] grad T . n g(S) )_{dw(S)-S}.
	+ ( Q g(S) )_w(S)

   Now, g(S) is zero on {dw(S)-S}, except maybe at some surface elements 
   adjacents to dS, but not in dS. 
   The second terme vanishes then if either S is closed, 
   or adjacent to a homogeneous Neumann boundary condition. 

   The third term also vanishes, except if a region 
   with a nonzero heatsource Q is in contact with the surface S.

   So we have nearly always the following practical formula
   to evaluate the heat flux across a surface S,
   in terms of a well-chosen auxiliary scalar function g(S).

   q(S) = ( -lambda[] grad T . grad g(S) )_{support of g(S)}


   Particular cases:

   For the heat flux through the boundary of a thermal electrode,
   one uses g(S) = BF(T_electrode).
   Note that this heat flux is equal to Q_electrode 
   in the stationary case.
   This is the case for the heat flux through the boundary of Window2 
   
   Auxiliary functions g(S) are also generated
   for the surfaces named "Surface_i", i=1,2,3 in the model. 
   Note that the flux computed through Surface_3 is incorrect
   because this surface is not adjacent to surfaces
   with homogeneous Neumann boundary conditions.
*/

Include "brick_common.pro";

QWindow1 = 
  DefineNumber[1e3, Name "Window1/Heat source [W]"];

/* The user is given the choice of setting either the global temperature
   or the global heat flux in Window2. 
   Check how the variable "Flag_ConstraintWin2" is used at different 
   places to alter the model according to that choice 
   and to manage the visibility of related input and output data.*/
QWindow2 = 
  DefineNumber[1e3, Name "Window2/Heat source [W]",
	       Visible Flag_ConstraintWin2];
TWindow2 = 
  DefineNumber[50, Name "Window2/Temperature [degC]",
	       Visible !Flag_ConstraintWin2];
outQWindow2 = 
  DefineNumber[0, Name "Output/Q window 2 [degC]",
	       Visible !Flag_ConstraintWin2, Highlight "Ivory"];
outTWindow2 = 
  DefineNumber[0, Name "Output/T window 2 [degC]",
	       Visible Flag_ConstraintWin2, Highlight "Ivory"];


ConvectionCoef = 
  DefineNumber[1000, Name"Surface2/hconv",
	       Label "Convection coefficient [W/(m^2K)]"];
T_Ambiance = 
  DefineNumber[20, Name"Surface2/Ambiance temperature [degC]"];
T_Dirichlet =
  DefineNumber[20, Name"Surface1/Imposed temperature [degC]"];

Group {
  /* Geometrical regions: */
  Brick = Region[100];
  LayerWindow1 = Region[115];
  Window1 = Region[111];
  Window2 = Region[112];
  SurfWindow1 = Region[211];
  SurfWindow2 = Region[212];

  Surface_1 = Region[201];
  Surface_2 = Region[202];
  Surface_3 = Region[203];
  NbSurface = 3;

  /* Abstract regions: 

     Vol_Lambda_Ele     : volume regions with a thermal conductivity lambda[]
     Sur_Dirichlet_The  : Dirichlet bondary condition surface
     Sur_Neu_Ele        : Neumann bondary condition surface
     Sur_Convection_The : convective surface q.n = h ( T-Tinf )
     Vol_Qsource_The    : volume heat source regions
     Tfloating_The      : thermal electrodes
  */

  Vol_Lambda_The = Region[ {Brick, Window1, Window2, LayerWindow1} ];
  Sur_Dirichlet_The = Region[ Surface_1 ];
  Sur_Neumann_The = Region[ Surface_3 ];
  Sur_Convection_The = Region[ Surface_2 ];

  Vol_Qsource_The = Region[ Window1 ];
  Tfloating_The = Region[ Window2 ];

  If( !Flag_Regularization )
    Vol_OneSide_The = Region[ Brick ];
    Sur_Tdisc_The = Region[ SurfWindow1 ];
  Else
    Vol_OneSide_The = Region[ {} ];
    Sur_Tdisc_The = Region[ {} ];
  EndIf
}


Function {
  lambda_Brick = 50.; // steel
  lambda[ Brick ] = lambda_Brick;
  lambda[ Region[ {Window1, Window2} ] ] = lambda_Brick ; // comment

  lambda_Layer = 1.0; // 50 time smaller than surrounding regions
  If( Flag_Regularization )
    lambda[ LayerWindow1 ] = lambda_Layer;
  EndIf
  lambda[ Sur_Tdisc_The ] = lambda_Layer;
  thickness[ Sur_Tdisc_The ] = e_layer;

  h[ Surface_2 ] = ConvectionCoef;
  Tinf[ Surface_2 ] = T_Ambiance;

  If( Flag_QFromFile )
    Qsource[ Window1 ] = ScalarField[XYZ[],0,1]{1};
  Else
    Qsource[ Window1 ] = 0*QWindow1/SurfaceArea[];
  EndIf
}

Constraint {
  { Name Dirichlet_The ;
    Case {
      { Region Sur_Dirichlet_The ; Value T_Dirichlet ; }
    }
  }
  { Name T_Discontinuity ;
    Case {
      { Region SurfWindow1 ; Value 20 ; }
    }
  }
  { Name T_electrode;
    Case {
      If( !Flag_ConstraintWin2 )
        { Region Window2 ; Value TWindow2 ; }
      EndIf
    }
  }
  { Name Q_electrode;
    Case {
      If( Flag_ConstraintWin2 )
        { Region Window2 ; Value QWindow2 ; }
      EndIf
    }
  }
  For i In {1:NbSurface}
  { Name FluxLayer~{i} ;
      Case {
        { Region Surface~{i} ; Value 1. ; }
      }
    }
  EndFor
}

Integration {
  { Name Int ;
    Case {
      {
        Type Gauss ;
        Case {
          { GeoElement Triangle   ; NumberOfPoints  4 ; }
          { GeoElement Quadrangle ; NumberOfPoints  4 ; }
          { GeoElement Line       ; NumberOfPoints  4 ; }
        }
      }
    }
  }
}

Jacobian {
  { Name Vol ;
    Case { 
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name Sur ;
    Case { 
      { Region All ; Jacobian Sur ; }
    }
  }
}

Group {
  Dom_Hgrad_T = Region[ {Vol_Lambda_The, 
			 Sur_Convection_The, 
			 Sur_Tdisc_The} ];
  DomainWithSurf_TL_The = 
    ElementsOf[ {Vol_OneSide_The, Sur_Tdisc_The}, 
		OnOneSideOf Sur_Tdisc_The ];
}

FunctionSpace {
  { Name Hgrad_T; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef Tn ; Function BF_Node ;
        Support Dom_Hgrad_T ; 
	Entity NodesOf[All, Not Tfloating_The] ; }
      { Name sf ; NameOfCoef Tf ; Function BF_GroupOfNodes ;
        Support Dom_Hgrad_T ; Entity GroupsOfNodesOf[ Tfloating_The ] ; }
      { Name sdn ; NameOfCoef Tdn ; Function BF_Node ;
        Support DomainWithSurf_TL_The ; Entity NodesOf[ Sur_Tdisc_The ] ; }
    }
    SubSpace {
      { Name Tcont ; NameOfBasisFunction { sn, sf } ; }
      { Name Tdisc ; NameOfBasisFunction { sdn } ; }
    }
    GlobalQuantity {
      { Name T_electrode ; Type AliasOf        ; NameOfCoef Tf ; }
      { Name Q_electrode ; Type AssociatedWith ; NameOfCoef Tf ; }
    }
    Constraint {
      { NameOfCoef Tn ; EntityType NodesOf ; 
        NameOfConstraint Dirichlet_The ; }
      // { NameOfCoef Tdn ; EntityType NodesOf ; //
      // 	NameOfConstraint T_Discontinuity ; }  //
      { NameOfCoef T_electrode ; EntityType GroupsOfNodesOf ; 
        NameOfConstraint T_electrode ; }
      { NameOfCoef Q_electrode ; EntityType GroupsOfNodesOf ; 
        NameOfConstraint Q_electrode ; }
    }
  }
  For i In {1:NbSurface}
  { Name FluxLayer~{i} ; Type Form0 ;
    BasisFunction {
      { Name gn ; NameOfCoef un ; Function BF_GroupOfNodes;
        Support Dom_Hgrad_T ; Entity GroupsOfNodesOf[ Surface~{i} ] ; }
    }
    Constraint {
      { NameOfCoef un ; EntityType GroupsOfNodesOf ; 
	NameOfConstraint FluxLayer~{i} ; }
    }
  }
  EndFor

}

Formulation {
  { Name Thermal_T ; Type FemEquation ;
    Quantity {
      { Name T ; Type Local ; NameOfSpace Hgrad_T ; }
      { Name Tcont ; Type Local ; NameOfSpace Hgrad_T[Tcont] ; }
      { Name Tdisc ; Type Local ; NameOfSpace Hgrad_T[Tdisc] ; }
      { Name Tglob ; Type Global ; NameOfSpace Hgrad_T[T_electrode] ; }
      { Name Qglob ; Type Global ; NameOfSpace Hgrad_T[Q_electrode] ; }
      For i In {1:NbSurface}
        { Name un~{i} ; Type Local ; NameOfSpace FluxLayer~{i} ; }
      EndFor

    }
    Equation {

      Galerkin { [ lambda[]  * Dof{d T} , {d T} ] ;
                 In Vol_Lambda_The; Integration Int ; Jacobian Vol ; }

      Galerkin { [ ( lambda[]/thickness[] ) * Dof{Tdisc} , {Tdisc} ] ;
	         In Sur_Tdisc_The; Integration Int ; Jacobian Sur ; }

      Galerkin { [ -Qsource[] , {T} ] ; 
                 In Vol_Qsource_The ; Integration Int ; Jacobian Vol ; }

      Galerkin { [ h[] * Dof{T} , {T} ] ;
                 In Sur_Convection_The ; Integration Int ; Jacobian Sur ; }

      Galerkin { [ -h[] * Tinf[] , {T} ] ;
                 In Sur_Convection_The ; Integration Int ; Jacobian Sur ; }

      GlobalTerm { [-Dof{Qglob} , {Tglob} ] ; In Tfloating_The ; }

      For i In {1:NbSurface}
        Galerkin { [ 0 * Dof{un~{i}} , {un~{i}} ] ; 
          In Vol_Lambda_The ; Integration Int ; Jacobian Vol ; }
      EndFor
    }
  }
}


Resolution {
  { Name Thermal_T ;
    System {
      { Name Sys_The ; NameOfFormulation Thermal_T ; }
    }
    Operation {
      If( Flag_QFromFile )
        GmshRead[ "Q.pos", 1];
      EndIf
      DeleteFile["output.txt"];
      Generate Sys_The ; Solve Sys_The ; SaveSolution Sys_The ;
    }
  }
}

PostProcessing {
  { Name Thermal_T ; NameOfFormulation Thermal_T ;
    PostQuantity {
      { Name T ; Value { Term { [ {T} ] ; 
	    In Dom_Hgrad_T ; Jacobian Vol ; } } }
      { Name q ; Value { Term { [ -lambda[] * {d T} ] ; 
	    In Dom_Hgrad_T ; Jacobian Vol ; } } }
      { Name Tcont ; Value { Term { [ {Tcont} ] ; 
	    In Dom_Hgrad_T ; Jacobian Vol ; } } }
      { Name Tdisc ; Value { Term { [ {Tdisc} ] ; 
	    In Dom_Hgrad_T ; Jacobian Vol ; } } }
      { Name Qglob; Value { Term { [ {Qglob} ]; In Tfloating_The; } } }
      { Name Tglob; Value { Term { [ {Tglob} ]; In Tfloating_The; } } }

      For i In {1:NbSurface}
      { Name un~{i} ; Value { Local { [ {un~{i}} ] ; 
	    In Vol_Lambda_The ; Jacobian Vol ; } } }
      { Name IFlux~{i} ; Value { Integral { [ -lambda[]*{d T} * {d un~{i}} ];
	    In Vol_Lambda_The ; Jacobian Vol ; Integration Int ; } } }
      EndFor

    }
  }
}

PostOperation Map_T UsingPost Thermal_T {
  If( !Flag_Regularization )
    Print[ Tcont, OnElementsOf Vol_Lambda_The, File "Tcont.pos"] ;
    Print[ Tdisc, OnElementsOf Vol_Lambda_The, File "Tdisc.pos"] ;
  EndIf

  If(Flag_ConstraintWin2) 
    Print[ Tglob, OnRegion Tfloating_The, File > "output.txt", Color "Ivory",
	   Format Table, SendToServer "Output/T window 2 [degC]" ];
  Else
    Print[ Qglob, OnRegion Tfloating_The, File > "output.txt", Color "Ivory",
	   Format Table, SendToServer "Output/Q window 2 [degC]" ];
  EndIf

  For i In {1:NbSurface}
    Print[un~{i}, OnElementsOf Vol_Lambda_The, 
    	  File Sprintf("FluxLayer_%g.pos",i)];
    Print[ IFlux~{i}[Vol_Lambda_The], OnGlobal, 
	   Format TimeTable, File > "Fluxes.dat", Color "Ivory",
	   SendToServer Sprintf("Output/Heat flux surface %g [W]", i)];
  EndFor
  Print[ T, OnElementsOf Vol_Lambda_The, File "T.pos" ] ;
  Echo[ StrCat["l=PostProcessing.NbViews-1;",
	       "View[l].IntervalsType = 3;",
	       "View[l].NbIso = 30;",
	       "View[l].NormalRaise = 0.0005;"],
	File "tmp.geo", LastTimeStepOnly] ;
  Print [ q, OnLine {{0.02,0.0001,0}{0.02,0.05,0}} {200}, 
	  Format SimpleTable, File "Cut.txt" ];
}


