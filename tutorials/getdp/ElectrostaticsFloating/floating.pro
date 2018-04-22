/* -------------------------------------------------------------------
   Tutorial 4 : floating potential of a microstrip electrode

   Features:
   - Global quantities and their special shape functions
   - Computation of the energy dual, i.e. of the armature charge of the electrode
   - More on ONELAB parameters (flags, model options, check boxes, menus, ...)

   To compute the solution interactively from the Gmsh GUI:
       File > Open > floating.pro
       Run (button at the bottom of the left panel)

   ------------------------------------------------------------------- */

/*
  A thing GetDP is pretty good at is the management of global (non-local) basis
  functions. Finite element expansions typically associate basis functions to
  individual nodes or edges in the mesh. But consider the situation
  where a scalar field is set to be uniform over a region of the problem 
  (Say a floating potential electrode in an Electrostatics problem, 
  to fix the idea). By factorizing the identical nodal value "v_electrode",
  a global (non-local) basis function "BF_electrode" is obtained as factor 
  which is the sum of the shape functions of all the nodes in the electrode
  region. This basis function "BF_electrode" 
  - is a continuous function, scalar in this case,
  - is equal to 1 at the nodes of the electrode region, and to 0 
    at all other nodes
  - decreases from 1 to 0 over the one element thick layer of outside 
    finite elements immediately in contact with the electrode region.
  One such glabal basis function can be associated with each electrode 
  in the system, so that the finite element expansion of the electric 
  scalar potential reads:

   v = Sum_k sn_k vn_k + Sum_electrode v_electrode BF_electrode

   We show in this tutorial how GetDP takes advantage of global quantities
   and the associated global basis functions
   - to reduce the number of unknowns
   - to compute efficiently the electrode charges "Q_electrode", 
     which are precisely the energy duals of the global "v_electrode" quantities
   - to deal with floating potentials, which are the computed electrode 
     potential when the electrode charge is imposed
   - to provide output quantities (charges, armature voltages, 
     capacitances, ...) that can be immediately used in a external circuit.
 */


Group { 
  /* Geometrical regions: */ 

  Air = Region[101]; 
  Diel1 = Region[111];

  Ground = Region[120];
  Microstrip = Region[121];
  SurfInf = Region[130];

  /* Abstract regions:
     Vol_Dielectric_Ele : dielectric volume regions where 
                          "Div ( epsr[] Grad v)" is solved
     Sur_Neu_Ele        : Neumann bondary condition ( epsr[] n.Grad v = 0 )
     Electrodes_Ele     : electrode regions
     
     No prefix (Vol_ or Sur_) for the region "Electrodes_Ele", 
     which may contain both surface or volume regions.
     There are two electrodes in this model: Ground and Microstrip
  */
 
  Vol_Dielectric_Ele = Region[ {Air, Diel1} ];
  Sur_Neu_Ele = Region[ {SurfInf} ];
  Electrodes_Ele = Region [ {Ground, Microstrip} ]; 
}

/* A number of ONELAB parameters are defined to define model parameters 
   or model options interactively. */
MicrostripTypeBC = 
  DefineNumber[0, Name "1Microstrip excitation/Type", 
	       Choices{ 0="Fixed voltage", 1="Fixed charge"} ] ;
MicrostripValueBC = 
  DefineNumber[1e-3, Name "1Microstrip excitation/Value"] ;
EpsilonRelDiel = 
  DefineNumber[9.8, Name "2Dielectric/Relative permittivity"] ;
DisplayGlobalBF = 
  DefineNumber[0, Name "3Options/Display global basis functions", Choices {0,1} ] ;
OverwriteOutput = 
  DefineNumber[1, Name "3Options/Overwrite output.txt file", Choices {0,1} ] ;

Function {
  epsr[Air] = 1.;
  epsr[Diel1] = EpsilonRelDiel;
  eps0 = 8.854187818e-12;  // permittivity of empty space
}


Constraint {
  /* Dirichlet boundary condition is no longer used. 
     The microstrip and the ground are now treated as electrodes,
     whose voltage is imposed with the "SetGlobalPotential" constraint below. */
    { Name Dirichlet_Ele; Type Assign;
      Case {}
    }

  { Name SetGlobalPotential; Type Assign;
    Case {
      /* Define the imposed potential regionwise on the different parts of
	 "Electrodes_Ele". No voltage imposed to the Microstrip electrode
	 when the "Fixed charge" option is enabled 
	 ( MicrostripTypeBC = true ). */
      { Region Ground; Value 0; }
      If( !MicrostripTypeBC )
	{ Region Microstrip; Value MicrostripValueBC; }
      EndIf
    }
  }
  { Name SetArmatureCharge; Type Assign;
    Case {
      If( MicrostripTypeBC )
	{ Region Microstrip; Value MicrostripValueBC; }
      EndIf
    }
  }
}

Group{
  /* The domain of definition lists all regions 
     on which the field "v" is defined.*/
  Dom_Hgrad_v_Ele =  Region[ {Vol_Dielectric_Ele,
			      Sur_Neu_Ele,
			      Electrodes_Ele } ];
}

FunctionSpace {
  /* The magic in the treatment of global quantitities by GetDP is in the fact 
     that nearly all the work is done at the level of the FunctionSpace
     definition. The finite element expansion of "v" is 

     v = Sum_k sn_k vn_k + Sum_electrode v_electrode BF_electrode

     with the the sum_k running over all nodes except those of the electrode
     regions. This is exactly what one finds in the FunctionSpace definotion 
     below with "sf" standing for "BF_electrode" and "vf" for "v_electrode".
     The global quantities are also be attributed a more explicit 
     and meaningful name. Moreover the "AssociatedWith" statement manifests 
     the fact that the global potential of an electrode is the (electrostatic)
     energy dual of the electric charge carried by that electrode. 
     By checking the "Display global basis functions" checkbox and running 
     the model, you can take a look on how the two "BF_electrode" basis 
     functions in this model look like. 
     Constraints can then be set on either component of the FunctionSpace.
     Besides the usual Dirichlet boundary condition conditions, which is 
     left here for the sake of completeness but is not used in this model, 
     there is the possibility to fix either the GlobalPotential or the
     ArmatureCharge of each indidual electrode (not both, of course). 
     When the ArmatureCharge is fixed, the computed GlobalPotential computed 
     for that electrode is the so-called floating potential. 
  */

  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Dom_Hgrad_v_Ele; Entity NodesOf[ All, Not Electrodes_Ele ]; }
      { Name sf; NameOfCoef vf; Function BF_GroupOfNodes; 
        Support Dom_Hgrad_v_Ele; Entity GroupsOfNodesOf[ Electrodes_Ele ]; }
    }
    GlobalQuantity {
      { Name GlobalPotential; Type AliasOf       ; NameOfCoef vf; }
      { Name ArmatureCharge ; Type AssociatedWith; NameOfCoef vf; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; 
        NameOfConstraint Dirichlet_Ele; } // unused in this model, left for completeness
      { NameOfCoef GlobalPotential; EntityType GroupsOfNodesOf;
	NameOfConstraint SetGlobalPotential; }
      { NameOfCoef ArmatureCharge; EntityType GroupsOfNodesOf; 
	NameOfConstraint SetArmatureCharge; }
    }
    // Subspace definition only needed to display BF_electrode in PostProcessing
    SubSpace { 
      { Name vf; NameOfBasisFunction sf; }
    }
  }
}

Jacobian {
  { Name Vol ;
    Case { 
      { Region All ; Jacobian Vol ; }
    }
  }
}

Integration {
  { Name Int ;
    Case { {Type Gauss ;
            Case { { GeoElement Triangle    ; NumberOfPoints  4 ; }
                   { GeoElement Quadrangle  ; NumberOfPoints  4 ; } }
      }
    }
  }
}

Formulation {
  /* Only minor changes in the formulation.
     The global quantities are declared in the "Quantity{}" section,
     and a "GlobalTerm" is added that triggers the assembly 
     of an additional equation per electrode in the system 
     to compute the charge Q_electrode accordint to:

     Q_electrode = (-epsr[] Grad v, Grad BF_electrode)_Vol_Dielectric_Ele
   */
  { Name Electrostatics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
      { Name U; Type Global; NameOfSpace Hgrad_v_Ele [GlobalPotential]; }
      { Name Q; Type Global; NameOfSpace Hgrad_v_Ele [ArmatureCharge]; }
      // next line only needed to display the BF_electrode in PostProcessing
      { Name vf; Type Local; NameOfSpace Hgrad_v_Ele [vf]; }
    }
    Equation {
      Galerkin { [ epsr[] * Dof{d v} , {d v} ]; In Vol_Dielectric_Ele; 
	Jacobian Vol; Integration Int; }
      GlobalTerm { [ -Dof{Q}/eps0 , {U} ]; In Electrodes_Ele; }
    }
  }
}

Resolution {
  { Name EleSta_v;
    System {
      { Name Sys_Ele; NameOfFormulation Electrostatics_v; }
    }
    Operation { 
      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
      If( OverwriteOutput )
	DeleteFile[ "output.txt" ];
      EndIf
    }
  }
}

PostProcessing {
 { Name EleSta_v; NameOfFormulation Electrostatics_v;
    Quantity {
      { Name v; Value { Term { [ {v} ]; In Dom_Hgrad_v_Ele; Jacobian Vol; } } }
      { Name e; Value { Term { [ -{d v} ]; In Dom_Hgrad_v_Ele; Jacobian Vol; } } }
      { Name d; Value { Term { [ -eps0*epsr[] * {d v} ]; 
	    In Dom_Hgrad_v_Ele; Jacobian Vol; } } }
      { Name Q; Value { Term { [ {Q} ]; In Electrodes_Ele; } } }
      { Name U; Value { Term { [ {U} ]; In Electrodes_Ele; } } }
      { Name C; Value { Term { [ {Q}/{U} ]; 
	    In Electrodes_Ele; } } }
      { Name energy;
	Value { Integral { Type Global;
	    [ eps0*epsr[] / 2. * SquNorm[{d v}] ];
	    In Vol_Dielectric_Ele; Jacobian Vol; Integration Int;
	  }
	}
      }
      // next lines only needed to display global BF in PostProcessing
      { Name BFGround; Value { Term { [ BF {vf} ]; In Dom_Hgrad_v_Ele; 
	    SubRegion Ground; Jacobian Vol; } } }
      { Name BFMicrostrip; Value { Term { [ BF {vf} ]; In Dom_Hgrad_v_Ele; 
	    SubRegion Microstrip; Jacobian Vol; } } }

    }
  }
}

/* Various output results are generated, which are both displayed 
   in the graphical user interface, and stored in disk files. 
   In particular, global quantities related results are stored 
   in the "output.txt" file. A user option allows to chose 
   to not overwrite the "output.txt" file when running a new simulation. */

PostOperation {
  { Name Map; NameOfPostProcessing EleSta_v;
     Operation {
       If( DisplayGlobalBF )
	 Print[ BFGround, OnElementsOf Dom_Hgrad_v_Ele, File "BFGround.pos" ];
         Echo[ StrCat["l=PostProcessing.NbViews-1;", 
		      "View[l].IntervalsType = 1;",
		      "View[l].NbIso = 40;",
		      "View[l].ShowElement = 1;"],
	       File "BFGround.opt", LastTimeStepOnly] ;

	 Print[ BFMicrostrip, OnElementsOf Dom_Hgrad_v_Ele, 
		File "BFMicrostrip.pos" ];
	 Echo[ StrCat["l=PostProcessing.NbViews-1;", 
		      "View[l].IntervalsType = 1;",
		      "View[l].NbIso = 40;",
		      "View[l].ShowElement = 1;"],
	       File "BFMicrostrip.opt", LastTimeStepOnly] ;
       EndIf

       Print[ v, OnElementsOf Dom_Hgrad_v_Ele, File "v.pos" ];
       Echo[ StrCat["l=PostProcessing.NbViews-1;", 
                    "View[l].IntervalsType = 3;",
                    "View[l].NbIso = 40;"],
             File "v.opt", LastTimeStepOnly] ;

       Echo[ "Microstrip charge [C]:", Format Table, File > "output.txt"] ;
       Print[ Q, OnRegion Microstrip, File > "output.txt", Color "AliceBlue",
	      Format Table, SendToServer "Output/Microstrip/Charge [C]" ];
       Echo[ "Microstrip potential [V]:", Format Table, File > "output.txt"] ;
       Print[ U, OnRegion Microstrip, File > "output.txt", Color "AliceBlue",
	      Format Table, SendToServer "Output/Microstrip/Potential [V]" ];
       Echo[ "Ground charge [C]:", Format Table, File > "output.txt"] ;
       Print[ Q, OnRegion Ground, File > "output.txt", Color "AliceBlue",
	      Format Table, SendToServer "Output/Ground/Charge [C]" ];
       Echo[ "Microstrip capacitance [F]:", Format Table, File > "output.txt"] ;
       Print[ C, OnRegion Microstrip, File > "output.txt", Color "AliceBlue",
	      Format Table, SendToServer "Output/Global/Capacitance [F]" ];
       Echo[ "Electrostatic energy [J]:", Format Table, File > "output.txt"] ;
       Print[ energy, OnRegion Vol_Dielectric_Ele, File > "output.txt", 
	      Color "AliceBlue",
	      Format Table, SendToServer "Output/Global/Energy [J]" ];
     }
  }
}


