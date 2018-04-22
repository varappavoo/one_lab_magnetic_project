/* -------------------------------------------------------------------
   Tutorial 1 : electrostatic field of a microstrip

   Features:
   - Physical and abstract regions
   - Scalar FunctionSpace with Dirichlet constraint
   - Galerkin term for stiffness matrix

   To compute the solution in a terminal:
       getdp microstrip -solve EleSta_v
       getdp microstrip -pos Map
       getdp microstrip -pos Cut

   To compute the solution interactively from the Gmsh GUI:
       File > Open > microstrip.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

Group {
  /* One starts by giving explicit meaningful names to
     the Physical regions defined in the "microstrip.msh" mesh file.
     There are 2 volume regions and 3 surface regions in this model. */

  Air = Region[101];
  Diel1 = Region[111];

  Ground = Region[120];
  Electrode = Region[121];
  SurfInf = Region[130];

  /* We now define abstract regions to be used below
     in the definition of the scalar electric potential formulation:

     Vol_Dielectric_Ele : dielectric volume regions where
                          "Div ( epsr[] Grad v)" is solved
     Sur_Dir_Ele        : Dirichlet boundary condition (v imposed)
     Sur_Neu_Ele        : Neumann bondary condition ( epsr[] n.Grad v = 0 )

     Vol_xxx groups contain only volume elements of the mesh (triangles here).
     Sur_xxx groups contain only surface elements of the mesh (lines here).
  */

  Vol_Dielectric_Ele = Region[ {Air, Diel1} ];
  Sur_Dir_Ele = Region[ {Ground, Electrode} ];
  Sur_Neu_Ele = Region[ {SurfInf} ];
}

Function {
  /* Material laws (here the relative permittivity)
     are defined piecewise in terms of the above defined physical regions */

  epsr[Air] = 1.;
  epsr[Diel1] = 9.8;
}

Constraint {
  /* As for material laws, the Dirichlet boundary condition
     is defined piecewise.
     The constraint "Dirichlet_Ele" is invoked in the FunctionSpace below */

  { Name Dirichlet_Ele; Type Assign;
    Case {
      { Region Ground; Value 0.; }
      { Region Electrode; Value 1.e-3; }
    }
  }
}

Group{
  /* The domain of definition of a FunctionSpace lists all regions
     on which a field is defined.

     Contrary to Vol_xxx and Sur_xxx regions,
     which contain only volume or surface regions, resp.,
     domains of definition Dom_xxx regions may contain
     both volume and surface regions.
     Hence the use of the prefixes Vol_, Sur_ and Dom_ to avoid confusions.*/

  Dom_Hgrad_v_Ele =  Region[ {Vol_Dielectric_Ele, Sur_Dir_Ele, Sur_Neu_Ele} ];
}

FunctionSpace {
  /* The function space in which we shall pick
     the electric scalar potential "v" solution
     is definied by
     - a domain of definition ("Dom_Hgrad_v_Ele")
     - a type ("Form0" means scalar field)
     - a set of scalar basis functions ("BF_Node" means nodal basis functions)
     - a constraint (here the Dirichlet boundary conditions)

     The FE expansion of the unknown field "v" reads

     v = Sum_k vn_k sn_k

     where the "vn_k" are the nodal values (connectors)
     and "sn_k" the basis functions.
     Not all connectors are unknowns of the FE problem,
     due to the "Constraint", which assigns particular values
     to the nodes of the region "Sur_Dir_Ele".
     GetDP deals with that automatically
     on basis of the definition of the FunctionSpace.
  */

  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Dom_Hgrad_v_Ele; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf;
        NameOfConstraint Dirichlet_Ele; }
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
  /* The syntax of the Formulation{} section is a hard nut to crack.
     So let's deal with it carefully.

     The weak form of the electrostatic problem is particularly simple
     in this model, as it has only one term:

     (epsr[] Grad v, Grad vp)_Vol_Dielectric_Ele = 0
     for all vp in S0(Vol_Dielectric_Ele).

     The corresponding Euler-Lagrange equations are:
     * Div ( epsr[] Grad v) = 0  on Vol_Dielectric_Ele
     * epsr[] n.Grad v = 0       on Sur_Neu_Ele

     The Galerkin{} statement is a symbolic representation
     of this weak formulation term.
     It has got 4 semicolon separated arguments:
     * the density [,] to be integrated,
     * the domain of integration,
     * the jacobian of the transformation reference element -> real element,
     * the integration method to be used.
     The symbol "d" represents the exterior derivative,
     and it is a synonym of "Grad" when applied to a scalar function.
     The expression "d v" stands thus for the gradient of the v field,
     i.e., for the electric field up to a minus sign.

     What is a bit confusing is that the two comma-separated terms
     of the bracket [,], in the first argument,
     are not interpreted the same way. Let us unravel this in detail.

     As the Galerkin method uses as trial (test) functions
     the basis functions "sn_k" of the unknown field "v",
     the density should be something like this

     [ epsr[] * {d v} , basis_functions_of {d v} ].

     Since the second argument is devoted to the trial functions,
     the operator "basis_functions_of" would always be there.
     It can therefore be made implicit, and,
     according to the GetDP syntax, it is omitted.
     So, one writes simply "{d v}".

     The first argument, on the other hand, can contain
     a much wider variety of expressions than the second one.
     In this problem, it contains

     epsr[] * {d v} = Sum_k   vn_k  epsr[]*{d sn_k}

     which is the eletric displacement vector, up to a sign.
     Here, we have two valid syntaxes, with very different meanings.

     The first argument can be expressed in terms of
     the FE expansion of "v" at the present system solution.
     This is indicated by invoking the Dof{} operator:

     [ epsr[] * Dof{d v} , {d v} ].

     Another option, which would not work here,
     is to evaluate the first argument
     with the last available already computed solution.
     For this the Dof{} operator is omitted:

     [ epsr[] * {d v} , {d v} ].

     Both choices are commonly used in different contexts,
     and we shall come back on this often in subsequent tutorials.

     To express the stiffness matrix of the electrostatic problem at hand,
     we have to take the first option.
     Hence the final expression of the density below.

   */
  { Name Electrostatics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
    }
    Equation {
      Galerkin { [ epsr[] * Dof{d v} , {d v} ];
	In Vol_Dielectric_Ele;
	Jacobian Vol; Integration Int; }
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
    }
  }
}


eps0 = 8.854187818e-12;  // permittivity of empty space

/* Post-processing is done in two parts.
   The first part defines, in terms of the Formulation,
   which itself refers to the FunctionSpace,
   a number of quantities that can be evaluated at the postprocessing stage.
   The three quantities defined here are :
   - the scalar vector potential,
   - the electric field,
   - the electric displacement.
   The second part consists in defining groups of post-processing operations,
   which can be invoked separately.
   The first group is invoked by default when Gmsh is run interactively.
   Each Operation specifies
   - a quantity to be eveluated,
   - the region on which the evaluation is done,
   - the name of the output file.
   The generated post-processing files are automatically displayed by Gmsh
   if the "Merge result automatically" option is enabled
   (which is the default). */

PostProcessing {
  { Name EleSta_v; NameOfFormulation Electrostatics_v;
    Quantity {
      { Name v;
        Value {
          Local { [ {v} ];
	    In Dom_Hgrad_v_Ele; Jacobian Vol; }
        }
      }
      { Name e;
        Value {
          Local { [ -{d v} ];
	    In Dom_Hgrad_v_Ele; Jacobian Vol; }
        }
      }
      { Name d;
        Value {
          Local { [ -eps0*epsr[] * {d v} ];
	    In Dom_Hgrad_v_Ele; Jacobian Vol; }
        }
      }
    }
  }
}

e = 1.e-7; // tolerance to ensure that the cut is inside the simulation domain
h = 2.e-3; // vertical position of the cut

PostOperation {
  { Name Map; NameOfPostProcessing EleSta_v;
     Operation {
       Print [ v, OnElementsOf Dom_Hgrad_v_Ele, File "mStrip_v.pos" ];
       Print [ e, OnElementsOf Dom_Hgrad_v_Ele, File "mStrip_e.pos" ];
       Echo[ StrCat["l=PostProcessing.NbViews-1;",
		    "View[l].IntervalsType = 3;",
		    "View[l].NbIso = 40;"],
	     File "tmp.geo", LastTimeStepOnly] ;
       Print [ e, OnLine {{e,h,0}{14.e-3,h,0}}{60}, File "Cut_e.pos" ];
     }
  }
  { Name Cut; NameOfPostProcessing EleSta_v;
    // same cut as above, with more points and exported in raw text format
    Operation {
      Print [ e, OnLine {{e,e,0}{14.e-3,e,0}} {500}, Format TimeTable, File "Cut_e.txt" ];
    }
  }
}
