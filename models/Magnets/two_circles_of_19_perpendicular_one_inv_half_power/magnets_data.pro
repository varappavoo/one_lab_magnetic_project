DefineConstant[
  NumMagnets = {38, Min 1, Max 55, Step 1, Name "Parameters/0Number of magnets"}
];
mm = 1.e-3;

//x = -1;
y = 0;
z = 0;

//radius = 0.0119281497082362427277580294513;
//radius = 0.0117;//19281497082362427277580294513;

// MY MAGNET DIAMETER AND HEIGHT ARE 9*3MM
radius = 0.0045; 
height = 0.003;

// *1
x = {-0.10280232,  0.10280232, -0.28086117,  0.28086117, -0.10280232,
         0.10280232, -0.38366349,  0.38366349, -0.20560465,  0.        ,
         0.20560465, -0.38366349,  0.38366349, -0.10280232,  0.10280232,
        -0.28086117,  0.28086117, -0.10280232,  0.10280232};

z = {-0.38366349, -0.38366349, -0.28086117, -0.28086117, -0.17805885,
        -0.17805885, -0.10280232, -0.10280232,  0.        ,  0.        ,
         0.        ,  0.10280232,  0.10280232,  0.17805885,  0.17805885,
         0.28086117,  0.28086117,  0.38366349,  0.38366349};

count = -1;
For i In {1:NumMagnets}
	//z=z+1;

	// If((i)%6 == 0)
	// 	x=x+1; 
	// 	z=0;
	// EndIf
	count = count +1;
	//If(count == 10)
	//	height  = 10 * mm;
	//EndIf

	If(count < 19)
		DefineConstant[
		      HC~{i} = {800e3,
        Name Sprintf("Parameters/Magnet %g/0Coercive magnetic field [Am^-1]", i)},

			X~{i} = {x(count)*45*mm, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0X position [m]", i) },
			Y~{i} = {0, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0Y position [m]", i) },
			Z~{i} = {z(count)*45*mm, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0Z position [m]", i) },

			M~{i} = {0, Choices{0="Cylinder",1="Cube"},
			  Name Sprintf("Parameters/Magnet %g/00Shape", i)},

			R~{i} = {radius, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1Radius [m]", i),
			  Visible (M~{i} == 0) },
			L~{i} = {height, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1Length [m]", i),
			  Visible (M~{i} == 0) },

			Lx~{i} = {50*mm, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1X length [m]", i),
			  Visible (M~{i} == 1) },
			Ly~{i} = {50*mm, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1Y length [m]", i),
			  Visible (M~{i} == 1) },
			Lz~{i} = {50*mm, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1Z length [m]", i),
			  Visible (M~{i} == 1) },

			//Rx~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			Rx~{i} = {Pi, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2X rotation [rad]", i) },
			Ry~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2Y rotation [rad]", i) },
			Rz~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2Z rotation [rad]", i) }
		];
	EndIf

	shift_y = 50*mm;
	shift_z = 50*mm;
	If(count > 18)
		DefineConstant[
		      HC~{i} = {400e3,
        Name Sprintf("Parameters/Magnet %g/0Coercive magnetic field [Am^-1]", i)},
			X~{i} = {x(count-19)*45*mm, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0X position [m]", i) },
			Y~{i} = {z(count-19)*45*mm + shift_y, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0Y position [m]", i) },
			Z~{i} = {shift_z, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0Z position [m]", i) },

			M~{i} = {0, Choices{0="Cylinder",1="Cube"},
			  Name Sprintf("Parameters/Magnet %g/00Shape", i)},

			R~{i} = {radius, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1Radius [m]", i),
			  Visible (M~{i} == 0) },
			L~{i} = {height, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1Length [m]", i),
			  Visible (M~{i} == 0) },

			Lx~{i} = {50*mm, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1X length [m]", i),
			  Visible (M~{i} == 1) },
			Ly~{i} = {50*mm, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1Y length [m]", i),
			  Visible (M~{i} == 1) },
			Lz~{i} = {50*mm, Min mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/1Z length [m]", i),
			  Visible (M~{i} == 1) },

			Rx~{i} = {-Pi/2, Min -Pi, Max Pi, Step Pi/180,
			//Rx~{i} = {-Pi/2, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2X rotation [rad]", i) },
			Ry~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2Y rotation [rad]", i) },
			Rz~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2Z rotation [rad]", i) }
		];
	EndIf
	
EndFor
