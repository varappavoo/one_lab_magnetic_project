DefineConstant[
  NumMagnets = {24, Min 1, Max 50, Step 1, Name "Parameters/0Number of magnets"}
];
mm = 1.e-3;

//x = -1;
y = 0;
z = 0;

x = {3,4,2,3,4,5,1,2,3,4,5,6,1,2,3,4,5,6,2,3,4,5,3,4};
z = {1,1,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,6,6};

count = -1;
For i In {1:NumMagnets}
	//z=z+1;

	// If((i)%6 == 0)
	// 	x=x+1; 
	// 	z=0;
	// EndIf
	count = count +1;


	DefineConstant[
		X~{i} = {x(count)*45*mm, Min -100*mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/0X position [m]", i) },
		Y~{i} = {0, Min -100*mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/0Y position [m]", i) },
		Z~{i} = {z(count)*45*mm, Min -100*mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/0Z position [m]", i) },

		M~{i} = {0, Choices{0="Cylinder",1="Cube"},
		  Name Sprintf("Parameters/Magnet %g/00Shape", i)},

		R~{i} = {20*mm, Min mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/1Radius [m]", i),
		  Visible (M~{i} == 0) },
		L~{i} = {10*mm, Min mm, Max 100*mm, Step mm,
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

		Rx~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
		  Name Sprintf("Parameters/Magnet %g/2X rotation [rad]", i) },
		Ry~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
		  Name Sprintf("Parameters/Magnet %g/2Y rotation [rad]", i) },
		Rz~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
		  Name Sprintf("Parameters/Magnet %g/2Z rotation [rad]", i) }
	];
	
EndFor
