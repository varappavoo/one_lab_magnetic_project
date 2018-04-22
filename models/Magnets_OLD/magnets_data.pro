DefineConstant[
  NumMagnets = {24, Min 1, Max 50, Step 1, Name "Parameters/0Number of magnets"}
];
mm = 1.e-3;

x = -1;
y = 0;
z = -1;

count = 0;
For i In {1:NumMagnets}
	z=z+1;

	If((i-1)%5 == 0)
		x=x+1; 
		z=0;
	EndIf
	count = count +1;
	If(count != 13)
		DefineConstant[
			X~{i} = {x*45*mm, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0X position [m]", i) },
			Y~{i} = {0, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0Y position [m]", i) },
			Z~{i} = {z*45*mm, Min -100*mm, Max 100*mm, Step mm,
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
	EndIf
EndFor
