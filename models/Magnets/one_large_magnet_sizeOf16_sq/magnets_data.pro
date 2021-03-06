DefineConstant[
  NumMagnets = {1, Min 1, Max 50, Step 1, Name "Parameters/0Number of magnets"}
];
mm = 1.e-3;

y = 0;


// MY MAGNET DIAMETER AND HEIGHT ARE 9*3MM
//radius = 0.0045; 
//height = 0.003;

// my square magnet dimensions
length_x = 0.034766762481839936;
length_y = 3*mm;
length_z = 0.034766762481839936;
spacing = 1*mm;

x = {0, 0};

z = {0, 0};


count = -1;
For i In {1:NumMagnets}

	count = count +1;


	DefineConstant[
		X~{i} = {x(count), Min -100*mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/0X position [m]", i) },
		Y~{i} = {0, Min -100*mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/0Y position [m]", i) },
		Z~{i} = {z(count), Min -100*mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/0Z position [m]", i) },

		M~{i} = {1, Choices{0="Cylinder",1="Cube"},
		  Name Sprintf("Parameters/Magnet %g/00Shape", i)},

/*
		R~{i} = {radius, Min mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/1Radius [m]", i),
		  Visible (M~{i} == 0) },
		L~{i} = {height, Min mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/1Length [m]", i),
		  Visible (M~{i} == 0) },
*/
		Lx~{i} = {length_x, Min mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/1X length [m]", i),
		  Visible (M~{i} == 1) },
		Ly~{i} = {length_y, Min mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/1Y length [m]", i),
		  Visible (M~{i} == 1) },
		Lz~{i} = {length_z, Min mm, Max 100*mm, Step mm,
		  Name Sprintf("Parameters/Magnet %g/1Z length [m]", i),
		  Visible (M~{i} == 1) },

		Rx~{i} = {Pi, Min -Pi, Max Pi, Step Pi/180,
		  Name Sprintf("Parameters/Magnet %g/2X rotation [rad]", i) },
		Ry~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
		  Name Sprintf("Parameters/Magnet %g/2Y rotation [rad]", i) },
		Rz~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
		  Name Sprintf("Parameters/Magnet %g/2Z rotation [rad]", i) }
	];
	
EndFor
