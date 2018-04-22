DefineConstant[
  NumMagnets = {112, Min 1, Max 112, Step 1, Name "Parameters/0Number of magnets"}
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
x = {-7.40861540e-02,  1.64158893e-01, -3.06895332e-01,
         3.90359230e-01, -5.17186905e-01,  5.87918027e-01,
        -2.66311085e-02,  2.11613938e-01, -2.59440287e-01,
        -7.02547593e-01,  4.09172735e-01, -4.69731860e-01,
         7.42339949e-01, -4.36624082e-02,  1.94582638e-01,
        -2.73793828e-01, -8.19315292e-01,  5.63594657e-01,
        -5.86499559e-01,  8.42294708e-01,  3.49004560e-01,
         1.02719085e-01, -1.34619045e-01, -3.90561527e-01,
         6.53353766e-01, -8.75967975e-01, -6.43043018e-01,
         8.80448399e-01,  4.38763668e-01,  2.02625622e-01,
        -1.82356266e-01, -4.20295651e-01, -8.68348905e-01,
         6.21234135e-01,  2.49853854e-02, -6.36388529e-01,
         8.54001599e-01,  3.45662912e-01, -2.88535051e-01,
        -7.93316819e-01,  5.32127305e-01,  1.67385771e-01,
        -5.04627930e-01, -8.01001155e-02,  7.64894768e-01,
         3.53797474e-01, -3.24327454e-01, -6.61556219e-01,
         6.19665878e-01,  1.31446021e-01, -1.05897163e-01,
        -4.81255744e-01,  4.28970701e-01, -2.65644445e-01,
         2.06800972e-01, -3.05422122e-02};
z = {-8.77596902e-01, -8.65284311e-01, -8.25518223e-01,
        -7.89483852e-01, -7.12869404e-01, -6.55757176e-01,
        -6.43801435e-01, -6.31488843e-01, -5.91722755e-01,
        -5.31123301e-01, -4.97762168e-01, -4.79073936e-01,
        -4.73916113e-01, -4.05847159e-01, -3.93534568e-01,
        -3.42983562e-01, -3.23090596e-01, -3.15921104e-01,
        -2.71041230e-01, -2.57302752e-01, -2.11693505e-01,
        -1.73367837e-01, -1.49224221e-01, -1.34950856e-01,
        -9.48880644e-02, -9.13520022e-02, -3.70961689e-02,
        -2.18105055e-02,  9.33953490e-03,  4.32677702e-02,
         8.45137943e-02,  1.01751877e-01,  1.47089294e-01,
         1.63017221e-01,  2.02504225e-01,  2.02826949e-01,
         2.15282027e-01,  2.35255894e-01,  3.00627333e-01,
         3.82509486e-01,  3.84314024e-01,  3.93905459e-01,
         4.01702404e-01,  4.16675511e-01,  4.36578830e-01,
         5.42777817e-01,  5.57920331e-01,  5.81384941e-01,
         6.25842856e-01,  6.29745737e-01,  6.53839620e-01,
         7.37602868e-01,  7.69187376e-01,  8.39701203e-01,
         8.56094878e-01,  8.80188761e-01};

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

	If(count < 56)
		DefineConstant[
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

			Rx~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2X rotation [rad]", i) },
			Ry~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2Y rotation [rad]", i) },
			Rz~{i} = {Pi, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2Z rotation [rad]", i) }
		];
	EndIf

	shift_y = 300*mm;
	//shift_z = 300*mm;
	If(count > 55)
		DefineConstant[
			X~{i} = {x(count-56)*45*mm, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0X position [m]", i) },
			Y~{i} = {shift_y, Min -100*mm, Max 100*mm, Step mm,
			  Name Sprintf("Parameters/Magnet %g/0Y position [m]", i) },
			Z~{i} = {z(count-56)*45*mm, Min -100*mm, Max 100*mm, Step mm,
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

			Rx~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			// Rx~{i} = {-Pi/2, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2X rotation [rad]", i) },
			Ry~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2Y rotation [rad]", i) },
			Rz~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
			  Name Sprintf("Parameters/Magnet %g/2Z rotation [rad]", i) }
		];
	EndIf
	
EndFor
