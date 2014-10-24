(* MSMModels Package *)

(* Messages *)
LoadDynamicalModel::usage = "LoadDynamicalModel['name'] return an object that contains all vital information about a dynamical model.";

(* 2D MODEL : Diffusion in a 3-well potential *)

LoadDynamicalModel["3Well2D"]:=
	Module[
		{
			potential2D,
			listOfEnergies2D,
			tMFull2D,
			eigenV2D,
			leftEV2D,
			rightEV2D,
			eqFull2D,
			plotBounds2D,
			rangeX2D,
			rangeY2D,
			listOfStatePositions2D,
			nearestFnc2D,
			minimumStates2D,
			pathToResources
		},
		pathToResources = "/Users/jan-hendrikprinz/Library/Mathematica/Applications/DynamicalModels/";
		potential2D=Import[FileNameJoin[{pathToResources,"resources","01_surface.energy"}],"Table"];
		listOfEnergies2D=Flatten[potential2D];
		tMFull2D = TransitionMatrix[ReadSparseMatrix[FileNameJoin[{pathToResources,"resources","02.t"}]]];
		{eigenV2D,leftEV2D,rightEV2D} = FullEigensystem[tMFull2D];
(*		clPCCA2D=PCCASoftClusterMatrix[tMFull2D,3];*)
		eqFull2D=(#/Total[#])&[leftEV2D[[1]]];
		plotBounds2D = {{0,0},{0,30},{30,30},{30,0}}+0.5;
		rangeX2D = {0,30}+0.5;
		rangeY2D = {0,30}+0.5;
		listOfStatePositions2D = Table[{Mod[ii,30]+1,Floor[ii/30]+1},{ii,0,30^2-1}];
		nearestFnc2D = Nearest[listOfStatePositions2D -> Automatic];
		minimumStates2D = {311,611,411};
		{
			"NumberOfStates"  -> Length[tMFull2D],
			"TMatrix"         -> tMFull2D,
			"RightEV"         -> rightEV2D,
			"LeftEV"          -> leftEV2D,
			"EValues"         -> eigenV2D,
			"Energies"        -> listOfEnergies2D,
			"Geometry"        -> {
				"ScaleUp"     -> Function[p, p * 30 + 0.5],
				"ScaleDown"   -> Function[x, (x - 0.5)/30],
				"Border"      -> plotBounds2D,
				"Ranges"      -> {rangeX2D, rangeY2D},
				"Positions"   -> listOfStatePositions2D,
				"GetPosition" -> Function[s, listOfStatePositions2D[[s]]],
				"GetState"    -> Function[p, nearestFnc2D[p][[1]]]
			},
			"MinimumStates"   -> minimumStates2D,
			"EQDist" -> eqFull2D,
			"PathToResources" -> pathToResources
		}
	];

LoadDynamicalModel["4Well1D",opts___]:=
	Module[
		{
			range,
			numberOfStates
		},
		range = {-1.0, 1.0};

		numberOfStates = States/.{opts}/.{States->100};
		beta = (Beta /.{opts}) /. {Beta->1.0};

		CreateDiffusionModel1D[
			Function[x,4(x^8-0.0x+0.5Exp[-40(x+0.5)^2]+0.8Exp[-80(x+0)^2]+0.2Exp[-80(x-0.5)^2])], 
			range,
			numberOfStates,
			Beta->beta
		]
	];

CreateSimulation[tMatrix_,length_] := Module[{},
		timeS1D = ToTimeSeries[traMMSM1D, 2, length];
		
		gaussTimeS1D = Map[RandomReal[Apply[NormalDistribution,gaussParameters1D[[#]]]]&,timeS1D[[1]]];
		gaussTimeS1D2 = Select[gaussTimeS1D,#>0&&#<1&];
		discTimeS1D = TimeSeries[Map[nearestFnc1D[#][[1]]&,gaussTimeS1D2]];
		discCountMatrix1D = CountMatrix[N[Normal[ToCountMatrix[discTimeS1D,Lagtime->1,Closed->False][[1]]]]];
		corMEstSim1D = Symmetrize[CorrelationMatrix[#/Total[Flatten[#]]&[discCountMatrix1D[[1]]]]];
		traMEstSim1D = ToTransitionMatrix[corMEstSim1D];

		{eigenVEst1D,leftEVEst1D,rightEVEst1D} = FullEigensystem[traMEstSim1D];
		{
			"TMatrix"     -> traMEstSim1D,
			"CMatrix"     -> corMEstSim1D,
			"RightEV"     -> rightEVEst1D,
			"LeftEV"      -> leftEVEst1D,
			"EValues"     -> eigenVEst1D,
			"EQDist"      -> Equilibrium[traMEstSim1D]
		}
]

LoadDynamicalModel["3StateGaussianNoise"] :=
	Module[
		{
			listOfStatePositions1D,
			nearestFnc1D,
			eqMSM1D,
			range1D,
			gaussParameters1D,
			numberOfStates,
			traMMSM1D,
			gaussFnc1D,
			unnormedMembership1D,
			perfectMembership1D,
			corMExact1D,
			traMExact1D,
			timeS1D,
			gaussTimeS1D,
			gaussTimeS1D2,
			discTimeS1D,
			discCountMatrix1D,
			corMEstSim1D,
			traMEstSim1D,
			eigenVEst1D, leftEVEst1D, rightEVEst1D,
			eigenVExact1D, leftEVExact1D, rightEVExact1D,
			eigenVMSM1D, leftEVMSM1D, rightEVMSM1D

		},
		gaussParameters1D = {{0.3,0.08},{0.4,0.08},{0.8,0.06}};
		range1D = {0.0, 1.0};
		numberOfStates = 100;

		listOfStatePositions1D = MeanRange[{range1D}, 1, numberOfStates];
		nearestFnc1D = Nearest[listOfStatePositions1D -> Automatic];

		traMMSM1D=GetExampleTransitionMatrix["3StateDBExample1"];
		eqMSM1D=Equilibrium[traMMSM1D];

		gaussFnc1D = Map[PDF[Apply[NormalDistribution,#]]&,gaussParameters1D];
		SetAttributes[gaussFnc1D,Listable];

		perfectMembership1D = ClusterMatrix[Transpose[Map[#/Apply[Plus,#]&,Transpose[Map[#[listOfStatePositions1D]&,gaussFnc1D]]]]];
		unnormedMembership1D = ClusterMatrix[Transpose[Map[#&,Transpose[Map[#[listOfStatePositions1D]&,gaussFnc1D]]]]/100];

		corMExact1D = Expand[unnormedMembership1D,ToCorrelationMatrix[traMMSM1D]];
		traMExact1D = ToTransitionMatrix[corMExact1D];

		timeS1D = ToTimeSeries[traMMSM1D, 2, 500000];
		
		gaussTimeS1D = Map[RandomReal[Apply[NormalDistribution,gaussParameters1D[[#]]]]&,timeS1D[[1]]];
		gaussTimeS1D2 = Select[gaussTimeS1D,#>0&&#<1&];
		discTimeS1D = TimeSeries[Map[nearestFnc1D[#][[1]]&,gaussTimeS1D2]];
		discCountMatrix1D = CountMatrix[N[Normal[ToCountMatrix[discTimeS1D,Lagtime->1,Closed->False][[1]]]]];
		corMEstSim1D = Symmetrize[CorrelationMatrix[#/Total[Flatten[#]]&[discCountMatrix1D[[1]]]]];
		traMEstSim1D = ToTransitionMatrix[corMEstSim1D];

		{eigenVEst1D,leftEVEst1D,rightEVEst1D} = FullEigensystem[traMEstSim1D];
		{eigenVExact1D,leftEVExact1D,rightEVExact1D} = FullEigensystem[traMExact1D];
		{eigenVMSM1D,leftEVMSM1D,rightEVMSM1D} = FullEigensystem[traMMSM1D];

		{
			"NumberOfStates"  -> Length[tMFull1D],
			"Exact"           -> {
				"TMatrix"     -> traMExact1D,
				"CMatrix"     -> corMExact1D,
				"RightEV"     -> rightEVExact1D,
				"LeftEV"      -> leftEVExact1D,
				"EValues"     -> eigenVExact1D,
				"EQDist"      -> Equilibrium[traMExact1D]
			},
			"Estimation"      -> {
				"TMatrix"     -> traMEstSim1D,
				"CMatrix"     -> corMEstSim1D,
				"RightEV"     -> rightEVEst1D,
				"LeftEV"      -> leftEVEst1D,
				"EValues"     -> eigenVEst1D,
				"EQDist"      -> Equilibrium[traMEstSim1D]
			},
			"MSM"             -> {
				"TMatrix"     -> traMMSM1D,
				"RightEV"     -> rightEVMSM1D,
				"LeftEV"      -> leftEVMSM1D,
				"EValues"     -> eigenVMSM1D,
				"EQDist"      -> Equilibrium[traMMSM1D],
				"Membership"  -> unnormedMembership1D,
				"GaussParameter" -> gaussParameters1D				
			},
			"Energies"        -> {},
			"Geometry"        -> {
				"ScaleUp"     -> Function[p, Rescale[p,{0,1},range1D[[1]]]],
				"ScaleDown"   -> Function[x, Rescale[x, range1D[[1]],{0,1}]],
				"Border"      -> {},
				"Ranges"      -> {range1D},
				"Positions"   -> listOfStatePositions1D,
				"GetPosition" -> Function[s, listOfStatePositions1D[[s]]],
				"GetState"    -> Function[p, nearestFnc1D[p][[1]]]
			},
			"TimeSeries"      -> timeS1D,
			"MinimumStates"   -> {Map[Function[p, nearestFnc1D[p][[1]]],gaussParameters1D[[All,1]]]}
		}
	];

CreateMultiStateGaussModel[] :=
	Module[
		{},
		Null
		];

CreateDiffusionModel1D[potential1D_, range1D_, numberOfStates1D_, 
  opts___] :=
    Module[ {diagonalEntries, upperDiagonalEntries, 
      lowerDiagonalEntries, rules, beta, partitionSum1D, 
      probabilityFunction1D, probabilityFunctionRelative1D, 
      listOfStatePositions1D, listOfStatePositions1DBins, nearestFnc1D, 
      tMFull1D, tMFull1DPlain, eigenV1D, leftEV1D, rightEV1D, 
      listOfEnergies1D, eq1D, potentialFunction},
        beta = (Beta /. List[opts]) /. {Beta -> 1.0};
        partitionSum1D = 
         NIntegrate[Exp[-potential1D[x]], {x, range1D[[1]], range1D[[2]]}];
        probabilityFunction1D[x_] = 1/partitionSum1D*Exp[-potential1D[x]];
        probabilityFunctionRelative1D[x_, y_] = 
         Min[1, Exp[-beta (potential1D[y] - potential1D[x])]];
        listOfStatePositions1D = 
         Positions /. 
           List[opts] /. {Positions -> 
            MeanRange[{range1D}, 1, numberOfStates1D]};
        listOfStatePositions1DBins = 
         MeanRangeBin[{range1D}, 1, listOfStatePositions1D];
        nearestFnc1D = First[Nearest[listOfStatePositions1D -> Automatic]];
        diagonalEntries = Table[{ii, ii} -> 1, {ii, 1, numberOfStates1D}];
        upperDiagonalEntries = 
         Table[{ii + 1, ii} -> 
           probabilityFunctionRelative1D[listOfStatePositions1D[[ii + 1]], 
            listOfStatePositions1D[[ii]]], {ii, 1, numberOfStates1D - 1}];
        lowerDiagonalEntries = 
         Table[{ii, ii + 1} -> 
           probabilityFunctionRelative1D[listOfStatePositions1D[[ii]], 
            listOfStatePositions1D[[ii + 1]]], {ii, 1, 
           numberOfStates1D - 1}];
        rules = Join[upperDiagonalEntries, lowerDiagonalEntries];
        array = SparseArray[rules];
        rowsums = Total /@ array;
        array = array/Max[rowsums]*0.5;
        array = array + DiagonalMatrix[1 - 0.5*rowsums/Max[rowsums]];
        tMFull1D = TransitionMatrix[Normal[array]];
        tMFull1DPlain = Normal[tMFull1D[[1]]];
        listOfEnergies1D = Map[potential1D, listOfStatePositions1D];
        {eigenV1D, leftEV1D, rightEV1D} = FullEigensystem[tMFull1D];
        leftEV1D = EquilizeLeftEVSet[leftEV1D];
        eq1D = (#/Total[#]) &[leftEV1D[[1]]];
        rightEV1D = EquilizeRightEVSet[rightEV1D, eq1D];
        potentialFunction = potential1D;
        {
        	"NumberOfStates" -> Length[tMFull1D], 
        	"TMatrix" -> tMFull1D, 
        	"RightEV" -> rightEV1D, 
        	"LeftEV" -> leftEV1D, 
        	"EValues" -> eigenV1D, 
        	"Energies" -> listOfEnergies1D, 
        	"Potential" -> potentialFunction, 
        	"Geometry" -> {
        		"ScaleUp" -> Function[p, Rescale[p, {0, 1}, range1D[[1]]]], 
	            "ScaleDown" -> Function[x, Rescale[x, range1D[[1]], {0, 1}]], 
                "Border" -> {}, 
                "Ranges" -> {range1D}, 
                "Positions" -> listOfStatePositions1D, 
                "GetPosition" -> Function[s, listOfStatePositions1D[[s]]], 
                "GetState" -> Function[p, nearestFnc1D[p][[1]]]
            }, 
            "MinimumStates" -> {}, 
            "EQDist" -> eq1D, 
         	"PathToResources" -> ""
         }
    ]



CreateDiffusionModel2D[potential2D_?VectorQ, opts___] :=
    Module[ {rules, normRules, rowRules, rowSum, beta, 
      tMFull2D, tMFull2DPlain, eigenV2D, leftEV2D, rightEV2D, 
      listOfEnergies2D, eq2D, potentialFunction, numberOfStates2D},
        beta = Beta /. List[opts] /. {Beta -> 1.0};
        numberOfStates2D = Sqrt@Length[potential2D];
        state2position = 
         Thread[Range[numberOfStates2D*numberOfStates2D] -> 
           Flatten[Outer[List, Range[numberOfStates2D], 
             Range[numberOfStates2D]], 1]];
        distList = 
         Outer[Norm[#1 - #2] &, state2position[[All, 2]], 
          state2position[[All, 2]], 1];
        neighborList = Table[Idx[dist, # == 1 &], {dist, distList}];
        matNeighbor = 
         Table[Table[{from, to} -> 
            Min[1, Exp[-beta (potential2D[[to]] - potential2D[[from]])]], {to, 
            neighborList[[from]]}], {from, numberOfStates2D^2}];
        matDiagonal = Table[{from, from} -> 1, {from, numberOfStates2D^2}];
        rules = Join[Flatten@matNeighbor, matDiagonal];
        normRules = 
         Flatten[Table[rowRules = Select[rules, #[[1, 1]] == st &];
                       rowSum = Apply[Plus, rowRules[[All, 2]]]*1.;
                       Map[#[[1]] -> #[[2]]/rowSum &, rowRules], {st, 1, 
            numberOfStates2D^2}], 1];
        tMFull2D = TransitionMatrix[Normal[SparseArray[normRules]]];
        tMFull2DPlain = Normal[tMFull2D[[1]]];
        listOfEnergies2D = potential2D;
        {eigenV2D, leftEV2D, rightEV2D} = FullEigensystem[tMFull2D];
        leftEV2D = EquilizeLeftEVSet[leftEV2D];
        eq2D = (#/Total[#]) &[leftEV2D[[1]]];
        rightEV2D = EquilizeRightEVSet[rightEV2D, eq2D];
        potentialFunction = potential2D;
        {
        	"NumberOfStates" -> Length[tMFull1D], 
        	"TMatrix" -> tMFull2D, 
         "RightEV" -> rightEV2D, 
         "LeftEV" -> leftEV2D, 
         "EValues" -> eigenV2D, 
         "Energies" -> listOfEnergies2D, 
         "Potential" -> potentialFunction, 
         "MinimumStates" -> {}, 
         "EQDist" -> eq2D, 
         "PathToResources" -> ""
         }
    ]

CreateDiscreteModel1D[potential1D_,topology_,opts___] := 
	Module[
		{
			diagonalEntries,
			offDiagonalEntries,
			rules,
			normRules,
			rowRules,
			rowSum,
			beta,
			partitionSum1D,
			probabilityFunction1D,
			probabilityFunctionRelative1D,
			listOfStatePositions1D,
			tMFull1D,
			tMFull1DPlain,
			eigenV1D,leftEV1D,rightEV1D,
			listOfEnergies1D,
			eq1D,
			numberOfStates1D
		},

		numberOfStates1D = Length[potential1D];
		beta = Beta /. List[opts] /. {Beta -> 1.0};
		timeRescale = TimeRescale /. List[opts] /. {TimeRescale -> 1};
		partitionSum1D=Total[Exp[-potential1D]];
		probabilityFunction1D[x_] = 1/partitionSum1D * Exp[-potential1D];
		probabilityFunctionRelative1D[x_,y_] := Min[1, Exp[-beta(potential1D[[y]] - potential1D[[x]])]];

		listOfStatePositions1D = Positions /. List[opts] /. {Positions->Range[numberOfStates1D]};

		diagonalEntries=Table[{ii,ii}->1,{ii,1,numberOfStates1D}];
		offDiagonalEntries=Map[If[#[[1,1]]!=#[[1,2]],#[[1]]->probabilityFunctionRelative1D[#[[1,1]],#[[1,2]]]]&,Select[ArrayRules[topology],#[[2]]>0&]];

		rules=Join[diagonalEntries,offDiagonalEntries];
		normRules=Flatten[Table[
			rowRules=Select[rules,#[[1,1]]==st&];
			rowSum=Apply[Plus,rowRules[[All,2]]]*1.;
			Map[#[[1]]->#[[2]]/rowSum&,rowRules]
		,{st,1,numberOfStates1D}],1];

		tMFull1D=TransitionMatrix[Normal[SparseArray[normRules]]];
		{eigenV1D,leftEV1D,rightEV1D} = FullEigensystem[tMFull1D];
		tMFull1D=TransitionMatrix[FromFullEigensystem[{Abs[eigenV1D]^timeRescale,leftEV1D,rightEV1D}]];
		
		tMFull1DPlain=Normal[tMFull1D[[1]]];

		listOfEnergies1D = potential1D;
		{eigenV1D,leftEV1D,rightEV1D} = FullEigensystem[tMFull1D];

		eq1D = (#/Total[#])&[leftEV1D[[1]]];
		{
			"NumberOfStates"  -> Length[tMFull1D],
			"TMatrix"         -> tMFull1D,
			"RightEV"         -> rightEV1D,
			"LeftEV"          -> leftEV1D,
			"EValues"         -> eigenV1D,
			"Energies"        -> listOfEnergies1D,
			"Potential"       -> potential1D,
			"Geometry"        -> {
				"ScaleUp"     -> Null,
				"ScaleDown"   -> Null,
				"Border"      -> {},
				"Ranges"      -> {{1,numberOfStates1D}},
				"Positions"   -> listOfStatePositions1D,
				"GetPosition" -> Function[s, listOfStatePositions1D[[s]]],
				"GetState"    -> Null
			},
			"MinimumStates"   -> {},
			"EQDist"          -> eq1D,
			"PathToResources" -> ""
		}		
]

CreateDiscreteBarrierModel1D[potential1D_,topology_,opts___] := 
	Module[
		{
			diagonalEntries,
			offDiagonalEntries,
			rules,
			normRules,
			rowRules,
			rowSum,
			beta,
			partitionSum1D,
			probabilityFunction1D,
			probabilityFunctionRelative1D,
			listOfStatePositions1D,
			tMFull1D,
			tMFull1DPlain,
			eigenV1D,leftEV1D,rightEV1D,
			listOfEnergies1D,
			eq1D,
			numberOfStates1D,
			barrier
		},

		numberOfStates1D = Length[potential1D];
		beta = Beta /. List[opts] /. {Beta -> 1.0};
		barrier = Barrier /. List[opts] /. {Barrier -> 4};
		partitionSum1D=Total[Exp[-potential1D]];
		probabilityFunction1D[x_] = 1/partitionSum1D * Exp[-potential1D];
		probabilityFunctionRelative1D[x_,y_] := Exp[-barrier]*Min[1, Exp[-beta(potential1D[[y]] - potential1D[[x]])]];

		listOfStatePositions1D = Positions /. List[opts] /. {Positions->Range[numberOfStates1D]};

		offDiagonalEntries=Map[If[#[[1,1]]!=#[[1,2]],#[[1]]->probabilityFunctionRelative1D[#[[1,1]],#[[1,2]]]]&,Select[ArrayRules[topology],#[[2]]>0&]];

		rules=offDiagonalEntries;
		diagonalEntries=Table[
			rowRules=Select[rules,#[[1,1]] == st&];
			rowSum=Apply[Plus,rowRules[[All,2]]]*1.;
			{st,st}->(1 - rowSum)
		,{st,1,numberOfStates1D}];

		normRules = Join[offDiagonalEntries,diagonalEntries];
		normRules = Map[#[[1]]->#[[2]]&,normRules];

		tMFull1D=TransitionMatrix[Normal[SparseArray[normRules]]];
		tMFull1DPlain=Normal[tMFull1D[[1]]];

		listOfEnergies1D = potential1D;
		{eigenV1D,leftEV1D,rightEV1D} = FullEigensystem[tMFull1D];

		eq1D = (#/Total[#])&[leftEV1D[[1]]];
		{
			"NumberOfStates"  -> Length[tMFull1D],
			"TMatrix"         -> tMFull1D,
			"RightEV"         -> rightEV1D,
			"LeftEV"          -> leftEV1D,
			"EValues"         -> eigenV1D,
			"Energies"        -> listOfEnergies1D,
			"Potential"       -> potential1D,
			"Geometry"        -> {
				"ScaleUp"     -> Null,
				"ScaleDown"   -> Null,
				"Border"      -> {},
				"Ranges"      -> {{1,numberOfStates1D}},
				"Positions"   -> listOfStatePositions1D,
				"GetPosition" -> Function[s, listOfStatePositions1D[[s]]],
				"GetState"    -> Null
			},
			"MinimumStates"   -> {},
			"EQDist"          -> eq1D,
			"PathToResources" -> ""
		}		
]

CreateDiscreteRateModel1D[potential1D_,topology_,opts___] := 
	Module[
		{
			diagonalEntries,
			offDiagonalEntries,
			rules,
			normRules,
			rowRules,
			rowSum,
			beta,
			partitionSum1D,
			probabilityFunction1D,
			probabilityFunctionRelative1D,
			listOfStatePositions1D,
			kMFull1D,
			timeRescale,
			tMFull1D,
			tMFull1DPlain,
			eigenV1D,leftEV1D,rightEV1D,
			listOfEnergies1D,
			eq1D,
			numberOfStates1D
		},

		numberOfStates1D = Length[potential1D];
		beta = Beta /. List[opts] /. {Beta -> 1.0};
		timeRescale = TimeRescale /. List[opts] /. {TimeRescale -> 1/100};
		partitionSum1D=Total[Exp[-potential1D]];
		probabilityFunction1D[x_] = 1/partitionSum1D * Exp[-potential1D];
		probabilityFunctionRelative1D[x_,y_] := Min[1, Exp[-beta(potential1D[[y]] - potential1D[[x]])]];

		listOfStatePositions1D = Positions /. List[opts] /. {Positions->Range[numberOfStates1D]};

		offDiagonalEntries=Map[If[#[[1,1]]!=#[[1,2]],#[[1]]->probabilityFunctionRelative1D[#[[1,1]],#[[1,2]]]]&,Select[ArrayRules[topology],#[[2]]>0&]];

		rules=offDiagonalEntries;
		diagonalEntries=Table[
			rowRules=Select[rules,#[[1,1]]==st&];
			rowSum=Apply[Plus,rowRules[[All,2]]]*1.;
			{st,st}->-rowSum
		,{st,1,numberOfStates1D}];

		normRules = Join[offDiagonalEntries,diagonalEntries];
		normRules = Map[#[[1]]->#[[2]]*timeRescale&,normRules];

		kMFull1D=RateMatrix[Normal[SparseArray[normRules]]];
		tMFull1D=ToTransitionMatrix[kMFull1D,1];

		tMFull1DPlain=Normal[tMFull1D[[1]]];

		listOfEnergies1D = potential1D;
		{eigenV1D,leftEV1D,rightEV1D} = FullEigensystem[tMFull1D];

		eq1D = (#/Total[#])&[leftEV1D[[1]]];
		{
			"NumberOfStates"  -> Length[tMFull1D],
			"KMatrix"         -> kMFull1D,
			"TMatrix"         -> tMFull1D,
			"RightEV"         -> rightEV1D,
			"LeftEV"          -> leftEV1D,
			"EValues"         -> eigenV1D,
			"Energies"        -> listOfEnergies1D,
			"Potential"       -> potential1D,
			"Geometry"        -> {
				"ScaleUp"     -> Null,
				"ScaleDown"   -> Null,
				"Border"      -> {},
				"Ranges"      -> {{1,numberOfStates1D}},
				"Positions"   -> listOfStatePositions1D,
				"GetPosition" -> Function[s, listOfStatePositions1D[[s]]],
				"GetState"    -> Null
			},
			"MinimumStates"   -> {},
			"EQDist"          -> eq1D,
			"PathToResources" -> ""
		}		
]

CreateTransitionMatrix[size_,opts___] := 
	Module[
		{
			countsPerRow,
			neighborRange,
			offdiagonalRange,
			matrix,
			rules,
			normRules,
			rowRules,
			rowSum,
			tmFull,
			r1,
			symmetric
		},

		symmetric = Symmetric /. List[opts] /. Symmetric->True;
		countsPerRow = RowCounts /. List[opts] /. {RowCounts -> 1000 * Table[1,{size}]};
		neighborRange = NeighborRange /. List[opts] /. NeighborRange -> {1,10};
		offdiagonalRange = OffDiagonalRange /. List[opts] /. OffDiagonalRange -> {0,10};
		If[offdiagonalRange[[2]] > 0,
			matrix = Table[
				Which[
					ii == jj, 0, 
					ii == jj+1, RandomInteger[neighborRange],
					ii == jj-1, RandomInteger[neighborRange],
					True, RandomInteger[offdiagonalRange]
				], {ii,1,size}, {jj,1,size}
			];
			Table[ matrix[[ii,ii]] = countsPerRow[[ii]] - Total[matrix[[ii]]], {ii,1,size} ];
			If[symmetric, matrix = 1/2(matrix + Transpose[matrix])];
			tmFull = Normalize[TransitionMatrix[matrix]];
		,
			rules = Flatten[Table[{{ii, ii+1}-> (r1 = RandomInteger[neighborRange]), {ii + 1, ii} -> If[symmetric, r1, RandomInteger[neighborRange]]}, {ii,1,size - 1} ],1];
			normRules=Flatten[Table[
				rowRules=Select[rules,#[[1,1]]==st&];
				rowSum=Apply[Plus,rowRules[[All,2]]];
				Join[{{st,st}->(countsPerRow[[st]] - rowSum)/countsPerRow[[st]]}, Map[#[[1]]->#[[2]]/countsPerRow[[st]]&, rowRules]]
			,{st,1,size}],1];
			tmFull = TransitionMatrix[Normal[SparseArray[normRules]]];
		];

		tmFull
	];

GetExampleTransitionMatrix["3StateDBExample1"] = Module[{transitionMatrixUSym1D,transitionMatrix1D},
	transitionMatrixUSym1D=TransitionMatrix[{{0.995,0.002,0.003},{0.002,0.99,0.008},{0.005,0.015,0.98}}];
	transitionMatrix1D=TransitionMatrix[Map[#/Total[#]&,0.5(#+Transpose[#])&[DiagonalMatrix[Equilibrium[transitionMatrixUSym1D]].transitionMatrixUSym1D[[1]]]]];
	transitionMatrix1D
];
GetExampleTransitionMatrix["3StateExample1"] = Module[{transitionMatrixUSym1D},
	transitionMatrixUSym1D=TransitionMatrix[{{0.995,0.002,0.003},{0.002,0.99,0.008},{0.005,0.015,0.98}}];
	transitionMatrixUSym1D
];

BeginPackage["DynamicalModels`"];

Begin["`Private`"];

End[ ];

EndPackage[ ];