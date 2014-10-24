(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *) 

RateMatrix::usage = "RateMatrix[matrix, options] creates a RateMatrix Object.";
RateMatrixQ::usage = "RateMatrixQ[system] returns true, if system is a valid RateMatrix object, false otherwise.";


Begin["`Private`"] (* Begin Private Context *) 

RateMatrix/:MakeBoxes[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]],StandardForm] := 
        If[ $DefaultShowMatrixSize>Length[m],
            InterpretationBox[#,RateMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalR]",Length[m]],"\[Rule]",ToBoxes[MatrixForm[m]]}]],
            InterpretationBox[#,RateMatrix[m,opts]]&[SubscriptBox["\[DoubleStruckCapitalR]",Length[m]]]
        ]
RateMatrix/:RateMatrix[RateMatrix[m_?MatrixQ,opts1:OptionsPattern[RateMatrix]],opts2:OptionsPattern[RateMatrix]] := 
        RateMatrix[m,UnionRules[opts1,opts2]]
RateMatrixQ[x_] :=
    If[ Head[x]==RateMatrix&&MatrixQ[x[[1]]]&&RowSumZeroQ[x[[1]]],
        True,
        False
    ];
RateMatrix/:NormalForm[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]]] := 
        {m,opts}
RateMatrix/:Normal[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]]] := 
        m
RateMatrix/:Matrix[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]]] := 
        m
RateMatrix/:Spectrum[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]],optsSP___] := 
        Eigenvalues[N[m],optsSP]
RateMatrix/:MatrixForm[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]]] := 
        MatrixForm[m]
RateMatrix/:Length[RateMatrix[m_?MatrixQ],opts:OptionsPattern[RateMatrix]] := 
        Length[m]
RateMatrix/:RowSumZeroQ[RateMatrix[m_?MatrixQ],opts:OptionsPattern[RateMatrix]] := 
        RowSumZeroQ[m]
RateMatrix/:MatrixPlot[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]],optsMP:OptionsPattern[MatrixPlot]] := 
        MatrixPlot[m,optsMP,$TMPkgDefaultMatrixPlotOptions]
RateMatrix/:ArrayPlot[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix],optsAP:OptionsPattern[ArrayPlot]]] := 
        ArrayPlot[m,optsAP,$TMPkgDefaultArrayPlotOptions]
RateMatrix/:N[RateMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        RateMatrix[N[m],opts]
RateMatrixQ[x_] :=
    If[ Head[x]==RateMatrix&&MatrixQ[x[[1]]],
        True,
        False
    ];
RateMatrix/:NormalizedQ[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]]] := 
        RowSumZeroQ[m]
RateMatrix/:Normalize[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]]] := 
        RateMatrix[Map[#-Apply[Plus,#]&,m],opts]
RateMatrix/:Equilibrium[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]]] := 
        Stochastize[Eigenvectors[Transpose[m]][[-1]]]
RateMatrix/:ToTransitionMatrix[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]],nLagtime_] := 
        TransitionMatrix[MatrixExp[m*nLagtime],Lagtime->nLagtime]

RateMatrix/:ImpliedTimeScales[RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]],optsITS___] := 
        Module[ {nEV},
            nEV = Eigenvalues/.List[UnionRules[optsITS,Eigenvalues->Length[m]]];
            Join[{Infinity},Reverse[-1.0/Eigenvalues[N[m]][[Range[1,Length[m]-1]]]]]
        ]


End[] (* End Private Context *)

EndPackage[]