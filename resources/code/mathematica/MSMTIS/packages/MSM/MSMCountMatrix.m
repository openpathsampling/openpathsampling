(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

CountMatrix::usage = "CountMatrix[matrix, options] creates a CountMatrix Object that contains number of occurrances of observations.";
CountMatrixQ::usage = "CountMatrixQ[system] returns true if system is a valid countmatrix object.";

Options[CountMatrix]                 = {Lagtime->1,Distribution->DirichletDistribution,Prior->1};

Begin["`Private`"] (* Begin Private Context *) 

(* CountMatrix functions *)

CountMatrix/:MakeBoxes[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]],StandardForm] := 
        If[ $DefaultShowMatrixSize>Length[m],
            InterpretationBox[#,CountMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalC]",Length[m]],"\[Rule]","[",Lagtime[CountMatrix[m,opts]],"]",ToBoxes[MatrixForm[m]]}]],
            InterpretationBox[#,CountMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalC]",Length[m]],"[",Lagtime[CountMatrix[m,opts]],"]"}]]
        ]
CountMatrix/:CountMatrix[CountMatrix[m_?MatrixQ,opts1:OptionsPattern[CountMatrix]],opts2:OptionsPattern[CountMatrix]] := 
        CountMatrix[m,UnionRules[opts1,opts2]]
CountMatrix/:NormalForm[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        {m,opts}
CountMatrix/:Normal[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        m
CountMatrix/:Matrix[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        m
CountMatrix/:Normalize[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        TransitionMatrix[Map[#/Apply[Plus,#]&,m]]
CountMatrix/:MatrixForm[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        MatrixForm[m]
CountMatrix/:Lagtime[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        Lagtime/.List[opts]/.{Lagtime->1}
CountMatrix/:Length[CountMatrix[m_?MatrixQ],opts:OptionsPattern[CountMatrix]] := 
        Length[m]
CountMatrix/:MatrixPlot[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]],optsMP:OptionsPattern[MatrixPlot]] := 
        MatrixPlot[m,optsMP,$TMPkgDefaultMatrixPlotOptions]
CountMatrixQ[x_] :=
    If[ Head[x]===CountMatrix&&MatrixQ[x[[1]]],
        True,
        False
    ];
CountMatrix/:VarianceMatrix[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        VarianceMatrix[OptionValue[CountMatrix,Distribution][CountMatrix[m,opts]]]
CountMatrix/:VarianceTensor[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]],states_] := 
        VarianceTensor[OptionValue[CountMatrix,Distribution][CountMatrix[m,opts]],states]
CountMatrix/:MeanMatrix[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        MeanMatrix[OptionValue[CountMatrix,Distribution][CountMatrix[m,opts]]]
CountMatrix/:MaximumLikelihoodMatrix[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]],optsA___] := 
        MaximumLikelihoodMatrix[OptionValue[CountMatrix,Distribution][CountMatrix[m,opts]],optsA]
CountMatrix/:DrawTransitionMatrix[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        DrawTransitionMatrix[OptionValue[CountMatrix,Distribution][CountMatrix[m,opts]]]
CountMatrix/:NewStateProbabilityVector[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        Normal[MeanMatrix[OptionValue[CountMatrix,Distribution][CountMatrix[ReplacePart[ArrayPad[m,{{0,1},{0,1}}],{Length[m]+1,Length[m]+1}->1],opts]]]][[Range[1,Length[m]],Length[m]+1]]
CountMatrixTotal[v_?VectorQ] :=
    Module[ {cMatrixQ, matchSizeQ},
        cMatrixQ = Apply[And,Select[v,CountMatrixQ]];
        If[ cMatrixQ,
            matchSizeQ = (Min[Map[Length,v]]==Max[Map[Length,v]]);
            If[ matchSizeQ,
                CountMatrix[Apply[Plus,v[[All,1]]]]
            ]
        ]
    ]

End[] (* End Private Context *)

EndPackage[]