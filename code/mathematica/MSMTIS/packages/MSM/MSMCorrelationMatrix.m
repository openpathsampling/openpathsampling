(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

CorrelationMatrix::usage = "CorrelationMatrix[matrix, options] creates a TransitionMatrix Object.";
CorrelationMatrixQ::usage = "CorrelationMatrixQ[system] return true if system is a valid correlation matrix, false otherwise.";

Options[CorrelationMatrix]           = {Lagtime->1};

Begin["`Private`"] (* Begin Private Context *) 

(* CorrelationMatrix functions *)

CorrelationMatrix/:MakeBoxes[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]],StandardForm] := 
        If[ $DefaultShowMatrixSize>Length[m],
            InterpretationBox[#,CorrelationMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalK]",Length[m]],"[",Lagtime[CorrelationMatrix[m,opts]],"]","\[Rule]",ToBoxes[MatrixForm[m]]}]],
            InterpretationBox[#,CorrelationMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalK]",Length[m]],"[",Lagtime[CorrelationMatrix[m,opts]],"]"}]]
        ]
CorrelationMatrix/:CorrelationMatrix[CorrelationMatrix[m_?MatrixQ,opts1:OptionsPattern[CorrelationMatrix]],opts2:OptionsPattern[CorrelationMatrix]] := 
        CorrelationMatrix[m,UnionRules[opts1,opts2]]
CorrelationMatrix/:NormalForm[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        {m,opts}
CorrelationMatrix/:Normal[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        m
CorrelationMatrix/:Matrix[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        m
CorrelationMatrix/:Lagtime[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        Lagtime/.List[opts]/.{Lagtime->1}
CorrelationMatrix/:MatrixForm[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        MatrixForm[m]
CorrelationMatrix/:Length[CorrelationMatrix[m_?MatrixQ],opts:OptionsPattern[CorrelationMatrix]] := 
        Length[m]
CorrelationMatrix/:SymmetricQ[CorrelationMatrix[m_?MatrixQ],opts:OptionsPattern[CorrelationMatrix]] := 
        Max[Abs[Chop[m-Transpose[m]]]]==0.
CorrelationMatrix/:MatrixPlot[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]],optsMP:OptionsPattern[MatrixPlot]] := 
        MatrixPlot[m,optsMP,$TMPkgDefaultMatrixPlotOptions]
CorrelationMatrix/:ArrayPlot[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix],optsAP:OptionsPattern[ArrayPlot]]] := 
        ArrayPlot[m,optsAP,$TMPkgDefaultArrayPlotOptions]
CorrelationMatrixQ[x_] :=
    If[ Head[x]===CorrelationMatrix&&MatrixQ[x[[1]]],
        True,
        False
    ];
CorrelationMatrix/:Chop[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        CorrelationMatrix[Chop[m],opts]
CorrelationMatrix/:NormalizedQ[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        TotalSumOneQ[m]
CorrelationMatrix/:Normalize[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        CorrelationMatrix[m/Total[Flatten[m]],opts]
CorrelationMatrix/:Equilibrium[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        #/Total[#]&[Map[Apply[Plus,#]&,m]]
CorrelationMatrix/:Symmetrize[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        Module[ {},
            CorrelationMatrix[(m+Transpose[m])/2, opts]
        ]
CorrelationMatrix/:LeftEigenvectors[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        ToLeftEigenvectors[m,All];
CorrelationMatrix/:RightEigenvectors[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        ToRightEigenvectors[m,All];
CorrelationMatrix/:FullEigensystem[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        ToFullEigensystem[m,All];
CorrelationMatrix/:Eigenvalues[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        Eigenvalues[m]

End[] (* End Private Context *)

EndPackage[]