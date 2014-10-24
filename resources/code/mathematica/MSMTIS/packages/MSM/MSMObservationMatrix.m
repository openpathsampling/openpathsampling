(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ObservationMatrix::usage = "";

Options[ObservationMatrix]           = {ShowSize->7, Lagtime->1};


Begin["`Private`"] (* Begin Private Context *) 

ObservationMatrix/:MakeBoxes[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]],StandardForm] := 
        If[ OptionValue[ObservationMatrix,ShowSize]>Length[m],
            InterpretationBox[#,ObservationMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalO]",Length[m]],"[",Lagtime[ObservationMatrix[m,opts]],"]","\[Rule]",ToBoxes[MatrixForm[m]]}]],
            InterpretationBox[#,ObservationMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalO]",Length[m]],"[",Lagtime[ObservationMatrix[m,opts]],"]"}]]
        ]

ObservationMatrix/:ObservationMatrix[ObservationMatrix[m_?MatrixQ,opts1:OptionsPattern[ObservationMatrix]],opts2:OptionsPattern[ObservationMatrix]] := 
        ObservationMatrix[m,UnionRules[opts1,opts2]]

ObservationMatrix/:NormalForm[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        {m,opts}

ObservationMatrix/:Normal[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        m

ObservationMatrix/:Matrix[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        m

ObservationMatrix/:Lagtime[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        Lagtime/.List[opts]/.{Lagtime->1}

ObservationMatrix/:MatrixForm[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        MatrixForm[m]

ObservationMatrix/:Length[ObservationMatrix[m_?MatrixQ, opts:OptionsPattern[ObservationMatrix]]] := 
        Length[m]

ObservationMatrix/:SymmetricQ[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        SymmetricQ[m]

ObservationMatrix/:MatrixPlot[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]],optsMP:OptionsPattern[MatrixPlot]] := 
        MatrixPlot[m,optsMP,$TMPkgDefaultMatrixPlotOptions]
ObservationMatrix/:ArrayPlot[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix],optsAP:OptionsPattern[ArrayPlot]]] := 
        ArrayPlot[m,optsAP,$TMPkgDefaultArrayPlotOptions]
ObservationMatrix[x_] :=
    If[ Head[x]===ObservationMatrix&&MatrixQ[x[[1]]],
        True,
        False
    ];
    
ObservationMatrix/:Chop[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        ObservationMatrix[Chop[m],opts]
        
ObservationMatrix/:Normalize[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        ObservationMatrix[m/Total[Flatten[m]],opts]
        
ObservationMatrix/:Symmetrize[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        Module[ {},
            ObservationMatrix[(m+Transpose[m])/2]
        ]

ObservationMatrix/:LeftEigenvectors[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        ToLeftEigenvectors[m,All];
        
ObservationMatrix/:RightEigenvectors[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        ToRightEigenvectors[m,All];
        
ObservationMatrix/:FullEigensystem[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        ToFullEigensystem[m,All];
        
ObservationMatrix/:Eigenvalues[ObservationMatrix[m_?MatrixQ,opts:OptionsPattern[ObservationMatrix]]] := 
        Eigenvalues[m]

End[] (* End Private Context *)

EndPackage[]