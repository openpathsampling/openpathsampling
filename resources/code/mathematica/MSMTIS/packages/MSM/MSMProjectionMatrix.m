(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ProjectionMatrix::usage = "ProjectionMatrix[matrix, options] creates a ProjectionMatrix Object.";
ProjectionMatrixQ::usage = "";

Options[ProjectionMatrix]               = {};

Begin["`Private`"] (* Begin Private Context *) 

ProjectionMatrix/:MakeBoxes[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]],StandardForm] := 
        InterpretationBox[#,ProjectionMatrix[m,opts]]&[SubscriptBox["\[DoubleStruckCapitalX]",RowBox[{Dimensions[m][[2]],"\[Rule]",Dimensions[m][[1]]}]]]
ProjectionMatrix/:Normal[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]]] := 
        m
ProjectionMatrixQ[x_] :=
    If[ Head[x]==ProjectionMatrix&&MatrixQ[x[[1]]],
        True,
        False
    ];
ProjectionMatrix/:NormalForm[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]]] := 
        {m,opts}
ProjectionMatrix/:Normal[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]]] := 
        m
ProjectionMatrix/:Matrix[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]]] := 
        m
ProjectionMatrix/:Reverse[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]]] := 
        ProjectionMatrix[Transpose[m],opts]
ProjectionMatrix/:MatrixForm[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]]] := 
        MatrixForm[m]
ProjectionMatrix/:Length[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]]] := 
        Reverse[Dimensions[m]]
ProjectionMatrix/:MatrixPlot[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]],optsMP:OptionsPattern[MatrixPlot]] := 
        MatrixPlot[m,optsMP,Frame->True,FrameTicks->{Automatic,Automatic,None,None},ImageSize->320,ColorFunction->"TemperatureMap"]
ProjectionMatrix/:InvolvedStates[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]]] :=
        Drop[ArrayRules[Transpose[Normal[m]]][[All,1,1]],-1]
ProjectionMatrix/:Reorder[ProjectionMatrix[m_?MatrixQ,opts:OptionsPattern[ProjectionMatrix]], or_?VectorQ] :=
        ProjectionMatrix[m[[or]],opts]

End[] (* End Private Context *)

EndPackage[]