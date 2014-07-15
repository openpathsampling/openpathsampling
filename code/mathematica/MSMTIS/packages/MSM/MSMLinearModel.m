(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

LinearModel::usage = "LinearModel[{eigenvalues, processvectors}, options] creates a LinearModel Object.";
LinearModelQ::usage = "";

EquilibriumMatrix::usage = "";

Options[LinearModel]            = {ShowSize->7, Lagtime -> 1, LagTimeUnit -> ""};

Begin["`Private`"] (* Begin Private Context *) 

LinearModel/:MakeBoxes[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]],StandardForm] := 
        If[ False,
            InterpretationBox[#,LinearModel[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalM]",RowBox[{Dimensions[m][[1]],"\[Rule]",Dimensions[m][[2]]}]],"(", SubscriptBox["\[Lambda]", "2"], "=", v[[2]],")"}]],
            InterpretationBox[#,LinearModel[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalM]",RowBox[{Dimensions[m][[1]],"\[Rule]",Dimensions[m][[2]]}]],"[",Lagtime[LinearModel[{v,m},opts]]*(LagTimeUnit/.{opts}),"]"}]]
        ]
        
LinearModel/:LinearModel[LinearModel[{v_?VectorQ,m_?MatrixQ}, opts1:OptionsPattern[LinearModel]], opts2:OptionsPattern[LinearModel]] := 
        LinearModel[{v,m}, UnionRules[opts1,opts2]]

LinearModel/:NormalForm[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        {{v,m}, opts}

LinearModel/:Normal[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        {v,m}

LinearModel/:Matrix[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        {v,m}

LinearModel/:Length[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        Length[m[[1]]]

LinearModel/:Dimensions[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        Dimensions[m]
        
LinearModelQ[x_] :=
    If[ Head[x]===LinearModel&&MatrixQ[x[[1,2]]]&&VectorQ[x[[1,2]]],
        True,
        False
    ];

LinearModel/:Chop[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        LinearModel[Chop[{v,m}],opts]

LinearModel/:Normalize[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        Module[ {},
            LinearModel[{v, EquilizeLeftEVSet[m]}, opts]
        ]
        
LinearModel/:EquilibriumMatrix[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        Transpose[m].m
        
LinearModel/:Eigenvalues[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        v
        
LinearModel/:Processes[LinearModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[LinearModel]]] := 
        m
        
LinearModel/:Lagtime[LinearModel[{v_?VectorQ, m_?MatrixQ},  opts : OptionsPattern[LinearModel]]] := Lagtime /. List[opts] /. {Lagtime -> 1}

LinearModel/:ChangeLagtime[LinearModel[{v_?VectorQ, m_?MatrixQ}, opts : OptionsPattern[LinearModel]], newLT_] :=
 		Module[ {},
     		LinearModel[{v^(newLT/Lagtime[LinearModel[{v,m},opts]]), m}, ReplaceRules[Lagtime -> newLT, opts]]
 		]
  
End[] (* End Private Context *)

EndPackage[]