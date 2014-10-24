(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

CorrelationModel::usage = "CorrelationModel[{eigenvalues, processvectors}] creates a CorrelationModel Object with eigenvalues and processvectors.";
CorrelationModelQ::usage = "CorrelationModelQ[system] returns true if system is a CorrelationModel.";

EquilibriumMatrix::usage = "EquilibriumMatrix[correlationmodel] returns the equilibrium matrix of the correlationmodel object.";

Options[CorrelationModel]            = {Lagtime -> 1};

Begin["`Private`"] (* Begin Private Context *) 

CorrelationModel/:MakeBoxes[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]],StandardForm] := 
        If[ False,
            InterpretationBox[#,CorrelationModel[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalM]",RowBox[{Dimensions[m][[1]],"\[Rule]",Dimensions[m][[2]]}]],"(", SubscriptBox["\[Lambda]", "2"], "=", v[[2]],")"}]],
            InterpretationBox[#,CorrelationModel[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalM]",RowBox[{Dimensions[m][[1]],"\[Rule]",Dimensions[m][[2]]}]],"[",Lagtime[CorrelationModel[{v,m},opts]],"]"}]]
        ]
        
CorrelationModel/:CorrelationModel[CorrelationModel[{v_?VectorQ,m_?MatrixQ}, opts1:OptionsPattern[CorrelationModel]], opts2:OptionsPattern[CorrelationModel]] := 
        CorrelationModel[{v,m}, UnionRules[opts1,opts2]]

CorrelationModel/:NormalForm[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        {{v,m}, opts}

CorrelationModel/:Normal[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        {v,m}

CorrelationModel/:Matrix[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        {v,m}

CorrelationModel/:Length[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        Length[m[[1]]]

CorrelationModel/:Dimensions[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        Dimensions[m]
        
CorrelationModelQ[x_] :=
    If[ Head[x]===CorrelationModel&&MatrixQ[x[[1,2]]]&&VectorQ[x[[1,2]]],
        True,
        False
    ];

CorrelationModel/:Chop[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        CorrelationModel[Chop[{v,m}],opts]

CorrelationModel/:Normalize[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        Module[ {},
            CorrelationModel[{v, EquilizeLeftEVSet[m]}, opts]
        ]
        
CorrelationModel/:EquilibriumMatrix[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        Transpose[m].m
        
CorrelationModel/:Eigenvalues[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        v
        
CorrelationModel/:Processes[CorrelationModel[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[CorrelationModel]]] := 
        m
        
CorrelationModel/:Lagtime[CorrelationModel[{v_?VectorQ, m_?MatrixQ},  opts : OptionsPattern[CorrelationModel]]] := Lagtime /. List[opts] /. {Lagtime -> 1}

CorrelationModel/:ChangeLagtime[CorrelationModel[{v_?VectorQ, m_?MatrixQ}, opts : OptionsPattern[CorrelationModel]], newLT_] :=
         Module[ {},
             CorrelationModel[{v^(newLT/Lagtime[CorrelationModel[{v,m},opts]]), m}, ReplaceRules[Lagtime -> newLT, opts]]
         ]
  
End[] (* End Private Context *)

EndPackage[]