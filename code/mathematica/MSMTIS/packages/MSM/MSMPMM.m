(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

PMM::usage = "PMM[{eigenvalues, processvectors}] creates a PMM Object with eigenvalues and processvectors.";
PMMQ::usage = "PMMQ[system] returns true if system is a PMM.";

EquilibriumMatrix::usage = "EquilibriumMatrix[pmm] returns the equilibrium matrix of the pmm object.";

Options[PMM]            = {Lagtime -> 1};

Begin["`Private`"] (* Begin Private Context *) 

PMM/:MakeBoxes[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]],StandardForm] := 
        If[ False,
            InterpretationBox[#,PMM[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalP]",RowBox[{Dimensions[m][[1]],"\[Rule]",Dimensions[m][[2]]}]],"(", SubscriptBox["\[Lambda]", "2"], "=", v[[2]],")"}]],
            InterpretationBox[#,PMM[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalP]",RowBox[{Dimensions[m][[1]],"\[Rule]",Dimensions[m][[2]]}]],"[",Lagtime[PMM[{v,m},opts]],"]"}]]
        ]
        
PMM/:PMM[PMM[{v_?VectorQ,m_?MatrixQ}, opts1:OptionsPattern[PMM]], opts2:OptionsPattern[PMM]] := 
        PMM[{v,m}, UnionRules[opts1,opts2]]

PMM/:NormalForm[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        {{v,m}, opts}

PMM/:Normal[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        {v,m}

PMM/:Matrix[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        {v,m}

PMM/:Length[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        Length[m[[1]]]

PMM/:Dimensions[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        Dimensions[m]
        
PMMQ[x_] :=
    If[ Head[x]===PMM&&MatrixQ[x[[1,2]]]&&VectorQ[x[[1,2]]],
        True,
        False
    ];

PMM/:Chop[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        PMM[Chop[{v,m}],opts]

PMM/:Normalize[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        Module[ {},
            PMM[{v, EquilizeLeftEVSet[m]}, opts]
        ]
        
PMM/:EquilibriumMatrix[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        Transpose[m].m
        
PMM/:Eigenvalues[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        v
        
PMM/:Processes[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := 
        m
        
PMM/:Lagtime[PMM[{v_?VectorQ, m_?MatrixQ},  opts : OptionsPattern[PMM]]] := Lagtime /. List[opts] /. {Lagtime -> 1}

PMM/:ChangeLagtime[PMM[{v_?VectorQ, m_?MatrixQ}, opts : OptionsPattern[PMM]], newLT_] :=
         Module[ {},
             PMM[{v^(newLT/Lagtime[PMM[{v,m},opts]]), m}, ReplaceRules[Lagtime -> newLT, opts]]
         ]
  
End[] (* End Private Context *)

EndPackage[]