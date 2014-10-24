(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

TransitionMatrix::usage = "TransitionMatrix[matrix, options] creates a TransitionMatrix Object.";
TransitionMatrixQ::usage = "TransitionMatrixQ[system] returns true if system is a valid TransitionMatrix object, false otherwise.";

ReversibleQ::usage = "Reversible[transitionmatrix] return true if transitionmatrix fulfills detailed balance, i.e. it is symmetric w.r.t. its stationary distribution.";

Options[TransitionMatrix]            = {Lagtime->1};

Begin["`Private`"] (* Begin Private Context *) 

(* TransitionMatrix functions *)

TransitionMatrix/:MakeBoxes[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],StandardForm] := 
        If[ $DefaultShowMatrixSize>Length[m],
            InterpretationBox[#,TransitionMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalT]",Length[m]],"[",Lagtime[TransitionMatrix[m,opts]],"]","\[Rule]",ToBoxes[MatrixForm[m]]}]],
            InterpretationBox[#,TransitionMatrix[m,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalT]",Length[m]],"[",Lagtime[TransitionMatrix[m,opts]],"]"}]]
        ]
TransitionMatrix/:TransitionMatrix[TransitionMatrix[m_?MatrixQ,opts1:OptionsPattern[TransitionMatrix]],opts2:OptionsPattern[TransitionMatrix]] := 
        TransitionMatrix[m,UnionRules[opts1,opts2]]
TransitionMatrixFromMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]] :=
    TransitionMatrix[Map[#/Apply[Plus,#]&,m]]
TransitionMatrix/:NormalForm[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        {m,opts}
TransitionMatrix/:Normal[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        m
TransitionMatrix/:Matrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        m
TransitionMatrix/:Lagtime[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Lagtime/.List[opts]/.{Lagtime->1}
TransitionMatrix/:MatrixForm[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        MatrixForm[m]
TransitionMatrix/:Length[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Length[m]
TransitionMatrix/:RowSumOneQ[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        RowSumOneQ[m]
TransitionMatrix/:MatrixPlot[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],optsMP:OptionsPattern[MatrixPlot]] := 
        MatrixPlot[m,optsMP,$TMPkgDefaultMatrixPlotOptions]
TransitionMatrix/:ArrayPlot[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix],optsAP:OptionsPattern[ArrayPlot]]] := 
        ArrayPlot[m,optsAP,$TMPkgDefaultArrayPlotOptions]
TransitionMatrixQ[x_] :=
    If[ Head[x]===TransitionMatrix&&MatrixQ[x[[1]]]&&RowSumOneQ[x[[1]]],
        True,
        False
    ];
TransitionMatrix/:Chop[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        TransitionMatrix[Chop[m],opts]
TransitionMatrix/:NormalizedQ[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        RowSumOneQ[m]
TransitionMatrix/:Normalize[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],optsN:OptionsPattern[Normalize]] := 
        Module[ {voidFnc},
            voidFnc = VoidFunction/.List[optsN]/.{VoidFunction->(#*0&)}/.{Automatic->Function[x,x*0],Zero->(#*0&),Infinity->(#/0&)};
            TransitionMatrix[Map[If[ Apply[Plus,#]>0,
                                     #/Apply[Plus,#],
                                     voidFnc[#]
                                 ]&, m], opts]
        ]
        
TransitionMatrix/:Symmetrize[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Module[ {eqDist,size,cM},
            size = Length[m];
            eqDist = Equilibrium[TransitionMatrix[m,opts]];
            cM = m;
            Table[cM[[ii]] = cM[[ii]]*eqDist[[ii]];,{ii,1,size}];
            TransitionMatrix[(cM+Transpose[cM])/2]
        ]

TransitionMatrix/:LeftEigenvectors[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        ToLeftEigenvectors[m,All];
TransitionMatrix/:RightEigenvectors[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        ToRightEigenvectors[m,All];
TransitionMatrix/:FullEigensystem[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        ToFullEigensystem[m,All];
TransitionMatrix/:Eigenvalues[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Eigenvalues[m]
        
TransitionMatrix/:ChangeLagtime[TransitionMatrix[m_?MatrixQ, opts : OptionsPattern[TransitionMatrix]], newLT_] :=
  		Module[{evSystem, lagtime},
  			evSystem = ToFullEigensystem[m, All];
  			lagtime = Lagtime[TransitionMatrix[m, opts]];
  			evSystem[[1]] = evSystem[[1]]^(newLT/lagtime);
  			TransitionMatrix[FromFullEigensystem[evSystem], ReplaceRules[Lagtime -> newLT, opts]]
  		]

TransitionMatrix/:Equilibrium[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        #/Apply[Plus,#]&[Eigenvectors[Transpose[m],Min[20,Length[m]]][[1]]]

TransitionMatrix /: 
 ReversibleQ[TransitionMatrix[m_?MatrixQ, opts___]] := Module[{eq},
  eq = Equilibrium[TransitionMatrix[m]];
  SymmetricQ@(DiagonalMatrix[eq].m)
  ]

End[] (* End Private Context *)

EndPackage[]