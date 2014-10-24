(* ::Package:: *)

(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

OOM::usage = "OOM[matrix, options] creates a OOM Object.";
OOMQ::usage = "OOMQ[system] return true if system is a valid OOM object, false otherwise.";

CountObservables::usage = "CountObservables[oom] returns the number of observables in the OOM.";
CountProcesses::usage = "CountProcesses[oom] returns the number of processes in the OOM, i.e. the size of the tau matrices.";

GetTauZero::usage = "GetTauZero[oom] returns the tauzero matrix, i.e. the sum of all tau matrices.";
GetListOfTauMatrices::usage = "GetListOfTauMatrices[oom] returns the list of tau matrices";
GetTauMatrix::usage = "GetTauMatrix[oom, index] returns the tau matrix for observation with the index given.";
GetInitialDistribution::usage = "GetInitialDistribution[oom] returns the initial distribution specified in the OOM/";

ComputeTraceProbability::usage = "ComputeTraceProbability[oom, {initialdistribution,} trace] returns the probability to observe the specified trace gien the OOM. If initialdistribution is not specified the native initial distribution of the OOM is used.";
GetRMatrix::usage = "GetRMatrix[oom] returns the R-matrix of the OOM, i.e. the matrix of eigenvectors of tauzero.";

SwitchToProbabilityRepresentation::usage = "SwitchToProbabilityRepresentation[oom] returns a similar OOM, where the tau matrices are given in the probability representation as e.g. estimated from an HMM.";

Similarity::usage = "Similarity[oom1, oom2] computes a similarity norm of two ooms in the space of tau matrices of oom1. The algorithm estimates a similarity matrix S that minimizes the difference norm and returns {norm, smatrix}.";

(*
Member Variables of OOMs

"StateSpace" -> ListOfSizeM
"TimeStep" -> Number
"TauMatries" -> ListOfKMatricesOfSizeMByM
"ObservableSpace" -> ListOfSizeK
"InitialDistribution" -> ListOfSizeM
*)



Options[OOM]            = {};

Begin["`Private`"] (* Begin Private Context *) 

(* OOM functions *)

OOM/:MakeBoxes[OOM[d_?ListQ,opts:OptionsPattern[OOM]],StandardForm] := 
            InterpretationBox[#,OOM[d,opts]]&[RowBox[{"\[ScriptCapitalM]","[",CountProcesses[OOM[d,opts]],"]","\[Rule]",SuperscriptBox["\[DoubleStruckCapitalR]",CountObservables[OOM[d,opts]]]}]];
OOM/:CountObservables[OOM[d_?ListQ,opts:OptionsPattern[OOM]]] := Length[GetListOfTauMatrices[OOM[d,opts]]];
OOM/:CountProcesses[OOM[d_?ListQ,opts:OptionsPattern[OOM]]] := Dimensions[GetListOfTauMatrices[OOM[d,opts]]][[2]];

OOMQ[x_] :=
    If[ Head[x]===OOM&&ListQ[x[[1]]],
        True,
        False
    ];

OOM/:GetTauZero[OOM[d_?ListQ,opts:OptionsPattern[OOM]]] := 
        Total["TauMatrices"/.d]

OOM/:GetListOfTauMatrices[OOM[d_?ListQ,opts:OptionsPattern[OOM]]] := 
        "TauMatrices"/.d

OOM/:GetTauMatrix[OOM[d_?ListQ,opts:OptionsPattern[OOM]],nn_] := 
        Take["TauMatrices"/.d,nn]

OOM/:GetRMatrix[OOM[d_?ListQ,opts:OptionsPattern[OOM]]] := 
        Eigenvectors[GetTauZero[OOM[d,opts]]]

OOM/:GetInitialDistribution[OOM[d_?ListQ,opts:OptionsPattern[OOM]]] := 
        "InitialDistribution"/.d

OOM/:ComputeTraceProbability[OOM[d_?ListQ,opts:OptionsPattern[OOM]],trace_] := 
        Module[ {initialDistribution},
            initialDistribution = GetInitialDistribution[OOM[d,opts]];
            ComputeTraceProbability[OOM[d,opts],initialDistribution,trace]
        ]

OOM/:ComputeTraceProbability[OOM[d_?ListQ,opts:OptionsPattern[OOM]],init_,trace_] := 
        Module[ {tauMatrices,tauZero,observationTypes,replaceMatrix,rMatrix,eigenvalues,tauInfinity,size,tauNull},
            tauMatrices = GetListOfTauMatrices[OOM[d,opts]];
            tauZero = GetTauZero[OOM[d,opts]];
            observationTypes = Union[trace];
            eigenvalues = Eigenvalues[tauZero];
            rMatrix = GetRMatrix[OOM[d,opts]];
            tauInfinity = Transpose[rMatrix][[All,{1}]].Inverse[Transpose[rMatrix]][[{1}]];
            size = CountProcesses[OOM[d,opts]];
            tauNull = IdentityMatrix[size];
            replaceMatrix = MapIndexed[#->
                    Switch[#,
                        All,
                        tauZero,
                        _?IntegerQ,
                        tauMatrices[[#]],
                        {All,0},
                        tauNull,
                        {All,_?IntegerQ},
                        MatrixPower[tauZero,#[[2]]],
                        {All,Infinity},
                        tauInfinity,
                        {{_?IntegerQ..},0},
                        tauNull,
                        {{_?IntegerQ..},_?IntegerQ},
                        MatrixPower[Total[tauMatrices[[#[[1]]]]],#[[2]]],
                        {{_?IntegerQ..},{_?IntegerQ,_?IntegerQ}},
                        MatrixRangeSum[Total[tauMatrices[[#[[1]]]]],#[[2]]],
                        {{_?IntegerQ..},{_?IntegerQ,Infinity}},
                        MatrixRangeSum[Total[tauMatrices[[#[[1]]]]],#[[2]]],
                        {_?IntegerQ..},
                        Total[tauMatrices[[#]]],
                        _,
                        0
                    ]&,
                observationTypes];
            AppendTo[tauMatrices, Total[tauMatrices]];
            Total[init.Apply[Dot,trace/.replaceMatrix]]
        ]
        
OOM/:Eigenvalues[OOM[d_?ListQ,opts:OptionsPattern[OOM]]] := 
        Eigenvalues[GetTauZero[OOM[d,opts]]]
        
OOM/:Equilibrium[OOM[d_?ListQ,opts:OptionsPattern[OOM]]] := 
        Module[ {},
            Table[
                ComputeTraceProbability[OOM[d,opts],{{All,Infinity},nn}],
                {nn,1,CountObservables[OOM[d,opts]]}
            ]
        ]
                
OOM/:Anything[OOM[d_?ListQ,opts:OptionsPattern[OOM]],st1_,st2_] := 
        Module[ {tauMatrices,tauZero,rMatrix,eigenvalues,tauInfinity,size,tauNull,init,sizeObservable,spaceComplement,matrixComplement},
            init = GetInitialDistribution[OOM[d,opts]];
            tauMatrices = GetListOfTauMatrices[OOM[d,opts]];
            tauZero = GetTauZero[OOM[d,opts]];
            eigenvalues = Eigenvalues[tauZero];
            rMatrix = GetRMatrix[OOM[d,opts]];
            tauInfinity = Transpose[rMatrix][[All,{1}]].Inverse[Transpose[rMatrix]][[{1}]];
            size = CountProcesses[OOM[d,opts]];
            tauNull = IdentityMatrix[size];
            sizeObservable = CountObservables[OOM[d,opts]];
            spaceComplement = Complement[Range[sizeObservable],{st2}];
            matrixComplement = Total[tauMatrices[[spaceComplement]]];
            Total[init.tauMatrices[[st1]].MatrixRangeMeanSum[matrixComplement,{0,Infinity}].tauMatrices[[st2]]]/ComputeTraceProbability[OOM[d,opts],{st1}]
        ]
       
OOM/:SwitchToProbabilityRepresentation[OOM[d_?ListQ,opts___]] :=
        Module[ {cTEx,cMEx},
            cTEx = Table[ComputeTraceProbability[OOM[d,opts],{sti,stk,stj}],
            {stk,1,3},{sti,1,3},{stj,1,3}];
            cMEx = Table[ComputeTraceProbability[OOM[d,opts],{sti,stj}],
            {sti,1,3},{stj,1,3}];
            OOM[{"TauMatrices"->Map[Inverse[cMEx].#&,cTEx],"InitialDistribution"->Map[Total,cMEx]}]
        ]

OOM/:Similarity[OOM[d1_?ListQ, opts1 : OptionsPattern[OOM]], OOM[d2_?ListQ, opts2 : OptionsPattern[OOM]]] :=
        Module[ {size1, size2, sizeP1, sizeP2, sizeP, sMatrix, tMats1, tMats2, sVars, func, totFunc, result},
            size1 = CountProcesses[OOM[d1, opts1]];
            size2 = CountProcesses[OOM[d2, opts2]];
            sizeP1 = CountObservables[OOM[d1, opts1]];
            sizeP2 = CountObservables[OOM[d2, opts2]];
            If[ sizeP1 == sizeP2,
                sizeP = sizeP1;
                sMatrix = 
                 Table[V[s, st1, st2, 3], {st1, 1, size1}, {st2, 1, size2}];
                Table[sMatrix[[1, st2]] += 1 - Total[sMatrix[[All, st2]]], {st2, 1,
                   size2}];
                sMatrix = Transpose[sMatrix];
                tMats1 = GetListOfTauMatrices[OOM[d1, opts1]];
                tMats2 = GetListOfTauMatrices[OOM[d2, opts2]];
                sVars = 
                 Select[Flatten[MapIndexed[{#2, #} &, sMatrix, {2}], 1], 
                   1 != #[[1, 2]] &][[All, 2]];
                func = 
                 Table[sMatrix.tMats1[[st]] - tMats2[[st]].sMatrix, {st, 1, sizeP}];
                totFunc = 
                 Total[Table[
                    Total[Chop[Flatten[func[[st]]*func[[st]]]]] // Expand, {st, 1, 
                     sizeP}]] // Expand;
                result = Chop[NMinimize[totFunc, sVars]];
                {result[[1]], sMatrix /. result[[2]]},
                Print["Error: Not equal number of observations!"];
            ]
        ]

End[] (* End Private Context *)

EndPackage[]
