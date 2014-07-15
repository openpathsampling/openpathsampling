(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ClusterMatrix::usage = "ClusterMatrix[matrix, options] creates a ClusterMatrix Object.";
ClusterMatrixQ::usage = "ClusterMatrixQ[system] returns true of system is a ClusterMatrix object.";

CrispClusterMatrix::usage = "CrispClusterMatrix[rules] creates a crisp ClusterMatrix from the rules given as a list of groups of states. If the the option Complete->True is given all not specified states are mapped one-to-one.";
RemoveClusterMatrix::usage = "RemoveClusterMatrix[states, size] creates a ClusterMatrix that removes all specified states of a total list of size states.";
RemoveEmptyStatesClusterMatrix::usage = "RemoveEmptyStatesClusterMatrix[system] returns a ClusterMatrix that removes all not visited states from system.";
IdentityClusterMatrix::usage = "IdentityClusterMatrix[size] returns a ClusterMatrix of size that does nothing.";
ToCrispRules::usage = "ToCrispRules[clustermatrix] returns a list of rules that assigns each state the most probable cluster given in ClusterMatrix.";
PivotClusterMatrix::usage = "PivotClusterMatrix[pivots] returns a ClusterMatrix that rearragnes states in the order specified in pivots.";

Crispify::usage = "Crispify[clustermatrix] returns a crisp version of the ClusterMatrix.";
InvolvedStates::usage = "InvolvedStates[clustermatrix] returns a list of initial states that are used in ClusterMatrix.";
ClosedEnd::usage = "ClosedEnd is an option for ToCountMatrix and ToCountTensor that ensures that so many states and the beginning or end of a timeseries are removed to ensure that the related graph is strongly connected.";

GetClusteringFromCenters::usage = "GetClusteringFromCenters[statecenters, clustercenters] creates a ClusterMatrix assigning each state in statecenters to the nearest cluster center in clustercenters.";

Reweight::usage = "Reweight[clustermatrix, equilibrium] returns the ClusterMatrix where each membership function is reweighted with the equilibrium specified.";
Reorder::usage = "Reorder[clustermatrix, pivots] returns a ClusterMatrix with reordered indices of the cluster centers. Cluster pivots[[1]] will be new cluster 1.";

Complete::usage = "Complete is an option for CrispClusterMatrix that specifies whether not mentioned states should be kept one-to-one or incomplete clusterings are allowed."

Options[ClusterMatrix]               = {};
Options[CrispClusterMatrix]          = {Complete->False,States->0};

Begin["`Private`"] (* Begin Private Context *) 

ClusterMatrix/:MakeBoxes[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]],StandardForm] := 
        InterpretationBox[#,ClusterMatrix[m,opts]]&[SubscriptBox["\[DoubleStruckCapitalL]",RowBox[{Dimensions[m][[2]],"\[Rule]",Dimensions[m][[1]]}]]]
ClusterMatrix/:Normal[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        m
ClusterMatrixQ[x_] :=
    If[ Head[x]==ClusterMatrix&&MatrixQ[x[[1]]],
        True,
        False
    ];
ClusterMatrix/:NormalForm[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        {m,opts}
ClusterMatrix/:Normal[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        m
ClusterMatrix/:Matrix[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        m
ClusterMatrix/:Inverse[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        ClusterMatrix[Transpose[m],opts]
ClusterMatrix/:MatrixForm[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        MatrixForm[m]
ClusterMatrix/:Length[ClusterMatrix[m_?MatrixQ],opts:OptionsPattern[ClusterMatrix]] := 
        Reverse[Dimensions[m]]
ClusterMatrix/:MatrixPlot[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]],optsMP:OptionsPattern[MatrixPlot]] := 
        MatrixPlot[m,optsMP,Frame->True,FrameTicks->{Automatic,Automatic,None,None},ImageSize->320,ColorFunction->"TemperatureMap"]
ClusterMatrix/:InvolvedStates[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] :=
        Drop[ArrayRules[Transpose[Normal[m]]][[All,1,1]],-1]
        
IdentityClusterMatrix[st_] :=
    ClusterMatrix[IdentityMatrix[st]]

RemoveClusterMatrix[states_,size_] :=
    ClusterMatrix[IdentityMatrix[size][[Complement[Range[size],states]]]]

ClusterMatrix/:ToCrispRules[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        Module[ {},
            Table[nn->#&[Select[Transpose[{Range[Length[m]],m[[All,nn]]}],#[[2]]==Max[m[[All,nn]]]&][[1,1]]],{nn,1,Length[m[[1]]]}]
        ]
ClusterMatrix/:Crispify[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        Module[ {},
            ClusterMatrix[Normal[SparseArray[Table[{#,nn}->1&[Select[Transpose[{Range[Length[m]],m[[All,nn]]}],#[[2]]==Max[m[[All,nn]]]&][[1,1]]],{nn,1,Length[m[[1]]]}]]]]
        ]
ClusterMatrix/:Reweight[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]],eq_?VectorQ] :=
        ClusterMatrix[Map[#*eq/(#.eq)&,m],opts]

ClusterMatrix/:Reorder[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]], or_?VectorQ] :=
        ClusterMatrix[m[[or]],opts]
PivotClusterMatrix[order_] :=
    Module[ {},
        ClusterMatrix[IdentityMatrix[Length[order]][[order]]]
    ]
GetClusteringFromCenters[stateCenterList_,clusteredCenterList_] :=
    Module[ {nearF,clustAssignments},
        nearF = Nearest[clusteredCenterList->Automatic];
        clustAssignments = Map[nearF[#][[1]]&,stateCenterList];
        Inverse[CrispClusterMatrix[clustAssignments,States->Length[clusteredCenterList]]]
    ]

CrispClusterMatrix[rules_,opts:OptionsPattern[CrispClusterMatrix]] :=
    Module[ {nStates,thisRules},
        nStates = Max[OptionValue[CrispClusterMatrix,States],Max[Flatten[rules]]];
        thisRules = If[ OptionValue[Complete],
                        Join[rules,Map[{#}&,Complement[Range[nStates], Flatten[rules]]]],
                        rules
                    ];
        ClusterMatrix[Map[Table[0,{nStates}]+Apply[Plus,IdentityMatrix[nStates][[Flatten[{#}]]]]&,thisRules]]
    ];


End[] (* End Private Context *)

EndPackage[]