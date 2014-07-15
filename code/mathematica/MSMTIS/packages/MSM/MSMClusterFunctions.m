(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

KMeansClustering::usage = "KMeansClustering[data, n] runs a K-Means Algorithn on data estimating n cluster centers.";

Cluster::usage = "Cluster[clustermatrix, system] clusters system using the ClusterMatrix specified.";
HasEmptyStatesQ::usage = "HasEmptyStatesQ[system] return true if system contains non-visited states, zero otherwise.";
PCCASoftClusterMatrix::usage = "PCCASoftClusterMatrix[transitionmatrix, n] returns a fuzzy ClusterMatrix based on the PCCA decomposition for n states using right eigenvectors of TransitionMatrix.";
PCCACrispClusterMatrix::usage = "PCCACrispClusterMatrix[transitionmatrix, n] returns a crisp ClusterMatrix based on the PCCA decomposition for n states using right eigenvectors of TransitionMatrix.";
Points::usage = "Points is an option for PCCAAlgorithm to specify the corners when Corners mode Points is chosen.";
RemoveVoid::usage = "RemoveVoid[timeseries, void] return changed timeseries with where all states in void are replace by the last allowed (not void) state.";
Projection::usage = "Projection[clustermatrix, equilibrium, vectors] projects a list of vectors onto the clustering given in ClusterMatrix assuming a euklidean scalar product.";
ToWeightedProjection::usage = "ToWeightedProjection[clustermatrix, equilibrium] returns a function that projects a vector onto the statespace given in the clustering assuming an equilibrium weighted scalar product. ";
ToProjection::usage = "ToProjection[clustermatrix] returns a function that projects a vector onto the statespace given in the clustering assuming a euklidian scalar product.";
WeightedProjection::usage = "WeightedProjection[clustermatrix, equilibrium, vectors] projects a list of vectors onto the clustering given in ClusterMatrix assuming an equilibriumweighted scalar product.";
ToMapping::usage = "ToMapping[clustermatrix] returns a list of most likely clustercenters for each initial state in ClusterMatrix.";
ToNormalizedBasis::usage = "ToNormalizedBasis[lustermatrix] returns a list of normalized (sum of one) membership vectors.";
Expand::usgae = "ExpandClustering[clustermatrix, system] expands system back into the full size using ClusterMatrix.";
PCCAAlgorithm::usage = "PCCAAlgorithm[righteigenvectors, n] runs a PCCA algorithm on righteigenvectors with n states and returns a clustering matrix.";

ClusterAtEquilibrium::usage = "ClusterAtEquilibrium[clustermatrix, vector, equilibrium] projects vector using clustermatrix assuming a scalarproduct based on equilibrium.";
RemoveEmptyStatesClusterMatrix::usage = "RemoveEmptyStatesClusterMatrix[system] creates a ClusterMatrix that removes all states that are not visited in system";

GetInverse::usage = "GetInverse is an option for PCCAAlgorithm to specify if the inverse transformation should be returned.";

Options[PCCASoftClusterMatrix]       = {Corners->Automatic, GetInverse->False};

Begin["`Private`"] (* Begin Private Context *) 

PMM/:Cluster[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]], PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := Module[{},
  PMM[{v, m.Transpose[c]}, opts]
  ]

PMM/:Expand[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]], PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]]] := Module[{},
  PMM[{v, m.c}, opts]
  ]
  
PMM/:Take[PMM[{v_?VectorQ,m_?MatrixQ},opts:OptionsPattern[PMM]], part_] := Module[{},
  PMM[{Take[v, part], Take[m, part]}, opts]
]

ClusterMatrix/:Projection[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]],v_?VectorQ] := 
        Module[ {fnc},
            fnc = ToProjection[ClusterMatrix[m,opts]];
            Map[fnc,v]
        ]

ClusterMatrix/:ToProjection[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := 
        Module[ {clNorm},
            clNorm = N[Map[Stochastize,Normal[m]]];
            Function[v, clNorm.v.m]
        ]
ClusterMatrix/:ToWeightedProjection[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]],eq_?VectorQ] := 
        Module[ {clNorm},
            clNorm = N[Map[Stochastize,Map[#*eq&,Normal[m]]]];
            Function[v, clNorm.v.m]
        ]

ClusterMatrix/:WeightedProjection[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]],eq_?VectorQ,v_?MatrixQ] := 
        Module[ {fnc},
            fnc = ToWeightedProjection[ClusterMatrix[m,opts],eq];
            Map[fnc,v]
        ]

ClusterMatrix/:ToMapping[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := Module[ {},
                                                                                              Table[
                                                                                                  Select[Transpose[{Range[Length[m]],m[[All,nn]]}],#[[2]]==Max[m[[All,nn]]]&][[1,1]]
                                                                                              ,
                                                                                                  {nn,1,Length[m[[1]]]}]
                                                                                          ]

ClusterMatrix/:ToNormalizedBasis[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]]] := Module[ {},
                                                                                                      N[Normal[Map[#/Total[#]&,m]]]
                                                                                                  ]

TransitionMatrix/:RemoveEmptyStatesClusterMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Module[ {cM},
            cM = Select[Table[{nn,Apply[Plus,m[[nn]]]},{nn,1,Length[m]}],#[[2]]>0&][[All,1]];
            ClusterMatrix[SparseArray[MapIndexed[{#2[[1]],#}->1&,cM],{Length[cM],Length[m]}]]
        ]

TransitionMatrix/:Expand[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]],TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],optsCL___] := 
        Module[ {eqDist,size},
            size = Length[m];
            eqDist = Equilibrium/.List[optsCL]/.{Equilibrium->{}};
            If[ Length[eqDist]!=size,
                eqDist = Equilibrium[TransitionMatrix[m,opts]];
            ];
            Print[DiagonalMatrix[eqDist].m];
            Normalize[TransitionMatrix[Transpose[c].DiagonalMatrix[eqDist].m.c],RemoveRules[States,opts]]
        ];
CorrelationMatrix/:Expand[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]],CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]],optsCL___] := 
        Module[ {lt},
        	lt = Lagtime[CorrelationMatrix[m,opts]];
            CorrelationMatrix[Transpose[c].m.c, Lagtime->lt]
        ];
CorrelationMatrix/:Cluster[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]],CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]],optsCL___] := 
        Module[ {lt},
        	lt = Lagtime[CorrelationMatrix[m,opts]];
            CorrelationMatrix[c.m.Transpose[c],Lagtime->lt]
        ];
RateMatrix/:Cluster[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]],RateMatrix[m_?MatrixQ,opts:OptionsPattern[RateMatrix]]] := 
        Module[ {size},
            size = Length[m];
            RateMatrix[c.m.Transpose[c],RemoveRules[States,opts]]
        ];
ClusterMatrix/:Cluster[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]],v_?VectorQ] := 
        Module[ {rules},
            rules = ArrayRules[c];
            Table[Map[{#[[1,1]],#[[2]]}&,Select[rules,#[[1,2]]==v[[nn]]&]],{nn,1,Length[v]}]
        ];
CountMatrix/:Cluster[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]],CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        CountMatrix[c.m.Transpose[c],RemoveRules[States,opts]]
TimeSeries/:Cluster[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]],TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]]] := 
        Module[ {rules},
            rules = ToCrispRules[ClusterMatrix[c,optsCM]];
            TimeSeries[v/.rules]
        ];
CountMatrix/:HasEmptyStatesQ[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        If[ Apply[Times,Map[Apply[Plus,#]&,m]]==0,
            True,
            False
        ]
TransitionMatrix/:HasEmptyStatesQ[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        If[ Apply[Times,Map[Apply[Plus,#]&,m]]==0,
            True,
            False
        ]
        
OOM /: Cluster[ClusterMatrix[cM_?MatrixQ, cMopts : OptionsPattern[ClusterMatrix]], OOM[d_?ListQ, opts : OptionsPattern[OOM]]] := 
         Module[ {clusteredTauMatrices, clusteredInitialDistribution, init},
             init = GetInitialDistribution[OOM[d, opts]];
             tauMatrices = GetListOfTauMatrices[OOM[d, opts]];
             clusteredTauMatrices = cM.tauMatrices;
             clusteredInitialDistribution = init;
             OOM[{"TauMatrices" -> clusteredTauMatrices, 
             "InitialDistribution" -> clusteredInitialDistribution}]
         ]        
        
TimeSeries/:RemoveVoid[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]],states_] := 
        Module[ {vR,state,vC},
(*            vVoidStates=(VoidStates/.List[UnionRules[optsRV,opts]])/.Options[TimeSeries]; *)
            (* Standard Method *)
            vR = v/.Thread[states->0];
            state = 0;
            vC = Table[
                If[ vR[[nn]]>0,
                    state = vR[[nn]],
                    state
                ]
            ,{nn,1,Length[v]}];
            TimeSeries[Select[vC,#!=0&]]
        ];
        
ClusterAtEquilibrium[ClusterMatrix[m_?MatrixQ,opts:OptionsPattern[ClusterMatrix]],v_?VectorQ,eq_?VectorQ] :=
    Module[ {},
        (m.(eq*v))/(m.eq)
    ]
    
CountMatrix/:RemoveEmptyStatesClusterMatrix[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]] := 
        Module[ {cM,mOld,last,emptyLine},
            mOld = Normal[m];
            last = Length[m]+1;
            cM = Range[Length[m]];
            emptyLine = Table[0,{Length[m]}];
            While[Length[cM]<last,
                last = Length[cM];
                Map[(mOld[[#,All]] = emptyLine;
                     mOld[[All,#]] = emptyLine)&,Complement[Range[Length[m]],cM]];
                cM = Select[Table[{nn,Apply[Plus,mOld[[nn]]],Apply[Plus,mOld[[All,nn]]]},{nn,1,Length[m]}],#[[2]]>0&&#[[3]]>0&][[All,1]];                
            ];
            ClusterMatrix[SparseArray[MapIndexed[{#2[[1]],#}->1&,cM],{Length[cM],Length[m]}]]
        ]



TransitionMatrix/:Cluster[ClusterMatrix[c_?MatrixQ,optsCM:OptionsPattern[ClusterMatrix]],TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],optsCL___] := 
        Module[ {eqDist,size,cM},
            size = Length[m];
            eqDist = Equilibrium/.List[optsCL]/.{Equilibrium->{}};
            If[ Length[eqDist]!=size,
                eqDist = Equilibrium[TransitionMatrix[m,opts]];
            ];
            cM = m;
            Table[cM[[ii]] = cM[[ii]]*eqDist[[ii]];,{ii,1,size}];
            TransitionMatrix[Map[If[ Apply[Plus,#]>0,
                                     #/Apply[Plus,#],
                                     Table[0,{Length[c]}]
                                 ]&,c.cM.Transpose[c]],opts]
        ];

(* PCCA+ *)

TransitionMatrix/:PCCASoftClusterMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],level_,optsPC___] := 
        Module[ {},
            PCCAAlgorithm[Map[#/AbsMax[#]&,Eigenvectors[N[m]][[Range[level]]]], level, optsPC]
        ];

TransitionMatrix/:PCCACrispClusterMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],level_] := 
        Module[ {cornerList,clustering},
            {cornerList,clustering} = PCCASoftClusterMatrix[TransitionMatrix[m,opts],level];
            {cornerList,Crispify[clustering]}
        ];
        
PCCAAlgorithm[pV_, level_, optsPC___] :=
    Module[ {size, maxDist, distTable, cornerList, eVecsNorm, list2, ll, iMatrix, mixingPoints, result, invMixingPointList, eVecs, ap},
        eVecs = Transpose[pV];
        size = Length[eVecs[[1]]];
        cornerList = {1};
        ap = Subset /. {optsPC} /. {Subset -> Range[Length[eVecs]]};
        Which[
        (Corners /. {optsPC} /. Options[PCCASoftClusterMatrix]) === Automatic,
	        distTable = Flatten[Table[{qq, pp, Abs[Norm[eVecs[[ap[[pp]], {2}]] - eVecs[[ap[[qq]], {2}]]]]}, {pp, 2, Length[ap]}, {qq, 1, pp - 1}], 1];
	        maxDist = Max[distTable[[All, 3]]];
    	    cornerList = Select[distTable, #[[3]] == maxDist &][[1, {1, 2}]];
	        For[ll = 3, ll <= level, ll++, 
    		    distTable = Table[{pp, Apply[Plus, Map[Norm[eVecs[[pp, ll]] - eVecs[[#, ll]]] &, cornerList]]}, {pp, Complement[ap, cornerList]}];
	        	maxDist = Max[distTable[[All, 2]]];
    		    AppendTo[cornerList, Select[distTable, #[[2]] == maxDist &][[1, 1]]];
	        ];
        	eVecsNorm = Map[# - eVecs[[cornerList[[1]]]] &, eVecs][[All, Range[2, level]]];
	        list2 = eVecsNorm.Inverse[Drop[eVecsNorm[[cornerList]], 1]];
    	    list2 = Chop[Join[Transpose[list2], {eVecs[[All, 1]] - Map[Apply[Plus, #] &, list2]}]];
        	If[ GetInverse /. {optsPC} /. {GetInverse -> False},
            	iMatrix = Join[-{eVecs[[cornerList[[1]], Range[2, level]]]}, Inverse[Drop[eVecsNorm[[cornerList]], 1]]];
	            iMatrix = Transpose[Append[Transpose[iMatrix], Map[(1 - Total[#]) &, iMatrix]]];
            	result = {cornerList, ClusterMatrix[list2], list2.PseudoInverse[Transpose[eVecs]], Transpose[eVecs]};,
	            result = {cornerList, ClusterMatrix[list2]};
    	    ];,
        Head[Corners /. {optsPC} /. Options[PCCASoftClusterMatrix]] === List,
            cornerList = Corners /. {optsPC} /. Options[PCCASoftClusterMatrix];
            eVecsNorm = Map[# - eVecs[[cornerList[[1]]]] &, eVecs][[All, Range[2, level]]];
            list2 = eVecsNorm.Inverse[Drop[eVecsNorm[[cornerList]], 1]];
            list2 = Chop[Join[Transpose[list2], {eVecs[[All, 1]] - Map[Apply[Plus, #] &, list2]}]];
            If[ GetInverse /. {optsPC} /. {GetInverse -> False},
                 iMatrix = Join[-{eVecs[[cornerList[[1]], Range[2, level]]]}, Inverse[Drop[eVecsNorm[[cornerList]], 1]]];
                 iMatrix = Transpose[Append[Transpose[iMatrix], Map[(1 - Total[#]) &, iMatrix]]];
                 result = {cornerList, ClusterMatrix[list2], list2.PseudoInverse[Transpose[eVecs]], Transpose[eVecs]};,
                 result = {cornerList, ClusterMatrix[list2]};
             ];,
         (Corners /. {optsPC} /. Options[PCCASoftClusterMatrix]) === Points,
             mixingPoints = Points /. {optsPC} /. Options[PCCASoftClusterMatrix];
             eVecsNorm = Transpose[Table[eVecs[[All, nn]] - eVecs[[All, 1]]*mixingPoints[[1, nn - 1]], {nn, 2, level}]];
             invMixingPointList = Map[# - mixingPoints[[1]] &, Drop[mixingPoints, 1]];
             list2 = eVecsNorm.Inverse[invMixingPointList];
             list2 = Chop[Join[Transpose[list2], {eVecs[[All, 1]] - Map[Apply[Plus, #] &, list2]}]];
             If[ GetInverse /. {optsPC} /. {GetInverse -> False},
                 iMatrix = Join[-{eVecs[[cornerList[[1]], Range[2, level]]]}, Inverse[Drop[eVecsNorm[[cornerList]], 1]]];
                 iMatrix = Transpose[Append[Transpose[iMatrix], Map[(1 - Total[#]) &, iMatrix]]];
                 result = {cornerList, ClusterMatrix[list2], list2.PseudoInverse[Transpose[eVecs]], Transpose[eVecs]};,
                 result = {cornerList, ClusterMatrix[list2]};
             ];
         ];
        result
    ];

End[] (* End Private Context *)

EndPackage[]