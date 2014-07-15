(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ToRateMatrix::usage = "ToRateMatrix[system] converts system to a RateMatrix if possible.";
ToTransitionMatrix::usage = "ToTransitionMatrix[system] converts system to a TransitionMatrix if possible.";
ToCountMatrix::usage = "ToCountMatrix[system] converts system to a CountMatrix if possible.";
ToTimeSeries::usage = "ToTimeSeries[system] converts system to a TimeSeries if possible.";
ToCorrelationMatrix::usage = "ToCorrelationMatrix[system] converts system to a CorrelationMatrix if possible.";
ToOOM::usage = "ToOOM[system] converts system to a OOM if possible.";
ToCountTensor::usage = "ToCountTensor[timeseries] computes the count tensor for transitiion i->j->k in timeseries.";
ToPMM::usage = "ToPMM[system] converts system to a PMM if possible.";

GetCorrelationMatrix::usage = "GetCorrelationMatrix[pmm, lagtime] returns a correlation matrix at lagtime from the PMM.";
GetXM::usage = "GetXM[pmm, lagtime] returns the extendend transitionmatrix at lagtime from the PMM.";
GetCorrelation::usage = "GetCorrelation[pmm, observation, lagtime] returns the autocorrelation value for the PMM and an observation given as a vector at lagtime.";
GetCovariance::usage = "GetCorrelation[pmm, observation, lagtime] returns the autocorrelation value for the PMM and an observation given as a vector at lagtime.";
GetTransitionMatrix::usage = "GetTransitionMatrix[pmm, lagtime] returns the transitionmatrix at lagtime from the PMM.";

Discretize::usage = "Discretize[data, range] discretizes data into timeseries specified by range. Range can be a NearestFunction, {min, max, bins} or a number of bins";

Options[ToRateMatrix]                = {Method->MatrixLog,Lagtime->"Inherit"};
Options[ToTimeSeries]                = {AssumeSparse->False};
Options[KMeansClustering]            = {MaxIterations->20};

Begin["`Private`"] (* Begin Private Context *) 

TransitionMatrix/:ToRateMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],optsTRM:OptionsPattern[ToRateMatrix]] := 
        Module[ {rM,eVecs,eVals,nLagtime},
            If[ OptionValue[ToRateMatrix,optsTRM,Method]===MatrixLog,
                {eVals,eVecs} = Eigensystem[m];
                nLagtime = If[ OptionValue[ToRateMatrix,optsTRM,Lagtime]=="Inherit",
                               OptionValue[TransitionMatrix,Lagtime],
                               OptionValue[ToRateMatrix,optsTRM,Lagtime]
                           ];
                rM = RateMatrix[Transpose[eVecs].DiagonalMatrix[Log[eVals]/nLagtime].Inverse[Transpose[eVecs]],Lagtime->nLagtime];
            ];
            rM
        ]


TimeSeries/:ToCountMatrix[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]],optsCM:OptionsPattern[Join[Options[CountMatrix],{Lagtime->1, ClosedEnd->True, States->0, SeriesLength->All, MatrixType}]]] := 
        Module[ {cM,options,nStates,nLagtime,cV,matrixType, normLagtime},
            options = UnionRules[opts,States->Max[v]];
            nStates = Max[OptionValue[TimeSeries,States],Max[v],States/.{optsCM}/.States->0];
            nLagtime = Lagtime/.List[optsCM]/.{Lagtime->1};
            normLagtime = Max[1, nLagtime];
            matrixType = MatrixType/.{optsCM}/.{MatrixType->"DenseNumeric"};
            cM = SparseArray[Map[{Mod[#[[1]],(nStates+1)],Floor[#[[1]]/(nStates+1)]}->#[[2]]/normLagtime&,Tally[Drop[Drop[Join[(nStates+1)*v,Table[0,{nLagtime}]]+Join[Table[0,{nLagtime}],v],-nLagtime],nLagtime]]],{nStates,nStates}];
(*
            If[SeriesLength/.optsCM/.{SeriesLength->All}===All
                ,
                cvPart=Partition[v,SeriesLength/.optsCM/.{SeriesLength->Length[v]},SeriesLength/.optsCM/.{SeriesLength->Length[v]},1,{}];
                cM=SparseArray[Table[Map[{Mod[#[[1]],(nStates+1)],Floor[#[[1]]/(nStates+1)]}->#[[2]]&,Tally[Drop[Drop[Join[(nStates+1)*cvPart[[nn]],Table[0,{nLagtime}]]+Join[Table[0,{nLagtime}],cvPart[[nn]]],-nLagtime],nLagtime]]],{nn,1,Length[cvPart]}],{nStates,nStates}];

            ]
*)
            cV = v;
            If[ ClosedEnd/.List[optsCM]/.{ClosedEnd->False},
                While[(Sort[Union[Drop[cV,nLagtime]]]!=Sort[Union[Drop[cV,-nLagtime]]])&&(Length[cV]>nLagtime),
                    cV = Drop[cV,-1];
                    cM = SparseArray[Map[{Mod[#[[1]],(nStates+1)],Floor[#[[1]]/(nStates+1)]}->#[[2]]&,Tally[Drop[Drop[Join[(nStates+1)*cV,Table[0,{nLagtime}]]+Join[Table[0,{nLagtime}],cV],-nLagtime],nLagtime]]],{nStates,nStates}];
                ];
                If[ Length[cV]>=nLagtime,
                    cM = SparseArray[Map[{Mod[#[[1]],(nStates+1)],Floor[#[[1]]/(nStates+1)]}->#[[2]]/normLagtime&,Tally[Drop[Drop[Join[(nStates+1)*cV,Table[0,{nLagtime}]]+Join[Table[0,{nLagtime}],cV],-nLagtime],nLagtime]]],{nStates,nStates}],
                    cM = SparseArray[{},{nStates,nStates}]
                ];
            ];
            If[ Reversible/.{optsCM}/.{Reversible->True},
                cM = 1 / 2 * (cM+Transpose[cM])
            ];
            Which[matrixType=="DenseNumeric",
                cM = N[Normal[cM]],
                matrixType=="Dense",
                cM = Normal[cM],
                matrixType=="SparseNumeric",
                cM = N[cM]
            ];
            CountMatrix[cM,UnionRules[RemoveRules[{AssumeSparse,Reversible,ClosedEnd,MatrixType},optsCM],Lagtime->nLagtime]]
        ]

CountMatrix/:ToCorrelationMatrix[CountMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]], optsCM___] := 
        Module[ {lagtime},
            symmetrize = Symmetrize/.{optsCM}/.{Symmetrize->True};
            lagtime = Lagtime/.{opts}/.{Lagtime->1};
            cM = CorrelationMatrix[#/Total[Flatten[#]]&[m],   Lagtime -> lagtime];
            If[ symmetrize,
                cM = Symmetrize[cM];
            ];
            cM
        ]

TransitionMatrix/:ToCorrelationMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Module[ {eqDist,size,lagtime},
            lagtime = Lagtime/.{opts}/.{Lagtime->1};
            size = Length[m];
            eqDist = Equilibrium[TransitionMatrix[m,opts]];
            CorrelationMatrix[DiagonalMatrix[eqDist].m,Lagtime->lagtime]
        ]

CorrelationMatrix/:ToTransitionMatrix[CorrelationMatrix[m_?MatrixQ,opts:OptionsPattern[CorrelationMatrix]]] := 
        Module[ {},
            lagtime = Lagtime/.{opts}/.{Lagtime->1};
            Normalize[TransitionMatrix[m,Lagtime->lagtime]]
        ]
        
OOM/:ToTransitionMatrix[OOM[d_?ListQ,opts:OptionsPattern[OOM]], lagtime_] :=
        Module[ {size},
            size = CountObservables[OOM[d,opts]];
            tm = Table[
                ComputeTraceProbability[OOM[d,opts],{stX,{All,lagtime - 1},stY}]/ComputeTraceProbability[OOM[d,opts],{stX}],
                {stX,size},{stY,size}
            ];
            TransitionMatrix[Chop[tm], Lagtime->lagtime]
        ]
        
TransitionMatrix/:ToOOM[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] :=
    Module[ {},
        OOM[{"InitialDistribution"->Equilibrium[TransitionMatrix[m,opts]],"TauMatrices"->Map[m.DiagonalMatrix[#]&,IdentityMatrix[Length[m]]]}]
    ]
    
TransitionMatrix/:ToOOM[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],processes_] :=
    Module[ {evals,evecL,evecR,equilibrium,lambdaHalf,oneVecM,oneVecK,mRotation,mRotationT,sMatrix,eqMHalfI},
        {evals,evecL,evecR} = FullEigensystem[TransitionMatrix[m,opts]][[All,processes]];
        equilibrium = evecL[[1]];
        eqMHalfI = DiagonalMatrix[1/Sqrt[equilibrium]];
        evecL = EquilizeLeftEVSet[evecL];
        evecR = evecL.Inverse[DiagonalMatrix[equilibrium]];
        lambdaHalf = DiagonalMatrix[Sqrt[evals]];
        oneVecM = Table[1,{Length[m]}];
        oneVecK = Table[1,{Length[processes]}];
        If[ Length[processes]>1,
            mRotation = RotationMatrix[{lambdaHalf.evecL.oneVecM,oneVecK}];,
            mRotation = {{1}};
        ];
        mRotationT = Transpose[mRotation];
        sMatrix = 1/Sqrt[Length[processes]]*Transpose[evecR].lambdaHalf.mRotationT;
        OOM[{"InitialDistribution"->equilibrium.sMatrix,"TauMatrices"->Map[mRotation.lambdaHalf.evecL.eqMHalfI.DiagonalMatrix[#].eqMHalfI.Transpose[evecL].lambdaHalf.mRotationT&,IdentityMatrix[Length[m]]]}]
    ]
    
    TimeSeries /: 
 ToCountTensor[
  TimeSeries[v_?VectorQ, opts : OptionsPattern[TimeSeries]], 
  optsCM : OptionsPattern[
    Join[Options[CountMatrix], {Lagtime -> 1, ClosedEnd -> True, 
      States -> 0, SeriesLength -> All, MatrixType, Clustering}]]] := 
 Module[ {cM, options, nStates, nLagtime, cV, matrixType, state3v},
     options = UnionRules[opts, States -> Max[v]];
     nStates = 
      Max[OptionValue[TimeSeries, States], Max[v], 
       States /. {optsCM} /. States -> 0];
     clustering = (Clustering /. {optsCM}) /. {Clustering -> 
         IdentityClusterMatrix[nStates]};
     nClStates = Length[clustering][[2]];
     nLagtime = Lagtime /. List[optsCM] /. {Lagtime -> 1};
     matrixType = 
      MatrixType /. {optsCM} /. {MatrixType -> "SparseNumeric"};
     clV = Cluster[clustering, TimeSeries[v, opts]][[1]];
     state3v = 
      Drop[Drop[
        Join[(nStates)^2*(clV - 1), Table[0, {nLagtime*2}]] + 
         Join[Table[0, {nLagtime}], (nStates)*(v - 1), 
          Table[0, {nLagtime}]] + 
         Join[Table[0, {nLagtime*2}], (clV - 1)], -2*nLagtime], 2*nLagtime];
     cM = SparseArray[
       Map[(1 + {Mod[#[[1]], nStates], 
             Mod[Floor[#[[1]]/(nStates)], nStates], 
             Mod[Floor[#[[1]]/(nStates^2)], nStates]}) -> #[[2]]/
           nLagtime &, Tally[state3v]], {nClStates, nStates, 
        nClStates}];
     cV = v;
     If[ Reversible /. {optsCM} /. {Reversible -> False},
         cM = 1/2*(cM + Transpose[cM])
     ];
     Which[matrixType == "DenseNumeric", cM = N[Normal[cM]], 
      matrixType == "Dense", cM = Normal[cM], 
      matrixType == "SparseNumeric", cM = N[cM]];
     cM
 ]

TimeSeries /: 
 ToOOM[TimeSeries[v_?VectorQ, opts : OptionsPattern[TimeSeries]], 
  ClusterMatrix[c_?MatrixQ, optsCM : OptionsPattern[ClusterMatrix]], 
  optsOOM___] := 
 Module[ {lagtime},
     lagtime = Lagtime /. {optsOOM} /. {Lagtime -> 1};
     cM = #/Total[Flatten[#]] &[
       ToCountMatrix[TimeSeries[Drop[v, lagtime], Lagtime -> lagtime], 
         Lagtime -> lagtime] // Normal];
     cM = c.cM.Transpose[c];
     cT = ToCountTensor[TimeSeries[v, opts], Lagtime -> lagtime, 
        Clustering -> ClusterMatrix[c, optsCM]] // Normal;
     cT = Table[
         1/2*(cT[[All, ii, All]] + Transpose[cT[[All, ii, All]]]), {ii, 
          1, Dimensions[cT][[2]]}]/(Length[v] - 2*lagtime)*lagtime;
     OOM[{"TauMatrices" -> Map[Inverse[cM].# &, cT], 
       "InitialDistribution" -> Map[Total, cM]}]
 ]
    
TimeSeries /: ToOOM[TimeSeries[v_?VectorQ, opts : OptionsPattern[TimeSeries]], optsOOM___] := 
 Module[ {lagtime},
     lagtime = Lagtime /. {optsOOM} /. {Lagtime -> 1};
     cM = #/Total[Flatten[#]] &[
       ToCountMatrix[TimeSeries[Drop[v, lagtime], Lagtime -> lagtime], 
         Lagtime -> lagtime] // Normal];
     cT = ToCountTensor[TimeSeries[v, opts], Lagtime -> lagtime] // 
    Normal;
     cT = Table[
         1/2*(cT[[All, ii, All]] + Transpose[cT[[All, ii, All]]]), {ii, 
          1, Dimensions[cT][[2]]}]/(Length[v] - 2*lagtime)*lagtime;
     OOM[{"TauMatrices" -> Map[Inverse[cM].# &, cT], 
       "InitialDistribution" -> Map[Total, cM]}]
 ]
        
TransitionMatrix/:ToTimeSeries[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],initialState_,num___:1,optsTTS:OptionsPattern[ToTimeSeries]] := 
        Module[ {endState,size,list,startState,listOfRules,state,eq},
            size = Length[m];
            state = If[ initialState===Equilibrium,
                        eq = Equilibrium[TransitionMatrix[m,opts]];
                        RandomChoice[eq -> Range[Length[eq]]],
                        initialState
                    ];
            If[ OptionValue[ToTimeSeries,optsTTS,AssumeSparse]==True,
                listOfRules = Table[
                    (#[[All,1]]->#[[All,2]])&[Map[#[[2]]->#[[1,1]]&,Drop[ArrayRules[m[[nn]]],-1]]]
                ,{nn,1,Length[m]}];
                list = NestList[RandomChoice[listOfRules[[#]]]&,state,num];,
                list = Table[0,{num+1}];
                endState = state;
                list[[1]] = state;
                Table[
                    startState = endState;
                    endState = RandomChoice[m[[startState]]->Range[size]];
                    list[[nn]] = endState;
                ,{nn,2,num+1}];
            ];
            TimeSeries[list,States->size]
        ]

KMeansClustering[m_?MatrixQ,nCenters_,optsKM:OptionsPattern[KMeansClustering]] :=
    Module[ {centerList,nearest,clusters,nClusters,nIterations,oCL,residuum,ii},
        centerList = If[ NumberQ[nCenters],
                         m[[Round[Length[m]/nCenters*Range[0.5,nCenters-0.5,1.0]]]],
                         nCenters
                     ];
        nIterations = OptionValue[KMeansClustering,optsKM,MaxIterations];
        For[ii = 1,ii<=nIterations,ii++,
            nClusters = Length[centerList];
            nearest = Nearest[centerList->Automatic];
            clusters = Map[nearest[#,1][[1]]&,m];
            oCL = centerList;
            centerList = Table[
                Mean[m[[Select[Transpose[{Range[Length[clusters]],clusters}],#[[2]]==nn&][[All,1]]]]]
            ,{nn,1,nClusters}];
            residuum = Apply[Plus,Map[Norm,centerList-oCL]];
            If[ residuum==0,
                ii = nIterations
            ];
        ];
        centerList
    ]

ToTimeSeries[m_?MatrixQ,centers_] :=
    Module[ {nearest},
        nearest = Nearest[centers->Automatic];
        TimeSeries[Map[nearest[#,1][[1]]&,m],States->Length[centers]]
    ]

TransitionMatrix/:ToPMM[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
Module[ {es},
    es = FullEigensystem[TransitionMatrix[m]];
    PMM[{es[[1]], EquilizeLeftEVSet[es[[2]]]}]
]
  
PMM/:GetCorrelationMatrix[PMM[{v_?VectorQ,m_?MatrixQ}, opts:OptionsPattern[PMM]], lag_] := 
 CorrelationMatrix[Transpose[m].DiagonalMatrix[Chop[N[v]^lag]].m, Lagtime->lag]
 
PMM/:GetXM[PMM[{v_?VectorQ,m_?MatrixQ}, opts:OptionsPattern[PMM]], lag_] := 
 CorrelationMatrix[Inverse[Transpose[m].m].Transpose[m].DiagonalMatrix[Chop[N[v]^lag]].m, Lagtime->lag]

PMM/:GetCorrelation[PMM[{v_?VectorQ,m_?MatrixQ}, opts:OptionsPattern[PMM]], obs_?VectorQ, lag_] := 
 Chop[obs.GetCorrelationMatrix[PMM[{v,m}],lag].Transpose[obs]]

PMM/:GetCovariance[PMM[{v_?VectorQ,m_?MatrixQ}, opts:OptionsPattern[PMM]], obs1_?VectorQ, obs2_?VectorQ, lag_] := 
 	Chop[obs1.GetCorrelationMatrix[PMM[{v,m}],lag].Transpose[obs2]]


PMM/:GetTransitionMatrix[PMM[{v_?VectorQ,m_?MatrixQ}, opts:OptionsPattern[PMM]], lag_] := 
 	ToTransitionMatrix[GetCorrelationMatrix[PMM[{v,m}], lag]]
 
Discretize[m_?ListQ, range_] :=
    Module[ {mi, ma, st},
        Which[
            Head[range] === NearestFunction,
               TimeSeries[
                   Map[First[range[#]] &, m],
                  States -> range[[2, 1]],
                "Centers" -> range[[5]],
                "Length" -> Length[m]
            ],
               VectorQ[range] && Length[range] == 3,
            {mi, ma, st} = range[[{1, 2, 3}]];
            TimeSeries[
                ((Floor[(m - mi)/(ma - mi)*st]) /. {st -> st - 1}) + 1,
                  States -> st,
                  "Length" -> Length[m],
                  "Centers" -> (Range[st] - 0.5)/st*(ma - mi) + mi,
                  "Maximum" -> ma,
                  "Minimum" -> mi
              ],
            NumberQ[range] && IntegerQ[range],
             {mi, ma} = {Min[#], Max[#]} &[m];
             st = range;
             TimeSeries[
                  ((Floor[(m - mi)/(ma - mi)*st]) /. {st -> st - 1}) + 1,
                  States -> st,
                   "Length" -> Length[m],
                   "Centers" -> (Range[st] - 0.5)/st*(ma - mi) + mi,
                   "Maximum" -> ma,
                "Minimum" -> mi
              ]
        ]
    ]
    

End[] (* End Private Context *)

EndPackage[]