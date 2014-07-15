(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

Spectrum::usage = "Spectrum[transitionmatrix] returns the eigenvalues of transitionmatrix.";
SpectrumSensitivityTensor::usage = "SpectrumSensitivityTensor[transitionmatrix] returns a tensor of sensitivities of the eigenvalues of transitionmatrix w.r.t. the entries in matrix.";
ImpliedTimeScales::usage = "ImpliedTimeScales[transitionmatrix] returns the implied timescales of transitionmatrix.";
EquilibriumSensitivityTensor::usage = "EquilibriumSensitivityTensor[transitionmatrix] returns a tensor of sensitivities of the equilibrium (first left eigenvector) of matrix w.r.t. the entries in matrix.";
MeanFirstPassageTime::usage = "MeanFirstPassageTime[transitionmatrix] returns the matri of mean first passage times of transitionmatrix.";
MeanFirstPassageTimeSensitivityTensor::usage = "MeanFirstPassageTimeSensitivityTensor[transitionmatrix] returns a tensor that contains the sensitivities of the meanfirstpassagetimes in terms of the entries of the transitionmatrix.";
Committor::usage = "Committor[transitionmatrix, states] returns the multi committor for dynamics in transitionmatrix assuming corestates given in states.";

VisitPredictionMatrix::usage = "VisitPredictionMatrix[transitionmatrix, steps] returns a matrix of estimated transitions under the dynamics in the given TransitionMatrix. If multiplied with an initial vector the returning vector will represent the average times a state is visited.";
SteadyStateEquilibrium::usage = "SteadyStateEquilibrium[transitionmatrix, {stA, stB}] returns the steady state equilibium assuming that stB is uniquely and only connected to stA";

ForwardCommittor::usage = "ForwardCommittor[transitionmatrix, states] returns the multi forward committor for dynamics in transitionmatrix assuming corestates given in states.";
BackwardCommittor::usage = "BackwardCommittor[transitionmatrix, states] returns the multi backward committor for dynamics in transitionmatrix assuming corestates given in states.";

ReactiveFluxMatrix::usage = "ReactiveFluxMatrix[transitionmatrix, stA, stB] returns the matrix of reactive fluxes connecting stA and stB.";
NetFluxMatrix::usage = "NetFluxMatrix[transitionamtrix, stA, stB] returns a matrix of net fluxes connecting stA and stB.";
ReactivityVector::usage = "ReactivityVector[transitionmatrix, stA, stB] returns a vector containing the reactivities, i.e. the probability to find a reactive trajectory at each state that connects stA and stB.";
StateExitFlowRate::usage = "StateExitFlowRate[transitionamtrix, stA, stB] returns the flow exit rate from stA into stB.";
StateExitFlowRateSensitivityMatrix::usage = "StateExitFlowRateSensitivityMatrix[transitionmatrix, stA, stB] returns a matrix of sensitivities of the exit rate of stA into stB depending on all entries of the TransitionMatrix.";

TransitionRate::usage = "TransitionRate[oom, {corestates,} stA, stB] returns the transition rate from stA to stB for the specified OOM. If coresets are given then these are to be excluded in the computation.";
ProbabilityCurrent::usage = "ProbabilityCurrent[oom, {corestates,} stA, stB] return the probability current between states stA and stB of the specified OOM. If coresets are given these are to be exluded in the computation.";

Reactivity::usage = "Reactivity[oom, {corestates,} stA, stB] returns the generalized reactivity, i.e. the probability that each state contains a reactive trajectory that connects stA and stB without hitting any of the corestates.";


Options[MeanFirstPassageTime]        = {Lagtime->1};

Begin["`Private`"] (* Begin Private Context *) 

OOM/:Committor[OOM[d_?ListQ,opts:OptionsPattern[OOM]],coreStates_] :=
        ForwardCommittor[OOM[d,opts],coreStates]
        
OOM/:ForwardCommittor[OOM[d_?ListQ,opts:OptionsPattern[OOM]],coreStates_] :=
        Module[ {size,stateSpace,voidStates},
            size = CountObservables[OOM[d,opts]];
            stateSpace = Range[size];
            voidStates = Complement[stateSpace,coreStates];
            Table[
                ComputeTraceProbability[OOM[d,opts],{stX,{voidStates,{0,Infinity}},stF}]/ComputeTraceProbability[OOM[d,opts],{stX}],
                {stF,coreStates},{stX,stateSpace}
            ]
        ]
        
OOM/:BackwardCommittor[OOM[d_?ListQ,opts:OptionsPattern[OOM]],coreStates_] :=
        Module[ {size,stateSpace,voidStates},
            size = CountObservables[OOM[d,opts]];
            stateSpace = Range[size];
            voidStates = Complement[stateSpace,coreStates];
            Table[
                ComputeTraceProbability[OOM[d,opts],{stF,{voidStates,{0,Infinity}},stX}]/ComputeTraceProbability[OOM[d,opts],{stX}],
                {stF,coreStates},{stX,stateSpace}
            ]
        ]



OOM/:MeanFirstPassageTime[OOM[d_?ListQ,opts:OptionsPattern[OOM]],st1_,st2_] := 
        Module[ {tauMatrices,init,sizeObservable,spaceComplement,matrixComplement},
            init = GetInitialDistribution[OOM[d,opts]];
            tauMatrices = GetListOfTauMatrices[OOM[d,opts]];
            sizeObservable = CountObservables[OOM[d,opts]];
            spaceComplement = Complement[Range[sizeObservable],{st2}];
            matrixComplement = Total[tauMatrices[[spaceComplement]]];
            Total[init.tauMatrices[[st1]].MatrixRangeMeanSum[matrixComplement,{0,Infinity}].tauMatrices[[st2]]]/ComputeTraceProbability[OOM[d,opts],{st1}]
        ]
        
OOM/:Reactivity[OOM[d_?ListQ,opts:OptionsPattern[OOM]],coreStates_?VectorQ,stA_?IntegerQ,stB_?IntegerQ] :=
        Module[ {size,stateSpace,voidStates},
            size = CountObservables[OOM[d,opts]];
            stateSpace = Range[size];
            voidStates = Complement[stateSpace,coreStates];
            Table[
                Which[
                    stX==stA,
                    ComputeTraceProbability[OOM[d,opts],{stX,{voidStates,{0,Infinity}},stB}],
                    stX==stB,
                    ComputeTraceProbability[OOM[d,opts],{stA,{voidStates,{0,Infinity}},stX}],
                    MemberQ[coreStates,stX],
                    0,
                    True,
                    ComputeTraceProbability[OOM[d,opts],{stA,{voidStates,{0,Infinity}},stX,{voidStates,{0,Infinity}},stB}]
                ],
                {stX,stateSpace}
            ]
        ]
        
OOM/:Reactivity[OOM[d_?ListQ,opts:OptionsPattern[OOM]],stA_?IntegerQ,stB_IntegerQ] :=
        Reactivity[OOM[d,opts],{stA,stB},stA,stB]
        
OOM/:ProbabilityCurrent[OOM[d_?ListQ,opts:OptionsPattern[OOM]],stA_?IntegerQ,stB_IntegerQ] :=
        ProbabilityCurrent[OOM[d,opts],{stA,stB},stA,stB]

OOM/:ProbabilityCurrent[OOM[d_?ListQ,opts:OptionsPattern[OOM]],coreStates_?VectorQ,stA_?IntegerQ,stB_IntegerQ] :=
        Module[ {size,stateSpace,voidStates},
            size = CountObservables[OOM[d,opts]];
            stateSpace = Range[size];
            voidStates = Complement[stateSpace,coreStates];
            Table[
                Which[
                    stX==stA&&stY==stB,
                    ComputeTraceProbability[OOM[d,opts],{stX,stY}],
                    stX==stA&&!MemberQ[coreStates,stY],
                    ComputeTraceProbability[OOM[d,opts],{stX,stY,{voidStates,{0,Infinity}},stB}],
                    !MemberQ[coreStates,stX]&&stY==stB,
                    ComputeTraceProbability[OOM[d,opts],{stA,{voidStates,{0,Infinity}},stX,stY}],
                    !MemberQ[coreStates,stX]&&!MemberQ[coreStates,stY],
                    ComputeTraceProbability[OOM[d,opts],{stA,{voidStates,{0,Infinity}},stX,stY,{voidStates,{0,Infinity}},stB}],
                    True,
                    0
                ],
                {stX,stateSpace},{stY,stateSpace}
            ]
        ]
        
OOM/:TransitionRate[OOM[d_?ListQ,opts:OptionsPattern[OOM]],stA_?IntegerQ,stB_IntegerQ] :=
        TransitionRate[OOM[d,opts],{stA,stB},stA,stB]

OOM/:TransitionRate[OOM[d_?ListQ,opts:OptionsPattern[OOM]],coreStates_?VectorQ,stA_?IntegerQ,stB_IntegerQ] :=
        Module[ {size,stateSpace,voidStates},
            size = CountObservables[OOM[d,opts]];
            stateSpace = Range[size];
            voidStates = Complement[stateSpace,coreStates];
            ComputeTraceProbability[OOM[d,opts],{stA,{voidStates,{0,Infinity}},stB}]
        ]

TransitionMatrix/:Spectrum[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Eigenvalues[N[m]]

TransitionMatrix/:SpectrumSensitivityTensor[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Module[ {evL,evR},
            evL = Eigenvectors[m];
            evR = Eigenvectors[Transpose[m]];
            Table[Outer[Times,evL[[nn]],evR[[nn]]]/(evL[[nn]].evR[[nn]]),{nn,1,Length[m]}]
        ]
        
TransitionMatrix/:ImpliedTimeScales[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],optsITS:OptionsPattern[{Lagtime->1}]] := 
        Module[ {nLagtime,nEV},
            nLagtime = Lagtime/.List[UnionRules[optsITS,opts,Lagtime->1]];
            nEV = Eigenvalues/.List[UnionRules[optsITS,Eigenvalues->Length[m]]];
            ToImpliedTimescales[Eigenvalues[N[m],nEV][[Range[1,nEV]]],nLagtime]
        ]

TransitionMatrix/:VisitPredictionMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],length_] := 
        Module[ {evals,evecs,nLagtime,nLength},
            nLagtime = Lagtime/.List[opts]/.{Lagtime->1};
            nLength = length / nLagtime;
            {evals,evecs} = Eigensystem[m];
            Transpose[evecs].DiagonalMatrix[
                Map[
                    If[ #<1,
                        (#^(nLength+1)-1)/(#-1),
                        (nLength+1)
                    ]&,evals
                ]
            ].Transpose[Inverse[evecs]]
        ]        
        
TransitionMatrix/:SteadyStateEquilibrium[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],states_] := 
        Module[ {cM},
            cM = m;
            cM[[states[[2]]]] = UnitVector[Length[m],states[[1]]];
            (#/Total[#])&[Eigenvectors[Transpose[cM],Min[Length[states]+50,Length[m]]][[1]]]
        ]
        
TransitionMatrix/:Committor[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],states_] := 
        Module[ {cM,cMEV},
            cM = N[Normal[m]];
            Map[(cM[[#]] = Table[If[ ii==#,
                                     1,
                                     0
                                 ],{ii,1,Length[cM]}])&,states];
            cMEV = Eigenvectors[cM][[Range[Length[states]]]];
            Inverse[cMEV[[All,states]]].cMEV
        ]

TransitionMatrix/:MeanFirstPassageTime[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],state_,optsMFPT:OptionsPattern[MeanFirstPassageTime]] := 
        Module[ {cM,nLagtime},
            nLagtime = Lagtime/.List[UnionRules[optsMFPT,opts,Lagtime->1]];
            cM = m-IdentityMatrix[Length[m]];
            cM[[state]] = Table[If[ ii==state,
                                    1,
                                    0
                                ],{ii,1,Length[cM]}];
            Inverse[cM].Table[If[ ii==state,
                                  0,
                                  -1
                              ],{ii,1,Length[cM]}]*nLagtime
        ]

TransitionMatrix/:MeanFirstPassageTimeSensitivityTensor[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],state_,optsMFPT:OptionsPattern[MeanFirstPassageTime]] := 
        Module[ {cM,nLagtime,mfpt,cMInv,size},
            size = Length[m];
            nLagtime = Lagtime/.List[UnionRules[optsMFPT,opts,Lagtime->1]];
            cM = m-IdentityMatrix[Length[m]];
            cM[[state]] = Table[If[ ii==state,
                                    1,
                                    0
                                ],{ii,1,Length[cM]}];
            cMInv = Inverse[cM];
            mfpt = cMInv.Table[If[ ii==state,
                                   0,
                                   -1
                               ],{ii,1,Length[cM]}]*nLagtime;
            Table[If[ ii==state,
                      0,
                      1
                  ]*cMInv[[kk,ii]]*mfpt[[jj]],{ii,1,size},{jj,1,size},{kk,1,size}]
        ]

TransitionMatrix/:EquilibriumSensitivityTensor[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Module[ {cM,cP,size},
            size = Length[m];
            cM = PseudoInverse[Join[Transpose[m]-IdentityMatrix[size],{Table[1,{size}]}]];
            cP = #/Apply[Plus,#]&[Eigenvectors[Transpose[m],1][[1]]];
            Table[cP[[aa]]*cM.Table[If[ ii==aa,
                                        1,
                                        0
                                    ]+If[ ii==bb,
                                          -1,
                                          0
                                      ],{ii,1,size+1}],{aa,1,size},{bb,1,size}]
        ]
        
(* TransitionPathTheory *)

TransitionMatrix/:ReactiveFluxMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],stateA_,stateB_] := 
        Module[ {eqDist,size,cM,tm},
            size = Length[m];
            tm = (m-IdentityMatrix[size]);
            cM = Committor[TransitionMatrix[m,opts],{stateA,stateB}][[2]];
            eqDist = Equilibrium[TransitionMatrix[m,opts]];
            Table[cM[[jj]]*(1-cM[[ii]])*eqDist[[ii]]*tm[[ii,jj]],{ii,1,size},{jj,1,size}]
        ];
TransitionMatrix/:NetFluxMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],stateA_,stateB_] := 
        Module[ {eqDist,size,cM,tm},
            size = Length[m];
            tm = (m-IdentityMatrix[size]);
            cM = Committor[TransitionMatrix[m,opts],{stateA,stateB}][[2]];
            eqDist = Equilibrium[TransitionMatrix[m,opts]];
            Map[Max[0,#]&,(#-Transpose[#])&[Table[cM[[jj]]*(1-cM[[ii]])*eqDist[[ii]]*tm[[ii,jj]],{ii,1,size},{jj,1,size}]],{2}]
        ];
TransitionMatrix/:ReactivityVector[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],stateA_,stateB_] := 
        Module[ {eqDist,size,cM},
            size = Length[m];
            cM = Committor[TransitionMatrix[m,opts],{stateA,stateB}][[2]];
            eqDist = Equilibrium[TransitionMatrix[m,opts]];
            Table[cM[[ii]]*(1-cM[[ii]])*eqDist[[ii]],{ii,1,size}]
        ];
TransitionMatrix/:StateExitFlowRate[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],stateA_,stateB_] := 
        Module[ {eqDist,size,cM},
            size = Length[m];
            cM = Committor[TransitionMatrix[m,opts],{stateA,stateB}][[2]];
            eqDist = Equilibrium[TransitionMatrix[m,opts]];
            Apply[Plus,Table[If[ jj==stateB,
                                 0,
                                 1
                             ]*cM[[jj]]*(1-cM[[stateA]])*eqDist[[stateA]]*m[[stateA,jj]],{jj,1,size}]]
        ];
TransitionMatrix/:StateExitFlowRateSensitivityMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],stateA_,stateB_] := 
        Module[ {committor, equilibrium, committorSensitivity, equilibriumSensitivity, eqInverse, size, dummySum1, dummySum2, matrix1, matrix2, matrix3},
            size = Length[m];
            committor = Committor[TransitionMatrix[m,opts],{stateA,stateB}][[2]];
            equilibrium = Equilibrium[TransitionMatrix[m,opts]];
            committorSensitivity = CommittorSensitivityMatrix[TransitionMatrix[m],{stateA,stateB}];
            eqInverse = PseudoInverse[Join[Transpose[m]-IdentityMatrix[size],{Table[1,{size}]}]];
            equilibriumSensitivity = Table[equilibrium[[kk]]*(eqInverse[[stateA,kk]]-eqInverse[[stateA,ll]]),{kk,1,size},{ll,1,size}];
            dummySum1 = Apply[Plus,committor*m[[stateA]]]-committor[[stateB]]*m[[stateA,stateB]];
            dummySum2 = Table[Apply[Plus,committorSensitivity[[nn]]*m[[stateA]]]-committorSensitivity[[nn,stateB]]*m[[stateA,stateB]],{nn,1,size}];
            matrix1 = SparseArray[{},{size,size}];
            matrix1[[stateA]] = committor-committor[[stateA]];
            matrix1[[stateA,stateB]] = -committor[[stateA]];
            matrix2 = Table[dummySum2[[kk]] * (committor[[kk]] - committor),{kk,1,size}];
            matrix3 = equilibrium[[stateA]]*(matrix1 + matrix2);
            equilibriumSensitivity * dummySum1 + Normal[matrix3]
        ];

End[] (* End Private Context *)

EndPackage[]