(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

Spectrum::usage = "";
SpectrumSensitivityMatrix::usage = "";
ImpliedTimeScales::usage = "";
EquilibriumSensitivityTensor::usage = "";
MeanFirstPassageTime::usage = "";
MeanFirstPassageTimeSensitivityTensor::usage = "";
Committor::usage = "";
CommittorSensitivityMatrix::usage = "";
PotentialOfMeanForce::usage = "";
PotentialOfMeanForceBins::usage = "";
PotentialOfMeanForceSensitivityMatrix::usage = "";
PotentialOfMeanForceError::usage = "";
CommittorReactionCoordinate = "";
VisitPredictionMatrix::usage = "";
SteadyStateEquilibrium::usage = "";
SteadyStateCommittor::usage = "";

Options[MeanFirstPassageTime]        = {Lagtime->1};


Begin["`Private`"] (* Begin Private Context *) 

TransitionMatrix/:Spectrum[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
        Eigenvalues[N[m]]

TransitionMatrix/:SpectrumSensitivityMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]]] := 
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

TransitionMatrix/:VisitPredictionMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],length_,optsMFPT:OptionsPattern[MeanFirstPassageTime]] := 
        Module[ {evals,evecs,nLagtime,nLength},
            nLagtime = Lagtime/.List[UnionRules[optsMFPT,opts,Lagtime->1]];
            nLength = length / nLagtime;
            {evals,evecs} = Eigensystem[m];
            Transpose[evecs].DiagonalMatrix[Map[If[ #<1,
                                                    (#^(nLength+1)-1)/(#-1),
                                                    (nLength+1)
                                                ]&,evals]].Transpose[Inverse[evecs]]
        ]        
        
TransitionMatrix/:SteadyStateEquilibrium[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],states_] := 
        Module[ {cM},
            cM = m;
            cM[[states[[2]]]] = Table[If[ ii==states[[1]],
                                          1,
                                          0
                                      ],{ii,1,Length[cM]}];
            (#/Total[#])&[Eigenvectors[Transpose[cM],Min[Length[states]+50,Length[m]]][[1]]]
        ]

TransitionMatrix/:SteadyStateCommittor[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],states_] := 
        Module[ {cM,eq},
            cM = m;
            eq = Equilibrium[TransitionMatrix[m]];
            cM[[All,states[[2]]]] = Table[If[ ii==states[[1]],
                                              eq[[states[[2]]]] / eq[[states[[1]]]],
                                              0
                                          ],{ii,1,Length[cM]}];
            {cM,(#/Total[#])&[Eigenvectors[cM,Min[Length[states]+50,Length[m]]][[1]]]}
        ]

TransitionMatrix/:PotentialOfMeanForce[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],states_] := 
        Module[ {committor, equilibrium, s2,s3,dd},
            committor = Committor[TransitionMatrix[m,opts],states][[2]];
            equilibrium = Equilibrium[TransitionMatrix[m,opts]];
            s2 = Sort[Transpose[{committor,equilibrium}],#1[[1]]<#2[[1]]&];
            dd = 0;
            s3 = Transpose[{s2[[All,1]],Table[dd = dd+s2[[ii,2]],{ii,1,Length[s2]}]}];
            Transpose[{s3[[All,1]],Identity[s3[[All,2]]]}]
        ]

TransitionMatrix/:PotentialOfMeanForceBins[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],states_,bins_] := 
        Module[ {committor, equilibrium, s2},
            committor = Committor[TransitionMatrix[m,opts],states][[2]];
            equilibrium = Equilibrium[TransitionMatrix[m,opts]];
            s2 = Sort[Transpose[{committor,equilibrium}],#1[[1]]<#2[[1]]&];
            Select[Table[{(bins[[nn]]+bins[[nn+1]])/2,-Log[Apply[Plus,Select[s2,(#[[1]]>=bins[[nn]]&&#[[1]]<=bins[[nn+1]])&][[All,2]]]/(bins[[nn+1]]-bins[[nn]])]},{nn,1,Length[bins]-1}],#[[2]]<\[Infinity]&]
        ]

TransitionMatrix/:CommittorReactionCoordinate[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],states_,expectations_] := 
        Module[ {committor,s2,s3,dd},
            committor = Committor[TransitionMatrix[m,opts],states][[2]];
            s2 = Sort[Transpose[{committor,expectations}],#1[[1]]<#2[[1]]&];
            dd = 0;
            s3 = Transpose[{s2[[All,1]],Table[dd = dd+s2[[ii,2]],{ii,1,Length[s2]}]}];
            Transpose[{s3[[All,1]],s3[[All,2]]/Apply[Plus,expectations]}]
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

TransitionMatrix/:PotentialOfMeanForceError[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],state1_,state2_] := 
        Module[ {cM,cP,size,cC,cL,cB},
            size = Length[m];
            cM = PseudoInverse[Join[Transpose[m]-IdentityMatrix[size],{Table[1,{size}]}]];
            cP = #/Apply[Plus,#]&[Eigenvectors[Transpose[m],1][[1]]];
            cC = Committor[TransitionMatrix[m,opts],{state1,state2}][[2]];
            cL = Sort[Transpose[{Range[Length[m]],cC}],#1[[2]]<#2[[2]]&];
            Table[
                cB = Map[Apply[Plus,#]&,cM[[All,cL[[Range[nn],1]]]]];
                {cL[[nn,2]],Apply[Plus,Abs[Flatten[Outer[Plus,cB,-cB]]]]}
            ,{nn,1,10}]
        ]

TransitionMatrix/:PotentialOfMeanForceSensitivityMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],state1_,state2_,c_] := 
        Module[ {cM,cP,size,cC,cL,cB},
            size = Length[m];
            cM = PseudoInverse[Join[Transpose[m]-IdentityMatrix[size],{Table[1,{size}]}]];
            cP = #/Apply[Plus,#]&[Eigenvectors[Transpose[m],1][[1]]];
            cC = Committor[TransitionMatrix[m,opts],{state1,state2}][[2]];
            cL = Select[Transpose[{Range[Length[m]],cC}],(#[[2]]<=c)&][[All,1]];
            cB = Map[Apply[Plus,#]&,cM[[All,cL]]];
            Outer[Plus,cB,-cB]
        (*Table[cP[[aa]]*(cB[[aa]]-cB[[bb]]),{aa,1,size},{bb,1,size}]*)
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

TransitionMatrix/:CommittorSensitivityMatrix[TransitionMatrix[m_?MatrixQ,opts:OptionsPattern[TransitionMatrix]],states_] := 
        Module[ {cM,cS,cA,size},
            size = Length[m];
            cS = IdentityMatrix[size][[Complement[Range[size],states]]];
            cM = Inverse[(m-IdentityMatrix[size])[[Complement[Range[size],states],Complement[Range[size],states]]]];
            cA = Transpose[cS].cM.cS;
            Transpose[cA]
        ]

End[] (* End Private Context *)

EndPackage[]