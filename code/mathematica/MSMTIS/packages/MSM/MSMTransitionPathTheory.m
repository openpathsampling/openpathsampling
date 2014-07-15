(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ReactiveFluxMatrix::usage = "";
NetFluxMatrix::usage = "";
ReactivityVector::usage = "";
StateExitFlowRate::usage = "";
StateExitFlowRateSensitivityMatrix::usage = "";

Begin["`Private`"] (* Begin Private Context *) 

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