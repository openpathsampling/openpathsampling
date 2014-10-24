(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ToRateMatrix::MethodNotImplemented = "The Method has not yet been implemented. Please try Method -> MatrixLog.";
Likelihood::usage = "Likelihood[system, countmatrix] computed the likelihood of system to produce the observation in countmatrix.";

DirichletDistribution::usage = "DirichletDistribution is an option for the distribution used in CountMatrix.";
DiracDistribution::usage = "DiracDistribution is an option for the distribution used in CountMatrix.";
DirichletWithDBDistribution::usage = "DirichletWithDBDistribution is an option for the distribution used in CountMatrix.";

Options[DirichletWithDBDistribution] = {MaxIterations->1000,MaximalResidue->0.001};

Begin["`Private`"] (* Begin Private Context *) 

Unprotect[DirichletDistribution];

(* Likelihood Function *)

TransitionMatrix/:Likelihood[TransitionMatrix[t_?MatrixQ,optsT:OptionsPattern[CountMatrix]],CountMatrix[m_?MatrixQ,optsM:OptionsPattern[CountMatrix]]] :=
        Module[ {matrixAbs,t2},
            matrixAbs = Abs[Sign[m]];
            t2 = t+(1-matrixAbs);
            -Re[Apply[Plus,Flatten[Log[t2]*m]]]
        ];
RateMatrix/:Likelihood[RateMatrix[r_?MatrixQ,optsT:OptionsPattern[CountMatrix]],CountMatrix[m_?MatrixQ,optsM:OptionsPattern[CountMatrix]]] :=
        Module[ {matrixAbs,t,t2},
            t = MatrixExp[r];
            matrixAbs = Abs[Sign[m]];
            t2 = t+(1-matrixAbs);
            -Re[Apply[Plus,Flatten[Log[t2]*m]]]
        ];

(* Distribution Implementations*)

DirichletDistribution/:MeanMatrix[DirichletDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        TransitionMatrixFromMatrix[Map[#/Apply[Plus,#]&,(m+OptionValue[CountMatrix,Prior])]]
DirichletDistribution/:MaximumLikelihoodMatrix[DirichletDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        TransitionMatrixFromMatrix[Map[#/Apply[Plus,#]&,(m)]]
        
DirichletDistribution/:VarianceMatrix[DirichletDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        Map[#*(Apply[Plus,#]-#)/(Apply[Plus,#]^2*(Apply[Plus,#]+1))&,(m+OptionValue[CountMatrix,Prior])]

DirichletDistribution/:VarianceTensor[DirichletDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]],states_] := 
        Map[((Apply[Plus,#]*DiagonalMatrix[#]-Outer[Times,#,#])/(Apply[Plus,#]^2*(Apply[Plus,#]+1)))&,(m[[states]]+OptionValue[CountMatrix,Prior])]

DirichletDistribution/:DrawTransitionMatrix[DirichletDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        TransitionMatrix[Map[#/Apply[Plus,#]&,Map[RandomReal[GammaDistribution[#,1],1][[1]]&,m+1,{2}]]]
DiracDistribution/:MeanMatrix[DiracDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        TransitionMatrixFromMatrix[Map[#/Apply[Plus,#]&,(m)]]
DiracDistribution/:MaximumLikelihoodMatrix[DiracDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        TransitionMatrixFromMatrix[Map[#/Apply[Plus,#]&,(m)]]
DiracDistribution/:VarianceMatrix[DiracDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        0*m
DiracDistribution/:DrawTransitionMatrix[DiracDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        TransitionMatrix[Map[#/Apply[Plus,#]&,m]]
DirichletWithDBDistribution/:MeanMatrix[DirichletWithDBDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        Message[MeanMatrix::NYI]
DirichletWithDBDistribution/:VarianceMatrix[DirichletWithDBDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        Message[VarianceMatrix::NYI]
DirichletWithDBDistribution/:DrawTransitionMatrix[DirichletWithDBDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]]] := 
        Message[DrawTransitionMatrix::NYI]
DirichletWithDBDistribution/:MaximumLikelihoodMatrix[DirichletWithDBDistribution[CountMatrix[m_?MatrixQ,opts:OptionsPattern[CountMatrix]]],optsML:OptionsPattern[DirichletWithDBDistribution]] := 
        Module[ {size,sumVec,symmMat,result,nIterations,sumResVec,nResidue,oRes,ii},
            size = Length[m];
            nIterations = OptionValue[DirichletWithDBDistribution,List[optsML],MaxIterations];
            sumVec = Map[Apply[Plus,#]&,m];
            symmMat = m+Transpose[m];
            result = Table[1./size,{size},{size}];
            nResidue = OptionValue[DirichletWithDBDistribution,List[optsML],MaximalResidue];
            For[ii = 1,ii<=nIterations,ii++,
                sumResVec = Map[Apply[Plus,#]&,result];
                oRes = result;
                result = Table[symmMat[[ii,jj]]/(sumVec[[ii]]/sumResVec[[ii]]+sumVec[[jj]]/sumResVec[[jj]]),{ii,1,size},{jj,1,size}];
                If[ Max[Abs[result-oRes]]<nResidue,
                    ii = nIterations
                ];
            ];
            TransitionMatrix[Map[#/Apply[Plus,#]&,result]]
        ]

Protect[DirichletDistribution];

End[] (* End Private Context *)

EndPackage[]