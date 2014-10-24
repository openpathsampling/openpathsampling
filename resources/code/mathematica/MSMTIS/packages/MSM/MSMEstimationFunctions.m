(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

EstimateTimeScales::usage = "EstimateTimeScales[correlationmatrix1, correlationmatrix2] returns the simple unique spectral estimation solution.";
EstimatePMM::usage = "EstimatePMM[cormatrix1, cormatrix2] estimates a PMM by a generalized eigendecompisition approach from the two specified CorrelationMatrix objects.";

FitTimeScales::usage = "FitTimeScales[correlationmatrixlist, size] estimates a correlation model by optimization for a specified list of CorrelationMatrix objects.";

EstimateTSFromTimeSeries::usage = "SpectralEstimation[timeseries] runs the spectral estimation algorithm on the timeseries to estimate the dominant timescales.";

EstimateHMM::usage = "EstimateHMM[timeseries] runs a EM algorithm to estimate the hidden markov model in the timeseries.";
Improvements::usage = "Improvements is an option for EstimateHMM that specifies the attempted number of improvement steps.";
Guesses::usage = "Guesses is an option for EstimateHMM that specifies the number of randomized initial guesses.";
HiddenStates::usage = "HiddenStates is an optin for EstimateHMM that specifies the number of assumed hidden states.";

Replacements::usage = "Replacements is an option for FitTimeScales that contains prior replacements that you want to set. Like fixed eigenvalues."
(*
FitPMM::usage = "FitPMM[...] estimates a PMM from the data specified.";
EstimateExpDecay::usage = "EstimateExpDecay[data, range] estimates the parameters of a single exponential decay by a linear to the logarithm of the data where";
*)

QLQDecomposition::usage = "QLQDecomposition[m1,m2] returns the QLQ Decomposition of two correlation matrices given two lagtimes t1 and t2. The times can be specified using the Lagtime option.";
SpectralEstimation::usage = "SpectralEstimation[n, c1, c2] computes the low rank QLQ approximation of two correlations matrices c1 and c2 of size n."

Begin["`Private`"] (* Begin Private Context *) 

QLQDecomposition[m1_, m2_, opts___] :=
    Module[ {eval2, evec2, t1, t2},
        {t1, t2} = (Lagtime /. {opts}) /. {Lagtime -> {1, 2}};        
        {eval2, evec2} = Eigensystem[{m1, m2}];
        {eval2, evec2} = 
         Chop@Transpose@
           Map[{If[ #[[1]] != 0 && NumericQ[#[[1]]],
                    #[[1]]^(1/(t1 - t2)),
                    0
                ], #[[2]]} &, Transpose@{eval2, evec2}];
        {eval2, evec2} = 
         Transpose@
          Map[{#[[1]], 
             Sign[#[[2, 1]]]*#[[2]]*Sqrt[(#[[1]]^t1)/(#[[2]].m1.#[[2]])]} &, 
           Transpose[{eval2, evec2}]];
        {eval2, evec2} = 
         Transpose@
          Sort[Transpose[{eval2, evec2}], Abs[#1[[1]]] > Abs[#2[[1]]] &];
        {eval2, Transpose@Inverse[evec2]}
    ]

SpectralEstimation[n_, CorrelationMatrix[m1_?MatrixQ,opts1:OptionsPattern[CorrelationMatrix]], CorrelationMatrix[m2_?MatrixQ,opts2:OptionsPattern[CorrelationMatrix]]] :=
    Module[ {
    		lEV, lED, svd1, svd2, svP11, svP21, svP12, svP22, err1, err2, svP13, svP23, lowRankC1,
	    	lowRankC2, erF1, erF2, mixM, mixM1, mixM2, pc21, pc22, eVecT, eValT, realT, wMD, wMV,
	    	appc21, appc22, lt1, lt2
    	},
    	lt1 = Lagtime /. {opts1} /. {Lagtime -> 1};
        lt2 = Lagtime /. {opts2} /. {Lagtime -> 2};
        lEV = Map[Normalize, LeftEigenvectors[ToTransitionMatrix@CorrelationMatrix@m2]][[Range[n]]];
        lED = Eigenvalues@ToTransitionMatrix@CorrelationMatrix@m2;
        svd1 = SingularValueDecomposition[m1];
        svd2 = SingularValueDecomposition[m2];
        svP11 = Transpose@svd1[[1, All, Range[n]]];
        svP21 = Transpose@svd2[[1, All, Range[n]]];
        svP12 = svd1[[2, Range[n], Range[n]]];
        svP22 = svd2[[2, Range[n], Range[n]]];
        If[n < Length[m1],
	        err1 = svd1[[2, n + 1, n + 1]]/svd1[[2, n, n]];
	        err2 = svd2[[2, n + 1, n + 1]]/svd2[[2, n, n]];
	        ,
	        err1 = 0;
	        err2 = 0;
        ];
        svP13 = Transpose@svd1[[3, All, Range[n]]];
        svP23 = Transpose@svd2[[3, All, Range[n]]];
        lowRankC1 = Transpose[svP11].svP12.svP13;
        lowRankC2 = Transpose[svP21].svP22.svP23;
        erF1 = FrobeniusNorm[lowRankC1 - m1]/FrobeniusNorm[m1];
        erF2 = FrobeniusNorm[lowRankC2 - m2]/FrobeniusNorm[m2];
        mixM = svP11.Transpose@svP21;
        mixM1 = svP11.Transpose@svP13;
        mixM2 = svP21.Transpose@svP23;
        pc21 = svP21.m1.Transpose[svP21];
        pc22 = svP21.m2.Transpose[svP21];
        eVecT = Eigenvectors[{svP21.m1.Transpose[svP21], svP21.m2.Transpose[svP21]}];
        eValT = Eigenvalues[{svP21.m1.Transpose[svP21], svP21.m2.Transpose[svP21]}];
        realT = Transpose[((svP21).Transpose@lEV)];
        {wMD, wMV} = QLQDecomposition[pc21, pc22, Lagtime->{lt1,lt2}];
        appc21 = Transpose[wMV].DiagonalMatrix[wMD^lt1].wMV;
        appc22 = Transpose[wMV].DiagonalMatrix[wMD^lt2].wMV;
        {
        	wMD, 
        	wMV.svP21
        }
    ]

EstimateTimeScales[CorrelationMatrix[m1_?MatrixQ,opts1:OptionsPattern[CorrelationMatrix]], CorrelationMatrix[m2_?MatrixQ,opts2:OptionsPattern[CorrelationMatrix]]] :=
    Module[ {lt1, lt2, tScales, sol, mu, muList, muListReal, muListComplex},
        lt1 = Lagtime /. {opts1} /. {Lagtime -> 1};
        lt2 = Lagtime /. {opts2} /. {Lagtime -> 2};
        sol = Solve[Det[m1 - mu*m2] == 0, mu];
        muList = mu /. sol;
        muListReal = Select[muList, Im[Chop[#]] == 0 &];
        muListComplex = Select[muList, Im[Chop[#]] != 0 &];
        tScales = muListReal^(1/(lt1 - lt2));
        tScales
    ]
  
EstimatePMM[CorrelationMatrix[m1_?MatrixQ,opts1:OptionsPattern[CorrelationMatrix]], CorrelationMatrix[m2_?MatrixQ,opts2:OptionsPattern[CorrelationMatrix]]] :=
    Module[ {lt1, lt2, tScales, orthoSpace, sol, mu, muList, oneVecSpace, 
    resultVectorListInv, muListReal, muListComplex, resultVectorListUN,
    resultVectorList, nullSpaceSize, pVec, imageSpaceSize},
        lt1 = Lagtime /. {opts1} /. {Lagtime -> 1};
        lt2 = Lagtime /. {opts2} /. {Lagtime -> 1};
        orthoSpace = NullSpace[m1];
        nullSpaceSize = Dimensions[orthoSpace][[1]];
        imageSpaceSize = Length[m1] - nullSpaceSize;
        pVec = Eigenvectors[m1][[Range[imageSpaceSize]]];
        sol = Solve[
          Det[pVec.m1.Transpose[pVec] - mu*pVec.m2.Transpose[pVec]] == 0, 
          mu
        ];
        muList = mu /. sol;
        muListReal = Select[muList, Head[#] === Real &];
        muListComplex = Select[muList, Head[#] =!= Real &];
        tScales = muListReal^(1/(lt1 - lt2));
        resultVectorListInv = Map[(
            oneVecSpace = NullSpace[m1 - #*m2];
            Total[
             oneVecSpace - oneVecSpace.(Transpose[orthoSpace].orthoSpace)]) &
          , muListReal
          ];
        resultVectorListUN = 
         MapIndexed[#/Sqrt[(#.m1.#/(tScales[[#2[[1]]]]^lt1))] &, 
          resultVectorListInv];
        resultVectorList = Transpose[PseudoInverse[resultVectorListUN]];
        PMM[{tScales, resultVectorList}, "ComplexEstimates" -> muListComplex]
    ]
    
EstimatePMM[m1_?MatrixQ, lt1_, m2_?MatrixQ, lt2_] :=
    EstimatePMM[CorrelationMatrix[m1, Lagtime -> lt1], CorrelationMatrix[m2, Lagtime -> lt2]]
(*
EstimateExpDecay[data_, range_] :=
    Module[ {},
        {1/Exp[c2], c1} /. 
         FindFit[Log[data[[range]]], c1 - x*c2, {c1, c2}, x]
    ]
    *)
  
(* Optimization Methods *)

qMatrix[size_, size2_] :=
    Table[V[q, ii, jj], {ii, size2}, {jj, size}]
qMatrixI[size_, size2_, stat_] := 
    Transpose@Join[{eq}, 
    	Transpose[Table[V[q, ii, jj], {ii, size2}, {jj, 2, size}]]]
qMatrix[size_] :=
    qMatrix[size, size]

lVectorF[size_] :=
    Table[V[l, ii], {ii, 0, size - 1}]
lVectorFInit[size_, values_] :=
    Table[{V[l, ii], values[[ii + 1]]}, {ii, 0, size - 1}]
lVector0[size_] :=
    Join[{1}, Table[V[l, ii], {ii, size - 1}]]

vListF[size_] :=
    vListF[size, size]
vListFInit[size_, values_] :=
    vListFInit[size, size, values]
vList0[size_] :=
    vList0[size, size]

vListF[size_, size2_] :=
    Flatten[{qMatrix[size, size2], lVectorF[size]}]
vListFInit[size_, size2_, values_] :=
    Flatten[{Flatten[qMatrix[size, size2]], lVectorFInit[size, values]}, 
     1]
vList0[size_, size2_] :=
    Flatten[{qMatrix[size, size2], Drop[lVectorF[size], 1]}]
vListI[size_, size2_] :=
    Flatten[{Transpose@Drop[Transpose@qMatrix[size, size2],1], Drop[lVectorF[size], 1]}]

mmF[size_] :=
    mmF[size, size]
mm0[size_] :=
    mm0[size, size]

mmF[size_, size2_] :=
    Function[k, 
     qMatrix[size, size2].DiagonalMatrix[lVectorF[size]^k].Transpose[
       qMatrix[size, size2]]]
mm0[size_, size2_] :=
    Function[k, 
     qMatrix[size, size2].DiagonalMatrix[lVector0[size]^k].Transpose[
       qMatrix[size, size2]]]

mmI[size_, size2_,stat_] :=
    Function[k, 
     qMatrixI[size, size2].DiagonalMatrix[lVector0[size]^k].Transpose[
       qMatrixI[size, size2]]]
       

FitTimeScales[observationList_, size_, opts___] :=
    Module[ {mode, fitOpts, result},
        mode = (Mode /. List[opts]) /. {Mode -> "AllTimeScales"};
        fitOpts = RemoveRules[Mode, opts];
        size2 = Length[observationList[[1, 1]]];
        data = Map[((Lagtime /. (Rest@(List @@ #))) /. 
              Lagtime -> 1) -> #[[1]] &, observationList];
        minLt = Min[data[[All, 1]]];
        replacements = (Replacements/.{opts})/.{Replacements->{}};
        equil = Map[Mean, Transpose@Map[Total, observationList[[All, 1]]]];
        result = Which[
          mode == "FixedStationary",
          eqnList = Evaluate[
            Apply[Plus,
             Map[
              Total[
                Flatten[(#[[2]] - 
                    mmI[size, size2, equil][#[[1]] - minLt + 1])^2]] &,
              data
              ]
             ]
            ];
          varList = Evaluate[vListI[size, size2]];
          lVec = lVector0[size];
          qMat = qMatrixI[size, size2, equil];
          consList = 
           Evaluate[
            Flatten[Join[Thread[Drop[lVec, 1] >= 0], 
              Thread[Drop[lVec, 1] <= 1]]]];
          ,
          mode == "FixedOne",
          eqnList = Evaluate[
            Apply[Plus,
             Map[
              Total[
                Flatten[(#[[2]] - 
                    mm0[size, size2][#[[1]] - minLt + 1])^2]] &,
              data
              ]
             ]
            ];
          varList = Evaluate[vList0[size, size2]];
          lVec = lVector0[size];
          qMat = qMatrix[size, size2];
          consList = 
           Evaluate[
            Flatten[Join[Thread[Drop[lVec, 1] >= 0], 
              Thread[Drop[lVec, 1] <= 1]]]];
          ,
          mode == "AllTimeScales",
          eqnList = Evaluate[
            Apply[Plus,
             Map[
              Total[
                Flatten[(#[[2]] - 
                    mmF[size, size2][#[[1]] - minLt + 1])^2]] &,
              data
              ]
             ]
            ];
          varList = Evaluate[vListF[size, size2]];
          lVec = lVectorF[size];
          qMat = qMatrix[size, size2];
          consList = 
           Evaluate[Flatten[Join[Thread[lVec >= 0], Thread[lVec <= 1]]]];
          ];
        result = {};
        If[ Length[varList] > 0,
        	replacements = Map[varList[[#[[1]]]] -> #[[2]]&, replacements];
        	eqnList = eqnList/.replacements;
            erg = NMinimize[
              {eqnList, consList},
              varList,
              fitOpts
              ];
            estimatedTimescales = (lVec /. erg[[2]]);
            pMatrix = 
             IdentityMatrix[
               size][[Sort[
                 Transpose[{Range[size], 
                   estimatedTimescales}], #1[[2]] > #2[[2]] &][[All, 1]]]];
            estimatedTimescales = pMatrix.estimatedTimescales;
            innerQMatrix = (qMat /. erg[[2]]).Transpose[pMatrix];
            result = {
              "Eigenvalues" -> estimatedTimescales,
              "QMatrix" -> innerQMatrix,
              "EstimationError" -> erg[[1]],
              "Lagtimes" -> data[[All, 1]],
              "Equilibrium" -> #/Total[#] &[innerQMatrix[[All, 1]]]
              };
        ];
        result
    ]

(*
FitPMM[size_, data_,projection_, {t1_, t2_, dt_}, {lt1_, lt2_, dlt_}, opts___] :=
    Module[ {
      estimate,
      innerqmatrix,
      projectedData,
      QMatrix = {{}},
      estimatedtimescales,
      eigenvalues = {},
      usedtimescalesList,
      realtimescalesList,
      correctedtimescales,
      realdt,
      realt1,
      correctionMatrix,
      permutationmatrix,
      qmatrix,
      maxiterations = 2000,
      realtimedist,
      simpleEstimation,
      equilibrium
      },
        QMatrix = qMatrix[size];
        eigenvalues = lVectorF[size];
        maxiterations =   MaxIterations /. {opts} /. {MaxIterations -> maxiterations};
        projectedData = Map[projection.#[[1]].Transpose[projection] &, data];
        usedtimescalesList = Range[t1, t2, dt];
        realtimescalesList = Range[lt1, lt2, dlt][[Range[t1, t2, dt]]];
        realdt = realtimescalesList[[2]] - realtimescalesList[[1]];
        realt1 = realtimescalesList[[1]];
        realtimedist = realtimescalesList[[-1]] - realtimescalesList[[1]];
        simpleEstimation = 
         Sort[Abs[mu /. 
             Solve[Det[projectedData[[t2]] - mu*projectedData[[t1]]] == 0, mu]]^(1/
             realtimedist)];
        Check[
         estimate = 
          FindMinimum[{Apply[Plus, 
             Table[Total[
               Flatten[(projectedData[[nn]] - mmF[size][(nn - t1)/dt + 1])^2]], {nn,
                t1, t2, dt}]]}, 
           vListFInit[size, Map[Min[#, 1] &, #] &[simpleEstimation]], 
           MaxIterations -> maxiterations];
         estimatedtimescales = eigenvalues /. estimate[[2]];
         correctedtimescales = Exp[Log[estimatedtimescales]/realdt];
         correctionMatrix = 
          DiagonalMatrix[1/Exp[realt1/(2*realdt) Log[estimatedtimescales]]];
         qmatrix = (QMatrix.correctionMatrix /. estimate[[2]]);
         innerqmatrix = qmatrix.Transpose[permutationmatrix];
         permutationmatrix = 
          IdentityMatrix[size][[
           Sort[Transpose[{Range[size], 
               estimatedtimescales}], #1[[2]] > #2[[2]] &][[All, 1]]]];
         {InnerQMatrix -> innerqmatrix,
          Accuracy -> estimate[[1]],
          Size -> size,
          OuterQMatrix -> Transpose[innerqmatrix].projection,
          TimeScales -> permutationmatrix.correctedtimescales,
          ListOfEstimationTimeFrames -> realtimescalesList,
          CorrectionMatrix -> correctionMatrix,
          LeftProcessVectors -> (Transpose[innerqmatrix].projection).DiagonalMatrix[
             equilibrium], PermutationMatrix -> permutationmatrix,
          Estimate -> estimate,
          Converged -> True,
          SimpleEstimation -> simpleEstimation,
          TimePoints -> realtimescalesList
          },
         estimatedtimescales = eigenvalues /. estimate[[2]];
         correctedtimescales = Exp[Log[estimatedtimescales]/realdt];
         correctionMatrix = 
          DiagonalMatrix[1/Exp[realt1/(2*realdt) Log[estimatedtimescales]]];
         qmatrix = (QMatrix.correctionMatrix /. estimate[[2]]);
         innerqmatrix = qmatrix.Transpose[permutationmatrix];
         permutationmatrix = 
          IdentityMatrix[size][[
           Sort[Transpose[{Range[size], 
               estimatedtimescales}], #1[[2]] > #2[[2]] &][[All, 1]]]];
         {InnerQMatrix -> innerqmatrix,
          Accuracy -> estimate[[1]],
          Size -> size,
          OuterQMatrix -> Transpose[innerqmatrix].projection,
          TimeScales -> permutationmatrix.correctedtimescales,
          ListOfEstimationTimeFrames -> realtimescalesList,
          CorrectionMatrix -> correctionMatrix,
          LeftProcessVectors -> (Transpose[innerqmatrix].projection).DiagonalMatrix[
             equilibrium], PermutationMatrix -> permutationmatrix,
          Estimate -> estimate,
          Converged -> False,
          SimpleEstimation -> simpleEstimation,
          TimePoints -> realtimescalesList
          }
         ]
    ]
    
    *)

TimeSeries /: EstimateTSFromTimeSeries[TimeSeries[v_?VectorQ, opts : OptionsPattern[TimeSeries]], optsSE___] := 
    Module[ {lt1, lt2, cM1, cM2},
        {lt1, lt2} = Lagtime /. {optsSE} /. {Lagtime -> {1, 2}};
        cM1 = #/Total[Flatten[#]] &[ToCountMatrix[TimeSeries[Drop[v, lt2 - lt1], Lagtime -> lt1], Lagtime -> lt1] // Normal];
        cM2 = #/Total[Flatten[#]] &[ToCountMatrix[TimeSeries[v, Lagtime -> lt2], Lagtime -> lt2] // Normal];
        EstimateTimeScales[cM1, cM2, Lagtime -> {lt1, lt2}]
    ]

EstimateTimeScales[m1_?MatrixQ, m2_?MatrixQ, opts___] :=
    Module[ {cl, size, m1Pos, m2Pos, lt1, lt2, tScales, evs, solPos, 
      solLim, sol, mu, muList, muListReal, muListComplex, 
      accuracy = 0.0001},
        size = Length[m1];
        If[ size == 1,
            {1},
            {lt1, lt2} = Lagtime /. {opts} /. {Lagtime -> {1, 2}};
            m1Pos = And @@ Thread[Eigenvalues[m1] > 0];
            m2Pos = And @@ Thread[Eigenvalues[m2] > 0];
            sol = Solve[Det[m1 - mu*m2] == 0, mu];
            muList = mu /. sol;
            muListReal = Select[muList, Im[Chop[#]] == 0 &];
            muListComplex = Select[muList, Im[Chop[#]] != 0 &];
            tScales = Reverse@Sort[Abs[muList]^(1/(lt1 - lt2))];
            solPos = And @@ Thread[tScales >= 0];
            solLim = And @@ Thread[tScales <= (1. + accuracy)];
            If[ m1Pos && m2Pos && solPos && solLim,
                Map[If[ Abs[1 - #] <= accuracy,
                        1,
                        If[ Abs[#] <= accuracy,
                            0,
                            #
                        ]
                    ] &,
                  tScales],
                cl = CrispClusterMatrix[
                   Join[Table[{ii}, {ii, 1, size - 2}], {{size - 1, 
                      size}}]][[1]];
                evs = FullEigensystem[TransitionMatrix[m1]];
                cl = evs[[3, Range[size - 1]]];
                Join[EstimateTimeScales[Chop[cl.m1.Transpose@cl], 
                  Chop[cl.m2.Transpose@cl], Lagtime -> {lt1, lt2}], {0}]
            ]
        ]
    ]
  
(* HMM Implementations *)

TimeSeries /: 
 EstimateHMM[
  TimeSeries[dtest_?VectorQ, opts : OptionsPattern[TimeSeries]], 
  optsHMM___] := 
 Module[ {n, nruns, mstat, kstates, BIG, BIGI, tend, 
   jlindex, alpha, beta, pstate, run, alphanorm, betanorm,
    A, b, powtab, irun, t, asum, sum, j, i, ii, bsum, lhood, lnorm, k, 
   bnew, denom, term, num, sumhold},
     n = Improvements /. {optsHMM} /. {Improvements -> 
         50}; (*number of levels of iterative improvement*)
     nruns = Guesses /. {optsHMM} /. {Guesses -> 
         1};(*number of runs from initial random guesses*)
     mstat = HiddenStates /. {optsHMM} /. {HiddenStates -> 
         3};(*number of states*)
     kstates = Max[dtest];(*number of outcomes*)
     
     (*data for the Hidden Markov Model should be stored in dtest*)
     BIG = 10^20;(*need this to handle the extremely low probabilities \
   encountered*)
     BIGI = N[1/BIG];
     SeedRandom[1];(*seed random numbers for reproducibility*)
     tend = Dimensions[dtest][[1]];
     jlindex = Table[Table[0, {i, 1, n}], {i, 1, nruns}];
     For[run = 1, run <= nruns, run++,
      alpha = Table[Table[0, {i, 1, mstat}], {i, 1, tend}];
      beta = Table[Table[0, {i, 1, mstat}], {i, 1, tend}];
      pstate = alpha;
      alphanorm = Table[0, {i, 1, tend}];
      betanorm = alphanorm;
      A = Table[Table[RandomReal[], {i, 1, mstat}], {i, 1, mstat}];
      b = Table[Table[RandomReal[], {i, 1, kstates}], {i, 1, mstat}];
      For[i = 1, i <= mstat, i++,
       A[[i, 1 ;; mstat]] = 
        A[[i, 1 ;; mstat]]/Total[A[[i, 1 ;; mstat]]];
       b[[i, 1 ;; kstates]] = 
        b[[i, 1 ;; kstates]]/Total[b[[i, 1 ;; kstates]]];
       ];
      
      (*Run the Hidden Markov Model*)
      powtab = Table[BIGI^(i - 6), {i, 0, 9}];
      For[irun = 1, irun <= n, irun++,
       alpha = 
        Table[Table[b[[ii, dtest[[i]]]], {ii, 1, mstat}], {i, 1, tend}];
       beta = Table[Table[1, {i, 1, mstat}], {i, 1, tend}];
       
       (*calculate alpha and beta and renormalize them-
       from Numerical Recipes*)
       For[t = 2, t <= tend, t++,
        asum = 0;
        For[j = 1, j <= mstat, j++,
         sum = 0;
         For[i = 1, i <= mstat, i++,
          sum += alpha[[t - 1, i]]*A[[i, j]]*b[[j, dtest[[t]]]];
       ];
         alpha[[t, j]] = sum;
         asum += sum;
         ];
        alphanorm[[t]] = alphanorm[[t - 1]];
        If[ asum < BIGI,
            ++alphanorm[[t]];
            For[j = 1, j <= mstat, j++,
             alpha[[t, j]] *= BIG;
       ];
        ];
        ];
       betanorm[[tend]] = 0;
       For[t = tend - 1, t >= 1, t--,
        bsum = 0;
        For[i = 1, i <= mstat, i++,
         sum = 0;
         For[j = 1, j <= mstat, j++,
          sum += A[[i, j]]*b[[j, dtest[[t + 1]]]]*beta[[t + 1, j]];
       ];
         beta[[t, i]] = sum;
         bsum += sum;
         ];
        betanorm[[t]] = betanorm[[t + 1]];
        If[ bsum < BIGI,
            ++betanorm[[t]];
            For[j = 1, j <= mstat, j++, beta[[t, j]] *= BIG;
       ];
        ];
        ];
       lhood = 0;
       For[i = 1, i <= mstat, i++,
        lhood += alpha[[1, i]]*beta[[1, i]];
     ];
       lnorm = alphanorm[[1]] + betanorm[[1]];
       While[lhood < BIGI,
        lhood *= BIG;
        lnorm++;
        ];
       jlindex[[run, irun]] = Log[lhood] + lnorm*Log[BIGI];
       For[t = 1, t <= tend, t++,
        sum = 0;
        For[i = 1, i <= mstat, i++,
         pstate[[t, i]] = alpha[[t, i]]*beta[[t, i]];
         sum += pstate[[t, i]];
         ];
        (*sum=lhood*BIGI^(lnorm-alphanorm[[t]]-betanorm[[t]]);*)
        For[i = 1, i <= mstat, i++,
        pstate[[t, i]] *= 1/sum;
      ];
        ];
       bnew = Table[Table[0, {ii, 1, kstates}], {i, 1, mstat}];
       For[i = 1, i <= mstat, i++,
        denom = 0;
        For[t = 1, t <= tend, t++,
         term = (alpha[[t, i]]*beta[[t, i]])/lhood*
           powtab[[alphanorm[[t]] + betanorm[[t]] - lnorm + 6]];
         (*term=(alpha[[t,i]]*beta[[t,i]])/lhood*
         BIGI^(alphanorm[[t]]+betanorm[[t]]-lnorm);*)
         denom += term;
         bnew[[i, dtest[[t]]]] += term;
         ];
        For[j = 1, j <= mstat, j++,
         num = 0;
         For[t = 1, t <= tend - 1, t++,
          num += alpha[[t, i]]*
            b[[j, dtest[[t + 1]]]]*beta[[t + 1, j]]*
            powtab[[alphanorm[[t]] + betanorm[[t + 1]] - lnorm + 
                6]]/lhood;
       (*num+=alpha[[t,i]]*b[[j,dtest[[t+1]]]]*beta[[t+1,j]]*
       BIGI^(alphanorm[[t]]+betanorm[[t+1]]-lnorm);*)
       ];
         num /= lhood;
         A[[i, j]] *= num/denom;
         ];
        sumhold = Total[A[[i, 1 ;; mstat]]];
        For[j = 1, j <= mstat, j++,
         A[[i, j]] /= sumhold;
      ];
        For[k = 1, k <= kstates, k++,
         bnew[[i, k]] *= 1/denom;
      ];
        ];
       b = bnew;
       ];
      ];
     {"MSM" -> A, "Output" -> b}
 ]

End[] (* End Private Context *)

EndPackage[]