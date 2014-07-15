(* Mathematica Package *)
BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

NormalForm::usage = "NormalForm[system] returns the embedded matrix of a TransitionMatrix, CountMatrix, Clustermatrix, etc...";
Normal::usage = "Normal[system] returns the embedded matrix of a TransitionMatrix, CountMatrix, Clustermatrix, etc...";
Normalize::usage = "Normalize[system] normalizes the system according to the type of system.";
NormalizedQ::usage = "NormalizedQ[system] checks, if the matrix is properly normlized according to the type of system.";

Matrix::usage = "Matrix[system] returns the embedded matrix of a TransitionMatrix, CountMatrix, Clustermatrix, etc...";

MatrixPlot::usage = "MatrixPlot[system] plots the embedded matrix of system.";
Length::usage = "Length[system] returns the number of states of system.";
Lagtime::usage = "Lagtime[system] returns the embedded lagtime of system.";

LeftEigenvectors::usage = "LeftEigenvectors[system] returns the normalized set of left eigenvectors of system assuming detailed balance.";
RightEigenvectors::usage = "RightEigenvectors[system] return the normalized set of right eigenvectors of system assuming detailed balance.";
FullEigensystem::usage = "FullEigensystem[system] returns the normalized set of {eigenvalues, left eigenvectors, right eigenvectors} of system assuming detailed balance.";

Equilibrium::usage = "Equilibrium[system] returns the stationary distribution of system.";

AssumeSparse::usage = "AssumeSparse is an option of ToTimeSeries to specify if the underlying transition matrix is of a sparse structure in which case the generation of a time series is much faster.";
Symmetrize::usage = "Symmetrize[system] symmetrizes the dynamics in system by averaging over forward and backward dynamics.";

VoidFunction::usage = "VoidFunction is an option for Normalize to specify the function that is used for states that have not been visited.";
Corners::usage = "Corners is an option for PCCA algorithms to specify which states are used as corners in the estimation.";
Reversible::usage = "Reversible is an option for ToCountMatrix to specify if counting should be done forward and backward.";

Processes::usage = "Processes[correlationmodel] returns the process vectors (the left projected eigenvectors) of the CorrelationModel.";

XPosition::usage = "XPosition is an option for GetHistogram to specify the type of the returned states.";
MatrixType::usage = "MatrixType is an option for ToCountMatrix and ToCountTensor to specify the number type of the output.";

$TMPkgDefaultMatrixPlotOptions       = Sequence[Frame->True,ImageSize->240,ColorFunction:>(PaperStyle`$DefaultColorFunction[2-#*2]&),ImagePadding->{{40,5},{5,40}},FrameTicks->{All,None,None,All}];
$TMPkgDefaultArrayPlotOptions        = Sequence[Frame->True,FrameTicks->All,DataReversed->True,ColorFunction:>(PaperStyle`$DefaultColorFunction[2-#*2]&)];
$DefaultShowMatrixSize = 5;
$DefaultShowVectorSize = 5;
VarianceMatrix::usage = "VarianceMatrix[countmatrix] returns the matrix of variances in all transition matrix elements around the mean value given the specified Distribution in the CountMatrix.";
VarianceTensor::usage = "VarianceTensor[countmatrix, states] returns the a list of correlation matrix between all transition matrix elements in the same row around the mean value given the specified Distribution in the CountMatrix for all states specied as a vector.";
MeanMatrix::usage = "MeanMatrix[countmatrix] the average TransitionMatrix from the countmatrix given the specified distribution";

MaximumLikelihoodMatrix::usage = "MaximumLikelihoodMatrix[countmatrix] returns the maxmimum likelihood estimator for the transitionmatrix that produced countmatrix.";
NewStateProbabilityVector::usage = "NewStateProbabilityVector[countmatrix] returns the probability to discover a new state per state assuming that there exists one hidden state that has not yet been visited from any state.";
DrawTransitionMatrix::usage = "DrawTransitionMatrix[countmatrix] draws an instantce of a transitionmatrix according to the distribution given by counts specified in countmatrix.";

Prior::usage = "Prior is an option for CountMatrix that specifies the number of prior counts to be used in the Bayesian computations. Prior can be a number or a matrix of the same size as the related CountMatrix."

MaximumLikelihoodMatrix::NYI = "Not Yet Implemented";
MeanMatrix::NYI = "Not Yet Implemented";
VarianceMatrix::NYI = "Not Yet Implemented";
DrawTransitionMatrix::NYI = "Not Yet Implemented";

States::usage = "States is an option for TimeSeries to specify the real number of states otherwise the number of states equals the largest found state in the TimeSeries.";
Distribution::usage = "Distribution is an option for CountMatrix an specifies the used Distribution of entries. Default is a DirichletDistribution.";

Centers::usage = "Centers is an options for TimeSeries to specify cluster centers for each discrete state."

MaximalResidue::usage = "MaximalResidue is an option for DirichletWithDBDistribution to specify the maximal residue for the optimization algorithm.";

Clustering::usage = "Clustering is an option for ToCountTensor to specify that clustered states are to be used for the outer states.";

Begin["`Private`"] (* Begin Private Context *) 

End[] (* End Private Context *)

EndPackage[]