(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

Trajectory::usage = "";
TrajectoryQ::usage = "";

PlotTrajectory1D::usage = "";

Options[Trajectory]                  = {};

Begin["`Private`"] (* Begin Private Context *) 

Trajectory/:MakeBoxes[Trajectory[v_?ListQ,opts:OptionsPattern[Trajectory]],StandardForm] := 
        InterpretationBox[#,Trajectory[v,opts]]&[RowBox[{"\[DoubleStruckCapitalT]race","\[LeftAngleBracket]",Dimensions[v][[1]],"\[RightAngleBracket]"}]]

TrajectoryQ[x_] :=
    If[ Head[x]==Trajectory&&ListQ[x[[1]]],
        True,
        False
    ];
Trajectory/:NormalForm[Trajectory[v_?ListQ,opts:OptionsPattern[Trajectory]]] := 
        {v,opts}
Trajectory/:Normal[Trajectory[v_?ListQ,opts:OptionsPattern[Trajectory]]] := 
        v
Trajectory/:MatrixForm[Trajectory[v_?ListQ,opts:OptionsPattern[Trajectory]]] := 
        MatrixForm[v]
Trajectory/:Length[Trajectory[v_?ListQ,opts:OptionsPattern[Trajectory]]] := 
        Length[v]
        
        
PlotTrajectory1D[trace_?TrajectoryQ, opts___] := 
 Module[{stepSize, timeScaling, hist},
  stepSize = StepSize /. {opts} /. {StepSize -> 20};
  timeScaling = TimeScaling /. {opts} /. {TimeScaling -> 1};
  hist = SmoothKernelDistribution[
    Take[pullingData[[1]], {1, Length[trace], stepSize}]];
  backdensity = 
   DensityPlot[PDF[hist, y], {x, 0, 60}, {y, 12, 15.5}, 
    ColorFunction -> (GrayLevel[Min[1, Max[1 - #, 0]]] &), 
    PlotPoints -> 40];
  Show[ListPlot[
    Take[Transpose[{Range[0, Length[trace] - 1]*timeScaling, 
       trace[[1]]}], {1, Length[trace], stepSize}], 
    RemoveRules[{TimeScaling, StepSize}, opts], 
    PlotStyle -> Directive[Black, Thin], AspectRatio -> 0.10, 
    ImageSize -> 800, FrameStyle -> Thickness[0.0005]
    ]]
  ]

End[] (* End Private Context *)

EndPackage[]