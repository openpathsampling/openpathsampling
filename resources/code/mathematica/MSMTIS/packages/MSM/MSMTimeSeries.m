(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

TimeSeries::usage = "TimeSeries[matrix, options] creates a TimeSeries Object.";
TimeSeriesQ::usage = "TimeSeriesQ[system] returns true if system is a TimeSeries, false otherwise.";

GetHistogram::usage = "GetHistogram[timeseries] returns a histogram of the states in timeseries.";

Options[TimeSeries] = {States->1,Lagtime->1,Centers->{}};

Begin["`Private`"] (* Begin Private Context *) 

(* TimeSeries & Trajectories *)

TimeSeries/:MakeBoxes[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]],StandardForm] := 
        If[ Length[v]<$DefaultShowVectorSize,
            InterpretationBox[#,TimeSeries[v,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalS]",RowBox[{Max[OptionValue[TimeSeries,States],Max[v]]}]],"[",Lagtime[TimeSeries[v,opts]],"]","\[Rule]","(",ToBoxes[TableForm[v,TableDirections->Row]],")"}]],
            InterpretationBox[#,TimeSeries[v,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalS]",RowBox[{Max[OptionValue[TimeSeries,States],Max[v]]}]],"[",Lagtime[TimeSeries[v,opts]],"]","\[Rule]","\[LeftAngleBracket]",Length[v],"\[RightAngleBracket]"}]]
        ]
TimeSeriesQ[x_] :=
    If[ Head[x]==TimeSeries&&VectorQ[x[[1]]],
        True,
        False
    ];
TimeSeries/:NormalForm[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]]] := 
        {v,opts}
TimeSeries/:Normal[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]]] := 
        v
TimeSeries/:MatrixForm[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]]] := 
        MatrixForm[v]
TimeSeries/:Length[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]]] := 
        Length[v]
TimeSeries/:Lagtime[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]]] := 
        Lagtime/.List[opts]/.{Lagtime->1}
TimeSeries/:States[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]]] := 
        Max[OptionValue[TimeSeries,States],Max[v]]
TimeSeries/:MatrixPlot[TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]],optsMP:OptionsPattern[MatrixPlot]] := 
        MatrixPlot[{v},optsMP,Frame->True,FrameTicks->{None,All,None,None},ImageSize->320,ColorFunction->"TemperatureMap"]
TimeSeries/:MakeBoxes[TimeSeries[v_?MatrixQ,opts:OptionsPattern[TimeSeries]],StandardForm] := 
            InterpretationBox[#,TimeSeries[v,opts]]&[RowBox[{SubscriptBox["\[DoubleStruckCapitalS]",RowBox[{Dimensions[v][[2]]}]],"\[LeftAngleBracket]",Length[v],"\[RightAngleBracket]"}]]

TimeSeries /: 
 GetCenters[
  TimeSeries[v_?VectorQ, opts : OptionsPattern[TimeSeries]]] := 
 Module[{},
  (Centers /. {opts}) /. {Centers -> 
     Range[States[TimeSeries[v, opts]]]}
  ]

TimeSeries /: 
 GetHistogram[
  TimeSeries[v_?VectorQ, opts : OptionsPattern[TimeSeries]], 
  optsGH___] := 
  Module[ {hist, histFull, st, centers, xPosition},
      xPosition = (XPosition /. {optsGH}) /. XPosition -> "State";
      st = States[TimeSeries[v, opts]];
      centers = GetCenters[TimeSeries[v, opts]];
      hist = Tally[v];
      histFull = 
       Join[hist, Map[{#, 0} &, Complement[Range[st], hist[[All, 1]]]]];
      histFull = Sort[histFull, #1[[1]] < #2[[1]] &];
      Which[
          xPosition == "State",
           histFull,
       xPosition == "Center",
       histFull /. Thread[Rule[Range[st], centers]],
       True,
       (* Wrong options *)
       histFull
       ]
  ]

End[] (* End Private Context *)

EndPackage[]