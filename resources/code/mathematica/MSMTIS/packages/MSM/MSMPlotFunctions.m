(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

PositionStreamPlot::usage = "PositionStreamPlot[positions, values] renders a streamline plot with values as a vector assumed at 2D points stored in positions. If positions is a positive integer n, a regular spaced grid of n x n points is assumed.";
PositionStreamDensityPlot::usage = "PositionStreamDensityPlot[positions, values] renders a streamline plus density plot with values as a vector assumed at 2D points stored in positions. If positions is a positive integer n, a regular spaced grid of n x n points is assumed.";
PositionDensityPlot::usage = "PositionDensityPlot[positions, values] renders a density plot with values as a vector assumed at 2D points stored in positions. If positions is a positive integer n, a regular spaced grid of n x n points is assumed.";
PositionListPlot::usage = "PositionListPlot[positions, timeseries] plots a timeseries in 2D with state centers given in postitions. If positions is a positive integer n, a regular spaced grid of n x n points is assumed.";
PhiPsiPlot::usage = "PhiPsiPlot[positions, values] renders a PhiPsiPlot for e.g. a Ramachandran plot which allows arbitrary positions of states to be set and the underlying Voronoi tesselation is periodic.";

FieldStyle::usage = "FieldStyle is an option for PositionStreamPlot and PositionStreamDensityPlot and can take \"Gradient\" or \"Rotation\" as values.";

ValueRange::usage = "ValueRange is an option for PhiPsiPlot and specifies the range of values that correspond to the colorrange [0,1].";
Periodicity::usage = "Periodicity is an options for PhiPsiPlot and specifies the periodic length in x and y dimension as {x,y}.";
PeriodicWindow::usage = "PeriodicWindow is an option for PhiPsiPlot and specifies the points of the polygon enclosing the periodic window. This is for visual purposes only.";
ColorRange::usage = "ColorRange is an option for PositionDensityPlot, PositionStreamDensityPlot which specifies the values that correspond to the colorrange [0,1].";

Vertices::usage = "Vertices is an option for PositionDensityPlot, PositionStreamDensityPlot which specifies if the positions of vertices should be plotted.";
VertexStyle::usage = "VertexStyle is an option for PositionDensityPlot, PositionStreamDensityPlot which specifies the style of the plotted vertices.";
VertexLabel::usage = "VertexLabel is an option for PositionDensityPlot, PositionStreamDensityPlot which specifies if to place labels instead of vertex points. If a list of strings is given, the length must equal the number of states and are then used as label names.";


Options[PositionDensityPlot]         = Join[Options[ListDensityPlot],{Vertices->False,VertexLabel->False,VertexStyle->{}}];

Begin["`Private`"] (* Begin Private Context *) 

(* Plot functions *)

PositionListPlot[positionList_,TimeSeries[v_?VectorQ,opts:OptionsPattern[TimeSeries]] ,opts:OptionsPattern[DensityPlot]] :=
    Module[ {pointList},
        pointList = If[ NumberQ[positionList],
                        Table[{Mod[ii,positionList]+1,Floor[ii/positionList]+1},{ii,0,positionList^2-1}],
                        positionList
                    ];
        ListPlot[Map[pointList[[#]]&,valueList],opts]
    ]

PositionListPlot[positionList_,valueList_,opts:OptionsPattern[DensityPlot]] :=
    Module[ {pointList},
        pointList = If[ NumberQ[positionList],
                        Table[{Mod[ii,positionList]+1,Floor[ii/positionList]+1},{ii,0,positionList^2-1}],
                        positionList
                    ];
        ListPlot[Map[pointList[[#]]&,valueList],opts]
    ]
PositionDensityPlot[positionList_, valueList_, 
  opts : OptionsPattern[PositionDensityPlot]] :=
    Module[ {arrayToPlot, verticesLabel, vertexList, epilog, colorRange, 
      data},
        data = valueList;
        If[ (ColorFunctionScaling /. {opts} /. {ColorFunctionScaling -> 
              False}) == False,
            colorRange = ColorRange /. {opts} /. {ColorRange -> {0, 1}};
            colorRange = colorRange /. {Min -> Min[data], Max -> Max[data]};
            data = Rescale[data, colorRange];
        ];
        arrayToPlot = 
         If[ NumberQ[positionList],
             Partition[data, positionList],
             Transpose[Join[Transpose[positionList], {data}]]
         ];
        vertexList = {};
        If[ (Vertices /. {opts} /. Options[PositionDensityPlot]) === True,
            verticesLabel = VertexLabel /. {opts} /. {VertexLabel -> False};
            vertexList = 
             If[ (verticesLabel) === False,
                 Map[Point, positionList],
                 If[ ! VectorQ[verticesLabel],
                     MapIndexed[
                      Text[Style[#2[[1]], 
                         KeepOptions[Style, 
                          Sequence @@ (VertexStyle /. {opts})]], #1] &, 
                      positionList],
                     MapIndexed[
                      Text[Style[verticesLabel[[#2[[1]]]], 
                         KeepOptions[Style, 
                          Sequence @@ (VertexStyle /. {opts})]], #1] &, 
                      positionList]
                 ]
             ];
        ];
        epilog = 
         Join[{VertexStyle} /. {opts} /. {VertexStyle -> {}}, 
          vertexList, {Epilog} /. {opts} /. {Epilog -> {}}];
        ListDensityPlot[arrayToPlot, Epilog -> epilog, 
         RemoveRules[{VertexStyle, ColorRange, VertexLabel, Epilog, 
           Vertices}, opts]]
    ]
PositionStreamDensityPlot[positionList_,valueList_,opts:OptionsPattern[Join[{InterpolationOrder->3, FieldStyle->{Gradient}},Options[StreamDensityPlot]]]] :=
    Module[ {arrayToPlot,interpolation,func,x,y,nOrder},
        arrayToPlot = If[ NumberQ[positionList],
                          Transpose[Join[Transpose[Table[#&[{Mod[ii,positionList]+1,Floor[ii/positionList]+1}],{ii,0,positionList^2-1}]],{valueList}]],
                          Transpose[Join[Transpose[positionList],{valueList}]]
                      ];
        nOrder = InterpolationOrder/.Flatten[Join[List[opts], {InterpolationOrder->3}]];
        interpolation = Interpolation[arrayToPlot, InterpolationOrder->nOrder];
        func = OptionValue[FieldStyle]/.{"Gradient"->{D[interpolation[x,y],{x}],D[interpolation[x,y],{y}]},"Rotation"->{D[interpolation[x,y],{y}],-D[interpolation[x,y],{x}]}};
        StreamDensityPlot[Evaluate@func,{x,interpolation[[1,1,1]],interpolation[[1,1,2]]},{y,interpolation[[1,2,1]],interpolation[[1,2,2]]},Evaluate[RemoveRules[InterpolationOrder,RemoveRules[FieldStyle,opts]]]]
    ]
PositionStreamPlot[positionList_,valueList_,opts:OptionsPattern[Join[{InterpolationOrder->3, FieldStyle->{Gradient}},Options[StreamPlot]]]] :=
    Module[ {arrayToPlot,interpolation,func,x,y,nOrder},
        arrayToPlot = If[ NumberQ[positionList],
                          Transpose[Join[Transpose[Table[#&[{Mod[ii,positionList]+1,Floor[ii/positionList]+1}],{ii,0,positionList^2-1}]],{valueList}]],
                          Transpose[Join[Transpose[positionList],{valueList}]]
                      ];
        nOrder = InterpolationOrder/.Flatten[Join[List[opts], {InterpolationOrder->3}]];
        interpolation = Interpolation[arrayToPlot, InterpolationOrder->nOrder];
        func = OptionValue[FieldStyle]/.{"Gradient"->{D[interpolation[x,y],{x}],D[interpolation[x,y],{y}]},"Rotation"->{D[interpolation[x,y],{y}],-D[interpolation[x,y],{x}]}};
        StreamPlot[Evaluate@func,{x,interpolation[[1,1,1]],interpolation[[1,1,2]]},{y,interpolation[[1,2,1]],interpolation[[1,2,2]]},Evaluate[RemoveRules[InterpolationOrder,RemoveRules[FieldStyle,opts]]]]
    ]
PhiPsiPlot[states_,data_,opts___] :=
    Module[ {
            pShift,pWindow,boundList,checkpoints,statePositionsPeriodic,
            nearestFnc,size,shiftTable,optimizedStateList,plotDataPeriodic,
            plotDataOptimized,statePositionsOptimized,vtLabelPeriodic,
            vtLabelOptimized,valueRange,vtLabel
        },
        size = Length[states];
        valueRange = ValueRange/.{opts}/.{ValueRange->{0,1}};
        pWindow = PeriodicWindow/.{opts}/.{PeriodicWindow->{{-180,-180},{180,-180},{180,180},{-180,180}}};
        pShift = Periodicity/.{opts}/.{Periodicity->{360,360}};
        shiftTable = Flatten[Table[{x,y},{x,{-pShift[[1]],0,pShift[[1]]}},{y,{-pShift[[2]],0,pShift[[2]]}}],1];
        boundList = Transpose[{pWindow,RotateLeft[pWindow]}];
        statePositionsPeriodic = Flatten[Table[Map[shift+#&,states],{shift,shiftTable}],1];
        checkpoints = Flatten[Map[Table[(#[[2]]-#[[1]])*scale+#[[1]],{scale,1/480,1,1/480}]&,boundList*1.01],1];
        nearestFnc = Nearest[statePositionsPeriodic->Automatic];
        optimizedStateList = Union[Map[nearestFnc[#][[1]]&,checkpoints],Range[size*4+1,size*5]];
        plotDataPeriodic = Flatten[Table[data,{shift,shiftTable}],1];
        plotDataOptimized = plotDataPeriodic[[optimizedStateList]];
        statePositionsOptimized = statePositionsPeriodic[[optimizedStateList]];
        vtLabel = VertexLabel/.{opts}/.{VertexLabel->Range[size]};
        vtLabelPeriodic = Flatten[Table[vtLabel[[Range[size]]],{shift,shiftTable}],1];
        vtLabelOptimized = vtLabelPeriodic[[optimizedStateList]];
        PositionDensityPlot[
        	statePositionsOptimized,
        	plotDataOptimized,
            RemoveRules[{ValueRange,Periodicity,PeriodicWindow,VertexLabel},opts],
            PlotRange->{{-192,192},{-192,192},valueRange},
            ColorFunctionScaling->False,
            VertexLabel->vtLabelOptimized,
            PlotRangePadding->5,
            VertexStyle->FontSize->5,
            BaseStyle->{FontFamily->"Palatino",FontSize->9},
            FrameLabel->{"\[Phi]-Angle [Degrees]","\[Psi]-Angle [Degrees]"},
            Epilog->{White,Opacity[0.8],Polygon[Map[Join[#*(193/180),Reverse[#]]&,boundList]]}
        ]
    ]


End[] (* End Private Context *)

EndPackage[]