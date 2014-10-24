(* ::Package:: *)

(* ProjectManagement Package *)
(* Version 0.4 *)

(* Last edited on Jul 26th, 2010 *)

BeginPackage["PaperStyle`"];


(* Messages *)

$Colors::usage = "A list of colors used in the package";
$DefaultColorFunction::usage = "";
$DefaultColorFunctionWhite::usage = "";
$DefaultColorFunctionOpac::usage = "";
$DefaultLineTypes::usage = "";
$DefaultPlotStyleDashed::usage = "";
$DefaultPlotStyleDotted::usage = "";
$DefaultPlotStyleSolid::usage = "";
$DefaultFillingStyle::usage = "";

$PaperColorBlue::usage = "";
$PaperColorGreen::usage = "";
$PaperColorRed::usage = "";
$PaperColorYellow::usage = "";
$PaperColorDarkRed::usage = "";
$PaperColorDarkBlue::usage = "";
$PaperColorViolet::usage = "";
$PaperColorLightBlue::usage = "";
$PaperColorBlack::usage = "";
$PaperFillingStyles::usage = "";

$ColorNameReplacements::usage = "";

LogTicks::usage = "";
LegendPlot::usage = "";
LegendPlotList::usage = "";
LegendTextPlot::usage = "";

Solid::usage = "";
Solid = Dashing[{}];

SetPlotStyle::usage = "";
ToFillingStyle::usage = "";
ToPlotStyle::usage = "";
ToImagePadding::usage = "";
ToImageSize::usage = "";
GetFillingStyle::usage = "";
GetPlotStyle::usage = "";

PLeft::usage = "";
PRight::usage = "";
PTop::usage = "";
PBottom::usage = "";
Spacing::usage = "";

MediumX::usage = "";
LargeX::usage = "";

Orientation::usage = "";
Height::usage = "";
Size::usage = "";
Step::usage = "";
RangeText::usage = "";
Title::usage = "";

cm::usage = "";
inch::usage = "";
points::usage = "";

MakeSubFigureGrid::usage = "";

$DefaultFillingOpacity::usage = "";


Begin["`Private`"];

(* Functions to handle projects *)


InvertColor[col_]:=Module[{},
	Which[Head[col]===RGBColor,
		RGBColor[1-col[[1]],1-col[[2]],1-col[[3]]],
		Head[col]===GrayLevel,
		GrayLevel[1-col[[1]]],
		True,
		White
	]
]

$DefaultFillingOpacity = 0.4;

$ColorNameReplacements = {"Red" -> 2, "Black" -> 1, "Yellow" -> 3, "Green" -> 4, "Blue" -> 5};
$ColorScheme = {
"Normal" -> {  
	"DarkRed" -> RGBColor[0.544549/3, 0, 0],
	"Red" -> RGBColor[0.544549, 0, 0],
	"Yellow" -> RGBColor[0.791638, 0.462715, 0.0744335],
	"Green" ->  RGBColor[0.525582, 0.607904, 0.164706],
	"LightBlue" -> RGBColor[0.197101, 0.559991, 0.726802],
	"Blue" -> RGBColor[0.0689708, 0.331441, 0.447944],
	"DarkBlue" -> RGBColor[0, 0.03, 0.405173],
	"Violet" -> RGBColor[0.236545, 0.107179, 0.390799],
	"Black" -> Black
	},
"BlackboardInverse" -> {
	"DarkRed" -> InvertColor[RGBColor[0.544549/3, 0, 0]],
	"Red" -> InvertColor[RGBColor[0.544549, 0, 0]],
	"Yellow" -> InvertColor[RGBColor[0.791638, 0.462715, 0.0744335]],
	"Green" ->  InvertColor[RGBColor[0.525582, 0.607904, 0.164706]],
	"LightBlue" -> InvertColor[RGBColor[0.197101, 0.559991, 0.726802]],
	"Blue" -> InvertColor[RGBColor[0.0689708, 0.331441, 0.447944]],
	"DarkBlue" -> InvertColor[RGBColor[0, 0.03, 0.405173]],
	"Violet" -> InvertColor[RGBColor[0.236545, 0.107179, 0.390799]],
	"Black" -> Black
	}
}

$DefaultPlotSettings = 
	{
		BaseStyle -> {FontSize->15, FontFamily->"Palatino"},
		FrameTicksStyle->Thickness[0.002],
		FrameStyle->Thickness[0.002],
		Frame->True,
		Axes->False
	};
$PlotFunctionsToSetToBaseStyle = 
	{
		Plot,
		ListPlot,
		MatrixPlot,
		Graphics,
		LogPlot,
		LogLogPlot,
		LogLinearPlot,
		ListLogPlot,
		ListLogLinearPlot,
		ListLogLogPlot
	};

$PaperColorDarkRed = RGBColor[0.544549/3, 0, 0];
$PaperColorRed = RGBColor[0.544549, 0, 0];
$PaperColorYellow = RGBColor[0.791638, 0.462715, 0.0744335];
$PaperColorGreen = RGBColor[0.525582, 0.607904, 0.164706];
$PaperColorLightBlue = RGBColor[0.197101, 0.559991, 0.726802];
$PaperColorBlue = RGBColor[0.0689708, 0.331441, 0.447944];
$PaperColorDarkBlue = RGBColor[0, 0.03, 0.405173];
$PaperColorViolet = RGBColor[0.236545, 0.107179, 0.390799];
$PaperColorBlack = Black;

$PaperStyleColors = 
	{
		$PaperColorBlack,
		$PaperColorRed,
		$PaperColorYellow,
		$PaperColorGreen,
		$PaperColorBlue,
		$PaperColorDarkBlue,
		$PaperColorViolet,
		$PaperColorBlack,
		$PaperColorRed,
		$PaperColorYellow,
		$PaperColorGreen,
		$PaperColorBlue,
		$PaperColorDarkBlue,
		$PaperColorViolet,
		$PaperColorBlack
	};

$DefaultColorFunctionFull = Function[x,Blend[{{0,$PaperColorDarkRed},{0.2,$PaperColorRed},{0.4,$PaperColorYellow},{0.6,$PaperColorGreen},{0.8,$PaperColorBlue},{1.0,$PaperColorDarkBlue}},x]];
$DefaultColorFunction = Function[x,Blend[{{0,$PaperColorDarkRed},{0.1,$PaperColorRed},{0.4,$PaperColorYellow},{0.55,$PaperColorGreen},{0.7,$PaperColorBlue},{1.0,$PaperColorDarkBlue}},x]];
$DefaultColorFunctionWhite = Function[x,Blend[{{0,$PaperColorDarkRed},{0.1,$PaperColorRed},{0.4,$PaperColorYellow},{0.55,$PaperColorGreen},{0.8,$PaperColorBlue},{1.0,White}},x]];
$DefaultColorFunctionOpac = Function[x,Blend[{{0,$PaperColorDarkRed},{0.1,$PaperColorRed},{0.4,$PaperColorYellow},{0.55,$PaperColorGreen},{0.8,$PaperColorBlue},{1.0,RGBColor[1,1,1,0]}},x]];
$DefaultColorFunctionPositive = $DefaultColorFunction;
$DefaultColorFunctionDiscrete = $DefaultColorFunction;

$DefaultLineTypes = {Dashing[{}], Dashed, Dotted};


$Colors = $PaperStyleColors;

ToPlotStyle[colIndices_,style_] := MapIndexed[Directive[#,Thick,If[VectorQ[style],style[[Mod[#2[[1]],Length[style],1]]],style]]&, $Colors[[colIndices]]/.$ColorNameReplacements];
ToFillingStyle[colIndices_] := MapIndexed[#2[[1]] -> Directive[Opacity[$DefaultFillingOpacity],#]&, $Colors[[colIndices]]];
ToImagePadding[opts___] := Module[{pReplacements,padding},
	pReplacements = {None->0, Small->5, Medium->20, MediumX->25, Large->40, LargeX-> 45};
	padding = ({{PLeft,PRight},{PBottom,PTop}} /. List[opts]) /. {PLeft ->Medium, PRight->Small, PBottom->Medium, PTop->Small}  /. pReplacements;
	padding
];
ToImageSize[size_] := Module[{sReplacements},
	sReplacements = {cm -> 72 / 2.54, inch -> 72, points -> 1};
	size /. sReplacements
]

$DefaultPlotStyleSolid = ToPlotStyle[Range[Length[$Colors]],Solid];
$DefaultPlotStyleDashed = ToPlotStyle[Range[Length[$Colors]],Dashed];
$DefaultPlotStyleDotted = ToPlotStyle[Range[Length[$Colors]],Dotted];
$DefaultFillingStyle = MapIndexed[#2[[1]] -> Directive[Opacity[0.2],#]&, $Colors];
$PaperFillingStyles = Map[Directive[Opacity[0.2],#]&, $Colors];

$ColorNameReplacements = {"Black"->1, "Red"->2, "Yellow"->3, "Green"->4, "Blue"->5};

Map[Table[SetOptions[#,opt],{opt,$DefaultPlotSettings}]&,$PlotFunctionsToSetToBaseStyle];

SetPlotStyle[colIndices_,style_] := Module[{col},
		col = colIndices /. $ColorNameReplacements;
		Sequence[{FillingStyle->ToFillingStyle[col],PlotStyle->ToPlotStyle[col,style]}]
];

GetFillingStyle[colIndex_] := Module[{col},
		col = colIndex /. $ColorNameReplacements;
		$PaperFillingStyles[[col]]
];
GetPlotStyle[colIndex_, style_] := Module[{col},
		col = colIndex /. $ColorNameReplacements;
		ToPlotStyle[style][[col]]
];


Map[
	SetOptions[#,PlotStyle->$DefaultPlotStyleSolid];
	SetOptions[#,FillingStyle->$DefaultFillingStyle];&
,
	{
		Plot,
		ListPlot,
		LogPlot,
		LogLogPlot,
		LogLinearPlot,
		ListLogPlot,
		ListLogLinearPlot,
		ListLogLogPlot
	}
];

Map[
	SetOptions[#,Joined->True];&
,
	{
		ListPlot,
		ListLogPlot,
		ListLogLinearPlot,
		ListLogLogPlot
	}
];


(* Utility Functions for Plotting *)


Options[LogTicks] = {Scaling->Identity,TickNumbers->{1,2,3,4,5,6,7,8,9}, MarkNumbers->{1,2,5},TickStyle->{},MarkStyle->{},MarkLength->{0.01,0},TickLength->{0.005,0},MarkTextFunction->ToString};
Options[LegendPlot] = {ColorFunction->$DefaultColorFunction, Height->0.1, Ticks->{}, RangeText->{"",""}, Title->None, Step->0.01, Orientation->Horizontal};


LogTicks[range_,OptionsPattern[LogTicks]] := Module[{},
	scaling = OptionValue[Scaling];
	logRange = {Floor[Log[10,range[[1]]]],Ceiling[Log[10,range[[2]]]]};
	ticksNumbers = OptionValue[TickNumbers];
	markNumbers = OptionValue[MarkNumbers];
	allTicks = Flatten[Table[ticksNumbers*10^ll,{ll, logRange[[1]], logRange[[2]]}]];
	allTicks = Select[allTicks,#>=range[[1]]&&#<=range[[2]]&];
	allMarks = Flatten[Table[markNumbers*10^ll,{ll, logRange[[1]], logRange[[2]]}]];
	allMarks = Select[allMarks,#>=range[[1]]&&#<=range[[2]]&];
	allMarks = Intersection[allMarks,allTicks];
	tickStyle = OptionValue[TickStyle];
	markStyle = OptionValue[MarkStyle];
	tickLength = OptionValue[TickLength];
	tickLength = Switch[Length[#],0,{#,0},1,{#,#},_,#]&[tickLength];
	markLength = OptionValue[MarkLength];
	markLength = Switch[Length[#],0,{#,0},1,{#,#},_,#]&[markLength];	
	markFunction = OptionValue[MarkTextFunction];
	Map[
		If[MemberQ[allMarks,#],
			{scaling[N[#]],markFunction[#],markLength,Directive[markStyle]},
			{scaling[N[#]],"",tickLength,Directive[tickStyle]}
		]&
	,
		allTicks
	]
]


LegendPlotList[opts___] := Module[{size,textSize,titleText, minText, maxText, step, height, cFunc, orientation, legend, additionalTicks},
	step=Step/.{opts}/.{Step->1/2};
	titleText=Title/.{opts}/.{Title->""};
	{minText,maxText} = RangeText/.{opts}/.{RangeText->{"Min","Max"}};
	additionalTicks = Ticks/.{opts}/.{Ticks->{}};
	height = Height/.{opts}/.{Height->0.15};
	cFunc=ColorFunction/.{opts}/.{ColorFunction->$DefaultColorFunction};
	orientation = (Orientation/.{opts})/.{Orientation->Horizontal};
	textSize=FontSize/.{opts}/.{FontSize->15};
	size = Size/.{opts}/.{Size->200};
	If[orientation===Horizontal,
		legend = Graphics[
			Table[
				{$Colors[[Round[x/step]+2]],Rectangle[{x-step/2,0},{x+step/2-0.0001,height}]}
			,{x,0,1,step}]
			,Frame->True,FrameTicks->{Join[{{0,minText},{1,maxText}},additionalTicks],None,None,None}
			,BaseStyle->{FontFamily->"Palatino",textSize},ImageSize->{size,Automatic},FrameLabel->{titleText,None,None,None},PlotRangePadding->0.02
		];
		,
		legend = Graphics[
			Table[
				{$Colors[[Round[x/step]+2]],Rectangle[{0,y-step/2},{height,y+step/2}]}
			,{y,0,1,step}]
			,Frame->True,FrameTicks->{None,None,None,Join[{{0,minText},{1,maxText}},additionalTicks]}
			,BaseStyle->{FontFamily->"Palatino",textSize},ImageSize->{Automatic,size},PlotRangePadding->0.02,FrameLabel->{None,None,None,titleText}
		];
	];
	legend
];

LegendPlot[opts___] := Module[{size,textSize,titleText, minText, maxText, step, height, cFunc, orientation, legend, additionalTicks},
	step=Step/.{opts}/.{Step->0.01};
	titleText=Title/.{opts}/.{Title->""};
	{minText,maxText} = RangeText/.{opts}/.{RangeText->{"Min","Max"}};
	additionalTicks = Ticks/.{opts}/.{Ticks->{}};
	height = Height/.{opts}/.{Height->0.15};
	cFunc=ColorFunction/.{opts}/.{ColorFunction->$DefaultColorFunction};
	orientation = (Orientation/.{opts})/.{Orientation->Horizontal};
	textSize=FontSize/.{opts}/.{FontSize->15};
	size = Size/.{opts}/.{Size->200};
	If[orientation===Horizontal,
		legend = Graphics[
			Table[
				{cFunc[x],Rectangle[{x-step/2,0},{x+step/2,height}]}
			,{x,0,1,step}]
			,Frame->True,FrameTicks->{Join[{{0,minText},{1,maxText}},additionalTicks],None,None,None}
			,BaseStyle->{FontFamily->"Palatino",textSize},ImageSize->{size,Automatic},FrameLabel->{titleText,None,None,None},PlotRangePadding->0.02
		];
		,
		legend = Graphics[
			Table[
				{cFunc[y],Rectangle[{0,y-step/2},{height,y+step/2}]}
			,{y,0,1,step}]
			,Frame->True,FrameTicks->{None,None,None,Join[{{0,minText},{1,maxText}},additionalTicks]}
			,BaseStyle->{FontFamily->"Palatino",textSize},ImageSize->{Automatic,size},PlotRangePadding->0.02,FrameLabel->{None,None,None,titleText}
		];
	];
	legend
];

LegendTextPlot[colIndices_,styles_,names_,opts___] := 
	Module[{position,spacing,textStyle,curves,col},
		position=Position/.List[opts]/.Position->{0.08,0.83};
		spacing=Spacing/.List[opts]/.Spacing->0.05;
		textStyle=TextStyle/.List[opts]/.TextStyle->{FontSize->14};
		col= colIndices /. $ColorNameReplacements;
		curves=Transpose[{col,styles,names}];
		MapIndexed[{$Colors[[#[[1]]]],Thick,#[[2]],
			Line[{Scaled[position+{0,-(#2[[1]]-1)*spacing}],Scaled[position+{0.08,-(#2[[1]]-1)*spacing}]}],
			Black,Text[Style[#[[3]],textStyle],Scaled[position+{0.10,-(#2[[1]]-1)*spacing}],{-1,0}]}&
		,curves]
	]


AddLeftText[grid_,texts_]:=Transpose[Join[
{Map[Rotate[Style[#,FontSize->20,FontFamily->"Palatino"],90Degree]&,texts]},Transpose[grid]]]
AddRightLegend[grid_,legend_]:=Transpose[Join[
Transpose[grid],{Table[If[nn==1,legend,SpanFromAbove],{nn,1,Length[grid]}]}]]


MakeSubFigureGrid[graphics_,width_]:=Module[{},
	g=Partition[
		PadRight[
			Flatten[
				MapIndexed[
					{Graphics[Text[Style["("<>FromCharacterCode[#2[[1]]+96]<>")",FontSize->18,FontFamily->"Palatino"],{0,0}],Frame->False,ImageSize->40],#1}&,
					Flatten[graphics]
				]
			],
			Ceiling[Length[Flatten[graphics]]/width]*width*2,Null],
		width*2];
	Grid[g]
]


End[ ];


EndPackage[ ];
