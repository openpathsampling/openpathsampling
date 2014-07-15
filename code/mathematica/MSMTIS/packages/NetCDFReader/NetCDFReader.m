(* ::Package:: *)

(* TransitionMatrix Package *)
(* Version 0.3 *)

(* Last edited on Apr 6th, 2009 *)

BeginPackage["NetCDFReader`"];


(* Messages *)

OpenNetCDFFile::usage = "Open a NetCDF File (.nc) and return an NetCDFFile-Object, that contains information about all containing datasets";
NetCDFFile::usage = "";
GetDataset::usage = "";
ListDatasets::usage = "";
ListDatasetsFull::usage = "";
WriteNetCDFFile::usage = "";

$ChkBox::usage = "";
$ListOfUnits::usage = "";
$DataDescription::usage = "";
OpenNetCDFExport::usage = "";
ExportAsNetCDF::usage = "";
GetDatasetNetCDF::usage = "";

View::usage = "";


(* All options *)

Options[OpenNetCDFFile]                  = {};
Options[NetCDFFile]                  = {};


Begin["`Private`"];


(* Functions to handle options *)

ReplaceRules[rule_,f___] := 
		Apply[Sequence,Join[{rule},Select[List[f],#[[1]]=!=rule[[1]]&]]]
UnionRules[rules1___] := 
		Apply[Sequence,Union[List[rules1][[All,1]]]/.Map[#[[1]]->#&,List[rules1]]]
RemoveRules[rule_,f___] := 
		Apply[Sequence,Select[List[f],#[[1]]=!=rule&]]


(* NetCDFFile functions *)

NetCDFFile/:MakeBoxes[NetCDFFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[NetCDFFile]],StandardForm] := 
		InterpretationBox[#,NetCDFFile[f,m,d,opts]]&[RowBox[{SubscriptBox["NetCDF",RowBox[{Length[d]}]],"(",StringCases[f,RegularExpression["[\\w._+-]+$"]][[1]],")"}]]


(* Open a dialog to save data into a NetCDFFile *)

OpenNetCDFExport[]:=Module[{listOfNumberObjects,arrayObjects,listOfObjects},
listOfObjects = Names["Global`*"];
listOfObjects=Select[listOfObjects,Length[Dimensions[Symbol[#]]]>0&];
listOfNumberObjects=Select[listOfObjects,DataArrayQ[Symbol[#]]&];
listOfNumberObjects=Select[listOfNumberObjects,StringFreeQ[#,"$"]&];
listOfNumberObjects=Select[listOfNumberObjects,!MemberQ[{"$IList","$ChkBox","listOfObjects","$ListOfUnits"},#]&];
arrayObjects=Map[#->ArrayInformation[Symbol[#]]&,listOfNumberObjects];
$IList=Map[GetChkBoxEntry,arrayObjects];
$ChkBox=Drop[$ChkBox,Complement[Range[Length[$ChkBox]],$IList]];
$IList=Map[GetChkBoxEntry,arrayObjects];
Grid[Join[{{Item["\[DownArrow]",ItemSize-> {Full,4}],
Item[Button["Export NetCDF",ExportAsNetCDF[SystemDialogInput["FileSave","data.nc"],$DataDescription],Alignment->{Center,Top},Appearance->"Palette",ImageSize->{100,15}, BaseStyle->{"Palatino",9},Method->"Queued"],Alignment->{Right,Top}],SpanFromLeft,
Item[InputField[Dynamic[$DataDescription],String,ImageSize->{500,12},ContentPadding->False,Enabled->True],Alignment->{Left,Top}],SpanFromLeft}
},
MapIndexed[(
{Checkbox[Dynamic[$ChkBox[[$IList[[#2[[1]]]],2]]]],
"\""<>#[[1]]<>"\"",
Item[("Type"/.#[[2]])<>ToShapeString[("Shape"/.#[[2]])],Alignment->Left],
Item[InputField[Dynamic[$ChkBox[[$IList[[#2[[1]]]],3]]],String,ImageSize->{100,12},ContentPadding->False]],
Item[InputField[Dynamic[$ChkBox[[$IList[[#2[[1]]]],4]]],String,ImageSize->{200,12},ContentPadding->False]],
Item[InputField[Dynamic[$ChkBox[[$IList[[#2[[1]]]],5]]],String,ImageSize->{70,12},ContentPadding->False]],Item[PopupMenu[Dynamic[$ChkBox[[$IList[[#2[[1]]]],5]]],$ListOfUnits,Style["[other]",FontSize->9],ContentPadding->False,ImageSize->{90,14},Alignment->{Center,Center},BaseStyle->{FontFamily->"Gill Sans",FontSize->10},Appearance->"Palette"],Alignment->{Right,Center}],
Item[ToString[("Memory"/.#[[2]])]<>" Bytes",Alignment->Right],Button["Remove",Remove[Evaluate[#[[1]]]],Alignment->{Center,Top},Appearance->"Palette",ImageSize->{50,14}, BaseStyle->{"Palatino",9},Method->"Queued",ContentPadding->False]})&,arrayObjects]],BaseStyle->{FontFamily->"Gill Sans",FontSize->10},Spacings->{0.5,0.2},Background->{None,{{LightGray,None}}},Frame->{None,All},FrameStyle->Thin,ItemSize->{Full,2},Alignment->{Center,Center}]
]


(* Access NetCDF Files *)

OpenNetCDFFile[filename_?StringQ]:=
		Module[{metadata,datasets,annotations,formats,dimensions},
			datasets = Import[filename,"NetCDF"];
			annotations = Import[filename,{"NetCDF","Annotations"}];
			metadata = Import[filename,{"NetCDF","Metadata"}];
			dimensions = Import[filename,{"NetCDF","Dimensions"}];
			formats = Import[filename,{"NetCDF","DataFormat"}];
			NetCDFFile[filename,metadata,Transpose[{datasets,dimensions,formats,annotations}]]
		]

WriteNetCDFFile[filename_?StringQ, datasets_]:=
		Module[{metadata,annotations,formats,dimensions},
			Export[filename,datasets[[All,2]],{"Datasets",datasets[[All,1]]}]
		]

NetCDFFile/:GetDataset[NetCDFFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[NetCDFFile]],set_?StringQ] := 
		Module[{},
			If[Apply[Or,StringMatchQ[d[[All,1]],set]],
				Import[f,{"Datasets",set}],
				Print["Error: Dataset not found"];
			]
		]

NetCDFFile/:GetDataset[NetCDFFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[NetCDFFile]],index_?IntegerQ] := 
		Module[{},
			If[index>0&&index<=Length[d],
				Import[f,{"Datasets",index}],
				Print["Error: Dataset not found"];
			]
		]

NetCDFFile/:ListDatasets[NetCDFFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[NetCDFFile]]] := 
		Module[{},
			d[[All,1]]
		]

NetCDFFile/:ListDatasetsFull[NetCDFFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[NetCDFFile]]] := 
		Module[{},
			d
		]


IntegerArrayQ[obj_]:=If[ArrayQ[obj],
And@@Flatten[Map[IntegerQ,obj,{-1}]],
False]

TypeArrayQ[type_,obj_]:=If[ArrayQ[obj],
And@@Flatten[Map[type,obj,{-1}]],
False]

RealQ[x_]:=(Im[x]==0)

ComplexQ[x_]:=(Im[x]!=0)

RealArrayQ[obj_]:=If[ArrayQ[obj],
And@@Flatten[Map[RealQ,obj,{-1}]],
False]

ComplexArrayQ[obj_]:=If[ArrayQ[obj],
And@@Flatten[Map[ComplexQ,obj,{-1}]],
False]

ArrayNumberType[obj_]:=
If[ArrayQ[obj],
If[
NumericArrayQ[obj],
Which[
IntegerArrayQ[obj],"Integer",
RealArrayQ[obj],"Real",
ComplexArrayQ[obj],"Complex",
True,"Numeric"
],
If[StringArrayQ[obj],"String",
"Symbolic"
]],
"NoArray"]
ArrayDimensions[obj_]:=If[ArrayQ[obj],
Dimensions[obj]//Length,
0]
ArrayShape[obj_]:=If[ArrayQ[obj],
Dimensions[obj],
{}]
ArraySize[obj_]:=If[ArrayQ[obj],
Times@@Dimensions[obj],
0]
ArrayInformation[obj_]:={
"Type"->ArrayNumberType[obj],
"Dimensions"->ArrayDimensions[obj],
"Shape"->ArrayShape[obj],
"Size"->ArraySize[obj],
"Memory"->ByteCount[obj]
}

NumericVectorQ[obj_]:=If[VectorQ[obj],
And@@Map[NumericQ,obj],
False]
DataArrayQ[obj_]:=If[ArrayQ[obj],
(And@@Flatten[Map[NumericQ,obj,{-1}]])
||(And@@Flatten[Map[StringQ,obj,{-1}]])
,
False]
StringArrayQ[obj_]:=If[ArrayQ[obj],
And@@Flatten[Map[StringQ,obj,{-1}]],
False]
NumericArrayQ[obj_]:=If[ArrayQ[obj],
And@@Flatten[Map[NumericQ,obj,{-1}]],
False]

ToShapeString[shape_]:="["<>StringJoin[Drop[Flatten[Map[{ToString[#],"\[Cross]"}&,shape]],-1]]<>"]"

NewChkBoxEntry:={#[[1]],False,#[[1]],"[description]",$ListOfUnits[[1,1]]}&

GetChkBoxEntry[obj_]:=Module[{chkPos},
chkPos=Position[$ChkBox[[All,1]],obj[[1]]];
If[Length[chkPos]>0,
chkPos[[1,1]],
AppendTo[$ChkBox,NewChkBoxEntry[obj]];
GetChkBoxEntry[obj]
]
]

ExportAsNetCDF[filename_,title_]:=Module[
{obList,dataList,annotationsList,metaList},
obList=Map[#[[1]]->Join[ArrayInformation[Symbol[#[[1]]]],{"name"->#[[3]],"desc"->#[[4]],"units"->#[[5]]}]&,(Select[$ChkBox,#[[2]]&])];
dataList=Map[("name"/.#[[2]])->Symbol[#[[1]]]&,obList];
annotationsList=Map[{"units"->("units"/.#[[2]]),"description"->(("desc"/.#[[2]])/."[description]"-> "[none]"),"mathematica_symbol"->#[[1]]}&,obList];
metaList={"title"->title,"notebook"-> NotebookFileName[],"application"->$Version,"date"->DateString[]};
Export[filename,{dataList,annotationsList,metaList},
{"NetCDF",{"Datasets","Annotations","Metadata"}}];
]

NetCDFFile/:RestoreDataset[NetCDFFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[NetCDFFile]],set_?StringQ,var_?StringQ] := 
		Module[{},
			If[Apply[Or,StringMatchQ[d[[All,1]],set]],
				data=GetDataset[NetCDFFile[f,m,d,opts],set];
			ToExpression[var<>"=data;"]
,
				Print["Error: Dataset not found"];
			]
		]

GetDatasetNetCDF[f_?StringQ,set_?StringQ] := 
		Module[{},
				Import[f,{"Datasets",set}]
		]

OpenNetCDFExport[]:=Module[{listOfNumberObjects,arrayObjects,listOfObjects},
listOfObjects = Names["Global`*"];
listOfObjects=Select[listOfObjects,Length[Dimensions[Symbol[#]]]>0&];
listOfNumberObjects=Select[listOfObjects,DataArrayQ[Symbol[#]]&];
listOfNumberObjects=Select[listOfNumberObjects,StringFreeQ[#,"$"]&];
listOfNumberObjects=Select[listOfNumberObjects,!MemberQ[{"$IList","$ChkBox","listOfObjects","$ListOfUnits"},#]&];
arrayObjects=Map[#->ArrayInformation[Symbol[#]]&,listOfNumberObjects];
$IList=Map[GetChkBoxEntry,arrayObjects];
$ChkBox=Drop[$ChkBox,Complement[Range[Length[$ChkBox]],$IList]];
$IList=Map[GetChkBoxEntry,arrayObjects];
Grid[Join[{{Item["\[DownArrow]",ItemSize-> {Full,4}],
Item[Button["Export NetCDF",ExportAsNetCDF[SystemDialogInput["FileSave","data.nc"],$DataDescription],Alignment->{Center,Top},Appearance->"Palette",ImageSize->{100,15}, BaseStyle->{"Palatino",9},Method->"Queued"],Alignment->{Right,Top}],SpanFromLeft,
Item[InputField[Dynamic[$DataDescription],String,ImageSize->{500,12},ContentPadding->False,Enabled->True],Alignment->{Left,Top}],SpanFromLeft}
},
MapIndexed[(
{Checkbox[Dynamic[$ChkBox[[$IList[[#2[[1]]]],2]]]],
"\""<>#[[1]]<>"\"",
Item[("Type"/.#[[2]])<>ToShapeString[("Shape"/.#[[2]])],Alignment->Left],
Item[InputField[Dynamic[$ChkBox[[$IList[[#2[[1]]]],3]]],String,ImageSize->{100,12},ContentPadding->False]],
Item[InputField[Dynamic[$ChkBox[[$IList[[#2[[1]]]],4]]],String,ImageSize->{200,12},ContentPadding->False]],
Item[InputField[Dynamic[$ChkBox[[$IList[[#2[[1]]]],5]]],String,ImageSize->{70,12},ContentPadding->False]],Item[PopupMenu[Dynamic[$ChkBox[[$IList[[#2[[1]]]],5]]],$ListOfUnits,Style["[other]",FontSize->9],ContentPadding->False,ImageSize->{90,14},Alignment->{Center,Center},BaseStyle->{FontFamily->"Gill Sans",FontSize->10},Appearance->"Palette"],Alignment->{Right,Center}],
Item[ToString[("Memory"/.#[[2]])]<>" Bytes",Alignment->Right],Button["Remove",Remove[Evaluate[#[[1]]]],Alignment->{Center,Top},Appearance->"Palette",ImageSize->{50,14}, BaseStyle->{"Palatino",9},Method->"Queued",ContentPadding->False]})&,arrayObjects]],BaseStyle->{FontFamily->"Gill Sans",FontSize->10},Spacings->{0.5,0.2},Background->{None,{{LightGray,None}}},Frame->{None,All},FrameStyle->Thin,ItemSize->{Full,2},Alignment->{Center,Center}]
]

NetCDFFile/:View[NetCDFFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[NetCDFFile]]]:=Module[{},
Grid[Join[{
{Item[Style["\[FilledSmallSquare] NetCDF ["<>f<>"]",FontWeight->Bold],ItemSize->{Full,2},Alignment->Left],SpanFromLeft},
{Item[Style["title"/.m,FontSlant->Italic],ItemSize->{Full,2},Alignment->Left],SpanFromLeft}
},
MapIndexed[(
{
Style[#[[1]],FontWeight->Bold],
Item[(#[[3]])<>ToShapeString[#[[2]]],Alignment->Left],
Item["description"/.#[[4]],Alignment->Left],
Item["units"/.#[[4]]],
Button["Restore [.]",
NotebookWrite[SelectedNotebook[],Cell[BoxData["\[Placeholder]=GetDatasetNetCDF[\""<>f<>"\",\""<>#[[1]]<>"\"];"],"Input"]]
,Alignment->{Center,Top},Appearance->"Palette",ImageSize->{50,16}, BaseStyle->{"Palatino",9},Method->"Queued",ContentPadding->False],
If[(!NameQ[("mathematica_symbol"/.#[[4]])/.{"mathematica_symbol"->""}])&&StringLength[("mathematica_symbol"/.#[[4]])/.{"mathematica_symbol"->""}]>0,
Button["Restore ["<>(("mathematica_symbol"/.#[[4]])/.{"mathematica_symbol"->""})<>"]",
RestoreDataset[NetCDFFile[f,m,d,opts],#[[1]],(("mathematica_symbol"/.#[[4]])/.{"mathematica_symbol"->""})]
,Alignment->{Center,Top},Appearance->"Palette",ImageSize->{150,16}, BaseStyle->{"Palatino",9},Method->"Queued",ContentPadding->False],""]
})&,d]],BaseStyle->{FontFamily->"Gill Sans",FontSize->10},Spacings->{0.5,0.2},Background->{None,{{LightGray,None}}},Frame->{None,All},FrameStyle->Thin,ItemSize->{Full,2.4},Alignment->{Center,Center}]
]
$ChkBox={};
$ListOfUnits=Map[#->Style[#,FontFamily->"Gill Sans",FontSize->9]&,{"[none]","K","m","J","N","kCal Mol^{-1}"}];
$DataDescription="[description]";


End[ ];


EndPackage[ ];
