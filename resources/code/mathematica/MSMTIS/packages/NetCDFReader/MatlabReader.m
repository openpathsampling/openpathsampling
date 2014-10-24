(* ::Package:: *)

(* TransitionMatrix Package *)
(* Version 0.3 *)

(* Last edited on Apr 6th, 2009 *)

BeginPackage["NetCDFReader`"];


(* Messages *)

OpenMatLabFile::usage = "Open a Matlab File (.mat) and return an MatlabFile-Object, that contains information about all containing datasets";
MatlabFile::usage = "";
GetDataset::usage = "";
ListDatasets::usage = "";
ListDatasetsFull::usage = "";
WriteMatlabFile::usage = "";


(* All options *)

Options[OpenMatlabFile]                  = {};
Options[MatlabFile]                  = {};


Begin["`Private`"];


(* Functions to handle options *)

ReplaceRules[rule_,f___] := 
		Apply[Sequence,Join[{rule},Select[List[f],#[[1]]=!=rule[[1]]&]]]
UnionRules[rules1___] := 
		Apply[Sequence,Union[List[rules1][[All,1]]]/.Map[#[[1]]->#&,List[rules1]]]
RemoveRules[rule_,f___] := 
		Apply[Sequence,Select[List[f],#[[1]]=!=rule&]]


(* MatlabFile functions *)

MatlabFile/:MakeBoxes[MatlabFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[MatlabFile]],StandardForm] := 
		InterpretationBox[#,MatlabFile[f,m,d,opts]]&[RowBox[{SubscriptBox["Matlab",RowBox[{Length[d]}]],"(",StringCases[f,RegularExpression["[\\w._+-]+$"]][[1]],")"}]]


(* Access Matlab Files *)

OpenMatlabFile[filename_?StringQ]:=
		Module[{metadata,datasets,annotations,formats,dimensions},
			datasets = Import[filename,"Matlab"];
			annotations = {};
			metadata = {"Comments"->Import[filename,"Comments"]};
			dimensions = {};
			formats = {};
			MatlabFile[filename,metadata,Transpose[{datasets,dimensions,formats,annotations}]]
		]

WriteMatlabFile[filename_?StringQ, datasets_]:=
		Module[{metadata,annotations,formats,dimensions},
			Export[filename,datasets[[All,2]],{"Datasets",datasets[[All,1]]}]
		]

MatlabFile/:GetDataset[MatlabFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[MatlabFile]],set_?StringQ] := 
		Module[{},
			If[Apply[Or,StringMatchQ[d[[All,1]],set]],
				Import[f,{"Datasets",set}],
				Print["Error: Dataset not found"];
			]
		]

MatlabFile/:GetDataset[MatlabFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[MatlabFile]],index_?IntegerQ] := 
		Module[{},
			If[index>0&&index<=Length[d],
				Import[f,{"Datasets",index}],
				Print["Error: Dataset not found"];
			]
		]

MatlabFile/:ListDatasets[MatlabFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[MatlabFile]]] := 
		Module[{},
			d[[All,1]]
		]

MatlabFile/:ListDatasetsFull[MatlabFile[f_?StringQ,m_?ListQ,d_?ListQ,opts:OptionsPattern[MatlabFile]]] := 
		Module[{},
			d
		]


End[ ];


EndPackage[ ];
