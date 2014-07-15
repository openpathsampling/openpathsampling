(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ShowDataStructure::usage = "ShowDataStructure[object] returns the data structure of object.";
ReduceDataStructure::usage = "ReduceDataStructure[data] reduces the datastructure recursively by aggregating equal objects into arrays.";

Begin["`Private`"] (* Begin Private Context *) 

StringDataFormat[str_] :=
    Module[ {
    locStr,
    style
    },
        locStr = str;
        style = {FontFamily -> "GillSans"};
        If[ StringFreeQ[locStr, "Private`"],
            style = Join[style, {}],
            style = Join[style, {Gray}]
        ];
        locStr = StringReplace[locStr, {"Private`" -> "", "MSM`" -> ""}];
        Style[locStr, Sequence@style]
    ]

ShowDataStructure[data_] :=
    Module[ {},
        Which[
         ArrayQ[data, _, NumericQ],
         StringDataFormat@("Array[" <> DimensionsToString@Dimensions[data] <>
             "]"),
         Head[data] === List,
         Map[ShowDataStructure, data],
         Head[data] === Rule,
         Rule[StringDataFormat[ToString[data[[1]]]], 
          ShowDataStructure[data[[2]]]],
         Head[data] === Symbol,
         Which[data == True || data == False,
          Boole,
          True,
          Symbol],
         True,
         StringDataFormat[ToString@Head[data]]
         ]
    ]

ReduceDataStructure[data_] :=
    Module[ {reduced, level, back},
        Which[
         ArrayQ[data],
         level = ArrayDepth[data];
         reduced = False;
         While[(level >= 0) && (reduced == False),
          If[ Apply[Equal, Flatten[data, level]] === True,
              reduced = True;
              back = 
               Superscript[
                Row[{"{", ReduceDataStructure[First[Flatten[data, level]]], 
                  "}"}], DimensionsToString[Dimensions[data]]];
              level = 0;
          ];
          level--;
          ];
         If[ reduced == True,
             back,
             Map[ReduceDataStructure, data]
         ],
         Head[data] === Rule,
         Rule[data[[1]], ReduceDataStructure[data[[2]]]],
         ListQ[data],
         Map[ReduceDataStructure, data],
         True,
         data
         ]
    ]

DimensionsToString[dim_] :=
    Module[ {},
        StringJoin@Drop[Flatten@Map[{ToString[#], "\[Times]"} &, dim], -1]
    ]
  
  End[] (* End Private Context *)

EndPackage[]