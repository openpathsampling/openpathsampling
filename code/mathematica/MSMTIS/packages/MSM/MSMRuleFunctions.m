(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ReplaceRules::usage = "ReplaceRules[rule, listofrules] replaces rule in listofrules.";
UnionRules::usage = "UnionRules[listofrules] flattens the list of rules of only the first appearance of each rules left side.";
RemoveRules::usage = "RemoveRules[list, listofrules] removes the leftsiderules in list from listofrules.";
KeepRules::usage = "KeepRules[list, listofrules] removes all but the rules listed in list.";
KeepOptions::usage = "KeepOptions[function, listofrules] keeps only rules in listofrules that are allowed for function";

Begin["`Private`"] (* Begin Private Context *) 

ReplaceRules[rule_,f___] :=
    Apply[Sequence,Join[{rule},Select[List[f],#[[1]]=!=rule[[1]]&]]]
UnionRules[rules1___] :=
    Apply[Sequence,Union[List[rules1][[All,1]]]/.Map[#[[1]]->#&,List[rules1]]]
RemoveRules[rules_,f___] :=
    Module[ {ruleList, result},
        ruleList = If[ ListQ[rules],
                       rules,
                       {rules}
                   ];
        result = List[f];
        Table[ result = Select[result,#[[1]]=!=rule&];, {rule, ruleList}];
        Apply[Sequence, result]
    ];
KeepRules[rules_,f___] :=
    Module[ {ruleList, result, r},
        ruleList = If[ ListQ[rules],
                       rules,
                       {rules}
                   ];
        result = {};
        Table[ 
            r = Select[List[f],#[[1]]===rule&];
            If[ Length[r]>0,
                AppendTo[result, r[[1]]];
            ];
        , {rule, ruleList}];
        Apply[Sequence, result]
    ]
KeepOptions[function_, f___] :=
 	Apply[Sequence, 
  		Select[List[f], MemberQ[Options[function][[All, 1]], #[[1]]] &]
  	]

End[] (* End Private Context *)

EndPackage[]