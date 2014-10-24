(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

ReadSparseMatrix::usage = "ReadSparseMatrix[filename] imports filename and assumes that each row contains idX idY value and constructs a sparse matrix object from it.";
WriteSparseMatrix::usage = "WriteSparseMatrix[filename, sparsematrix] exports a sparematrix to filename and writes each non-zero element as idX idY value.";

Begin["`Private`"] (* Begin Private Context *) 

ReadSparseMatrix[fileName_,opts:OptionsPattern[Import]] :=
    Module[ {readImport},
        readImport = Import[fileName,"Table",opts];
        SparseArray[Map[1+#[[{1,2}]]->#[[3]]&,readImport]]
    ]
WriteSparseMatrix[fileName_,m_?MatrixQ,optsE:OptionsPattern[Export]] :=
    Module[ {arrayToExport},
        arrayToExport = Map[{#[[1,1]]-1,#[[1,2]]-1,N[#[[2]]]}&,Drop[ArrayRules[SparseArray[m]],-1]];
        Export[fileName,arrayToExport,"Table",optsE]
    ]


End[] (* End Private Context *)

EndPackage[]