(* Mathematica Package *)

BeginPackage["MSM`"]
(* Exported symbols added here with SymbolName::usage *)  

RowSumOneQ::usage = "RowSumOneQ[matrix] returns true if the row sum of all rows of matrix equals 1, false otherwise.";
RowSumZeroQ::usage = "RowSumZeroQ[matrix] returns true if the row sum of all rows of matrix equals 0, false otherwise.";
TotalSumOneQ::usage = "TotalSumOneQ[matrix] returns true if the total sum of all entries of matrix is one, false otherwise.";

SymmetricQ::usage = "SymmetricQ[matrix] return true if the matrix is symmetric, false otherwise.";

AbsMax::usage = "AbsMax[list] returns the largest in magnitude value in list including its sign.";
OrientVector::usage = "OrientVector[vector] returns vector oriented(flipped) so that the value of maximal magnitude is positive.";
Stochastize::usage = "Stochastize[vector] returns a stochastized (row sum 1) version of vector.";

ToImpliedTimescales::usage = "ToImpliedTimescales[list, lagtime] converts eigenvalues in list to implied timescales at time steps lagtime.";
Unphysical::usage = "Unphysical is a symbol that represents unphysical timescales, i.e. timescales that correspond to eigenvalues larger than 1 or smaller than 0.";

AutoCorrelation::usage = "AutoCorrelation[list, range] return a vector of autocorrelations of list for a set of time lags specified in range.";
AutoCovariance::usage = "AutoCovariance[matrix, range] return a list of correlation matrices computed from matrix that contains a list of timeseries for a set of time lags specified in range.";

ToLeftEigenvectors::usage = "ToLeftEigenvectors[matrix, n] Computes the n left eigenvectors of matrix normalized to a symmetry in the stationary vector.";
ToRightEigenvectors::usage = "ToRightEigenvectors[matrix, n] Computes the n right eigenvectors of matrix normalized to a symmetry in the stationary vector.";
ToFullEigensystem::usage = "ToFullEigensystem[matrix, n] Returns a vector of the first n eigenvalues, the first n left eigenvectors and the first right eigenvectors of matrix normalized to a symmetry in the stationary vector.";
FromFullEigensystem::usage = "FromFullEigensystem[{evalues, leftev, rightev}] reconstruct a matrix from its eigenvalues evalues and left and right eigenvectors given in leftev and rightev.";

EquilizeLeftEVSet::usage = "EquilizeLeftEVSet[matrix] normalizes a set of left eigenvectors given in matrix assuming symmetry in the stationary distribution (the first left eigenvectors)";
EquilizeRightEVSet::usage = "EquilizeRightEVSet[matrix, equilibrium] normalizes a set of right eigenvectors given in matrix assuming symmetry in the stationary distribution equilibrium";

EVSet::usage = "EVSet is an option that can be used in FromFullEigensystem to specifiy the set of eigenvectors to be used in the reconstruction.";
EquilizeEVSet::usage = "EquilizeEVSet[{evalues, leftev, rightev}] return the normalized versions of the specified full eigensystem according to the stationary distribution in the first left eigenvector.";

MatrixMaxNorm::usage = "MatrixMaxNorm[matrix] returns the maximum norm of matrix.";
FrobeniusNorm::usage = "FrobeniusNorm[matrix] returns the Frobenius norm of matrix.";
SpectralNorm::usage = "SpectralNorm[matrix] returns the spectral norm of matrix.";
MatrixEqual::usage = "MatrixEqual[matrix1, matrix2] returns true of matrix1 equals matrix2 numerically, false otherwise.";

V::usage = "V[string, states, ...] constructs a mathematica variable that starts with string and adds numbers given by states.";
Idx::usage = "Idx[list, criterion] returns a list of indices where criterion evaluates to true in list. If criterion is not a function this is similar to Position[list, criterion].";

GetOptions::usage = "GetOptions[object] returns a list of options of object. This should become obsolete. Try not to use it.";

MatrixRangeSum::usage = "MatrixRangeSum[matrix, range] returns the sum of the powers of matrix in range.";
MatrixRangeMeanSum::usage = "MatrixRangeMeanSum[matrix, range] returns the sum of the powers of matrix in range each multiplied by the power.";

PadInteger::usage = "PadInteger[integer, padding] returns a string version of integer with zeros padded on the left to a total length of padding.";
ExpandList::usage = "ExpandList[list] threads over all elements in list.";

SteadyEigenvector::usage = "SteadyEigenvector[matrix] returns the normalized (total sum one) eigenvector whose eigenvalue closest matched one.";

MeanRange::usage = "MeanRange[range, dimension, bins] returns a list of bin centers of the rangeobject for the specified dimension.";
MeanRangeBin::usage = "MeanRange[range, dimension, bins] returns a list of bins of the rangeobject for the specified dimension.";

Begin["`Private`"] (* Begin Private Context *) 

MeanRange[range_,dim_,steps_] := 
	Module[{dStep},
		dStep=(range[[dim,2]]-range[[dim,1]])/(steps);
		Range[range[[dim,1]]+dStep/2,range[[dim,2]]-dStep/2,dStep]
	]

MeanRangeBin[range_,dim_,steps_] := 
	Module[{dStep, bins},
		bins = {};
		If[NumberQ[steps],
			dStep=(range[[dim,2]]-range[[dim,1]])/(steps);
			bins = Transpose[{Range[range[[dim,1]],range[[dim,2]]-dStep/2,dStep], Range[range[[dim,1]]+dStep/2,range[[dim,2]],dStep]}]
		];
		If[VectorQ[steps],
			bins = Partition[Flatten[Join[range[[dim,1]],Table[{1,1}*(steps[[nn]] + steps[[nn + 1]])/2,{nn,1,Length[steps] - 1}],range[[dim,2]]]],2];
		];
		bins
	]

V[h_, ii_] :=
    ToExpression[ToString[h] <> ToString[ii]]
V[h_, ii_, jj_] :=
    ToExpression[ToString[h] <> ToString[ii] <> ToString[jj]]

V[var_, st1_, st2_, size_] :=
    Module[ {},
        ToExpression[
         ToString[var] <> 
          PadInteger[st1, size] <> 
          PadInteger[st2, size]
        ]
    ]
    
V[var_, st_?ListQ, opts___] :=
    Module[ {size},
    	size = Padding/.{opts}/.{Padding->2};
        ToExpression[
         ToString[var] <> StringJoin[Map[PadInteger[#, size]&,st]]
        ]
    ]

ExpandList[list_] :=
    Map[Thread, list]

PadInteger[int_, digits_] :=
    StringTake["00000000" <> ToString[int], -digits]

RowSumOneQ[m_?MatrixQ] :=
    Apply[And,Map[#==1&,Map[Apply[Plus,#]&,m]]]
RowSumZeroQ[m_?MatrixQ] :=
    Apply[And,Map[#==0&,Map[Apply[Plus,#]&,m]]]
TotalSumOneQ[m_?MatrixQ] :=
    Total[Flatten[m]]==1
SymmetricQ[m_?MatrixQ] :=
    m==Transpose[m];
    
    
AbsMax[v_?VectorQ] :=
    Module[ {a1,a2},
        a1 = Max[v,0];
        a2 = -Max[-v,0];
        If[ a1>-a2,
            a1,
            a2
        ]
    ]
OrientVector[v_] :=
    Module[ {},
        Which[
            VectorQ[v], v/AbsMax[v],
            MatrixQ[v], Map[#/AbsMax[#]&, v],
            True, v
        ]
    ]
Stochastize[v_?VectorQ,opts___] :=
    Module[ {voidFnc},
        voidFnc = VoidFunction/.List[opts]/.{VoidFunction->(#*0&)}/.{Automatic->Function[x,x*0],Zero->(#*0&),Infinity->(#/0&)};
        If[ Apply[Plus,v]>0,
            v/Total[v],
            voidFnc[v]
        ]
    ];
ToImpliedTimescales[v_?VectorQ, lagtime_] :=
    Module[ {},
        Reverse[Sort[
            Map[Which[
                (1.-Chop[Abs[#]])<10^-8,Infinity,
                Im[#]!=0,Unphysical,
                Abs[#]<0,Unphysical,
                Abs[#]>1,Unphysical,
                #==0,0,
                True,-lagtime / Log[#]]&,
                v
            ]
        ,#1<#2&]]
    ];
AutoCorrelation[v_?VectorQ, range_] :=
    Module[ {len},
        len = If[ range>0,
                  range,
                  Length[v]+range+1
              ];
        ListCorrelate[v,v,{1,-range},0]/Range[Length[v],Length[v]-len+1,-1]
    ]
AutoCovariance[m_?MatrixQ, range_] :=
    Module[ {length, size, rangeLength,ergTab,singleResult},
        length = Dimensions[m][[1]];
        size = Dimensions[m][[2]];
        rangeLength = If[ range>0,
                          range,
                          length+range+1
                      ];
        ergTab = Table[{},{size},{size}];
        Table[
            singleResult = ListCorrelate[m[[All,ii]],m[[All,jj]],{1,-range},0]/Range[length,length-rangeLength+1,-1];
            Which[ii==jj,
                ergTab[[ii,jj]] = singleResult;,
            jj>ii,
                ergTab[[ii,jj]] = singleResult;
                ergTab[[jj,ii]] = singleResult;
            ];
        ,{ii,1,size}, {jj,1,size}];
        Transpose[ergTab,{3,2,1}]
    ]
ToLeftEigenvectors[m_?MatrixQ,n_] :=
    ToFullEigensystem[m,n][[2]];
ToRightEigenvectors[m_?MatrixQ,n_] :=
    ToFullEigensystem[m,n][[3]];
ToFullEigensystem[m_?MatrixQ,n_] :=
    Module[ {l,r,w,s,d,size},
        size = Min[(n/.{All->Length[m]}),Length[m]];
        {d,l} = Eigensystem[Transpose[m],size];
        l[[1]] = #/Total[#]&[l[[1]]];
        If[ SymmetricQ[m],
            r = l;,
            r = Eigenvectors[m];
        ];
        r[[1]] = #/Max[#]&[r[[1]]];
        w = Diagonal[Outer[Dot,l,r,1]];
        s = Sign[w];
        {
        	d,
        	EquilizeLeftEVSet[l],
        	EquilizeRightEVSet[r,l[[1]]]
        }
    ];
FromFullEigensystem[{d_?VectorQ,l_?MatrixQ,r_?MatrixQ},opts___] :=
    Module[ {evs},
        evs = (EVSet/.{opts})/.{EVSet->Range[Length[d]], All->Range[Length[d]]};
        Transpose[r[[evs]]].DiagonalMatrix[d[[evs]]].l[[evs]]
    ]
    
Idx[list_, criterion_] := Module[{},
	If[Head[criterion]===Function,
	  Select[Transpose[{Range[Length[list]], list}], criterion[#[[2]]] &][[All, 1]],
 	  Flatten[Position[list,criterion]]
	]
  ]
    
SteadyEigenvector[m_?MatrixQ] := Module[{eval, evec, eq},
	{eval,evec} = Eigensystem[m];
	eq = Total@evec[[Idx[eval,#==1.0&]]];
	eq/Total[eq]	
]
    
EquilizeLeftEVSet[evSet_] :=
    Module[ {firstEV, fac},
        firstEV = evSet[[1]]/Total[evSet[[1]]];
        Join[{firstEV}, Table[
          fac = evSet[[nn]].(evSet[[nn]]/firstEV);
          evSet[[nn]]/Sqrt[Abs[fac]]
          ,{nn, 2, Length[evSet]}]
        ]
    ]
  
EquilizeRightEVSet[evSet_,equilibrium_] :=
    Module[ {firstEV, fac},
        firstEV = equilibrium/Total[equilibrium];
        Table[
          fac = evSet[[nn]].(evSet[[nn]]*firstEV);
          evSet[[nn]]/Sqrt[Abs[fac]]
          ,{nn, 1, Length[evSet]}
        ]
    ]
  
EquilizeEVSet[{eVal_?VectorQ,leftEV_?MatrixQ,rightEV_?MatrixQ}] :=
    Module[ {eq},
        eq = leftEV[[1,1]]/Total[leftEV[[1]]];
        {eVal,
           EquilizeLeftEVSet[leftEV],
           EquilizeRightEVSet[rightEV, eq]
        }
    ]
  
MatrixMaxNorm[m_?MatrixQ] :=
    Length[m]*Max[Abs[Chop[m]]]
MatrixEqual[m1_?MatrixQ, m2_?MatrixQ] :=
    MatrixMaxNorm[m1-m2] == 0
FrobeniusNorm[m_?MatrixQ] :=
    Sqrt@Tr[Transpose[m].m]
SpectralNorm[m_?MatrixQ] :=
    Norm[m]

GetOptions[m_] :=
    Drop[Apply[List, m], 1]

MatrixRangeSum[m_?MatrixQ,range_] :=
    Module[ {size,laplace,inv},
        size = Length[m];
        laplace = IdentityMatrix[size]-m;
        inv = Inverse[laplace];
        inv.(If[ range[[1]]>0,
                 (MatrixPower[m,range[[1]]]),
                 (IdentityMatrix[size])
             ] - If[ range[[2]] =!= Infinity,
                     MatrixPower[m,range[[2]]+1],
                     0
                 ])
    ]
MatrixRangeMeanSum[m_?MatrixQ,range_] :=
    Module[ {size,inv,laplace},
        size = Length[m];
        laplace = IdentityMatrix[size]-m;
        inv = Inverse[laplace];
        inv.inv.(If[ range[[1]]>0,
                     (MatrixPower[m,range[[1]]].(IdentityMatrix[size]+range[[1]]*laplace)),
                     IdentityMatrix[size]
                 ] - If[ range[[2]]=!=Infinity,
                         MatrixPower[m,range[[2]]+1].(IdentityMatrix[size]+(range[[2]]+1)*laplace),
                         0
                     ])
    ]

End[] (* End Private Context *)

EndPackage[]