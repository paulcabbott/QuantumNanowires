(* ::Package:: *)

(* ::Title:: *)
(*Wigner-Eisenbud Functions*)


(* ::Section:: *)
(*Package Information*)


(* ::Text:: *)
(*This package computes the Wigner-Eisenbud Functions (WEFs) for arbitrary potentials in 1D, and 2D in cylindrical coordinates using the Fourier discrete cosine transform (DCT).*)


(* ::Section:: *)
(*Package Beginnings*)


BeginPackage["QuantumDynamics`WignerEisenbud`"];


(* ::Text:: *)
(*Clear pre-existing definitions.*)


Unprotect["QuantumDynamics`WignerEisenbud`*"];


ClearAll["QuantumDynamics`WignerEisenbud`*","QuantumDynamics`WignerEisenbud`Private`*"];


(* ::Subsection:: *)
(*Usage statements*)


SampledFunction::usage="SampledFunction[f,n,{a,b}] samples the 1D function f uniformly n times over the interval [a,b]. SampledFunction[f,{m,n}] samples the 2D function f[r,z] uniformly, m\[Times]n times over [-1,1]\[Times][0,1]";


DataPlot::usage="DataPlot[f,{a,b}] plots a uniformly sampled 1D function f over [a,b]";


DataPlot3D::usage="DataPlot3D[f] plots a uniformly sampled 2D function f over [-1,1]\[Times][0,1]";


HamiltonianMatrix::usage="HamiltonianMatrix[v] gives the Hamiltonian matrix for a sampled potential v"; 


HamiltonianMatrix::usage="HamiltonianMatrix[m,v,s] gives the Hamiltonian matrix for magnetic quantum number m, sampled potential v, and parameter s = R/d."; 


WignerEisenbudEnergies::usage="WignerEisenbudEnergies[m,v,s] returns the energies for sampled potential v, and parameter s = R/d."; 


WignerEisenbudSystem::usage="WignerEisenbudSystem[m,v,s] gives the set of Wigner-Eisenbud energies and eigenvectors for sampled potential v, and parameter s = R/d.";


EnergyOrdered::usage="EnergyOrdered[es] orders the eigensystem in increasing eigenvalue.";


AxialBasis::usage="AxialBasis[n] uses the Fourier DCT to generate the sampled axial basis.";


RadialBasis::usage="RadialBasis[n] uses the Fourier DCT to generate the sampled radial basis.";


AxialMatrix::usage="AxialMatrix[f] computes the matrix elements for a sampled function f of z over [-1,1].";


RadialMatrix::usage="RadialMatrix[f] computes the matrix elements for a sampled function f of r over [0,1].";


AxialKineticMatrix::usage="AxialKineticMatrix[n] gives the n\[Times]n axial Kinetic energy matrix";


RadialKineticMatrix::usage="RadialKineticMatrix[m,n] gives the n\[Times]n radial Kinetic energy matrix for magnetic quantum number m.";


AxialOverlapMatrix::usage="AxialOverlapMatrix[n] computes the n\[Times]n axial overlap matrix.";


RadialOverlapMatrix::usage="RadialOverlapMatrix[n] Compute the n\[Times]n radial overlap matrix.";


MatrixElements::usage="MatrixElements[f] gives the matrix elements for a sampled function f in 1D over [-1,1] and 2D over over [-1,1]\[Times][0,1].";


WignerEisenbudFunctions::usage="WignerEisenbudFunctions[v] returns the Wigner-Eisenbud functions for sampled potential v. WignerEisenbudFunctions[m,v,l,s] returns the l-th set of 2D Wigner-Eisenbud functions for magnetic quantum number m, sampled potential v, and parameter s = d/R."; 


(* ::Subsection:: *)
(*Begin Private context*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Kronecker product*)


CircleTimes=KroneckerProduct;


(* ::Subsection:: *)
(*Sampled function*)


SampledFunction[f_,n_,{a_,b_}]:= With[{d=1/n}, f/@(a+(b-a)Range[d/2,1-d/2,d])]


SampledFunction[f_,{m_,n_}]:= Module[{d=1/n,e=1/m}, Table[f[r,z],{z,-1+e,1-e,2e},{r,d/2,1-d/2,d}]]


(* ::Subsection:: *)
(*Bases*)


AxialBasis[n_]:=AxialBasis[n]=(z\[Function]FourierDCT[z,3])/@IdentityMatrix[n]


RadialBasis[n_]:=RadialBasis[n]=(z\[Function]FourierDCT[z,4])/@IdentityMatrix[n]


(* ::Section:: *)
(*Matrices*)


(* ::Subsection:: *)
(*Kinetic matrix*)


AxialKineticMatrix[n_]:=AxialKineticMatrix[n]=\[Pi]^2/4.0 DiagonalMatrix[Range[0,n-1]^2]


(* ::Subsubsection:: *)
(*m=0*)


RadialKineticMatrix[n_]:=RadialKineticMatrix[n]=
SparseArray[{
	{k_,k_}->1/12 \[Pi]^2 (1-2 k)^2,
	{k_,l_}/;k!=l->((2 k-1) (2 l-1) (-1)^(k+l) ((k-1) k+3 (l-1) l+1))/(2 (k-l)^2 (k+l-1)^2)},
	{n,n}
]


(* ::Subsubsection:: *)
(*General m*)


RadialKineticMatrix[m_,n_]:=RadialKineticMatrix[m,n]=RadialKineticMatrix[n]+m^2 IdentityMatrix[n]


(* ::Subsection:: *)
(*Axial and radial matrix elements*)


AxialMatrix[f_List]:=FourierDCT/@ Transpose[Transpose[AxialBasis[Length[f]]] f]


RadialMatrix[f_List]:=(z\[Function]FourierDCT[z,4])/@ Transpose[Transpose[RadialBasis[Length[f]]] f]


(* ::Subsection:: *)
(*Axial and radial overlap matrices*)


AxialOverlapMatrix[n_]:=AxialOverlapMatrix[n]=IdentityMatrix[n]


RadialOverlapMatrix[n_]:=RadialOverlapMatrix[n]=RadialMatrix[SampledFunction[r\[Function]r^2,n,{0,1}]]


(* ::Subsection:: *)
(*2D matrix elements*)


MatrixElements[f_?MatrixQ]:=Chop[ArrayFlatten[Map[AxialMatrix, Transpose[RadialMatrix/@f,{3,2,1}],{2}]]]


MatrixElements[f_?MatrixQ,{p_,k_}]:=
	Chop[ArrayFlatten[Map[AxialMatrix[#][[1;;p,1;;p]]&,Transpose[RadialMatrix/@f,{3,2,1}][[1;;k,1;;k]],{2}]]]


(* ::Subsection:: *)
(*Hamiltonian matrix*)


HamiltonianMatrix[v_List]:=AxialKineticMatrix[Length[v]]+AxialMatrix[v]


HamiltonianMatrix[m_,v_?MatrixQ,{p_,k_},s_:1]:=
	s^2 RadialKineticMatrix[m,k]\[CircleTimes]AxialOverlapMatrix[p]+
	RadialOverlapMatrix[k]\[CircleTimes]AxialKineticMatrix[p]+
	MatrixElements[v SampledFunction[{r,z}\[Function]r^2,Dimensions[v]],{p,k}]


(* ::Section:: *)
(*Energies and Eigenfunctions*)


(* ::Subsection:: *)
(*Energy-ordered eigensystems*)


EnergyOrdered[es_]:= Transpose[Sort[Transpose[Chop[es]]]]


(* ::Subsection:: *)
(*Wigner-Eisenbud eigensystem*)


(* ::Subsubsection:: *)
(*1D*)


WignerEisenbudSystem[v_List]:= WignerEisenbudSystem[v] = EnergyOrdered[Eigensystem[N[HamiltonianMatrix[v]]]]


(* ::Subsubsection:: *)
(*2D*)


WignerEisenbudSystem[m_,v_?MatrixQ,{p_,k_},s_:1]:= WignerEisenbudSystem[m,v,{p,k},s]=
	EnergyOrdered[
		Eigensystem[{HamiltonianMatrix[m,v,{p,k},s], RadialOverlapMatrix[k]\[CircleTimes]AxialOverlapMatrix[p]}]
	]


(* ::Subsection:: *)
(*Wigner-Eisenbud energies*)


(* ::Subsubsection:: *)
(*1D*)


WignerEisenbudEnergies[v_List]:=First[WignerEisenbudSystem[v]]


(* ::Subsubsection:: *)
(*2D*)


WignerEisenbudEnergies[m_,v_?MatrixQ,{p_,k_},s_:1]:=First[WignerEisenbudSystem[m,v,{p,k},s]]


(* ::Subsection:: *)
(*Wigner-Eisenbud functions*)


(* ::Subsubsection:: *)
(*1D*)


WignerEisenbudFunctions[v_List]:= 
Module[
	{n=Length[v], evecs=Last[WignerEisenbudSystem[v]]},
	evecs[[All,1]]=2 evecs[[All,1]];
	Sqrt[n/2] Normalize/@(evecs.AxialBasis[n])
]


(* ::Subsubsection:: *)
(*2D*)


\[LeftAngleBracket]f_,g_\[RightAngleBracket]:=Total[SampledFunction[r\[Function]r,Length[f],{0,1}].(f g)] 2/(Times@@Dimensions[f]) /;
	Dimensions[f]==Dimensions[g]


WignerEisenbudFunctions[m_,v_,{p_,k_},l_,s_:1]:=
Module[
	{evecs=Partition[WignerEisenbudSystem[m,v,{p,k},s][[2,l]],p]},
	Normalize[RadialBasis[k].evecs.AxialBasis[p],f\[Function]Sqrt[\[LeftAngleBracket]f,f\[RightAngleBracket]]]
]


(* ::Section:: *)
(*Plotting functions*)


DataPlot[f_List,{a_,b_},opts___]:= Module[{e=(b-a)/(2 Last[Dimensions[f]])},
	ListPlot[f,
		opts,
		DataRange->{a+e,b-e},
		PlotRange->{{a,b},Automatic},
		Joined->True
	]
]


DataPlot3D[f_?MatrixQ, opts___]:=Module[{d,e},
	{d,e}=1/Dimensions[f];
	ListPlot3D[f,
		opts,
		DataRange->{{-1+e,1-e},{d/2,1-d/2}},
		Mesh->Full,
		PlotRange->{{-1,1},{0,1},Automatic},
		AxesLabel->{Style["z",Italic],Style["r",Italic]}
	]
]


(* ::Section:: *)
(*Package Ending*)


End[]; 


(* ::Subsection:: *)
(*Setting Attributes*)


(*SetAttributes[{...}, {Listable, Protected, ReadProtected, NumericFunction}]; *)


(* ::Subsection:: *)
(*End package*)


EndPackage[];
