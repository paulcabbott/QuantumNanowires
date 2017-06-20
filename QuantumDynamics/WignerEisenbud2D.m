(* ::Package:: *)

(* ::Title:: *)
(*2D Wigner-Eisenbud Functions*)


(* ::Section:: *)
(*Package Information*)


(* ::Text:: *)
(*This package computes 2D Wigner-Eisenbud Functions (WEFs) in cylindrical coordinates using the Fourier-Bessel functions as the radial basis.*)


(* ::Section:: *)
(*Package Beginnings*)


(* ::Text:: *)
(*Begin package and load associated pacakge:*)


BeginPackage["QuantumDynamics`WignerEisenbud2D`",{"QuantumDynamics`WignerEisenbud`","QuantumDynamics`TransportCoefficients`"}];


(* ::Text:: *)
(*Clear pre-existing definitions.*)


Unprotect["QuantumDynamics`WignerEisenbud2D`*"];


ClearAll["QuantumDynamics`WignerEisenbud2D`*","QuantumDynamics`WignerEisenbud2D`Private`*"];


(* ::Subsection:: *)
(*Usage statements*)


JZeros::usage="JZeros[m] computes and saves the first 500 zeros of BesselJ[m,z]";


FourierBesselFunction::usage="FourierBesselFunction[m,k,r] is the k-th Fourier-Bessel function of r for quantum number m.";


TransversalEnergy::usage="TransversalEnergy[m,k,s] is the transversal eigenenergy of the k-th state for quantum number m and scale factor s.";


TotalEnergy::usage="TotalEnergy[{l,m,s},v] gives the total energy of the system.";


OrderedPairs::usage="OrderedPairs[list1,list2,n,f] returns the first n indexed elements of the sorted list obtained by applying f to pairs of elements taken from list1 and list2.";


kp::usage="kp[{m,s},v] gives the lowest-order Wigner-Eisenbud energies with associated k and p indices for quantum number m, scale factor s, and sampled potential v."; 


WignerEisenbudPlot::usage="WignerEisenbudPlot[v,n,scale,opts] plots the first n scaled Wigner-Eisenbud functions for sampled potential v, displaced by the Wigner-Eisenbud energies.";


WignerEisenbud2DPlot::usage="WignerEisenbud2DPlot[{l,m,s},v] plots the 2D Wigner-Eisenbud functions for quantum number m, potential v, index l, and scale factor s.";


WignerEisenbud2DGrid::usage="WignerEisenbud2DGrid[{lMax,m,s},v,n] plots a grid of plots of 2D Wigner-Eisenbud functions, with n per row, for l \[LessEqual] lMax.";


WignerEisenbud2DContourPlot::usage="WignerEisenbudS2DContourPlot[{l,m,s},v] produces a contour plot of the 2D Wigner-Eisenbud functions for quantum number m, potential v, index l, and scale factor s.";


(* ::Subsection:: *)
(*Begin Private context*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Kronecker product*)


CircleTimes=KroneckerProduct;


(* ::Section:: *)
(*Fourier-Bessel basis*)


JZeros/:JZeros[m_]:=JZeros/:JZeros[m]=N[BesselJZero[m,Range[500]]]


FourierBesselFunction[m_,k_,r_]:= (Sqrt[2] BesselJ[m,JZeros[m][[k]] r])/BesselJ[m+1,JZeros[m][[k]]]


(* ::Section:: *)
(*Energies Relation*)


TransversalEnergy[m_,k_,s_]:=s^2 JZeros[m][[k]]^2


OrderedPairs[list1_List,list2_List,n_Integer?Positive,f_:Plus]:=
	Take[Sort[Flatten[MapIndexed[Flatten[{##}]&,Outer[f,Sort[list1],Sort[list2]],{2}],1]],n]


kp[{m_,s_},v_List]:=kp[{m,s},v]=
	OrderedPairs[s^2 JZeros[m]^2, WignerEisenbudEnergies[v], 3 Length[JZeros[m]]]


TotalEnergy[{l_,m_,s_},v_List]:=kp[{m,s},v][[l,1]]


(* ::Section:: *)
(*Plotting functions*)


WignerEisenbudPlot[v_List,n_:3,scale_:1,opts___]:=
	DataPlot[scale WignerEisenbudFunctions[v][[1;;n]]^2+WignerEisenbudEnergies[v][[1;;n]],{-1,1},
		Filling->Thread[Range[n]->WignerEisenbudEnergies[v][[1;;n]]]
]


WignerEisenbud2DPlot[{l_,m_,s_:1},v_List]:=
	DataPlot3D[
		FourierBesselFunction[m,kp[{m,s},v][[l,2]],Range[0,1,1/20]] \[CircleTimes] WignerEisenbudFunctions[v][[kp[{m,s},v][[l,3]]]],
		Mesh->False,
		PlotRange->All,
		PlotLabel->{l,TotalEnergy[{l,m,s},v]}
	]


WignerEisenbud2DGrid[{lMax_,m_,s_:1},v_List,n_:3]:=
	GraphicsGrid[
		Partition[
			Table[WignerEisenbud2DPlot[{l,m,s},v],{l,1,lMax}],
		n]
	]


WignerEisenbud2DContourPlot[{l_,m_,s_:1},v_List]:= Module[{e=1/Length[v]},
	ListContourPlot[
		FourierBesselFunction[m,kp[{m,s},v][[l,2]],Range[0,1,1/20]] \[CircleTimes] WignerEisenbudFunctions[v][[kp[{m,s},v][[l,3]]]],
		DataRange->{{-1+e,1-e},{0,1}},
		FrameLabel->{Style["z",Italic],Style["r",Italic]}
	]
]


(* ::Section:: *)
(*Package Ending*)


End[]; 


(* ::Subsection:: *)
(*End package*)


EndPackage[];
