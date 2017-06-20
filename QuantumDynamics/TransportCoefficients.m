(* ::Package:: *)

(* ::Title:: *)
(*Transport Coefficients via 2D WEFs*)


(* ::Section:: *)
(*Package Beginnings*)


(* ::Text:: *)
(*Begin package and load associated pacakges:*)


BeginPackage["QuantumDynamics`TransportCoefficients`",
	{"QuantumDynamics`WignerEisenbud`","QuantumDynamics`WignerEisenbud2D`"}];


(* ::Text:: *)
(*Clear pre-existing definitions:*)


Unprotect["QuantumDynamics`TransportCoefficients`*"];


ClearAll["QuantumDynamics`TransportCoefficients`*","QuantumDynamics`TransportCoefficients`Private`*"];


(* ::Subsection:: *)
(*Usage statements*)


kMax::usage="kMax[{m,s},v] gives the maximum value of k for quantum number m, scale factor s, and sampled potential v.";


RMatrix::usage="RMatrix[{m,s},v] gives the R-matrix for quantum number m, scale factor s, and sampled potential v as a function of energy.";


StepMatrix::usage="StepMatrix[{m,s,n}] gives the n\[Times]n step matrix for quantum number m and scale factor s as a function of energy.";


KMatrix::usage="KMatrix[{m,s,n}] gives the n\[Times]n K-matrix for quantum number m and scale factor s as a function of energy.";


OmegaMatrix::usage="OmegaMatrix[{m,s},v] gives the \[CapitalOmega]-matrix for quantum number m, scale factor s, and sampled potential v as a function of energy.";


CurrentScatteringMatrix::usage="CurrentScatteringMatrix[{m,s},v] gives the current scattering matrix for quantum number m, scale factor s, and sampled potential v as a function of energy.";


TransmissionCoefficient::usage="TransmissionCoefficient[{m,s},v] gives the transmission coefficient for quantum number m, scale factor s, and sampled potential v as a function of energy.";


TransmissionProbability::usage="TransmissionProbability[{m,s},v] gives the transmission probability for quantum number m, scale factor s, and sampled potential v as a function of energy.";


ReflectionProbability::usage="ReflectionProbability[{m,s},v] gives the reflection probability for quantum number m, scale factor s, and sampled potential v as a function of energy.";


ReflectionCoefficient::usage="ReflectionCoefficient[{m,s},v] gives the reflection coefficient for quantum number m, scale factor s, and sampled potential v as a function of energy.";


TransmissionPlot::usage="TransmissionPlot[{m,s},v,r] plots the transmission coefficient T for quantum number m, scale factor s, and sampled potential v as a function of energy with energy resolution r.";


ReflectionPlot::usage="ReflectionPlot[{m,s},v,r] plots the reflection coefficient R for quantum number m, scale factor s, and sampled potential v as a function of energy with energy resolution r.";


(* ::Subsection:: *)
(*Begin Private context*)


Begin["`Private`"];


(* ::Section:: *)
(*R-matrix*)


kMax[{m_,s_},v_List]:=Max[kp[{m,s},v][[All,2]]]


dyadic[v_List,p_]:=dyadic[v,p]=WignerEisenbudFunctions[v][[p,{1,-1}]] \[CircleTimes] WignerEisenbudFunctions[v][[p,{1,-1}]]


RBlock[{m_,s_},v_][energy_][{e_,k_,p_}]:=SparseArray[{{k,k}->\[Pi]/2 1/(energy-e)},{kMax[{m,s},v],kMax[{m,s},v]}] \[CircleTimes] dyadic[v,p]


RMatrix[{m_,s_},v_List]:= RMatrix[{m,s},v] = Function[energy, Evaluate[Total[RBlock[{m,s},v][energy] /@ kp[{m,s},v]]]]


(* ::Section:: *)
(*\[CapitalTheta]-matrix and K-matrix*)


StepMatrix/:StepMatrix[{m_,s_,n_}]:=StepMatrix/:StepMatrix[{m,s,n}]=
	Function[energy, SparseArray[{{k_,k_}:>1/; Length[Cases[Thread[s^2 JZeros[m]^2<energy],True]]>=k},{n,n}] \[CircleTimes] IdentityMatrix[2]]


KMatrix[{m_,s_,n_}]:=KMatrix[{m,s,n}]=
	Function[energy, SparseArray[{{k_,k_}:>2/\[Pi] Sqrt[energy-TransversalEnergy[m,k,s]]},{n,n}] \[CircleTimes] IdentityMatrix[2]]


(* ::Section:: *)
(*Current scattering matrix *)


OmegaMatrix[{m_,s_},v_List][energy_]:=Module[{k=kMax[{m,s},v],K12},
	K12=MatrixPower[KMatrix[{m,s,k}][energy],0.5];
	K12.RMatrix[{m,s},v][energy].K12
]


CurrentScatteringMatrix[{m_,s_},v_List][energy_]:=Module[{k=kMax[{m,s},v],\[CapitalTheta],Id,\[CapitalOmega]},
	\[CapitalTheta] = StepMatrix[{m,s,k}][energy];
	Id = IdentityMatrix[2k];
	\[CapitalOmega] = OmegaMatrix[{m,s},v][energy];
	\[CapitalTheta].(Id-2 Inverse[Id+I \[CapitalOmega]]).\[CapitalTheta]
]


(* ::Section:: *)
(*Reflection and transmission probabilities*)


TransmissionProbability[{m_,s_},v_List][energy_]:=
	Total[Abs[Diagonal[CurrentScatteringMatrix[{m,s},v][energy],1]]^2]


ReflectionProbability[{m_,s_},v_List][energy_]:=
	Total[Abs[Diagonal[CurrentScatteringMatrix[{m,s},v][energy]]]^2]/2


(* ::Section:: *)
(*Transport plots*)


TransportPlot[f_,{m_,s_},v_List,resolution_,opts___]:=
ListPlot[
	Table[{energy,N[f[{m,s},v][energy]]},{energy,0,1500,resolution}], 
	opts,
	Mesh->All,
	PlotRange->All
]


TransmissionPlot[{m_,s_},v_List,resolution_,opts___]:=
	TransportPlot[TransmissionProbability,{m,s},v,resolution,AxesLabel->{Style["E",Italic],Style["T(E)",Italic]},opts]


ReflectionPlot[{m_,s_},v_List,resolution_,opts___]:=
	TransportPlot[ReflectionProbability,{m,s},v,resolution,AxesLabel->{Style["E",Italic],Style["R(E)",Italic]},opts]


(* ::Section:: *)
(*Package Ending*)


End[]; 


(* ::Subsection:: *)
(*Setting Attributes*)


(*SetAttributes[{...}, {Listable, Protected, ReadProtected, NumericFunction}]; *)


(* ::Subsection:: *)
(*End package*)


EndPackage[];
