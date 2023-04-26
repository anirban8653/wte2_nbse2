(* ::Package:: *)

SetDirectory["~/Desktop/wte2/case_a"]

kernum=80;
ldos  = 0;
\[Eta] = 0.001;
wmax = 0.2; 
nw = 400;
it=0;
Nx=10000;
Ny=200;
Nband=8;
mu=0.00;
delta=0;
For[iy= 1, iy<= Ny,iy++,it++;
    ypos[it] = iy;];
      Ny=it;


gso[kxs_,rxyp_]:={{0.295217 DiracDelta[rxyp]-0.327133 Cos[1. kxs] DiracDelta[rxyp]-0.00638308 DiracDelta[rxyp] Sin[1. kxs],
(0.-0.00478731 I) DiracDelta[0.61-rxyp]-(0.+0.00478731 I) DiracDelta[0.39+rxyp]+0.0558519 E^((0.-1. I) kxs) DiracDelta[0.39+rxyp]-
0.0558519 E^((0.+1. I) kxs) DiracDelta[0.39+rxyp],0.203461 E^((0.-0.5 I) kxs) DiracDelta[0.36-rxyp]+0.203461 E^((0.+0.5 I) kxs) DiracDelta[0.36-rxyp],
0.155587 E^((0.-0.5 I) kxs) DiracDelta[0.25+rxyp]-0.155587 E^((0.+0.5 I) kxs) DiracDelta[0.25+rxyp]+
0.115693 E^((0.-1.5 I) kxs) DiracDelta[0.25+rxyp]-0.115693 E^((0.+1.5 I) kxs) DiracDelta[0.25+rxyp],
(0.+0.0247344 I) DiracDelta[rxyp] Sin[1. kxs],-0.0199471 DiracDelta[0.61-rxyp]-0.0203461 DiracDelta[0.39+rxyp],
0,-0.00438837 E^((0.-0.5 I) kxs) DiracDelta[0.25+rxyp]-0.00438837 E^((0.+0.5 I) kxs) DiracDelta[0.25+rxyp]},
{(0.+0.00478731 I) DiracDelta[0.39-rxyp]-0.0558519 E^((0.-1. I) kxs) DiracDelta[0.39-rxyp]+0.0558519 E^((0.+1. I) kxs) DiracDelta[0.39-rxyp]+
(0.+0.00478731 I) DiracDelta[0.61+rxyp],0.0518625 DiracDelta[1.-rxyp]-0.698149 DiracDelta[rxyp]+
(0.450805-0.00398942 I) E^((0.-1. I) kxs) DiracDelta[rxyp]+(0.450805+0.00398942 I) E^((0.+1. I) kxs) DiracDelta[rxyp]+
0.0518625 DiracDelta[1.+rxyp],-0.155587 E^((0.-0.5 I) kxs) DiracDelta[0.25+rxyp]+0.155587 E^((0.+0.5 I) kxs) DiracDelta[0.25+rxyp]-
0.115693 E^((0.-1.5 I) kxs) DiracDelta[0.25+rxyp]+0.115693 E^((0.+1.5 I) kxs) DiracDelta[0.25+rxyp],0.159577 E^((0.-0.5 I) kxs) DiracDelta[0.14-rxyp]+
0.159577 E^((0.+0.5 I) kxs) DiracDelta[0.14-rxyp],0.0203461 DiracDelta[0.39-rxyp]+0.0199471 DiracDelta[0.61+rxyp],(0.+0.0319154 I) DiracDelta[rxyp] Sin[1. kxs],
-0.00438837 E^((0.-0.5 I) kxs) DiracDelta[0.25+rxyp]-0.00438837 E^((0.+0.5 I) kxs) DiracDelta[0.25+rxyp],0},{0.203461 E^((0.-0.5 I) kxs) DiracDelta[0.36+rxyp]+
0.203461 E^((0.+0.5 I) kxs) DiracDelta[0.36+rxyp],0.155587 E^((0.-0.5 I) kxs) DiracDelta[0.25-rxyp]-0.155587 E^((0.+0.5 I) kxs) DiracDelta[0.25-rxyp]+
0.115693 E^((0.-1.5 I) kxs) DiracDelta[0.25-rxyp]-0.115693 E^((0.+1.5 I) kxs) DiracDelta[0.25-rxyp],0.295217 DiracDelta[rxyp]-
0.327133 Cos[1. kxs] DiracDelta[rxyp]+0.00638308 DiracDelta[rxyp] Sin[1. kxs],(0.+0.00478731 I) DiracDelta[0.39-rxyp]+
0.0558519 E^((0.-1. I) kxs) DiracDelta[0.39-rxyp]-0.0558519 E^((0.+1. I) kxs) DiracDelta[0.39-rxyp]+(0.+0.00478731 I) DiracDelta[0.61+rxyp],
0,0.00438837 E^((0.-0.5 I) kxs) DiracDelta[0.25-rxyp]+0.00438837 E^((0.+0.5 I) kxs) DiracDelta[0.25-rxyp],(0.-0.0247344 I) DiracDelta[rxyp] Sin[1. kxs],
0.0203461 DiracDelta[0.39-rxyp]+0.0199471 DiracDelta[0.61+rxyp]},{-0.155587 E^((0.-0.5 I) kxs) DiracDelta[0.25-rxyp]+
0.155587 E^((0.+0.5 I) kxs) DiracDelta[0.25-rxyp]-0.115693 E^((0.-1.5 I) kxs) DiracDelta[0.25-rxyp]+0.115693 E^((0.+1.5 I) kxs) DiracDelta[0.25-rxyp],
0.159577 E^((0.-0.5 I) kxs) DiracDelta[0.14+rxyp]+0.159577 E^((0.+0.5 I) kxs) DiracDelta[0.14+rxyp],(0.-0.00478731 I) DiracDelta[0.61-rxyp]-
(0.+0.00478731 I) DiracDelta[0.39+rxyp]-0.0558519 E^((0.-1. I) kxs) DiracDelta[0.39+rxyp]+0.0558519 E^((0.+1. I) kxs) DiracDelta[0.39+rxyp],
0.0518625 DiracDelta[1.-rxyp]-0.698149 DiracDelta[rxyp]+(0.450805+0.00398942 I) E^((0.-1. I) kxs) DiracDelta[rxyp]+
(0.450805-0.00398942 I) E^((0.+1. I) kxs) DiracDelta[rxyp]+0.0518625 DiracDelta[1.+rxyp],0.00438837 E^((0.-0.5 I) kxs) DiracDelta[0.25-rxyp]+
0.00438837 E^((0.+0.5 I) kxs) DiracDelta[0.25-rxyp],0,-0.0199471 DiracDelta[0.61-rxyp]-0.0203461 DiracDelta[0.39+rxyp],
(0.-0.0319154 I) DiracDelta[rxyp] Sin[1. kxs]},{(0.-0.0247344 I) DiracDelta[rxyp] Sin[1. kxs],0.0199471 DiracDelta[0.61-rxyp]+
0.0203461 DiracDelta[0.39+rxyp],0,0.00438837 E^((0.-0.5 I) kxs) DiracDelta[0.25+rxyp]+0.00438837 E^((0.+0.5 I) kxs) DiracDelta[0.25+rxyp],
0.295217 DiracDelta[rxyp]-0.327133 Cos[1. kxs] DiracDelta[rxyp]+0.00638308 DiracDelta[rxyp] Sin[1. kxs],(0.+0.00478731 I) DiracDelta[0.61-rxyp]+
(0.+0.00478731 I) DiracDelta[0.39+rxyp]+0.0558519 E^((0.-1. I) kxs) DiracDelta[0.39+rxyp]-0.0558519 E^((0.+1. I) kxs) DiracDelta[0.39+rxyp],
0.203461 E^((0.-0.5 I) kxs) DiracDelta[0.36-rxyp]+0.203461 E^((0.+0.5 I) kxs) DiracDelta[0.36-rxyp],0.155587 E^((0.-0.5 I) kxs) DiracDelta[0.25+rxyp]-
0.155587 E^((0.+0.5 I) kxs) DiracDelta[0.25+rxyp]+0.115693 E^((0.-1.5 I) kxs) DiracDelta[0.25+rxyp]-0.115693 E^((0.+1.5 I) kxs) DiracDelta[0.25+rxyp]},
{-0.0203461 DiracDelta[0.39-rxyp]-0.0199471 DiracDelta[0.61+rxyp],(0.-0.0319154 I) DiracDelta[rxyp] Sin[1. kxs],0.00438837 E^((0.-0.5 I) kxs) DiracDelta[0.25+rxyp]+
0.00438837 E^((0.+0.5 I) kxs) DiracDelta[0.25+rxyp],0,(0.-0.00478731 I) DiracDelta[0.39-rxyp]-0.0558519 E^((0.-1. I) kxs) DiracDelta[0.39-rxyp]+
0.0558519 E^((0.+1. I) kxs) DiracDelta[0.39-rxyp]-(0.+0.00478731 I) DiracDelta[0.61+rxyp],0.0518625 DiracDelta[1.-rxyp]-0.698149 DiracDelta[rxyp]+
(0.450805+0.00398942 I) E^((0.-1. I) kxs) DiracDelta[rxyp]+(0.450805-0.00398942 I) E^((0.+1. I) kxs) DiracDelta[rxyp]+0.0518625 DiracDelta[1.+rxyp],
-0.155587 E^((0.-0.5 I) kxs) DiracDelta[0.25+rxyp]+0.155587 E^((0.+0.5 I) kxs) DiracDelta[0.25+rxyp]-0.115693 E^((0.-1.5 I) kxs) DiracDelta[0.25+rxyp]+
0.115693 E^((0.+1.5 I) kxs) DiracDelta[0.25+rxyp],0.159577 E^((0.-0.5 I) kxs) DiracDelta[0.14-rxyp]+0.159577 E^((0.+0.5 I) kxs) DiracDelta[0.14-rxyp]},
{0,-0.00438837 E^((0.-0.5 I) kxs) DiracDelta[0.25-rxyp]-0.00438837 E^((0.+0.5 I) kxs) DiracDelta[0.25-rxyp],(0.+0.0247344 I) DiracDelta[rxyp] Sin[1. kxs],
-0.0203461 DiracDelta[0.39-rxyp]-0.0199471 DiracDelta[0.61+rxyp],0.203461 E^((0.-0.5 I) kxs) DiracDelta[0.36+rxyp]+0.203461 E^((0.+0.5 I) kxs) DiracDelta[0.36+rxyp],
0.155587 E^((0.-0.5 I) kxs) DiracDelta[0.25-rxyp]-0.155587 E^((0.+0.5 I) kxs) DiracDelta[0.25-rxyp]+0.115693 E^((0.-1.5 I) kxs) DiracDelta[0.25-rxyp]-
0.115693 E^((0.+1.5 I) kxs) DiracDelta[0.25-rxyp],0.295217 DiracDelta[rxyp]-0.327133 Cos[1. kxs] DiracDelta[rxyp]-0.00638308 DiracDelta[rxyp] Sin[1. kxs],
(0.-0.00478731 I) DiracDelta[0.39-rxyp]+0.0558519 E^((0.-1. I) kxs) DiracDelta[0.39-rxyp]-0.0558519 E^((0.+1. I) kxs) DiracDelta[0.39-rxyp]-
(0.+0.00478731 I) DiracDelta[0.61+rxyp]},{-0.00438837 E^((0.-0.5 I) kxs) DiracDelta[0.25-rxyp]-0.00438837 E^((0.+0.5 I) kxs) DiracDelta[0.25-rxyp],
0,0.0199471 DiracDelta[0.61-rxyp]+0.0203461 DiracDelta[0.39+rxyp],(0.+0.0319154 I) DiracDelta[rxyp] Sin[1. kxs],-0.155587 E^((0.-0.5 I) kxs) DiracDelta[0.25-rxyp]+
0.155587 E^((0.+0.5 I) kxs) DiracDelta[0.25-rxyp]-0.115693 E^((0.-1.5 I) kxs) DiracDelta[0.25-rxyp]+0.115693 E^((0.+1.5 I) kxs) DiracDelta[0.25-rxyp],
0.159577 E^((0.-0.5 I) kxs) DiracDelta[0.14+rxyp]+0.159577 E^((0.+0.5 I) kxs) DiracDelta[0.14+rxyp],(0.+0.00478731 I) DiracDelta[0.61-rxyp]+
(0.+0.00478731 I) DiracDelta[0.39+rxyp]-0.0558519 E^((0.-1. I) kxs) DiracDelta[0.39+rxyp]+0.0558519 E^((0.+1. I) kxs) DiracDelta[0.39+rxyp],
0.0518625 DiracDelta[1.-rxyp]-0.698149 DiracDelta[rxyp]+(0.450805-0.00398942 I) E^((0.-1. I) kxs) DiracDelta[rxyp]+(0.450805+0.00398942 I) E^((0.+1. I) kxs) DiracDelta[rxyp]+0.0518625 DiracDelta[1.+rxyp]}};


Clear[t,tijmn];
chempot=SparseArray[{x_,x_}->-mu,{Ny*Nband,Ny*Nband}]



(*Superconductivity Block defined*)
(*empty=SparseArray[{i_,i_}\[Rule]0,{Nx*Nband/2,Nx*Nband/2}];*)
(*Hs1=SparseArray[{ila_,ila_}->delta,{Ny*Nband/2,Ny*Nband/2}]
Hs2=SparseArray[{ila_,jla_}/;NumericQ[abracadabra]->delta,{Ny*Nband/2,Ny*Nband/2}]
Hs=ArrayFlatten[{{Hs2,Hs1},{-Hs1,Hs2}}];
Clear[Hs1];
Clear[Hs2];*)


r[1]={0.,0.};
r[2]={0.,-0.39};
r[3]={0.5,-0.64};
r[4]={0.5,-0.25};
r[5]=r[1];r[6]=r[2];r[7]=r[3];r[8]=r[4];
(*tijmn=Table[0,
{i,1,Nx},{j,1,Nx},{m,1,Nband},{n,1,Nband}];*)
(*ival=Flatten[Table[{i,i+1},{i,  1,(nsite - 1)*Nband + Nband,800}]];*)
(*ival={3,4,11,12,19,20,27,28,35,36,43,44,51,52,59,60,67,68,75,76,83,84,91,92,99,100,107,108,115,116,123,124,131,132,139,140,147,148,803,804}*)
ival={1,2,3,4,97,98,99,100,197,198,199,200};


For[lx=1,lx<=Ny*Nband,lx++,
If[lx<= Ny*Nband/2,
m=Mod[lx-1/10,Nband/2]+1/10;i=Quotient[lx-0.1,Nband/2]+1;,
m=Mod[lx-1/10,Nband/2]+1/10+4;i=Quotient[lx-0.1,Nband/2]+1-Ny;];
(*If[ly\[LessEqual] Ny*Nband/2,
n=Mod[ly-1/10,Nband/2]+1/10;j=Quotient[ly-0.1,Nband/2]+1;,
n=Mod[ly-1/10,Nband/2]+1/10+4;j=Quotient[ly-0.1,Nband/2]+1-Ny;];*)
iv[lx]=ypos[i];mv[lx]=m;
];


CloseKernals[];
LaunchKernals[kernum];
Print["start calculation"];

ldos=Total[ParallelTable[
kxs=(2*jx/Nx)*Pi;
For[m=1,m<=8,m++,
For[n=1,n<=8,n++,
For[cy=-2,cy<=2,cy++,
rys=cy-r[m][[2]]+r[n][[2]];
t[jx][-cy][m][n]=Chop[Sqrt[2*Pi]*gso[kxs,rys][[m,n]]/DiracDelta[0]];
]]];

Hop1=SparseArray[{lxs_,lys_}/;NumericQ[t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]]]->t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]],{Ny*Nband,Ny*Nband}];(*
Hop2=SparseArray[{lxs_,lys_}/;NumericQ[t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]]]->If[(lxs<=Ny*Nband/2 && lys<=Ny*Nband/2 ||lxs>= Ny*Nband/2 && lys>= Ny*Nband/2),
 -Conjugate[t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]]],Conjugate[t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]]]],{Ny*Nband,Ny*Nband}];*)

(*H=ArrayFlatten[{{(Hop1+chempot),Hs},{ConjugateTranspose[Hs],(Hop2-chempot)}}];
*)

ES=Eigensystem[Hop1];
	
elist = Table[ES[[1,l]], {l, 1, Nband*Ny}];
ulist = Table[Abs[ES[[2,l,ival]]]^2, {l, 1, Nband*Ny}];
(*vlist = Table[Abs[ES[[2,l,ival + Nband*Ny]]]^2, {l, 1, 2*Nband*Ny}]; *) 
      
ldossc =Table[Im[Total[Table[ulist[[l,1 ;; All]]*(1/(-wmax + wmax*2*(w/nw) - elist[[l]] + I*\[Eta]))(*+ 
      vlist[[l,1 ;; All]]*(1/(-wmax + 2*wmax*(w/nw) + elist[[l]] + I*\[Eta]))*), {l, 1, Nband*Ny}]]], {w, 0, nw}]; 

ldossc,{jx,0,Nx}]];


fl = StringJoin["x-edge","wmax",ToString[wmax],"Ny",ToString[Ny],"Nx",ToString[Nx],
"eta",ToString[\[Eta]],"Ns",ToString[Ns],"mu",ToString[mu]];
DeleteFile[StringJoin[fl,"Ldos.dat"]];
str1 = OpenAppend[StringJoin[fl,"Ldos.dat"], FormatType -> StandardForm];
For[n=1,n<=nw+1,n++,
For[i=1,i<=Length[ival],i++,
Write[str1,-wmax + wmax*2*((n-1)/nw),"  ", (-Pi^(-1))*ldos[[n,i]]/(Nx+1)];
]
];
Close[str1];

Exit[];
