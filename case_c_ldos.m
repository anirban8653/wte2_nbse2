(* ::Package:: *)

SetDirectory["~/Desktop/NTU/WTe2/case_a"]


ldos  = 0;
\[Eta] = 0.004;
wmax =0.2; 
nw = 400;
it=0;
Nx=200;
Ny=5000;
Nband=8;
mu=0.0;
delta=0;
For[ix = 1, ix <= Nx, ix++,it++;
    xpos[it] = ix;];
      Nx=it;


gso[kys_,rxsp_]:={{(-0.163566-0.00319154 I) DiracDelta[1.-rxsp]+0.295217 DiracDelta[rxsp]-(0.163566-0.00319154 I) DiracDelta[1.+rxsp],0.0558519 E^((0.+0.39 I) kys) DiracDelta[1.-rxsp]-
(0.+0.00478731 I) E^((0.+0.39 I) kys) DiracDelta[rxsp]-(0.+0.00478731 I) E^((0.-0.61 I) kys) DiracDelta[rxsp]-0.0558519 E^((0.+0.39 I) kys) DiracDelta[1.+rxsp],
0.203461 E^((0.-0.36 I) kys) DiracDelta[0.5-rxsp]+0.203461 E^((0.-0.36 I) kys) DiracDelta[0.5+rxsp],0.155587 E^((0.+0.25 I) kys) DiracDelta[0.5-rxsp]+
0.115693 E^((0.+0.25 I) kys) DiracDelta[1.5-rxsp]-0.155587 E^((0.+0.25 I) kys) DiracDelta[0.5+rxsp]-0.115693 E^((0.+0.25 I) kys) DiracDelta[1.5+rxsp],
-0.0123672 DiracDelta[1.-rxsp]+0.0123672 DiracDelta[1.+rxsp],-0.0203461 E^((0.+0.39 I) kys) DiracDelta[rxsp]-0.0199471 E^((0.-0.61 I) kys) DiracDelta[rxsp],
0,-0.00438837 E^((0.+0.25 I) kys) DiracDelta[0.5-rxsp]-0.00438837 E^((0.+0.25 I) kys) DiracDelta[0.5+rxsp]},{-0.0558519 E^((0.-0.39 I) kys) DiracDelta[1.-rxsp]+
(0.+0.00478731 I) E^((0.-0.39 I) kys) DiracDelta[rxsp]+(0.+0.00478731 I) E^((0.+0.61 I) kys) DiracDelta[rxsp]+0.0558519 E^((0.-0.39 I) kys) DiracDelta[1.+rxsp],
(0.450805-0.00398942 I) DiracDelta[1.-rxsp]-0.698149 DiracDelta[rxsp]+0.0518625 E^((0.-1. I) kys) DiracDelta[rxsp]+0.0518625 E^((0.+1. I) kys) DiracDelta[rxsp]+
(0.450805+0.00398942 I) DiracDelta[1.+rxsp],-0.155587 E^((0.+0.25 I) kys) DiracDelta[0.5-rxsp]-0.115693 E^((0.+0.25 I) kys) DiracDelta[1.5-rxsp]+
0.155587 E^((0.+0.25 I) kys) DiracDelta[0.5+rxsp]+0.115693 E^((0.+0.25 I) kys) DiracDelta[1.5+rxsp],0.159577 E^((0.-0.14 I) kys) DiracDelta[0.5-rxsp]+
0.159577 E^((0.-0.14 I) kys) DiracDelta[0.5+rxsp],0.0203461 E^((0.-0.39 I) kys) DiracDelta[rxsp]+0.0199471 E^((0.+0.61 I) kys) DiracDelta[rxsp],
-0.0159577 DiracDelta[1.-rxsp]+0.0159577 DiracDelta[1.+rxsp],-0.00438837 E^((0.+0.25 I) kys) DiracDelta[0.5-rxsp]-0.00438837 E^((0.+0.25 I) kys) DiracDelta[0.5+rxsp],0}
,{0.203461 E^((0.+0.36 I) kys) DiracDelta[0.5-rxsp]+0.203461 E^((0.+0.36 I) kys) DiracDelta[0.5+rxsp],0.155587 E^((0.-0.25 I) kys) DiracDelta[0.5-rxsp]+
0.115693 E^((0.-0.25 I) kys) DiracDelta[1.5-rxsp]-0.155587 E^((0.-0.25 I) kys) DiracDelta[0.5+rxsp]-0.115693 E^((0.-0.25 I) kys) DiracDelta[1.5+rxsp],
(-0.163566+0.00319154 I) DiracDelta[1.-rxsp]+0.295217 DiracDelta[rxsp]-(0.163566+0.00319154 I) DiracDelta[1.+rxsp],0.0558519 E^((0.-0.39 I) kys) DiracDelta[1.-rxsp]
+(0.+0.00478731 I) E^((0.-0.39 I) kys) DiracDelta[rxsp]+(0.+0.00478731 I) E^((0.+0.61 I) kys) DiracDelta[rxsp]-0.0558519 E^((0.-0.39 I) kys) DiracDelta[1.+rxsp]
,0,0.00438837 E^((0.-0.25 I) kys) DiracDelta[0.5-rxsp]+0.00438837 E^((0.-0.25 I) kys) DiracDelta[0.5+rxsp],0.0123672 DiracDelta[1.-rxsp]-0.0123672 DiracDelta[1.+rxsp],
0.0203461 E^((0.-0.39 I) kys) DiracDelta[rxsp]+0.0199471 E^((0.+0.61 I) kys) DiracDelta[rxsp]},{-0.155587 E^((0.-0.25 I) kys) DiracDelta[0.5-rxsp]-
0.115693 E^((0.-0.25 I) kys) DiracDelta[1.5-rxsp]+0.155587 E^((0.-0.25 I) kys) DiracDelta[0.5+rxsp]+0.115693 E^((0.-0.25 I) kys) DiracDelta[1.5+rxsp],
0.159577 E^((0.+0.14 I) kys) DiracDelta[0.5-rxsp]+0.159577 E^((0.+0.14 I) kys) DiracDelta[0.5+rxsp],-0.0558519 E^((0.+0.39 I) kys) DiracDelta[1.-rxsp]-
(0.+0.00478731 I) E^((0.+0.39 I) kys) DiracDelta[rxsp]-(0.+0.00478731 I) E^((0.-0.61 I) kys) DiracDelta[rxsp]+0.0558519 E^((0.+0.39 I) kys) DiracDelta[1.+rxsp],
(0.450805+0.00398942 I) DiracDelta[1.-rxsp]-0.698149 DiracDelta[rxsp]+0.0518625 E^((0.-1. I) kys) DiracDelta[rxsp]+0.0518625 E^((0.+1. I) kys) DiracDelta[rxsp]+
(0.450805-0.00398942 I) DiracDelta[1.+rxsp],0.00438837 E^((0.-0.25 I) kys) DiracDelta[0.5-rxsp]+0.00438837 E^((0.-0.25 I) kys) DiracDelta[0.5+rxsp],0,
-0.0203461 E^((0.+0.39 I) kys) DiracDelta[rxsp]-0.0199471 E^((0.-0.61 I) kys) DiracDelta[rxsp],0.0159577 DiracDelta[1.-rxsp]-0.0159577 DiracDelta[1.+rxsp]},
{0.0123672 DiracDelta[1.-rxsp]-0.0123672 DiracDelta[1.+rxsp],0.0203461 E^((0.+0.39 I) kys) DiracDelta[rxsp]+0.0199471 E^((0.-0.61 I) kys) DiracDelta[rxsp],
0,0.00438837 E^((0.+0.25 I) kys) DiracDelta[0.5-rxsp]+0.00438837 E^((0.+0.25 I) kys) DiracDelta[0.5+rxsp],(-0.163566+0.00319154 I) DiracDelta[1.-rxsp]+
0.295217 DiracDelta[rxsp]-(0.163566+0.00319154 I) DiracDelta[1.+rxsp],0.0558519 E^((0.+0.39 I) kys) DiracDelta[1.-rxsp]+(0.+0.00478731 I) E^((0.+0.39 I) kys) DiracDelta[rxsp]+
(0.+0.00478731 I) E^((0.-0.61 I) kys) DiracDelta[rxsp]-0.0558519 E^((0.+0.39 I) kys) DiracDelta[1.+rxsp],0.203461 E^((0.-0.36 I) kys) DiracDelta[0.5-rxsp]+
0.203461 E^((0.-0.36 I) kys) DiracDelta[0.5+rxsp],0.155587 E^((0.+0.25 I) kys) DiracDelta[0.5-rxsp]+0.115693 E^((0.+0.25 I) kys) DiracDelta[1.5-rxsp]-
0.155587 E^((0.+0.25 I) kys) DiracDelta[0.5+rxsp]-0.115693 E^((0.+0.25 I) kys) DiracDelta[1.5+rxsp]},{-0.0203461 E^((0.-0.39 I) kys) DiracDelta[rxsp]-
0.0199471 E^((0.+0.61 I) kys) DiracDelta[rxsp],0.0159577 DiracDelta[1.-rxsp]-0.0159577 DiracDelta[1.+rxsp],0.00438837 E^((0.+0.25 I) kys) DiracDelta[0.5-rxsp]+
0.00438837 E^((0.+0.25 I) kys) DiracDelta[0.5+rxsp],0,-0.0558519 E^((0.-0.39 I) kys) DiracDelta[1.-rxsp]-(0.+0.00478731 I) E^((0.-0.39 I) kys) DiracDelta[rxsp]-
(0.+0.00478731 I) E^((0.+0.61 I) kys) DiracDelta[rxsp]+0.0558519 E^((0.-0.39 I) kys) DiracDelta[1.+rxsp],(0.450805+0.00398942 I) DiracDelta[1.-rxsp]-
0.698149 DiracDelta[rxsp]+0.0518625 E^((0.-1. I) kys) DiracDelta[rxsp]+0.0518625 E^((0.+1. I) kys) DiracDelta[rxsp]+(0.450805-0.00398942 I) DiracDelta[1.+rxsp],
-0.155587 E^((0.+0.25 I) kys) DiracDelta[0.5-rxsp]-0.115693 E^((0.+0.25 I) kys) DiracDelta[1.5-rxsp]+0.155587 E^((0.+0.25 I) kys) DiracDelta[0.5+rxsp]+
0.115693 E^((0.+0.25 I) kys) DiracDelta[1.5+rxsp],0.159577 E^((0.-0.14 I) kys) DiracDelta[0.5-rxsp]+0.159577 E^((0.-0.14 I) kys) DiracDelta[0.5+rxsp]},
{0,-0.00438837 E^((0.-0.25 I) kys) DiracDelta[0.5-rxsp]-0.00438837 E^((0.-0.25 I) kys) DiracDelta[0.5+rxsp],-0.0123672 DiracDelta[1.-rxsp]+
0.0123672 DiracDelta[1.+rxsp],-0.0203461 E^((0.-0.39 I) kys) DiracDelta[rxsp]-0.0199471 E^((0.+0.61 I) kys) DiracDelta[rxsp],0.203461 E^((0.+0.36 I) kys) DiracDelta[0.5-rxsp]
+0.203461 E^((0.+0.36 I) kys) DiracDelta[0.5+rxsp],0.155587 E^((0.-0.25 I) kys) DiracDelta[0.5-rxsp]+0.115693 E^((0.-0.25 I) kys) DiracDelta[1.5-rxsp]-
0.155587 E^((0.-0.25 I) kys) DiracDelta[0.5+rxsp]-0.115693 E^((0.-0.25 I) kys) DiracDelta[1.5+rxsp],(-0.163566-0.00319154 I) DiracDelta[1.-rxsp]+
0.295217 DiracDelta[rxsp]-(0.163566-0.00319154 I) DiracDelta[1.+rxsp],0.0558519 E^((0.-0.39 I) kys) DiracDelta[1.-rxsp]-(0.+0.00478731 I) E^((0.-0.39 I) kys) DiracDelta[rxsp]-
(0.+0.00478731 I) E^((0.+0.61 I) kys) DiracDelta[rxsp]-0.0558519 E^((0.-0.39 I) kys) DiracDelta[1.+rxsp]},{-0.00438837 E^((0.-0.25 I) kys) DiracDelta[0.5-rxsp]-
0.00438837 E^((0.-0.25 I) kys) DiracDelta[0.5+rxsp],0,0.0203461 E^((0.+0.39 I) kys) DiracDelta[rxsp]+0.0199471 E^((0.-0.61 I) kys) DiracDelta[rxsp],
-0.0159577 DiracDelta[1.-rxsp]+0.0159577 DiracDelta[1.+rxsp],-0.155587 E^((0.-0.25 I) kys) DiracDelta[0.5-rxsp]-0.115693 E^((0.-0.25 I) kys) DiracDelta[1.5-rxsp]+
0.155587 E^((0.-0.25 I) kys) DiracDelta[0.5+rxsp]+0.115693 E^((0.-0.25 I) kys) DiracDelta[1.5+rxsp],0.159577 E^((0.+0.14 I) kys) DiracDelta[0.5-rxsp]+
0.159577 E^((0.+0.14 I) kys) DiracDelta[0.5+rxsp],-0.0558519 E^((0.+0.39 I) kys) DiracDelta[1.-rxsp]+(0.+0.00478731 I) E^((0.+0.39 I) kys) DiracDelta[rxsp]+
(0.+0.00478731 I) E^((0.-0.61 I) kys) DiracDelta[rxsp]+0.0558519 E^((0.+0.39 I) kys) DiracDelta[1.+rxsp],(0.450805-0.00398942 I) DiracDelta[1.-rxsp]-
0.698149 DiracDelta[rxsp]+0.0518625 E^((0.-1. I) kys) DiracDelta[rxsp]+0.0518625 E^((0.+1. I) kys) DiracDelta[rxsp]+(0.450805+0.00398942 I) DiracDelta[1.+rxsp]}};


Clear[t,tijmn];
chempot=SparseArray[{x_,x_}->-mu,{Nx*Nband,Nx*Nband}];


(*
(*Superconductivity Block defined*)
(*empty=SparseArray[{i_,i_}\[Rule]0,{Nx*Nband/2,Nx*Nband/2}];*)
Hs1=SparseArray[{ila_,ila_}\[Rule]delta,{Nx*Nband/2,Nx*Nband/2}]
Hs2=SparseArray[{ila_,jla_}/;NumericQ[abracadabra]\[Rule]delta,{Nx*Nband/2,Nx*Nband/2}]
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
(*ival=Join[{1,2,3,4},{Nx*Nband/2-3,Nx*Nband/2-2,Nx*Nband/2-1,Nx*Nband/2},{Nx*Nband-3,Nx*Nband-2,Nx*Nband-1,Nx*Nband}];*)
ival ={1,2,3,4,401,402,403,404,797,798,799,800}


For[lx=1,lx<=Nx*Nband,lx++,
If[lx<= Nx*Nband/2,
m=Mod[lx-1/10,Nband/2]+1/10;i=Quotient[lx-0.1,Nband/2]+1;,
m=Mod[lx-1/10,Nband/2]+1/10+4;i=Quotient[lx-0.1,Nband/2]+1-Nx;];
(*If[ly\[LessEqual] Nx*Nband/2,
n=Mod[ly-1/10,Nband/2]+1/10;j=Quotient[ly-0.1,Nband/2]+1;,
n=Mod[ly-1/10,Nband/2]+1/10+4;j=Quotient[ly-0.1,Nband/2]+1-Nx;];*)
iv[lx]=xpos[i];mv[lx]=m;
];



AbsoluteTiming[
ldos=Total[ParallelTable[
ky=(2*jy/(Ny)-1)*(Pi);
For[m=1,m<=8,m++,
For[n=1,n<=8,n++,
For[cx=-2,cx<=2,cx++,
rxs=cx+r[m][[1]]-r[n][[1]];(*rys=cy-r[m][[2]]+r[n][[2]];*)
t[-cx][m][n]=Chop[Sqrt[2 Pi]*gso[ky,rxs][[m,n]]/DiracDelta[0]];
]]];
Hop1=SparseArray[{lxs_,lys_}/;NumericQ[t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]]]->t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]],{Nx*Nband,Nx*Nband}]+chempot;
(*Hop2=SparseArray[{lxs_,lys_}/;NumericQ[t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]]]\[Rule]If[(lxs\[LessEqual]Nx*Nband/2 && lys\[LessEqual]Nx*Nband/2 ||lxs\[GreaterEqual] Nx*Nband/2 && lys\[GreaterEqual] Nx*Nband/2),
 -Conjugate[t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]]],Conjugate[t[iv[lxs]-iv[lys]][mv[lxs]][mv[lys]]]],{Nx*Nband,Nx*Nband}];*)

(*H=ArrayFlatten[{{(Hop1+chempot),Hs},{ConjugateTranspose[Hs],(Hop2-chempot)}}];*)


ES=Eigensystem[Hop1];

elist = Table[ES[[1,l]], {l, 1, Nband*Nx}];
ulist = Table[Abs[ES[[2,l,ival]]]^2, {l, 1, Nband*Nx}];
(*vlist = Table[Abs[ES[[2,l,ival + Nband*Nx]]]^2, {l, 1, 2*Nband*Nx}]; *) 
      
ldossc =Table[Im[Total[Table[ulist[[l,1 ;; All]]*(1/(-wmax + wmax*2*(w/nw) - elist[[l]] + I*\[Eta]))(*+ 
      vlist[[l,1 ;; All]]*(1/(-wmax + 2*wmax*(w/nw) + elist[[l]] + I*\[Eta]))*), {l, 1, Nband*Nx}]]], {w, 0, nw}]; 


 (*   ldossc = Table[Im[Total[Table[Abs[ES[[2,l,ival]]]^2*(1/(-wmax + wmax*2*(w/nw) - ES[[1,l]] + I*\[Eta]))+ 
      Abs[ES[[2,l,Nband*Nx+ival]]]^2*(1/(-wmax + 2*wmax*(w/nw) + ES[[1,l]] + I*\[Eta])), {l, 1, 2*Nband*Nx}]]], {w, 0, nw}];   *)

ldossc,{jy,0,Ny}]];]

  



Export["ldos_1_2_xstrip_sup_0.0.dat", Table[{(-wmax + wmax*2*((n-1)/nw))*1000,-ldos[[n,i]]/(Ny+1)},{n,1,nw+1},{i,1,Length[ival]}]];
Print["Done"]


(* ::Input::Initialization:: *)
fl = StringJoin["wmax",ToString[wmax],"Ny",ToString[Ny],"Nx",ToString[Nx],
"eta",ToString[\[Eta]]];
DeleteFile[StringJoin[fl,"Ldos_wte2_case_a.dat"]];
str1 = OpenAppend[StringJoin[fl,"Ldos_wte2_case_a.dat"], FormatType -> StandardForm];
For[n=1,n<=nw+1,n++,
For[i=1,i<=Length[ival],i++,
Write[str1,-wmax + wmax*2*((n-1)/nw),"  ", -1/(Pi*(Ny+1))*ldos[[n,i]]];
]
];
Close[str1];
Print["Done"]


Exit[];
