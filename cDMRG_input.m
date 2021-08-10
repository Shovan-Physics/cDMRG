(* ::Package:: *)

(* ::Input::Initialization:: *)
(*Monomials*)
allmonomialsfixedorder[n_,d_]:=Block[{part=IntegerPartitions[d+n,{n}]},Flatten[Permutations/@(part-1),1]];
allmonomials[n_,dmax_]:=Join@@Table[allmonomialsfixedorder[n,d],{d,0,dmax}];
monint=Compile[{{mon,_Real,1}},1./Times@@Accumulate[mon+1.],"RuntimeOptions"->"Speed","CompilationTarget"->"C"];
fmonint[mon_]:=fmonint[mon]=monint[mon];
makesym[lower_]:=MapThread[Join,{lower,Rest/@Flatten[lower,{{2},{1}}]}];
symOuter[f_,list_]:=makesym@Table[f[list[[i]],list[[j]]],{i,Length[list]},{j,i}];
(*Basis construction*)
legendreKE=
Compile[{{mon1,_Real,1},{mon2,_Real,1}},
Block[{sum=mon1+mon2,prod=mon1*mon2,res=0.},
Do[If[prod[[i]]!=0,res+=prod[[i]]*monint@MapAt[#-1&,sum,i]],{i,Length[prod]}];
res-Total[prod]*monint[sum]],
"CompilationTarget"->"C","RuntimeOptions"->"Speed"];
merge=Compile[{{mon,_Integer,1},{i,_Integer}},Join[Take[mon,i-1],{mon[[i]]+mon[[i+1]]},Drop[mon,i+1]]];
legendreIE=
Compile[{{sum,_Integer,1}},
Block[{merged,res=0.},
Do[
merged=merge[sum,i];
res+=monint@MapAt[#+1&,merged,i]-monint@MapAt[#+2&,merged,i];,
{i,Length[sum]-1}];
res],
"CompilationTarget"->"C","RuntimeOptions"->"Speed"];
flegendreIE[sum_]:=flegendreIE[sum]=legendreIE[sum];
legendrebasis[n_,dmax_,c_,keep_]:=
If[n==0,{{1.}},
Block[{monlist=allmonomials[n,dmax],ov,kin,int,akeep,eigenvecs,scale=n!},
ov=symOuter[fmonint[#1+#2]&,monlist];
kin=symOuter[legendreKE,monlist];
int=symOuter[flegendreIE[#1+#2]&,monlist];
akeep=Min[Length[monlist],keep];
eigenvecs=Reverse@Eigenvectors[{kin+c*int ,ov},-akeep];
#/Sqrt[scale*(# . ov . #)]&/@eigenvecs]];
(*Local operators*)
KE=
Compile[{{mon1,_Real,1},{mon2,_Real,1}},
Block[{sum=mon1+mon2,prod=mon1*mon2,res=0.},
Do[If[prod[[i]]!=0,res+=prod[[i]]*monint@MapAt[#-2&,sum,i]],{i,Length[prod]}];
res],
"CompilationTarget"->"C","RuntimeOptions"->"Speed"];
kinoverlaps[n_,dmax_]:=Block[{monlist=allmonomials[n,dmax]},(n!/2)*symOuter[KE,monlist]];
IE=
Compile[{{sum,_Integer,1}},
Block[{merged,res=0.},
Do[
merged=merge[sum,i];
res+=monint@merged;,
{i,Length[sum]-1}];
res],
"CompilationTarget"->"C","RuntimeOptions"->"Speed"];
fIE[sum_]:=fIE[sum]=IE[sum];
intoverlaps[n_,dmax_]:=Block[{monlist=allmonomials[n,dmax]},(n!/2)*symOuter[fIE[#1+#2]&,monlist]];
denscoefs=
Compile[{{sum,_Integer,1}},
Block[{n=Length[sum],powers,coefs,phases,leftint,rightint,middleint,partialsums={-1.},result,scale},
powers=Accumulate[Join[{sum[[1]]},Rest[sum]+1]];
result=ConstantArray[0.,1+powers[[-1]]];
phases=(-1.)^Range[n];
leftint=phases*Join[{1.},Table[monint@Take[sum,i],{i,n-1}]];
rightint=phases*Join[Table[monint@Drop[sum,j],{j,n-1}],{1.}];
Do[
AppendTo[partialsums,Sum[leftint[[i]]*monint[Take[sum,{j,i+1,-1}]],{i,j-1}]+leftint[[j]]];,
{j,2,n}];
coefs=partialsums*rightint;
scale=N@n!;
Do[result[[1+powers[[i]]]]=scale*coefs[[i]],{i,n}];
result],
"RuntimeOptions"->"Speed","CompilationTarget"->"C"];
fdenscoefs[sum_]:=fdenscoefs[sum]=denscoefs[sum];
densoverlaps[n_,dmax_]:=If[n==0,{{{}}},Block[{monlist=allmonomials[n,dmax]},symOuter[fdenscoefs[#1+#2]&,monlist]]];
poly[coefs_?VectorQ,x_]:=coefs . x^Range[0,Length@coefs-1];
poly[coefs_?VectorQ,0]:=If[coefs!={},coefs[[1]],0];
poly[table_/;MatrixQ[table,VectorQ],x_]:=Map[poly[#,x]&,table,{2}];
psicoefs=
Compile[{{smallmon,_Integer,1},{bigmon,_Integer,1}},
Block[{suml=smallmon+Most[bigmon],sumr=smallmon+Rest[bigmon],n=Length[bigmon],powers,coefs,phases,leftint,rightint,middleint,partialsums={-1.},result,scale},
powers=Accumulate[Join[{bigmon[[1]]},sumr+1]];
result=ConstantArray[0.,1+powers[[-1]]];
phases=(-1.)^Range[n];
leftint=phases*Join[{1.},Table[monint@Take[suml,i],{i,n-1}]];
rightint=phases*Join[Table[monint@Drop[sumr,j],{j,0,n-2}],{1.}];
Do[
AppendTo[partialsums,Sum[leftint[[i]]*monint[Take[sumr,{j,i,-1}]],{i,j}]+leftint[[j+1]]];,
{j,n-1}];
coefs=partialsums*rightint;
scale=Sqrt[N@n]*Gamma[N@n];
Do[result[[1+powers[[i]]]]=scale*coefs[[i]],{i,n}];
result],
"RuntimeOptions"->"Speed","CompilationTarget"->"C"];
psioverlaps[n_,dmaxSmall_,dmaxBig_]:=
Block[{smallmons=allmonomials[n-1,dmaxSmall],bigmons=allmonomials[n,dmaxBig]},Outer[psicoefs,smallmons,bigmons,1]];
incgamma=
Compile[{{p,_Integer},{z,_Complex}},
Block[{sum=1.+0.*I,residue,k=1.},
residue=sum;
While[Abs[residue/sum]>10^-15.,
residue*=z/(p+k);
sum+=residue;
k++];
Exp[-z]*z^p*sum/p],
"RuntimeOptions"->"Speed","CompilationTarget"->"C"];
PE=
Compile[{{ncoefs,_Real,1},{a,_Real},{b,_Real}},
Block[{len=Length[ncoefs],res=0.},
Do[
If[ncoefs[[p]]!=0,res+=ncoefs[[p]]*(1./p+Re[I^p*Exp[I*b]*incgamma[p,-I*a]]/a^p)],
{p,len}];
res],
"RuntimeOptions"->"Speed","CompilationTarget"->"C"];
symPE[matncoefs_,{a_,b_}]:=makesym@Table[PE[matncoefs[[i,j]],a,b],{i,Length[matncoefs]},{j,i}];
blockmatrix[blocks_]:=SparseArray[{Band[{1,1}]->blocks}];
sblockmatrix[blocks_]:=SparseArray[{Band[{1,2}]->blocks},1+Plus@@(Last[Dimensions[#]]&/@blocks)];
makebasis[basisparams_?MatrixQ,numsegments_Integer,g_?NumericQ,numwells_?NumericQ,V0_?NumericQ,maxpower_:Infinity]:=
Block[{dmaxlist,keeps,\[Delta]x,nmax,bases,txmat,kinblocks,kinmonomial,kin,intblocks,intmonomial,int,Hsitebox,psiblocks,psimaxlen,psiblockspadded,psicoefmonomial,densblocks,densmaxlen,densblockspadded,denscoefmonomial,boundaries,widths,arglist,pots},
{dmaxlist,keeps}=Transpose@Join[{{0,1}},basisparams];
\[Delta]x=1./numsegments;
nmax=Length[basisparams];
bases=legendrebasis[#-1,dmaxlist[[#]],g*\[Delta]x,keeps[[#]]]&/@Range[nmax+1];
txmat=blockmatrix[bases];
basisdims=Length/@bases;
localqn=Flatten@MapIndexed[ConstantArray[#2-1,#1]&,basisdims];
kinblocks=kinoverlaps[#-1,dmaxlist[[#]]]&/@Range[nmax+1];
kinmonomial=blockmatrix[kinblocks];
kin=(txmat/\[Delta]x^2) . kinmonomial . Transpose[txmat];
intblocks=intoverlaps[#-1,dmaxlist[[#]]]&/@Range[nmax+1];
intmonomial=blockmatrix[intblocks];
int=((g/\[Delta]x)*txmat) . intmonomial . Transpose[txmat];
Hsitebox=kin+int;
psiblocks=psioverlaps[#,dmaxlist[[#]],dmaxlist[[#+1]]]&/@Range[nmax];
psimaxlen=Max@Map[Length,psiblocks,{3}];
apsimaxlen=Min[psimaxlen,maxpower+1];
psiblockspadded=Map[PadRight[#,apsimaxlen]&,psiblocks,{3}];
psicoefmonomial=SparseArray[sblockmatrix/@Transpose[TensorTranspose[#,{2,3,1}]&/@psiblockspadded]];
psicoef=TensorTranspose[(txmat/Sqrt[\[Delta]x]) . Transpose[psicoefmonomial] . Transpose[txmat],{1,3,2}];
psiL=SparseArray@poly[psicoef,0];
psiR=SparseArray@poly[psicoef,1];
densblocks=densoverlaps[#-1,dmaxlist[[#]]]&/@Range[nmax+1];
densmaxlen=Max@Map[Length,densblocks,{3}];
adensmaxlen=Min[densmaxlen,maxpower+1];
densblockspadded=Map[PadRight[#,adensmaxlen]&,densblocks,{3}];
denscoefmonomial=SparseArray[blockmatrix/@Transpose[TensorTranspose[#,{2,3,1}]&/@densblockspadded]];
denscoef=TensorTranspose[(txmat/\[Delta]x) . Transpose[denscoefmonomial] . Transpose[txmat],{1,3,2}];
boundaries=N@Subdivide[numsegments];
widths=ConstantArray[\[Delta]x,numsegments];
arglist=(2*numwells*\[Pi])*Transpose[{widths,Most[boundaries]}];
pots=(V0*\[Delta]x/2)*SparseArray@symPE[Normal@denscoef,#]&/@arglist;
Do[Hsite[j]=Hsitebox+pots[[j]],{j,numsegments}];
];
(*Initial state*)
equispacedpos[particles_]:=If[particles==0,{},Ceiling[(Range@particles-1/2)*numsegments/particles]];
equispacedvec:=
Block[{filling,remainder,occupations,posrules},
{filling,remainder}=QuotientRemainder[particles,numsegments];
occupations=MapAt[#+1&,ConstantArray[filling,numsegments],Transpose[{equispacedpos[remainder]}]];
posrules=Join[{0->1},MapIndexed[#2[[1]]->#1&,Most[Accumulate[basisdims]+1]]];
occupations/.posrules];
(*Creating input file for DMRG*)
cfloat[q_?NumericQ]:=ToString[CForm[N@q]];
cformat[val_]:=
Block[{},
If[Head[N@val]==Real && Not[Head[val]===Integer],Return@cfloat[N@val]];
If[val==True,Return["true"]];
If[val==False,Return["false"]];
ToString[val]
];
outputParams[file_OutputStream,name_String,parameters_List]:=
Block[{},
WriteLine[file,name];
WriteLine[file,"{"];
writeParam[file,#]&/@parameters;
WriteLine[file,"}"];
WriteLine[file,""];
];
writeParam[file_OutputStream,prule_Rule]:=
Block[{param=prule[[1]],val=prule[[2]]},
WriteString[file,param<>"="];
WriteLine[file,cformat[val]];
];
outputVec[file_OutputStream,name_String,vec_List]:=
Block[{},
WriteLine[file,name];
WriteLine[file,"{"];
WriteString[file,cformat[#]," "]&/@vec;
WriteLine[file,""];
WriteLine[file,"}"];
WriteLine[file,""];
];
outputMat[file_OutputStream,name_String,rules_List]:=
Block[{numrules=Length[rules]},
WriteLine[file,name];
WriteLine[file,"{"];
WriteLine[file,ToString[numrules]];
writeRule[file,#]&/@rules;
WriteLine[file,"}"];
WriteLine[file,""];
];
writeRule[file_OutputStream,line_Rule]:=
Block[{dat=Flatten[List@@line]},
WriteString[file,cformat[#]," "]&/@dat;
WriteLine[file,""];
];
createInputFile[filename_String,paramset_List]:=
Block[{file,parameters="parameters"/.paramset,vecparams="vecparams"/.paramset,localmats="localmats"/.paramset},
file=OpenWrite[filename];
outputParams[file,"parameters",parameters];
outputVec[file,#1,#2]&@@@vecparams;
outputMat[file,#1,Most@ArrayRules@Transpose[#2]]&@@@localmats;
Close[file];
];
(*Command-line inputs*)
particles=ToExpression[$CommandLine[[5]]];
gamma=ToExpression[$CommandLine[[6]]];
numwells=ToExpression[$CommandLine[[7]]];
V0byEr=ToExpression[$CommandLine[[8]]];
numsegments=ToExpression[$CommandLine[[9]]];
basisid=ToExpression[$CommandLine[[10]]];
epsid=ToExpression[$CommandLine[[11]]];
saveid=ToExpression[$CommandLine[[12]]];
(*Basis, DMRG, and save parameters*)
basistable=Import["Parameters/basistable.m"];
If[MemberQ[Range@Length[basistable],basisid],basisparams=basistable[[basisid,2]],Print["Basis ID not found"];Exit[]];
Clear[basistable];
epstable=Import["Parameters/epstable.m"];
If[MemberQ[Range@Length[epstable],epsid],sweepparams=epstable[[epsid,2]],Print["Eps ID not found"];Exit[]];
targetdiscontinuity="targetdisc"/.sweepparams;epsparams=Rest[sweepparams];
Clear[epstable,sweepparams];
savetable=Import["Parameters/savetable.m"];
If[MemberQ[Range@Length[savetable],saveid],fullsaveparams=savetable[[saveid,2]],Print["Save ID not found"];Exit[]];
maxpower="maxpower"/.fullsaveparams;saveparams=Rest[fullsaveparams];
Clear[savetable,fullsaveparams];
(*Calculate basis and local operators*)
g=gamma*particles;
Er=(numwells*\[Pi])^2./2;
V0=V0byEr*Er;
makebasis[basisparams,numsegments,g,numwells,V0,maxpower];
Clear[fmonint,flegendreIE,fIE,fdenscoefs];
(*Set initial state*)
initialvec=equispacedvec;
loadfromfile=False;
(*Set file name*)
namestr[x_?StringQ,y_?NumericQ]:="_"<>x<>If[0<y<0.001||y>=10^4,"_1E"<>ToString@Rationalize[Log[10,y],0],"_"<>ToString[y]];
namestr[x_?StringQ,y_?StringQ]:="_"<>x<>"_"<>y;
namestr[list_?MatrixQ]:=StringDrop[StringJoin[namestr@@@list],1];
id={{"N",particles},{"gamma",gamma},{"Nwells",numwells},{"V0",V0byEr},{"M",numsegments},{"basisparams",basisid},{"sweepparams",epsid}};
name="Runs/"<>namestr[id]<>"_ITensor";
(*Write parameters to input file*)
paramlist=
{
"numsites"->numsegments,
"targetdiscontinuity"->targetdiscontinuity,
"numeps"->Length["eps"/.epsparams],
"numstates"->Total[basisdims],
"maxorder"->Max[apsimaxlen,adensmaxlen],
"loadfromfile"->loadfromfile
};
paramset=
{"parameters"->Join[paramlist,saveparams,If[loadfromfile,{"loadfromfilename"->loadfromfilename},{}]],
"vecparams"->Join[epsparams,{"qns"->localqn,"initialvec"->initialvec}],
"localmats"->Join[{"psiL"->psiL,"psiR"->psiR,"psi"->psicoef,"n"->denscoef},"H"<>ToString[#]->Hsite[#]&/@Range[numsegments]]};
createInputFile[name<>".txt",paramset];
(*Run DMRG*)
Run["./cDMRG "<>name<>".txt >"<>name<>".out 2>"<>name<>".error"];
(*Post processing*)
If[("savelocaldms"/.saveparams)&&("savefinalmeas"/.saveparams),
localdms=SparseArray@Import[name<>"_localdms.m"];
varPartition[l_,p_]:=MapThread[l[[#;;#2]]&,{{0}~Join~Most@#+1,#}&@Accumulate@p];
avgoccs=varPartition[Normal@Mean[Diagonal/@localdms],basisdims];
finalmeas=Join[Import[name<>"_finalmeas.m"],{"avgoccupations"->avgoccs}];
Export[name<>"_finalmeas.m",finalmeas];
Clear[localdms,avgoccs,finalmeas]
];
If["savewfMMA"/.saveparams,
psi=Import[name<>"_psi_MMAformat.m"];
DeleteFile[name<>"_psi_MMAformat.m"];
Export[name<>"_psi_MMAformat.dat",Compress@psi,"String"];
Clear[psi]
];
Exit[];
