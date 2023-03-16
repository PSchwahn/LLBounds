#last missing Fourier mode 2e2 of D6/T6

cartantype="D6"
print("G="+cartantype)
print("Initializing...")
famG=cartantype[0]
rankG=int(cartantype[1:])
LiEG=lie(cartantype)
LiEH=lie("T"+cartantype[1:])
dimG=lie.dim(LiEG).sage()
G=WeylCharacterRing(cartantype,style="coroots")
f=open(str(LiEG)+str(LiEH)+".txt","w")
f.write("====== Sym2LowerBounds ("+version+") for G="+str(LiEG)+", H="+str(LiEH)+"======\n\n")
print("Computing isotropy representation...")
LiEtrivial=lie.null(rankG)
LiEAdG=lie.adjoint(LiEG)
LiEAdBranched=LiEAdG.Demazure(LiEG)
branchdict={LiEAdG:LiEAdBranched}
LiEIsotropy=LiEAdBranched-lie(rankG)*lie.poly_one(rankG)
li=LiEIsotropy.length().sage()
print("LiE isotropy representation: "+str(LiEIsotropy))
f.write("LiE isotropy representation: "+str(LiEIsotropy)+"\n")
print("Computing symmetric powers: Sym^2_0m...")
LiESym20m=lie.sym_tensor(2,LiEIsotropy,LiEH)-lie.poly_one(rankG)
lsm=LiESym20m.length().sage()
print("Computing symmetric powers: Sym^2_0g...")
LiESym20g=lie.sym_tensor(2,LiEAdG,LiEG)-lie.poly_one(rankG)
lsg=LiESym20g.length().sage()
print("Computing Einstein constant...         ")
C=CartanMatrix([famG,rankG])^(-1)
CasIsotropyUnnorm=[]
CasIsotropy=rankG/(dimG-rankG)
for i in range(1,li+1):
    v=vector(LiEIsotropy.expon(i).sage())
    CasIsotropyUnnorm+=[v*C*v]
if CasIsotropyUnnorm!=[CasIsotropyUnnorm[0] for i in range(0,li)]:
    print(CasIsotropyUnnorm)
    raise ValueError("The Casimir on the isotropy representation is not constant.")
normconst=CasIsotropy/CasIsotropyUnnorm[0]
#now torus Casimir is normconst*v*C*v
E=(CasIsotropy+1/2)/2
print("                                       ")
print("Einstein constant: "+str(E))
f.write("Einstein constant: E="+str(E)+"\n")
print("Branching isotypes of Sym^2_0g...      ")
CasSym20g=[]
LiEGisotypesbranched=[]
for j in range(1,lsg+1):
    LiEGisotype=lie.expon(LiESym20g,j)
    Gisotype=SageRepfromList(G,LiEGisotype.sage())
    CasSym20g+=[CasIrred(Gisotype,[1],[G])]
    if LiEGisotype in branchdict:
        LiEGisotypesbranched+=[branchdict[LiEGisotype]]
    else:
        branched=LiEGisotype.Demazure(LiEG)
        branchdict[LiEGisotype]=branched
        LiEGisotypesbranched+=[branched]
#Now let's iterate over the isotypical components of Sym20m.
spacing=max([25,2*rankG+6])
#format_string1="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}"
#print(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min"))
PotInstab=[]
CasHDict={}
CasGMinDict={}
AAMinList=[]
AAMaxList=[]
for i in range(0,lsm):
    LiEisotype=LiESym20m.expon(i+1)
    print("Scanning Sym20m, "+str(i+1)+" of "+str(lsm)+"        ",end="\r")
    v=vector(LiEisotype.sage())
    casGofIsotype=[]
    for j in range(0,lsg):
        LiEGisotype=LiESym20g.expon(j+1)
        branched=LiEGisotypesbranched[j]
        if lie(branched._name+'|'+LiEisotype._name).sage()!=0:
            gisotype=SageRepfromList(G,LiEGisotype.sage())
            casGofIsotype+=[CasIrred(gisotype,[1],[G])]
    CasGMin=min(casGofIsotype)
    CasGMinDict[LiEisotype]=CasGMin
    CasGMax=max(casGofIsotype)
    CasH=normconst*v*C*v
    CasHDict[LiEisotype]=CasH
    #formula: A*A=prCasG-CasH-2CasIsotropy
    AAMin=CasGMin-CasH-2*CasIsotropy
    AAMinList+=[AAMin]
    AAMax=CasGMax-CasH-2*CasIsotropy
    AAMaxList+=[AAMax]
    #formula: 4q(R)=A*A+4CasH
    qRMin=1/4*AAMin+CasH
    #print(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
    #                             str(AAMin),str(AAMax),str(qRMin)))
    if qRMin<=E:
        PotInstab+=[LiEisotype]
print("PotInstab: "+str(PotInstab)+"     ")
#Single out those Fourier modes that have homs into PotInstab
format_string2="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}"
print(format_string2.format("Fourier mode","CasG","LL >=","LL >= 2E?"))
f.write(format_string2.format("Fourier mode","CasG","LL >=","LL >= 2E?")+"\n")
LiERepG=lie([0,2,0,0,0,0])
repstring="2e2"
CasG=11/5
print(repstring+": Branching...")
if LiERepG in branchdict:
    LiERepGbranched=branchdict[LiERepG]
else:
    LiERepGbranched=LiERepG.Demazure(LiEG)
    branchdict[LiERepG]=LiERepGbranched
FibrewiseTermMin=9/10
print("FibrewiseTermMin: "+str(FibrewiseTermMin))
#only consider relevant part of Sym^2_0g:
RelevantSym20g=lie("1X[0,0,0,1,0,0] +1X[0,2,0,0,0,0] +1X[2,0,0,0,0,0]")
print("RelevantSym20g: "+str(RelevantSym20g))
#only consider relevant part of repG*Sym^2_0g:
#to be precise it would be Hom(repG,Sym^2_0g), but Casimir stays the same
#this is without multiplicities
TensorProd=lie("1X[0,0,0,0,0,0] +2X[0,0,0,1,0,0] +1X[0,0,0,2,0,0] +2X[0,0,2,0,0,0] + 1X[0,1,0,0,0,0] +1X[0,1,0,0,0,2] +1X[0,1,0,0,2,0] +2X[0,1,0,1,0,0] + 4X[0,2,0,0,0,0] +2X[0,2,0,1,0,0] +1X[0,3,0,0,0,0] +1X[0,4,0,0,0,0] + 1X[1,0,0,0,1,1] +3X[1,0,1,0,0,0] +1X[1,0,1,1,0,0] +1X[1,1,0,0,1,1] + 4X[1,1,1,0,0,0] +1X[1,2,1,0,0,0] +2X[2,0,0,0,0,0] +2X[2,0,0,1,0,0] + 1X[2,0,2,0,0,0] +2X[2,1,0,0,0,0] +2X[2,2,0,0,0,0] +1X[3,0,1,0,0,0] + 1X[4,0,0,0,0,0]")
print("TensorProd: "+str(TensorProd))
print(repstring+": Branching terms of RepG*Sym20g...               ")
tplength=TensorProd.length().sage()
skipped=0
RelevantTensorProd=[]
for i in range(0,tplength):
    print(repstring+": Branching terms of RepG*Sym20g... ("\
          +str(i+1)+" of "+str(tplength)+")         ")
    TensorProdTerm=TensorProd.expon(i+1)
    print("TensorProdTerm: "+str(TensorProdTerm))
    if TensorProdTerm in branchdict:
        branched=branchdict[TensorProdTerm]
        skipped+=1
    else:
        branched=TensorProdTerm.Demazure(LiEG)
        branchdict[TensorProdTerm]=branched
    if lie(branched._name+'|'+LiEtrivial._name).sage()!=0:
        RelevantTensorProd+=[TensorProdTerm]
#MixedTerm:
print(repstring+": Computing Casimirs...                           ")
CasMixed=[CasKilling(SageRepfromList(G,x.sage())) for x in RelevantTensorProd]
MixedTermMax=max(CasMixed)
#Lower bound for LL:
#formula: LL=3/2*CasG-1/2*MixedTerm+FibrewiseTerm-CasIsotropy
#thus LL>=3/2*CasG-1/2*max(MixedTerm)+min(FibrewiseTerm)-CasIsotropy
LLlower=3/2*CasG-1/2*MixedTermMax+FibrewiseTermMin-CasIsotropy
goodbound=""
if LLlower>2*E: goodbound=">"
if LLlower==2*E: goodbound=">="
print(repstring+": Checking for Killing tensors...                 ")
print(format_string2.format(repstring,str(CasG),str(LLlower),goodbound)\
      +" (tplength: "+str(tplength)+", skipped: "+str(skipped)+")")
f.write(format_string2.format(repstring,str(CasG),str(LLlower),goodbound)+"\n")
print("\nDone!")
f.write("\nDone!")
f.close()
