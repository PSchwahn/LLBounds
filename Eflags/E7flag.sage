#scan for E7T7
cartantype="E7"
print("G="+cartantype)
print("Initializing...",end="\r")
famG=cartantype[0]
rankG=int(cartantype[1:])
LiEG=lie(cartantype)
LiEH=lie("T"+cartantype[1:])
dimG=lie.dim(LiEG).sage()
G=WeylCharacterRing(cartantype,style="coroots")
f=open(str(LiEG)+str(LiEH)+".txt","w")
f.write("====== Sym2LowerBounds ("+version+") for G="+str(LiEG)+", H="+str(LiEH)+"======\n\n")
print("Computing isotropy representation...",end="\r")
LiEtrivial=lie.null(rankG)
LiEAdG=lie.adjoint(LiEG)
LiEAdBranched=LiEAdG.Demazure(LiEG)
branchdict={LiEAdG:LiEAdBranched}
LiEIsotropy=LiEAdBranched-lie(rankG)*lie.poly_one(rankG)
li=LiEIsotropy.length().sage()
print("LiE isotropy representation: "+str(LiEIsotropy))
f.write("LiE isotropy representation: "+str(LiEIsotropy)+"\n")
print("Computing symmetric powers: Sym^2_0m...",end="\r")
LiESym20m=lie.sym_tensor(2,LiEIsotropy,LiEH)-lie.poly_one(rankG)
lsm=LiESym20m.length().sage()
print("Computing symmetric powers: Sym^2_0g...",end="\r")
LiESym20g=lie.sym_tensor(2,LiEAdG,LiEG)-lie.poly_one(rankG)
lsg=LiESym20g.length().sage()
print("Computing Einstein constant...         ",end="\r")
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
print("                                       ",end="\r")
print("Einstein constant: "+str(E))
f.write("Einstein constant: E="+str(E)+"\n")
print("Branching isotypes of Sym^2_0g...      ",end="\r")
CasSym20g=[]
LiEGisotypesbranched=[]
for j in range(1,lsg+1):
    #hier gibts manchmal nen hash table overflow? warum??
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
for i in range(0,1000):
    LiEisotype=LiESym20m.expon(i+1)
    isotype=tuple(LiEisotype.sage())
    print("Scanning Sym20m, "+str(i+1)+" of "+str(lsm),end="\r")
    v=vector(LiEisotype.sage())
    casGofIsotype=[]
    for j in range(0,lsg):
        LiEGisotype=LiESym20g.expon(j+1)
        branched=LiEGisotypesbranched[j]
        if lie(branched._name+'|'+LiEisotype._name).sage()!=0:
            gisotype=SageRepfromList(G,LiEGisotype.sage())
            casGofIsotype+=[CasIrred(gisotype,[1],[G])]
    CasGMin=min(casGofIsotype)
    CasGMinDict[isotype]=CasGMin
    CasGMax=max(casGofIsotype)
    CasH=normconst*v*C*v
    CasHDict[isotype]=CasH
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
if PotInstab==[]:
    print("q(R)>E holds globally, implying stability.")
    f.write("q(R)>E holds globally, implying stability.\n")
else:
    print("q(R)>E is potentially violated on a subspace of dimension: "+str(len(PotInstab)))
    f.write("q(R)>E is potentially violated on a subspace of dimension: "+str(len(PotInstab))+"\n")
f.write("\nDone!")
f.close() 
