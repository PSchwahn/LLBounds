def Casimir(rep):
    hweight=rep.highest_weight()
    #might need to be careful with rho() for groups of type A_n (GL/SL)
    delta=rep.parent().space().rho()
    return hweight.inner_product(hweight+2*delta)
def CasimirEff(rep,delta):
    hweight=rep.highest_weight()
    return hweight.inner_product(hweight+2*delta)
#scheint richtige Werte zu liefern. Immer mit adjoint normalisieren!
#careful with non-simple groups
def CasAdjoint(G):
    return Casimir(G.adjoint_representation())
def dimhom(rep1,rep2):
    #monom1=rep1.monomials()
    #coef1=rep1.coefficients()
    #return sum([coef1[i]*rep2.multiplicity(monom1[i]) for i in range(0,len(monom1))])
    return rep1.inner_product(rep2)
def CasKilling(rep):
    return Casimir(rep)/CasAdjoint(rep.parent())
def ScanCas(G,q):
    upperbound=q*CasAdjoint(G)
    rank=G.rank()
    #tuple is immutable/non-subscriptable, but list is not hashable..
    nullweight=tuple(0 for a in range(0,rank))
    out={nullweight : 0}
    for i in range(0,rank):
        tempout=out.copy()
        for rep in tempout:
            replist=list(rep)
            replist[i]+=1
            c=Casimir(SageRepfromList(G,replist))
            while c<=upperbound:
                out[tuple(replist)]=c
                replist[i]+=1
                c=Casimir(SageRepfromList(G,replist))
    return dict(sorted(out.items(), key=lambda item: item[1]))
def ScanCasKillingList(G,q):
    CasAd=CasAdjoint(G)
    rank=G.rank()
    #tuple is immutable/non-subscriptable, but list is not hashable..
    nullweight=[0 for a in range(0,rank)]
    out=[[nullweight,0]]
    for i in range(0,rank):
        tempout=list(out)
        for pair in tempout:
            rep=list(pair[0])
            rep[i]+=1
            c=Casimir(SageRepfromList(G,rep))/CasAd
            while c<=q:
                out+=[[list(rep),c]]
                rep[i]+=1
                c=Casimir(SageRepfromList(G,rep))/CasAd
    return sorted(out,key=lambda pair: pair[1])
def GroupFactors(H):
    cartanstr=str(H.cartan_type())
    if cartanstr.count("x")==0:
        return [H]
    else:
        return [WeylCharacterRing(s,style="coroots") for s in cartanstr.split("x")]
def reptolist(rep):
    return [int(x) for x in str(rep).split("(")[1].removesuffix(")").split(",")]
def SemisimpleAdjointRep(Hlist,H):
    ranklist=[Hi.rank() for Hi in Hlist]
    adjointlist=[Hi.adjoint_representation() for Hi in Hlist]
    out=0
    rankleft=0
    rankright=sum(ranklist)
    for i in range(0,len(Hlist)):
        rankright-=ranklist[i]
        out+=SageRepfromList(H,[0 for k in range(0,rankleft)]+reptolist(adjointlist[i])+[0 for k in range(0,rankright)])
        rankleft+=ranklist[i]
    return out
def RepFactors(rep,Hlist):
    ranklist=[Hi.rank() for Hi in Hlist]
    repindex=reptolist(rep)
    cumrank=0
    out=[]
    for i in range(0,len(Hlist)):
        out+=[Hlist[i](tuple(repindex[cumrank:cumrank+ranklist[i]]))]
        cumrank+=ranklist[i]
    return out
#Let H = H_1 x ... x H_k.
#Input: H, isotropy representation of H
#Output: normalization constants for Casimir operator
def CasNormConstants(isorep,Hlist):
    monomials=isorep.monomials()
    out=[]
    for i in range(0,len(Hlist)):
        s=0
        for monom in monomials:
            coef=isorep.monomial_coefficients()[monom.highest_weight()]
            s+=coef*monom.degree()*Casimir(RepFactors(monom,Hlist)[i])/CasAdjoint(Hlist[i])
        h=Hlist[i].adjoint_representation().degree()
        out+=[h/(h+s)]
    return out
def CasIrred(rep,normconstants,Hlist):
    factors=RepFactors(rep,Hlist)
    return sum([normconstants[i]*Casimir(factors[i])/CasAdjoint(Hlist[i]) for i in range(0,len(Hlist))])
def CasReducible(rep,normconstants,Hlist):
    factoredcomponents=[RepFactors(x,Hlist) for x in rep.monomials()]
    out=[]
    for j in range(0,len(factoredcomponents)):
        out+=[sum([normconstants[i]*Casimir(factoredcomponents[j][i])/CasAdjoint(Hlist[i])\
                   for i in range(0,len(Hlist))])]
    return out
def EinsteinConstant(IsotropyRep,normconstants,Hlist):
    cas=CasReducible(IsotropyRep,normconstants,Hlist)
    c1=cas[0]
    if cas==[c1 for i in range(0,len(cas))]:
        return (c1+1/2)/2
    else: raise ValueError("The Casimir operator has multiple eigenvalues on the isotropy representation. This normal homogeneous space is not Einstein.")
#LiE parser + branching in LiE:
def LiEgroupfromWCR(G):
    if G.space().root_system.is_irreducible():
        ct=G.cartan_type()
        ctstring=ct[0]+str(ct[1])
    else:
        ct=G.cartan_type()
        ctstring=""
        for x in ct.component_types():
            ctstring+=x[0]+str(x[1])
    return lie(ctstring)
#G simple, H semisimple, b branching rule
def resmatrix(G,H,b):
    rkG=G.rank()
    fw=G.fundamental_weights()
    HWR=WeightRing(H)
    return [[int(x) for x in str(HWR(list(b(fw[i])))).split("(")[1].removesuffix(")").split(",")]\
               for i in range(1,rkG+1)]
def fastLiEbranch(Grep,G,H,rm):
    LiEG=LiEgroupfromWCR(G)
    LiEH=LiEgroupfromWCR(H)
    LiErm=lie(str(rm))
    LiEGRep=lie(str(reptolist(Grep)))
    return lie.branch(LiEGRep,LiEH,LiErm,LiEG)
#reverse LiE parser:
def LiEreptomonomialslist(LiERep):
    return [[int(y) for y in x.split("]")[0].split("[")[1].split(",")] for x in str(LiERep).split("+")]
def LiErepHasInvariantPart(LiERep,rank):
    return str(LiERep).find(",".join(["0" for i in range(0,rank)]))!=-1
def omega(rep):
    repstring=""
    for a in range(0,len(rep)):
        if rep[a]!=0:
            repstring+=str(rep[a])+"e_"+str(a+1)+"+"
    repstring=repstring.removesuffix("+")
    if repstring=="": repstring="0"
    return repstring
#Sage is weird with A1 representations
def SageRepfromList(G,rep):
    if str(G.cartan_type()).count("x")==0:
        if G.cartan_type()[1]==1:
            return G(*rep,0)
    return G(*rep)
#PROBLEM: lie interface cannot read in strings/data that is too long.
#-> data must be generated in LiE, not fed into LiE.
#GOAL: Large data generated in LiE itself.
def ScanCasKillingLiE(G,q):
    #output: [[LiErep,CasG],...]
    CasAd=CasAdjoint(G)
    rank=G.rank()
    nullweight=lie.null(rank)
    out=[[nullweight,0]]
    for i in range(1,rank+1):
        tempout=list(out)
        for pair in tempout:
            LiErep=lie(pair[0]._name)
            lie.eval(LiErep._name+"["+str(i)+"]+=1")
            c=Casimir(SageRepfromList(G,LiErep.sage()))/CasAd
            while c<=q:
                out+=[[lie(LiErep._name),c]]
                lie.eval(LiErep._name+"["+str(i)+"]+=1")
                c=Casimir(SageRepfromList(G,LiErep.sage()))/CasAd
    return sorted(out,key=lambda pair: pair[1])
#here, restriction matrix is not directly read in, but via a .lie file.
#Problem: Sometimes LiE doesn't want to read files. Why?
def Sym2LowerBoundsRead(G,H,b,startscan=0,endscan=-1):
    print("Initializing groups...",end="\r")
    rankG=G.rank()
    rankH=H.rank()
    HFactors=GroupFactors(H)
    LiEG=LiEgroupfromWCR(G)
    LiEH=LiEgroupfromWCR(H)
    LiEHtrivial=lie.null(rankH)
    if type(b) is list:
        rm=b
    else:
        print("Computing restriction matrix...",end="\r")
        rm=resmatrix(G,H,b)
    f=open(str(LiEG)+str(LiEH)+".txt","w")
    f.write("====== Sym2LowerBoundsRead ("+version+") for G="+str(LiEG)+", H="+str(LiEH)+"======\n\n")
    print("Computing isotropy representation...",end="\r")
    LiEAdG=lie.adjoint(LiEG)
    LiEAdH=lie.adjoint(LiEH)
    #we don't want the interface to read in a possibly huge matrix, so we go via LiE's read command:
    g=open("temprm.lie","w")
    g.write("#G="+str(LiEG)+"\n")
    g.write("#H="+str(LiEH)+"\n")
    g.write("rm="+str(rm)+";")
    g.close()
    read=lie.eval("read temprm.lie")
    if read!="": raise ValueError(read)
    LiErm=lie("rm")
    LiEAdBranched=lie.branch(LiEAdG,LiEH,LiErm,LiEG)
    branchdict={LiEAdG:LiEAdBranched}
    LiEIsotropy=LiEAdBranched-LiEAdH
    IsotropyRep=0
    for i in range(0,LiEIsotropy.length().sage()):
        IsotropyRep+=LiEIsotropy.coef(i+1).sage()*SageRepfromList(H,LiEIsotropy.expon(i+1).sage())
    print("LiE isotropy representation: "+str(LiEIsotropy))
    f.write("LiE isotropy representation: "+str(LiEIsotropy)+"\n")
    print("Computing Einstein constant...         ",end="\r")
    normconstants=CasNormConstants(IsotropyRep,HFactors)
    E=EinsteinConstant(IsotropyRep,normconstants,HFactors)
    CasIsotropy=2*E-1/2
    print("                                ",end="\r")
    print("Einstein constant: E="+str(E))
    f.write("Einstein constant: E="+str(E)+"\n")
    print("Computing symmetric powers: Sym^2_0g...",end="\r")
    LiESym20g=lie.sym_tensor(2,LiEAdG,LiEG)-lie.poly_one(rankG)
    lsg=LiESym20g.length().sage()
    print("Computing symmetric powers: Sym^2_0m...",end="\r")
    LiESym20m=lie.sym_tensor(2,LiEIsotropy,LiEH)-lie.poly_one(rankH)
    lsm=LiESym20m.length().sage()
    print("Computing symmetric powers: Sym^3m...  ",end="\r")
    LiESym3m=lie.sym_tensor(3,LiEIsotropy,LiEH)
    print("Branching isotypes of Sym^2_0g...      ",end="\r")
    CasSym20g=[]
    LiEGisotypesbranched=[]
    for j in range(1,lsg+1):
        LiEGisotype=lie.expon(LiESym20g,j)
        Gisotype=SageRepfromList(G,LiEGisotype.sage())
        CasSym20g+=[CasIrred(Gisotype,[1],[G])]
        if LiEGisotype in branchdict:
            LiEGisotypesbranched+=[branchdict[LiEGisotype]]
        else:
            branched=lie.branch(LiEGisotype,LiEH,LiErm,LiEG)
            branchdict[LiEGisotype]=branched
            LiEGisotypesbranched+=[branched]
    #Now let's iterate over the isotypical components of Sym20m.
    spacing=max([25,2*rankH+6])
    format_string1="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}"
    print(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min"))
    f.write(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min")+"\n")
    PotInstab=[]
    CasHDict={}
    CasGMinDict={}
    AAMinList=[]
    AAMaxList=[]
    for i in range(0,lsm):
        LiEisotype=LiESym20m.expon(i+1)
        isotype=tuple(LiEisotype.sage())
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
        CasH=CasIrred(SageRepfromList(H,LiEisotype.sage()),normconstants,HFactors)
        CasHDict[isotype]=CasH
        #formula: A*A=prCasG-CasH-2CasIsotropy
        AAMin=CasGMin-CasH-2*CasIsotropy
        AAMinList+=[AAMin]
        AAMax=CasGMax-CasH-2*CasIsotropy
        AAMaxList+=[AAMax]
        #formula: 4q(R)=A*A+4CasH
        qRMin=1/4*AAMin+CasH
        print(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
                                     str(AAMin),str(AAMax),str(qRMin)))
        f.write(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
                                     str(AAMin),str(AAMax),str(qRMin))+"\n")
        if qRMin<=E:
            PotInstab+=[LiEisotype]
    if PotInstab==[]:
        print("q(R)>E holds globally, implying stability.")
        f.write("q(R)>E holds globally, implying stability.\n")
        return 0
    else:
        print("q(R)>E is potentially violated on: "+str(PotInstab))
        f.write("q(R)>E is potentially violated on: "+str(PotInstab)+"\n")
    #Now: give C such that LL>2E if CasG>C.
    if endscan==-1:
        c=var("c")
        eqn=(c+min(AAMinList)/2-2*E)^2==max(AAMaxList)*(c-min(CasHDict.values()))
        C=max([x.right() for x in solve(eqn,c)])
        print("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5)))
        f.write("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5))+"\n")
        print("Scanning for "+str(startscan)+" <= CasG <= C.")
        f.write("Scanning for "+str(startscan)+" <= CasG <= C.\n")
        LiEscancas=[x for x in ScanCasKillingLiE(G,C) if x[1]>=startscan]
    else:
        print("Scanning for "+str(startscan)+" <= CasG <= "+str(endscan)+".")
        f.write("Scanning for "+str(startscan)+" <= CasG <= "+str(endscan)+".\n")
        LiEscancas=[x for x in ScanCasKillingLiE(G,endscan) if x[1]>=startscan]
    if LiEscancas==[]:
        print("Nothing to scan.")
        f.write("Nothing to scan.\n")
        f.close()
        return 0
    #Single out those Fourier modes that have homs into PotInstab
    format_string2="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}"
    print(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?"))
    f.write(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?")+"\n")
    for FourierMode in LiEscancas:
        LiERepG=FourierMode[0]
        repstring=omega(LiERepG.sage())
        CasG=FourierMode[1]
        print(repstring+": Branching...",end="\r")
        if LiERepG in branchdict:
            LiERepGbranched=branchdict[LiERepG]
        else:
            LiERepGbranched=lie.branch(LiERepG,LiEH,LiErm,LiEG)
            branchdict[LiERepG]=LiERepGbranched
        #we only need those modes where we don't know q(R)>E:
        print(repstring+": Checking whether it has homs to PotInstab...",end="\r")
        if any(lie(LiERepGbranched._name+'|'+x._name).sage()!=0 for x in PotInstab):
            #first, consider relevant part of Sym^2_0m:
            print(repstring+": Tensoring with Sym20m...                        ",end="\r")
            RelevantSym20m=[]
            FibrewiseTermList=[]
            for i in range(0,lsm):
                LiEisotype=LiESym20m.expon(i+1)
                SageRepfromList(H,LiEisotype.sage())
                branchedhom=lie.tensor(LiERepGbranched,lie.X(LiEisotype),LiEH)
                if lie(branchedhom._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantSym20m+=[LiEisotype]
                    FibrewiseTermList+=[CasGMinDict[isotype]-3/2*CasHDict[isotype]]
            FibrewiseTermMin=min(FibrewiseTermList)
            #only consider relevant part of Sym^2_0g:
            RelevantSym20g=lie.poly_null(rankG)
            for j in range(0,lsg):
                LiEGisotype=LiESym20g.expon(j+1)
                branched=LiEGisotypesbranched[j]
                if any(lie(branched._name+'|'+x._name).sage()!=0 for x in RelevantSym20m):
                    RelevantSym20g+=lie.X(LiEGisotype)
            #only consider relevant part of repG*Sym^2_0g:
            #to be precise it would be Hom(repG,Sym^2_0g), but Casimir stays the same
            print(repstring+": Tensoring with Sym20g...                        ",end="\r")
            #this is without multiplicities
            TensorProd=lie.tensor(lie.X(LiERepG),RelevantSym20g,LiEG)
            print(repstring+": Branching terms of RepG*Sym20g...               ",end="\r")
            tplength=TensorProd.length().sage()
            skipped=0
            RelevantTensorProd=[]
            for i in range(0,tplength):
                print(repstring+": Branching terms of RepG*Sym20g... ("\
                      +str(i+1)+" of "+str(tplength)+")         ",end="\r")
                TensorProdTerm=TensorProd.expon(i+1)
                if TensorProdTerm in branchdict:
                    branched=branchdict[TensorProdTerm]
                    skipped+=1
                else:
                    branched=lie.branch(TensorProdTerm,LiEH,LiErm,LiEG)
                    branchdict[TensorProdTerm]=branched
                if lie(branched._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantTensorProd+=[TensorProdTerm]
            #MixedTerm:
            print(repstring+": Computing Casimirs...                           ",end="\r")
            CasMixed=[CasIrred(SageRepfromList(G,x.sage()),[1],[G]) for x in RelevantTensorProd]
            MixedTermMax=max(CasMixed)
            #Lower bound for LL:
            #formula: LL=3/2*CasG-1/2*MixedTerm+FibrewiseTerm-CasIsotropy
            #thus LL>=3/2*CasG-1/2*max(MixedTerm)+min(FibrewiseTerm)-CasIsotropy
            LLlower=3/2*CasG-1/2*MixedTermMax+FibrewiseTermMin-CasIsotropy
            goodbound=""
            if LLlower>2*E: goodbound=">"
            if LLlower==2*E: goodbound=">="
            print(repstring+": Checking for Killing tensors...                 ",end="\r")
            #check for Killing tensors
            killing=""
            dimHomSym3=0
            dimHomSym20=0
            for i in range(0,LiERepGbranched.length().sage()):
                summand=LiERepGbranched.expon(i+1)
                coef=LiERepGbranched.coef(i+1).sage()
                dimHomSym20+=coef*lie(LiESym20m._name+'|'+summand._name).sage()
                dimHomSym3+=coef*lie(LiESym3m._name+'|'+summand._name).sage()
            if dimHomSym3<dimHomSym20: killing="!"
            print(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound)\
                  +" (tplength: "+str(tplength)+", skipped: "+str(skipped)+")")
            f.write(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound)+"\n")
        else: print((spacing+50)*" ",end="\r")
    print("\nDone!")
    f.write("\nDone!")
    f.close()
    return 0
def Sym2LowerBounds(G,H,b,startscan=0,endscan=-1):
    rankG=G.rank()
    rankH=H.rank()
    HFactors=GroupFactors(H)
    LiEG=LiEgroupfromWCR(G)
    LiEH=LiEgroupfromWCR(H)
    print("G="+str(LiEG)+", H="+str(LiEH))
    LiEHtrivial=lie.null(rankH)
    if type(b) is list:
        #input restriction matrix as list of lists
        rm=b
    else:
        print("Computing restriction matrix...",end="\r")
        rm=resmatrix(G,H,b)
    f=open(str(LiEG)+str(LiEH)+".txt","w")
    f.write("====== Sym2LowerBounds ("+version+") for G="+str(LiEG)+", H="+str(LiEH)+"======\n\n")
    print("Computing isotropy representation...",end="\r")
    LiEAdG=lie.adjoint(LiEG)
    LiEAdH=lie.adjoint(LiEH)
    #this could cause problems
    LiErm=lie(rm)
    LiEAdBranched=lie.branch(LiEAdG,LiEH,LiErm,LiEG)
    branchdict={LiEAdG:LiEAdBranched}
    LiEIsotropy=LiEAdBranched-LiEAdH
    IsotropyRep=0
    for i in range(0,LiEIsotropy.length().sage()):
        IsotropyRep+=LiEIsotropy.coef(i+1).sage()*SageRepfromList(H,LiEIsotropy.expon(i+1).sage())
    print("LiE isotropy representation: "+str(LiEIsotropy))
    f.write("LiE isotropy representation: "+str(LiEIsotropy)+"\n")
    print("Computing Einstein constant...         ",end="\r")
    normconstants=CasNormConstants(IsotropyRep,HFactors)
    E=EinsteinConstant(IsotropyRep,normconstants,HFactors)
    CasIsotropy=2*E-1/2
    print("                                ",end="\r")
    print("Einstein constant: E="+str(E))
    f.write("Einstein constant: E="+str(E)+"\n")
    print("Computing symmetric powers: Sym^2_0g...",end="\r")
    LiESym20g=lie.sym_tensor(2,LiEAdG,LiEG)-lie.poly_one(rankG)
    lsg=LiESym20g.length().sage()
    print("Computing symmetric powers: Sym^2_0m...",end="\r")
    LiESym20m=lie.sym_tensor(2,LiEIsotropy,LiEH)-lie.poly_one(rankH)
    lsm=LiESym20m.length().sage()
    print("Computing symmetric powers: Sym^3m...  ",end="\r")
    LiESym3m=lie.sym_tensor(3,LiEIsotropy,LiEH)
    print("Branching isotypes of Sym^2_0g...      ",end="\r")
    CasSym20g=[]
    LiEGisotypesbranched=[]
    for j in range(1,lsg+1):
        LiEGisotype=lie.expon(LiESym20g,j)
        Gisotype=SageRepfromList(G,LiEGisotype.sage())
        CasSym20g+=[CasIrred(Gisotype,[1],[G])]
        if LiEGisotype in branchdict:
            LiEGisotypesbranched+=[branchdict[LiEGisotype]]
        else:
            branched=lie.branch(LiEGisotype,LiEH,LiErm,LiEG)
            branchdict[LiEGisotype]=branched
            LiEGisotypesbranched+=[branched]
    #Now let's iterate over the isotypical components of Sym20m.
    spacing=max([25,2*rankH+6])
    format_string1="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}"
    print(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min"))
    f.write(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min")+"\n")
    PotInstab=[]
    CasHDict={}
    CasGMinDict={}
    AAMinList=[]
    AAMaxList=[]
    for i in range(0,lsm):
        LiEisotype=LiESym20m.expon(i+1)
        isotype=tuple(LiEisotype.sage())
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
        CasH=CasIrred(SageRepfromList(H,LiEisotype.sage()),normconstants,HFactors)
        CasHDict[isotype]=CasH
        #formula: A*A=prCasG-CasH-2CasIsotropy
        AAMin=CasGMin-CasH-2*CasIsotropy
        AAMinList+=[AAMin]
        AAMax=CasGMax-CasH-2*CasIsotropy
        AAMaxList+=[AAMax]
        #formula: 4q(R)=A*A+4CasH
        qRMin=1/4*AAMin+CasH
        print(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
                                     str(AAMin),str(AAMax),str(qRMin)))
        f.write(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
                                     str(AAMin),str(AAMax),str(qRMin))+"\n")
        if qRMin<=E:
            PotInstab+=[LiEisotype]
    if PotInstab==[]:
        print("q(R)>E holds globally, implying stability.")
        f.write("q(R)>E holds globally, implying stability.\n")
        return 0
    else:
        print("q(R)>E is potentially violated on: "+str(PotInstab))
        f.write("q(R)>E is potentially violated on: "+str(PotInstab)+"\n")
    #Now: give C such that LL>2E if CasG>C.
    if endscan==-1:
        c=var("c")
        eqn=(c+min(AAMinList)/2-2*E)^2==max(AAMaxList)*(c-min(CasHDict.values()))
        C=max([x.right() for x in solve(eqn,c)])
        print("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5)))
        f.write("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5))+"\n")
        print("Scanning for "+str(startscan)+" <= CasG <= C.")
        f.write("Scanning for "+str(startscan)+" <= CasG <= C.\n")
        LiEscancas=[x for x in ScanCasKillingLiE(G,C) if x[1]>=startscan]
    else:
        print("Scanning for "+str(startscan)+" <= CasG <= "+str(endscan)+".")
        f.write("Scanning for "+str(startscan)+" <= CasG <= "+str(endscan)+".\n")
        LiEscancas=[x for x in ScanCasKillingLiE(G,endscan) if x[1]>=startscan]
    if LiEscancas==[]:
        print("Nothing to scan.")
        f.write("Nothing to scan.\n")
        f.close()
        return 0
    #Single out those Fourier modes that have homs into PotInstab
    format_string2="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}"
    print(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?"))
    f.write(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?")+"\n")
    for FourierMode in LiEscancas:
        LiERepG=FourierMode[0]
        repstring=omega(LiERepG.sage())
        CasG=FourierMode[1]
        print(repstring+": Branching...",end="\r")
        if LiERepG in branchdict:
            LiERepGbranched=branchdict[LiERepG]
        else:
            LiERepGbranched=lie.branch(LiERepG,LiEH,LiErm,LiEG)
            branchdict[LiERepG]=LiERepGbranched
        #we only need those modes where we don't know q(R)>E:
        print(repstring+": Checking whether it has homs to PotInstab...",end="\r")
        if any(lie(LiERepGbranched._name+'|'+x._name).sage()!=0 for x in PotInstab):
            #first, consider relevant part of Sym^2_0m:
            print(repstring+": Tensoring with Sym20m...                        ",end="\r")
            RelevantSym20m=[]
            FibrewiseTermList=[]
            for i in range(0,lsm):
                LiEisotype=LiESym20m.expon(i+1)
                isotype=tuple(LiEisotype.sage())
                branchedhom=lie.tensor(LiERepGbranched,lie.X(LiEisotype),LiEH)
                if lie(branchedhom._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantSym20m+=[LiEisotype]
                    FibrewiseTermList+=[CasGMinDict[isotype]-3/2*CasHDict[isotype]]
            FibrewiseTermMin=min(FibrewiseTermList)
            #only consider relevant part of Sym^2_0g:
            RelevantSym20g=lie.poly_null(rankG)
            for j in range(0,lsg):
                LiEGisotype=LiESym20g.expon(j+1)
                branched=LiEGisotypesbranched[j]
                if any(lie(branched._name+'|'+x._name).sage()!=0 for x in RelevantSym20m):
                    RelevantSym20g+=lie.X(LiEGisotype)
            #only consider relevant part of repG*Sym^2_0g:
            #to be precise it would be Hom(repG,Sym^2_0g), but Casimir stays the same
            print(repstring+": Tensoring with Sym20g...                        ",end="\r")
            #this is without multiplicities
            TensorProd=lie.tensor(lie.X(LiERepG),RelevantSym20g,LiEG)
            print(repstring+": Branching terms of RepG*Sym20g...               ",end="\r")
            tplength=TensorProd.length().sage()
            skipped=0
            RelevantTensorProd=[]
            for i in range(0,tplength):
                print(repstring+": Branching terms of RepG*Sym20g... ("\
                      +str(i+1)+" of "+str(tplength)+")         ",end="\r")
                TensorProdTerm=TensorProd.expon(i+1)
                if TensorProdTerm in branchdict:
                    branched=branchdict[TensorProdTerm]
                    skipped+=1
                else:
                    branched=lie.branch(TensorProdTerm,LiEH,LiErm,LiEG)
                    branchdict[TensorProdTerm]=branched
                if lie(branched._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantTensorProd+=[TensorProdTerm]
            #MixedTerm:
            print(repstring+": Computing Casimirs...                           ",end="\r")
            CasMixed=[CasIrred(SageRepfromList(G,x.sage()),[1],[G]) for x in RelevantTensorProd]
            MixedTermMax=max(CasMixed)
            #Lower bound for LL:
            #formula: LL=3/2*CasG-1/2*MixedTerm+FibrewiseTerm-CasIsotropy
            #thus LL>=3/2*CasG-1/2*max(MixedTerm)+min(FibrewiseTerm)-CasIsotropy
            LLlower=3/2*CasG-1/2*MixedTermMax+FibrewiseTermMin-CasIsotropy
            goodbound=""
            if LLlower>2*E: goodbound=">"
            if LLlower==2*E: goodbound=">="
            print(repstring+": Checking for Killing tensors...                 ",end="\r")
            #check for Killing tensors
            killing=""
            dimHomSym3=0
            dimHomSym20=0
            for i in range(0,LiERepGbranched.length().sage()):
                summand=LiERepGbranched.expon(i+1)
                coef=LiERepGbranched.coef(i+1).sage()
                dimHomSym20+=coef*lie(LiESym20m._name+'|'+summand._name).sage()
                dimHomSym3+=coef*lie(LiESym3m._name+'|'+summand._name).sage()
            if dimHomSym3<dimHomSym20: killing="!"
            print(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound)\
                  +" (tplength: "+str(tplength)+", skipped: "+str(skipped)+")")
            f.write(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound)+"\n")
        else: print((spacing+50)*" ",end="\r")
    print("\nDone!")
    f.write("\nDone!")
    f.close()
    return 0
def NormConstantsSS(LiEIsotropy,HssFactors,LiEH):
    ranklist=[Hi.rank() for Hi in HssFactors]
    cumrank=0
    out=[]
    li=LiEIsotropy.length().sage()
    for i in range(0,len(HssFactors)):
        s=0
        for j in range(1,li+1):
            monom=LiEIsotropy.expon(j).sage()
            monomssi=monom[cumrank:cumrank+ranklist[i]]
            dim=LiEIsotropy.coef(j).sage()*LiEIsotropy.expon(j).dim(LiEH).sage()
            s+=dim*Casimir(SageRepfromList(HssFactors[i],monomssi))/CasAdjoint(HssFactors[i])
        h=HssFactors[i].adjoint_representation().degree()
        out+=[h/(h+s)]
        cumrank+=ranklist[i]
    return out
#Now for non-semisimple H. Hss is the semisimple part.
#be careful to use the right restriction matrix!
def Sym2LowerBoundsWithTorus(G,Hss,rm,startscan=0,endscan=-1):
    rankG=G.rank()
    rankH=len(rm[0])
    rankHss=Hss.rank()
    rankT=rankH-rankHss
    HssFactors=GroupFactors(Hss)
    LiEG=LiEgroupfromWCR(G)
    LiEHss=LiEgroupfromWCR(Hss)
    LiEH=LiEHss*lie("T"+str(rankT))
    print("G="+str(LiEG)+", H="+str(LiEH))
    f=open(str(LiEG)+str(LiEH)+".txt","w")
    f.write("====== Sym2LowerBounds ("+version+") for G="+str(LiEG)+", H="+str(LiEH)+"======\n\n")
    LiEHtrivial=lie.null(rankH)
    LiErm=lie(rm)
    LiEAdG=lie.adjoint(LiEG)
    LiEAdH=lie.adjoint(LiEH)
    LiEAdBranched=lie.branch(LiEAdG,LiEH,LiErm,LiEG)
    branchdict={LiEAdG:LiEAdBranched}
    print("Computing isotropy representation...",end="\r")
    LiEIsotropy=LiEAdBranched-LiEAdH
    li=LiEIsotropy.length().sage()
    print("LiE isotropy representation: "+str(LiEIsotropy))
    f.write("LiE isotropy representation: "+str(LiEIsotropy)+"\n")
    print("Computing symmetric powers: Sym^2_0m...",end="\r")
    LiESym20m=lie.sym_tensor(2,LiEIsotropy,LiEH)-lie.poly_one(rankH)
    lsm=LiESym20m.length().sage()
    print("Computing symmetric powers: Sym^2_0g...",end="\r")
    LiESym20g=lie.sym_tensor(2,LiEAdG,LiEG)-lie.poly_one(rankG)
    lsg=LiESym20g.length().sage()
    print("Computing symmetric powers: Sym^3m...  ",end="\r")
    LiESym3m=lie.sym_tensor(3,LiEIsotropy,LiEH)
    print("Computing Einstein constant...         ",end="\r")
    normconstants=NormConstantsSS(LiEIsotropy,HssFactors,LiEH)
    cm=CartanMatrix(G.cartan_type())
    r=Matrix(rm)
    torusinprod=((transpose(r)*cm*r)[rankHss:rankH,rankHss:rankH])^(-1)
    CasHssList=[]
    CasTUnnormList=[]
    CasTUnnormTraceList=[]
    for i in range (1,li+1):
        LiEisotype=LiEIsotropy.expon(i)
        isotypess=SageRepfromList(Hss,LiEisotype.sage()[0:rankHss])
        isotypet=vector(LiEisotype.sage()[rankHss:rankH])
        dim=lie.dim(LiEisotype,LiEH).sage()
        CasHss=CasIrred(isotypess,normconstants,HssFactors)
        CasTUnnorm=isotypet*torusinprod*isotypet
        CasHssList+=[CasHss]
        CasTUnnormList+=[CasTUnnorm]
        CasTUnnormTraceList+=[dim*CasTUnnorm]
    tconst=rankT/sum(CasTUnnormTraceList)
    #from now: CasT=tconst*isotypet*torusinprod*isotypet
    CasIsotropyList=[CasHssList[i]+tconst*CasTUnnormList[i] for i in range(0,li)]
    CasIsotropy=CasIsotropyList[0]
    if CasIsotropyList!=[CasIsotropy for i in range(0,li)]:
        raise ValueError("The Casimir on the isotropy representation is not constant.")
    E=(CasIsotropy+1/2)/2
    print("                                       ",end="\r")
    print("Einstein constant: "+str(E))
    f.write("Einstein constant: E="+str(E)+"\n")
    print("Branching isotypes of Sym^2_0g...      ",end="\r")
    CasSym20g=[]
    LiEGisotypesbranched=[]
    for j in range(1,lsg+1):
        LiEGisotype=lie.expon(LiESym20g,j)
        Gisotype=SageRepfromList(G,LiEGisotype.sage())
        CasSym20g+=[CasIrred(Gisotype,[1],[G])]
        if LiEGisotype in branchdict:
            LiEGisotypesbranched+=[branchdict[LiEGisotype]]
        else:
            branched=lie.branch(LiEGisotype,LiEH,LiErm,LiEG)
            branchdict[LiEGisotype]=branched
            LiEGisotypesbranched+=[branched]
    #Now let's iterate over the isotypical components of Sym20m.
    spacing=max([25,2*rankH+6])
    format_string1="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}"
    print(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min"))
    f.write(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min")+"\n")
    PotInstab=[]
    CasHDict={}
    CasGMinDict={}
    AAMinList=[]
    AAMaxList=[]
    for i in range(0,lsm):
        LiEisotype=LiESym20m.expon(i+1)
        isotype=tuple(LiEisotype.sage())
        isotypess=SageRepfromList(Hss,LiEisotype.sage()[0:rankHss])
        isotypet=vector(LiEisotype.sage()[rankHss:rankH])
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
        CasHss=CasIrred(isotypess,normconstants,HssFactors)
        CasT=tconst*(isotypet*torusinprod*isotypet)
        CasH=CasHss+CasT
        CasHDict[isotype]=CasH
        #formula: A*A=prCasG-CasH-2CasIsotropy
        AAMin=CasGMin-CasH-2*CasIsotropy
        AAMinList+=[AAMin]
        AAMax=CasGMax-CasH-2*CasIsotropy
        AAMaxList+=[AAMax]
        #formula: 4q(R)=A*A+4CasH
        qRMin=1/4*AAMin+CasH
        print(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
                                     str(AAMin),str(AAMax),str(qRMin)))
        f.write(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
                                     str(AAMin),str(AAMax),str(qRMin))+"\n")
        if qRMin<=E:
            PotInstab+=[LiEisotype]
    if PotInstab==[]:
        print("q(R)>E holds globally, implying stability.")
        f.write("q(R)>E holds globally, implying stability.\n")
        return 0
    else:
        print("q(R)>E is potentially violated on: "+str(PotInstab))
        f.write("q(R)>E is potentially violated on: "+str(PotInstab)+"\n")
    #Now: give C such that LL>2E if CasG>C.
    if endscan==-1:
        c=var("c")
        eqn=(c+min(AAMinList)/2-2*E)^2==max(AAMaxList)*(c-min(CasHDict.values()))
        C=max([x.right() for x in solve(eqn,c)])
        print("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5)))
        f.write("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5))+"\n")
        print("Scanning for "+str(startscan)+" <= CasG <= C.")
        f.write("Scanning for "+str(startscan)+" <= CasG <= C.\n")
        LiEscancas=[x for x in ScanCasKillingLiE(G,C) if x[1]>=startscan]
    else:
        print("Scanning for "+str(startscan)+" <= CasG <= "+str(endscan)+".")
        f.write("Scanning for "+str(startscan)+" <= CasG <= "+str(endscan)+".\n")
        LiEscancas=[x for x in ScanCasKillingLiE(G,endscan) if x[1]>=startscan]
    if LiEscancas==[]:
        print("Nothing to scan.")
        f.write("Nothing to scan.\n")
        f.close()
        return 0
    #Single out those Fourier modes that have homs into PotInstab
    format_string2="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}"
    print(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?"))
    f.write(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?")+"\n")
    for FourierMode in LiEscancas:
        LiERepG=FourierMode[0]
        repstring=omega(LiERepG.sage())
        CasG=FourierMode[1]
        print(repstring+": Branching...",end="\r")
        if LiERepG in branchdict:
            LiERepGbranched=branchdict[LiERepG]
        else:
            LiERepGbranched=lie.branch(LiERepG,LiEH,LiErm,LiEG)
            branchdict[LiERepG]=LiERepGbranched
        #we only need those modes where we don't know q(R)>E:
        print(repstring+": Checking whether it has homs to PotInstab...",end="\r")
        if any(lie(LiERepGbranched._name+'|'+x._name).sage()!=0 for x in PotInstab):
            #first, consider relevant part of Sym^2_0m:
            print(repstring+": Tensoring with Sym20m...                        ",end="\r")
            RelevantSym20m=[]
            FibrewiseTermList=[]
            for i in range(0,lsm):
                LiEisotype=LiESym20m.expon(i+1)
                isotype=tuple(LiEisotype.sage())
                branchedhom=lie.tensor(LiERepGbranched,lie.X(LiEisotype),LiEH)
                if lie(branchedhom._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantSym20m+=[LiEisotype]
                    FibrewiseTermList+=[CasGMinDict[isotype]-3/2*CasHDict[isotype]]
            FibrewiseTermMin=min(FibrewiseTermList)
            #only consider relevant part of Sym^2_0g:
            RelevantSym20g=lie.poly_null(rankG)
            for j in range(0,lsg):
                LiEGisotype=LiESym20g.expon(j+1)
                branched=LiEGisotypesbranched[j]
                if any(lie(branched._name+'|'+x._name).sage()!=0 for x in RelevantSym20m):
                    RelevantSym20g+=lie.X(LiEGisotype)
            #only consider relevant part of repG*Sym^2_0g:
            #to be precise it would be Hom(repG,Sym^2_0g), but Casimir stays the same
            print(repstring+": Tensoring with Sym20g...                        ",end="\r")
            #this is without multiplicities
            TensorProd=lie.tensor(lie.X(LiERepG),RelevantSym20g,LiEG)
            print(repstring+": Branching terms of RepG*Sym20g...               ",end="\r")
            tplength=TensorProd.length().sage()
            skipped=0
            RelevantTensorProd=[]
            for i in range(0,tplength):
                print(repstring+": Branching terms of RepG*Sym20g... ("\
                      +str(i+1)+" of "+str(tplength)+")         ",end="\r")
                TensorProdTerm=TensorProd.expon(i+1)
                if TensorProdTerm in branchdict:
                    branched=branchdict[TensorProdTerm]
                    skipped+=1
                else:
                    branched=lie.branch(TensorProdTerm,LiEH,LiErm,LiEG)
                    branchdict[TensorProdTerm]=branched
                if lie(branched._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantTensorProd+=[TensorProdTerm]
            #MixedTerm:
            print(repstring+": Computing Casimirs...                           ",end="\r")
            CasMixed=[CasKilling(SageRepfromList(G,x.sage())) for x in RelevantTensorProd]
            MixedTermMax=max(CasMixed)
            #Lower bound for LL:
            #formula: LL=3/2*CasG-1/2*MixedTerm+FibrewiseTerm-CasIsotropy
            #thus LL>=3/2*CasG-1/2*max(MixedTerm)+min(FibrewiseTerm)-CasIsotropy
            LLlower=3/2*CasG-1/2*MixedTermMax+FibrewiseTermMin-CasIsotropy
            goodbound=""
            if LLlower>2*E: goodbound=">"
            if LLlower==2*E: goodbound=">="
            print(repstring+": Checking for Killing tensors...                 ",end="\r")
            #check for Killing tensors
            killing=""
            dimHomSym3=0
            dimHomSym20=0
            for i in range(0,LiERepGbranched.length().sage()):
                summand=LiERepGbranched.expon(i+1)
                coef=LiERepGbranched.coef(i+1).sage()
                dimHomSym20+=coef*lie(LiESym20m._name+'|'+summand._name).sage()
                dimHomSym3+=coef*lie(LiESym3m._name+'|'+summand._name).sage()
            if dimHomSym3<dimHomSym20: killing="!"
            print(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound)\
                  +" (tplength: "+str(tplength)+", skipped: "+str(skipped)+")")
            f.write(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound)+"\n")
        else: print((75)*" ",end="\r")
    print("\nDone!")
    f.write("\nDone!")
    f.close()
#Sym2LowerBounds for flag manifolds:
#input: just the cartan type string of G
def Sym2LowerBoundsFullFlag(cartantype,startscan=0,endscan=-1):
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
    for i in range(0,lsm):
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
        return 0
    else:
        print("q(R)>E is potentially violated on a subspace of dimension: "+str(len(PotInstab)))
        f.write("q(R)>E is potentially violated on a subspace of dimension: "+str(len(PotInstab))+"\n")
    #Now: give C such that LL>2E if CasG>C.
    if endscan==-1:
        c=var("c")
        eqn=(c+min(AAMinList)/2-2*E)^2==max(AAMaxList)*(c-min(CasHDict.values()))
        C=max([x.right() for x in solve(eqn,c)])
        print("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5)))
        f.write("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5))+"\n")
        print("Scanning for "+str(startscan)+" <= CasG <= C.")
        f.write("Scanning for "+str(startscan)+" <= CasG <= C.\n")
        LiEscancas=[x for x in ScanCasKillingLiE(G,C) if x[1]>=startscan]
    else:
        print("Scanning for "+str(startscan)+" <= CasG <= "+str(endscan)+".")
        f.write("Scanning for "+str(startscan)+" <= CasG <= "+str(endscan)+".\n")
        LiEscancas=[x for x in ScanCasKillingLiE(G,endscan) if x[1]>=startscan]
    if LiEscancas==[]:
        print("Nothing to scan.")
        f.write("Nothing to scan.\n")
        f.close()
        return 0
    #Single out those Fourier modes that have homs into PotInstab
    format_string2="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}"
    print(format_string2.format("Fourier mode","CasG","LL >=","LL >= 2E?"))
    f.write(format_string2.format("Fourier mode","CasG","LL >=","LL >= 2E?")+"\n")
    for FourierMode in LiEscancas:
        LiERepG=FourierMode[0]
        repstring=omega(LiERepG.sage())
        CasG=FourierMode[1]
        print(repstring+": Branching...",end="\r")
        if LiERepG in branchdict:
            LiERepGbranched=branchdict[LiERepG]
        else:
            LiERepGbranched=LiERepG.Demazure(LiEG)
            branchdict[LiERepG]=LiERepGbranched
        #we only need those modes where we don't know q(R)>E:
        print(repstring+": Checking whether it has homs to PotInstab...",end="\r")
        if any(lie(LiERepGbranched._name+'|'+x._name).sage()!=0 for x in PotInstab):
            #first, consider relevant part of Sym^2_0m:
            print(repstring+": Tensoring with Sym20m...                        ",end="\r")
            RelevantSym20m=[]
            FibrewiseTermList=[]
            for i in range(0,lsm):
                LiEisotype=LiESym20m.expon(i+1)
                isotype=tuple(LiEisotype.sage())
                branchedhom=lie.tensor(LiERepGbranched,lie.X(LiEisotype),LiEH)
                if lie(branchedhom._name+'|'+LiEtrivial._name).sage()!=0:
                    RelevantSym20m+=[LiEisotype]
                    FibrewiseTermList+=[CasGMinDict[isotype]-3/2*CasHDict[isotype]]
            FibrewiseTermMin=min(FibrewiseTermList)
            #only consider relevant part of Sym^2_0g:
            RelevantSym20g=lie.poly_null(rankG)
            for j in range(0,lsg):
                LiEGisotype=LiESym20g.expon(j+1)
                branched=LiEGisotypesbranched[j]
                lr=len(RelevantSym20m)
                i=0
                while i<=lr:
                    x=RelevantSym20m[i]
                    if lie(branched._name+'|'+x._name).sage()!=0:
                        i=lr+1
                        RelevantSym20g+=lie.X(LiEGisotype)
                    else:
                        i+=1
                #for some reason, this doesn't work anymore:
                #if any(lie(branched._name+'|'+x._name).sage()!=0 for x in RelevantSym20m):
                #    RelevantSym20g+=lie.X(LiEGisotype)
            #only consider relevant part of repG*Sym^2_0g:
            #to be precise it would be Hom(repG,Sym^2_0g), but Casimir stays the same
            print(repstring+": Tensoring with Sym20g...                        ",end="\r")
            #this is without multiplicities
            TensorProd=lie.tensor(lie.X(LiERepG),RelevantSym20g,LiEG)
            print(repstring+": Branching terms of RepG*Sym20g...               ",end="\r")
            tplength=TensorProd.length().sage()
            skipped=0
            RelevantTensorProd=[]
            for i in range(0,tplength):
                print(repstring+": Branching terms of RepG*Sym20g... ("\
                      +str(i+1)+" of "+str(tplength)+")         ",end="\r")
                TensorProdTerm=TensorProd.expon(i+1)
                if TensorProdTerm in branchdict:
                    branched=branchdict[TensorProdTerm]
                    skipped+=1
                else:
                    branched=TensorProdTerm.Demazure(LiEG)
                    branchdict[TensorProdTerm]=branched
                if lie(branched._name+'|'+LiEtrivial._name).sage()!=0:
                    RelevantTensorProd+=[TensorProdTerm]
            #MixedTerm:
            print(repstring+": Computing Casimirs...                           ",end="\r")
            CasMixed=[CasKilling(SageRepfromList(G,x.sage())) for x in RelevantTensorProd]
            MixedTermMax=max(CasMixed)
            #Lower bound for LL:
            #formula: LL=3/2*CasG-1/2*MixedTerm+FibrewiseTerm-CasIsotropy
            #thus LL>=3/2*CasG-1/2*max(MixedTerm)+min(FibrewiseTerm)-CasIsotropy
            LLlower=3/2*CasG-1/2*MixedTermMax+FibrewiseTermMin-CasIsotropy
            goodbound=""
            if LLlower>2*E: goodbound=">"
            if LLlower==2*E: goodbound=">="
            print(repstring+": Checking for Killing tensors...                 ",end="\r")
            print(format_string2.format(repstring,str(CasG),str(LLlower),goodbound)\
                  +" (tplength: "+str(tplength)+", skipped: "+str(skipped)+")")
            f.write(format_string2.format(repstring,str(CasG),str(LLlower),goodbound)+"\n")
        else: print((75)*" ",end="\r")
    print("\nDone!")
    f.write("\nDone!")
    f.close()
    
#added startscan and endcan parameters for non-semisimple isotropy / full flags
#changed key type for CasHDict and CasGMinDict from LiEobject to tuple.
version="2023/03/04 19:15"
print("Sym2Bounds version "+version+".")
