#Input G (e.g. "D64"), H (e.g. "D8"), branched standard rep. (e.g. "1X[0,0,0,0,0,0,1,0]")
def Sym2LowerBoundsBig(Gstr,Hstr,o1,endscan=-1,verbose=False):
    print("Initializing...",end="\r")
    global LiEG
    LiEG=lie(Gstr)
    global dynkintype
    dynkintype=lie.Lie_code(LiEG).sage()[0]
    #dynkintype 1..4 should be classical
    if not(dynkintype in [1,2,3,4]): raise ValueError("Procedure only works for classical G.")
    H=WeylCharacterRing(Hstr,style="coroots")
    global LiEH
    LiEH=lie(Hstr.replace("x",""))
    f=open(Gstr+Hstr+"BigBranch.txt","w")
    fprint("====== Sym2BoundsBigBranch for G="+Gstr+", H="+Hstr+"======\n",f)
    if verbose: fprint("(verbose)",f)
    global rankG
    rankG=lie.Lie_rank(LiEG).sage()
    global rankH
    rankH=lie.Lie_rank(LiEH).sage()
    HFactors=GroupFactors(H)
    LiEHtrivial=lie.null(rankH)
    global LiEo1
    LiEo1=lie(o1)
    global branchdict
    branchdict={(1,):LiEo1}
    print("Computing isotropy representation...",end="\r")
    #AdG=G.adjoint_representation()
    #AdH=SemisimpleAdjointRep(HFactors,H)
    #IsotropyRep=b.branch(AdG)-AdH
    LiEAdG=lie.adjoint(LiEG)
    AdGCasRaw=LiECasSimpleRaw(LiEAdG.expon(1),LiEG)
    LiEAdH=lie.adjoint(LiEH)
    LiEAdBranched=RecBranch(LiEAdG.expon(1).sage())
    branchdict={LiEAdG.expon(1):LiEAdBranched}
    LiEIsotropy=LiEAdBranched-LiEAdH
    IsotropyRep=0
    for i in range(0,LiEIsotropy.length().sage()):
        IsotropyRep+=LiEIsotropy.coef(i+1).sage()*SageRepfromList(H,LiEIsotropy.expon(i+1).sage())
    fprint("Isotropy representation: "+str(LiEIsotropy),f)
    print("Computing Einstein constant...         ",end="\r")
    normconstants=CasNormConstants(IsotropyRep,HFactors)
    E=EinsteinConstant(IsotropyRep,normconstants,HFactors)
    CasIsotropy=2*E-1/2
    fprint("Einstein constant: E="+str(E)+"                ",f)
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
        CasSym20g+=[LiECasSimpleRaw(LiEGisotype,LiEG)/AdGCasRaw]
        if LiEGisotype in branchdict:
            LiEGisotypesbranched+=[branchdict[LiEGisotype]]
        else:
            #branched=lie.branch(LiEGisotype,LiEH,LiErm,LiEG)
            branched=RecBranch(LiEGisotype.sage())
            branchdict[LiEGisotype]=branched
            LiEGisotypesbranched+=[branched]
    #Now let's iterate over the isotypical components of Sym20m.
    spacing=max([25,2*rankH+6])#increase for composite H
    format_string1="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}"
    fprint(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min"),f)
    PotInstab=[]
    CasHDict={}
    CasGMinDict={}
    AAMinList=[]
    AAMaxList=[]
    for i in range(0,lsm):
        LiEisotype=LiESym20m.expon(i+1)
        isotype=SageRepfromList(H,LiEisotype.sage())
        casGofIsotype=[]
        for j in range(0,lsg):
            LiEGisotype=LiESym20g.expon(j+1)
            branched=LiEGisotypesbranched[j]
            if lie(branched._name+'|'+LiEisotype._name).sage()!=0:
                casGofIsotype+=[LiECasSimpleRaw(LiEGisotype,LiEG)/AdGCasRaw]
        CasGMin=min(casGofIsotype)
        CasGMinDict[LiEisotype]=CasGMin
        CasGMax=max(casGofIsotype)
        CasH=CasIrred(isotype,normconstants,HFactors)
        CasHDict[LiEisotype]=CasH
        #formula: A*A=prCasG-CasH-2CasIsotropy
        AAMin=CasGMin-CasH-2*CasIsotropy
        AAMinList+=[AAMin]
        AAMax=CasGMax-CasH-2*CasIsotropy
        AAMaxList+=[AAMax]
        #formula: 4q(R)=A*A+4CasH
        qRMin=1/4*AAMin+CasH
        fprint(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
                                     str(AAMin),str(AAMax),str(qRMin)),f)
        if qRMin<=E:
            PotInstab+=[LiEisotype]
    if PotInstab==[]:
        fprint("q(R)>E holds globally, implying stability.",f)
        f.close()
        return 0
    fprint("q(R)>E is potentially violated on: "+str(PotInstab),f)
    #Now: give C such that LL>2E if CasG>C.
    if endscan==-1:
        c=var("c")
        eqn=(c+min(AAMinList)/2-2*E)^2==max(AAMaxList)*(c-min(CasHDict.values()))
        C=max([x.right() for x in solve(eqn,c)])
        fprint("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5)),f)
    else: C=endscan
    print("Scanning Fourier modes...                        ",end="\r")
    LiEscancas=ScanCasJustLiE(LiEG,C)
    #Single out those Fourier modes that have homs into PotInstab
    format_string2="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}"
    if not(verbose): fprint(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?"),f)
    for FourierMode in LiEscancas:
        LiERepG=FourierMode[0]
        repstring=omega(LiERepG.sage())
        if verbose: fprint("LiERepG: "+repstring,f)
        CasG=FourierMode[1]
        LiERepGbranched=RecBranch(LiERepG.sage())
        if verbose: fprint("LiERepGbranched: "+str(LiERepGbranched),f)
        #we only need those modes where we don't know q(R)>E:
        if verbose: fprint(repstring+": Checking whether it has homs to PotInstab...",f)
        if any(lie(LiERepGbranched._name+'|'+x._name).sage()!=0 for x in PotInstab):
            #first, consider relevant part of Sym^2_0m:
            if verbose: fprint(repstring+": Tensoring with Sym20m...",f)
            RelevantSym20m=[]
            FibrewiseTermList=[]
            for i in range(0,lsm):
                LiEisotype=LiESym20m.expon(i+1)
                if verbose: fprint("LiEisotype: "+str(LiEisotype),f)
                branchedhom=lie.tensor(LiERepGbranched,lie.X(LiEisotype),LiEH)
                if lie(branchedhom._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantSym20m+=[LiEisotype]
                    FibrewiseTermList+=[CasGMinDict[LiEisotype]-3/2*CasHDict[LiEisotype]]
            FibrewiseTermMin=min(FibrewiseTermList)
            #only consider relevant part of Sym^2_0g:
            if verbose: fprint(repstring+": Collecting relevant part of Sym20g...",f)
            RelevantSym20g=lie.poly_null(rankG)
            for j in range(0,lsg):
                LiEGisotype=LiESym20g.expon(j+1)
                branched=LiEGisotypesbranched[j]
                if any(lie(branched._name+'|'+x._name).sage()!=0 for x in RelevantSym20m):
                    if verbose: fprint("Adding: "+omega(LiEGisotype.sage()),f)
                    RelevantSym20g+=lie.X(LiEGisotype)
            #only consider relevant part of repG*Sym^2_0g:
            #to be precise it would be Hom(repG,Sym^2_0g), but Casimir stays the same
            if verbose: fprint(repstring+": Tensoring with RelevantSym20g...",f)
            #this is without multiplicities
            TensorProd=lie.tensor(lie.X(LiERepG),RelevantSym20g,LiEG)
            if verbose: fprint(repstring+": Branching terms of RepG*Sym20g...",f)
            tplength=TensorProd.length().sage()
            RelevantTensorProd=[]
            for i in range(0,tplength):
                if verbose: fprint(repstring+": Branching terms of RepG*Sym20g... ("+str(i+1)+" of "+str(tplength)+")",f)
                TensorProdTerm=TensorProd.expon(i+1)
                branched=RecBranch(TensorProdTerm.sage())
                if verbose: fprint("Checking whether TensorProdTerm is relevant...",f)
                if lie(branched._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantTensorProd+=[TensorProdTerm]
                    if verbose: fprint("Yes.",f)
                else:
                    if verbose: fprint("No.",f)
            #MixedTerm:
            if verbose: fprint(repstring+": Computing Casimirs...",f)
            CasMixed=[LiECasSimpleRaw(x,LiEG)/AdGCasRaw for x in RelevantTensorProd]
            MixedTermMax=max(CasMixed)
            #Lower bound for LL:
            #formula: LL=3/2*CasG-1/2*MixedTerm+FibrewiseTerm-CasIsotropy
            #thus LL>=3/2*CasG-1/2*max(MixedTerm)+min(FibrewiseTerm)-CasIsotropy
            LLlower=3/2*CasG-1/2*MixedTermMax+FibrewiseTermMin-CasIsotropy
            goodbound=""
            if LLlower>2*E: goodbound=">"
            if LLlower==2*E: goodbound=">="
            if verbose: fprint(repstring+": Checking for Killing tensors...",f)
            #check for Killing tensors
            killing=""
            dimHomSym3=0
            dimHomSym20=0
            for i in range(0,LiERepGbranched.length().sage()):
                summand=LiERepGbranched.expon(i+1)
                if verbose: fprint("Term of LiERepGbranched: "+str(summand),f)
                coef=LiERepGbranched.coef(i+1).sage()
                dimHomSym20+=coef*lie(LiESym20m._name+'|'+summand._name).sage()
                dimHomSym3+=coef*lie(LiESym3m._name+'|'+summand._name).sage()
            if dimHomSym3<dimHomSym20:
                killing="!"
                if verbose: fprint("Killing tensors found!",f)
            if verbose: fprint(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?"),f)
            fprint(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound),f)
        else: print((spacing+50)*" ",end="\r")
    fprint("\nDone!",f)
    f.close()

#Input G (e.g. "C20"), semisimple part of H (e.g. "A12xC7"), branched standard rep., restriction matrix (as in familyXIV)
def Sym2LowerBoundsBigWithTorus(Gstr,Hssstr,o1,rm,endscan=-1,verbose=False):
    print("Initializing G...",end="\r")
    global LiEG
    LiEG=lie(Gstr)
    global rankG
    rankG=lie.Lie_rank(LiEG).sage()
    global dynkintype
    dynkintype=lie.Lie_code(LiEG).sage()[0]
    if not(dynkintype in [1,2,3,4]): raise ValueError("Procedure only works for classical G.")
    print("Initializing H...",end="\r")
    global rankH
    global LiEH
    rankH=len(rm[0])
    Hss=WeylCharacterRing(Hssstr,style="coroots")
    rankHss=Hss.rank()
    rankT=rankH-rankHss
    HssFactors=GroupFactors(Hss)
    LiEHss=lie(Hssstr.replace("x",""))
    LiEH=LiEHss*lie("T"+str(rankT))
    LiEHtrivial=lie.null(rankH)
    print("Initializing rm and o1...",end="\r")
    LiErm=listtoliebyhalf(rm)
    global LiEo1
    LiEo1=lie(o1)
    global branchdict
    branchdict={(1,):LiEo1}
    f=open(Gstr+str(LiEH)+"BigBranch.txt","w")
    fprint("====== Sym2BoundsBigBranch for G="+Gstr+", H="+str(LiEH)+"======\n",f)
    if verbose: fprint("(verbose)",f)
    print("Computing isotropy representation...",end="\r")
    #AdG=G.adjoint_representation()
    #AdH=SemisimpleAdjointRep(HFactors,H)
    #IsotropyRep=b.branch(AdG)-AdH
    LiEAdG=lie.adjoint(LiEG)
    AdGCasRaw=LiECasSimpleRaw(LiEAdG.expon(1),LiEG)
    LiEAdH=lie.adjoint(LiEH)
    LiEAdBranched=RecBranch(LiEAdG.expon(1).sage())
    branchdict={LiEAdG.expon(1):LiEAdBranched}
    LiEIsotropy=LiEAdBranched-LiEAdH
    li=LiEIsotropy.length().sage()
    fprint("Isotropy representation: "+str(LiEIsotropy),f)
    print("Computing Einstein constant...         ",end="\r")
    normconstants=NormConstantsSS(LiEIsotropy,HssFactors,LiEH)
    cm=Matrix(lie.Cartan(LiEG))
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
    fprint("Einstein constant: E="+str(E)+"                ",f)
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
        CasSym20g+=[LiECasSimpleRaw(LiEGisotype,LiEG)/AdGCasRaw]
        if LiEGisotype in branchdict:
            LiEGisotypesbranched+=[branchdict[LiEGisotype]]
        else:
            #branched=lie.branch(LiEGisotype,LiEH,LiErm,LiEG)
            branched=RecBranch(LiEGisotype.sage())
            branchdict[LiEGisotype]=branched
            LiEGisotypesbranched+=[branched]
    #Now let's iterate over the isotypical components of Sym20m.
    spacing=max([25,2*rankH+6])#increase for composite H
    format_string1="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}"
    fprint(format_string1.format("Isotype in Sym^2_0","CasGMin","CasGMax","CasH","A*AMin","A*AMax","q(R)Min"),f)
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
                casGofIsotype+=[LiECasSimpleRaw(LiEGisotype,LiEG)/AdGCasRaw]
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
        fprint(format_string1.format(str(LiEisotype),str(CasGMin),str(CasGMax),str(CasH),\
                                     str(AAMin),str(AAMax),str(qRMin)),f)
        if qRMin<=E:
            PotInstab+=[LiEisotype]
    if PotInstab==[]:
        fprint("q(R)>E holds globally, implying stability.",f)
        f.close()
        return 0
    fprint("q(R)>E is potentially violated on: "+str(PotInstab),f)
    #Now: give C such that LL>2E if CasG>C.
    if endscan==-1:
        c=var("c")
        eqn=(c+min(AAMinList)/2-2*E)^2==max(AAMaxList)*(c-min(CasHDict.values()))
        C=max([x.right() for x in solve(eqn,c)])
        fprint("Crude estimate: LL > 2E if CasG > "+str(C)+", approx. "+str(numerical_approx(C,digits=5)),f)
    else: C=endscan
    print("Scanning Fourier modes...                        ",end="\r")
    LiEscancas=ScanCasJustLiE(LiEG,C)
    #Single out those Fourier modes that have homs into PotInstab
    format_string2="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}"
    if not(verbose): fprint(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?"),f)
    for FourierMode in LiEscancas:
        LiERepG=FourierMode[0]
        repstring=omega(LiERepG.sage())
        if verbose: fprint("LiERepG: "+repstring,f)
        CasG=FourierMode[1]
        LiERepGbranched=RecBranch(LiERepG.sage())
        if verbose: fprint("LiERepGbranched: "+str(LiERepGbranched),f)
        #we only need those modes where we don't know q(R)>E:
        if verbose: fprint(repstring+": Checking whether it has homs to PotInstab...",f)
        if any(lie(LiERepGbranched._name+'|'+x._name).sage()!=0 for x in PotInstab):
            #first, consider relevant part of Sym^2_0m:
            if verbose: fprint(repstring+": Tensoring with Sym20m...",f)
            RelevantSym20m=[]
            FibrewiseTermList=[]
            for i in range(0,lsm):
                LiEisotype=LiESym20m.expon(i+1)
                isotype=tuple(LiEisotype.sage())
                if verbose: fprint("LiEisotype: "+str(LiEisotype),f)
                branchedhom=lie.tensor(LiERepGbranched,lie.X(LiEisotype),LiEH)
                if lie(branchedhom._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantSym20m+=[LiEisotype]
                    FibrewiseTermList+=[CasGMinDict[isotype]-3/2*CasHDict[isotype]]
            FibrewiseTermMin=min(FibrewiseTermList)
            #only consider relevant part of Sym^2_0g:
            if verbose: fprint(repstring+": Collecting relevant part of Sym20g...",f)
            RelevantSym20g=lie.poly_null(rankG)
            for j in range(0,lsg):
                LiEGisotype=LiESym20g.expon(j+1)
                branched=LiEGisotypesbranched[j]
                if any(lie(branched._name+'|'+x._name).sage()!=0 for x in RelevantSym20m):
                    if verbose: fprint("Adding: "+omega(LiEGisotype.sage()),f)
                    RelevantSym20g+=lie.X(LiEGisotype)
            #only consider relevant part of repG*Sym^2_0g:
            #to be precise it would be Hom(repG,Sym^2_0g), but Casimir stays the same
            if verbose: fprint(repstring+": Tensoring with RelevantSym20g...",f)
            #this is without multiplicities
            TensorProd=lie.tensor(lie.X(LiERepG),RelevantSym20g,LiEG)
            if verbose: fprint(repstring+": Branching terms of RepG*Sym20g...",f)
            tplength=TensorProd.length().sage()
            RelevantTensorProd=[]
            for i in range(0,tplength):
                if verbose: fprint(repstring+": Branching terms of RepG*Sym20g... ("+str(i+1)+" of "+str(tplength)+")",f)
                TensorProdTerm=TensorProd.expon(i+1)
                branched=RecBranch(TensorProdTerm.sage())
                if verbose: fprint("Checking whether TensorProdTerm is relevant...",f)
                if lie(branched._name+'|'+LiEHtrivial._name).sage()!=0:
                    RelevantTensorProd+=[TensorProdTerm]
                    if verbose: fprint("Yes.",f)
                else:
                    if verbose: fprint("No.",f)
            #MixedTerm:
            if verbose: fprint(repstring+": Computing Casimirs...",f)
            CasMixed=[LiECasSimpleRaw(x,LiEG)/AdGCasRaw for x in RelevantTensorProd]
            MixedTermMax=max(CasMixed)
            #Lower bound for LL:
            #formula: LL=3/2*CasG-1/2*MixedTerm+FibrewiseTerm-CasIsotropy
            #thus LL>=3/2*CasG-1/2*max(MixedTerm)+min(FibrewiseTerm)-CasIsotropy
            LLlower=3/2*CasG-1/2*MixedTermMax+FibrewiseTermMin-CasIsotropy
            goodbound=""
            if LLlower>2*E: goodbound=">"
            if LLlower==2*E: goodbound=">="
            if verbose: fprint(repstring+": Checking for Killing tensors...",f)
            #check for Killing tensors
            killing=""
            dimHomSym3=0
            dimHomSym20=0
            for i in range(0,LiERepGbranched.length().sage()):
                summand=LiERepGbranched.expon(i+1)
                if verbose: fprint("Term of LiERepGbranched: "+str(summand),f)
                coef=LiERepGbranched.coef(i+1).sage()
                dimHomSym20+=coef*lie(LiESym20m._name+'|'+summand._name).sage()
                dimHomSym3+=coef*lie(LiESym3m._name+'|'+summand._name).sage()
            if dimHomSym3<dimHomSym20:
                killing="!"
                if verbose: fprint("Killing tensors found!",f)
            if verbose: fprint(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?"),f)
            fprint(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound),f)
        else: print((spacing+50)*" ",end="\r")
    fprint("\nDone!",f)
    f.close()
    
########################
    
def fprint(string,file):
    print(string)
    file.write(string+"\n")
def listtoliebyhalf(inputlist):
    l=len(inputlist)
    trunc=int(l/2)
    if trunc==0: return lie(inputlist)
    auxlist1=lie(inputlist[0:trunc])
    auxlist2=lie(inputlist[trunc:l])
    return lie(auxlist1._name+"^"+auxlist2._name)
def poltoliebychunks(inputstr,maxlength):
    termlist=inputstr.split("+")
    l=len(termlist)
    splitlist=[termlist[i:i+maxlength] for i in range(0,l,maxlength)]
    auxpol=[lie("+".join(x)) for x in splitlist]
    result=auxpol[0]
    for x in auxpol[1:]:
        result=result+x
    return result   
def RecBranch(rep):
    if len(rep)>rankG-2:
        if (dynkintype==2&rep[-1]>0)|(dynkintype==4&rep[-2]+rep[-1]>0):
            raise ValueError("This is a spin representation. Please use a different branching procedure.")
    s=sum(rep)
    #trivial rep.:
    if s==0: return lie.poly_one(rankH)
    repc=rep.copy()
    while repc[-1]==0:
        repc.pop(-1)
    lookupindex=tuple(repc)
    #print("Looking up: "+str(lookupindex))
    if tuple(repc) in branchdict:
        LiEresult=branchdict[lookupindex]
        #print("Found: "+str(LiEresult))
        return LiEresult
    #print("Not found.")
    #print("Branching: "+str(repc))
    #global steps
    #steps+=1
    #print("Step: "+str(steps),"\r")
    #multiples of a fund. weight:
    nonzero=[x for x in repc if x!=0]
    if len([x for x in repc if x!=0])==1:
        mult=nonzero[0]
        pos=rep.index(mult)
        if mult==1:
            LiEresult=RecBranchFund(pos)
        else:
            LiEresult=RecBranchSym(pos,mult)
        branchdict[tuple(repc)]=LiEresult
        return LiEresult
    #cartan product
    l=len(repc)
    rep2=(l-1)*[0]+[repc[-1]]
    repc.pop(-1)
    LiEresult=RecBranchCartan(repc,rep2)
    #print("Saving "+str(LiEresult)+" for "+str(lookupindex))
    branchdict[lookupindex]=LiEresult
    return LiEresult
def RecBranchSym(pos,mult):
    #print("RecBranchSym. pos="+str(pos)+", mult="+str(mult))
    if mult==1: return RecBranchAlt(pos)
    n=(pos+1)*mult+3
    if n>rankG: raise ValueError("Rank of G is too small for this procedure.")
    LiED=lie.Lie_group(dynkintype,n)
    DLambda=n*[0]
    DLambda[pos]=1
    DGoal=n*[0]
    DGoal[pos]=mult
    LiEDLambda=lie("1X"+str(DLambda))
    LiEDGoal=lie("1X"+str(DGoal))
    LiEDSym=lie.sym_tensor(mult,LiEDLambda,LiED)
    LiEDRemainder=LiEDSym-LiEDGoal
    LiEresult=lie.sym_tensor(mult,RecBranch(DLambda),LiEH)
    for i in range(0,LiEDRemainder.length().sage()):
        LiEresult=LiEresult-LiEDRemainder.coef(i+1)*RecBranch(LiEDRemainder.expon(i+1).sage())
    return LiEresult
def RecBranchFund(pos):
    #print("RecBranchFund. pos="+str(pos))
    if dynkintype==3:
        if pos==0:
            return LiEo1
        else:
            return lie.alt_tensor(pos+1,LiEo1,LiEH)-lie.alt_tensor(pos-1,LiEo1,LiEH)
    else:
        return lie.alt_tensor(pos+1,LiEo1,LiEH)
def RecBranchCartan(rep1,rep2):
    #print("RecBranchCartan. rep1="+str(rep1)+", rep2="+str(rep2))
    n=len(rep1)+len(rep2)+2
    if n>rankG: raise ValueError("Rank of G is too small for this procedure.")
    LiED=lie.Lie_group(dynkintype,n)
    #Drep1=rep1+(n-len(rep1))*[0]
    #Drep2=rep2+(n-len(rep2))*[0]
    LiEDrep1=lie(rep1+(n-len(rep1))*[0])
    LiEDrep2=lie(rep2+(n-len(rep2))*[0])
    LiEDGoal=lie("X("+LiEDrep1._name+"+"+LiEDrep2._name+")")
    LiEDtensor=lie.tensor(LiEDrep1,LiEDrep2,LiED)
    LiEDRemainder=LiEDtensor-LiEDGoal
    #print("Remainder: "+str(LiEDRemainder))
    LiErep1br=RecBranch(rep1)
    LiErep2br=RecBranch(rep2)
    #print("Now tensoring over E8: rep1: "+str(rep1)+", branched: "+str(LiErep1br)+", and rep2: "\
    #      +str(rep2)+", branched: "+str(LiErep2br))
    LiEresult=lie.tensor(LiErep1br,LiErep2br,LiEH)
    #print("Branched tensor product of "+str(rep1)+" and "+str(rep2)+": "+str(LiEresult))
    for i in range(0,LiEDRemainder.length().sage()):
        subtr=LiEDRemainder.coef(i+1)*RecBranch(LiEDRemainder.expon(i+1).sage())
        #print("Subtracting branched "+str(LiEDRemainder.expon(i+1))+", namely "+str(subtr))
        LiEresult=LiEresult-subtr
    #print("Cartan product of "+str(rep1)+" and "+str(rep2)+": "+str(LiEresult))
    return LiEresult
#unnormalized Casimir through LiE (for simple G)
def LiECasSimpleRaw(LiErep,LiEG):
    a=lie(LiErep._name+"*i_Cartan("+LiEG._name+")")
    b=lie("("+LiErep._name+"+2*all_one(size("+LiErep._name+")))*i_Cartan("+LiEG._name+")")
    return lie.inprod(a,b,LiEG).sage()
def ScanCasJustLiE(LiEG,q):
    #output: [[LiErep,CasG],...]
    #CasAd=CasAdjoint(G)
    ad=lie.adjoint(LiEG).expon(1)
    casadraw=LiECasSimpleRaw(ad,LiEG)
    rank=lie.Lie_rank(LiEG).sage()
    nullweight=lie.null(rank)
    out=[[nullweight,0]]
    for i in range(1,rank+1):
        tempout=list(out)
        for pair in tempout:
            LiErep=lie(pair[0]._name)
            lie.eval(LiErep._name+"["+str(i)+"]+=1")
            c=LiECasSimpleRaw(LiErep,LiEG)/casadraw
            while c<=q:
                out+=[[lie(LiErep._name),c]]
                lie.eval(LiErep._name+"["+str(i)+"]+=1")
                c=LiECasSimpleRaw(LiErep,LiEG)/casadraw
    return sorted(out,key=lambda pair: pair[1])
    
version="2023/05/04 20:00"
print("Sym2BoundsBigBranch version "+version+".")
