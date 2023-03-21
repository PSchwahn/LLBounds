#the second part of the algorithm in more verbose. The first part has to be performed beforehand.


def fprint(string,file):
    print(string)
    file.write(string+"\n")

f=open(str(LiEG)+str(LiEH)+"verbose.txt","w")
#Single out those Fourier modes that have homs into PotInstab
format_string2="{:<"+str(spacing)+"}{:<15}{:<15}{:<15}{:<15}"
for FourierMode in LiEscancas:
    LiERepG=FourierMode[0]
    repstring=omega(LiERepG.sage())
    fprint("LiERepG: "+repstring,f)
    CasG=FourierMode[1]
    fprint(repstring+": Branching...",f)
    if LiERepG in branchdict:
        LiERepGbranched=branchdict[LiERepG]
    else:
        LiERepGbranched=lie.branch(LiERepG,LiEH,LiErm,LiEG)
        branchdict[LiERepG]=LiERepGbranched
    fprint("LiERepGbranched: "+str(LiERepGbranched),f)
    #we only need those modes where we don't know q(R)>E:
    fprint(repstring+": Checking whether it has homs to PotInstab...",f)
    if any(lie(LiERepGbranched._name+'|'+x._name).sage()!=0 for x in PotInstab):
        #first, consider relevant part of Sym^2_0m:
        fprint(repstring+": Tensoring with Sym20m...",f)
        RelevantSym20m=[]
        FibrewiseTermList=[]
        for i in range(0,lsm):
            LiEisotype=LiESym20m.expon(i+1)
            fprint("LiEisotype: "+str(LiEisotype),f)
            branchedhom=lie.tensor(LiERepGbranched,lie.X(LiEisotype),LiEH)
            if lie(branchedhom._name+'|'+LiEHtrivial._name).sage()!=0:
                RelevantSym20m+=[LiEisotype]
                FibrewiseTermList+=[CasGMinDict[LiEisotype]-3/2*CasHDict[LiEisotype]]
        FibrewiseTermMin=min(FibrewiseTermList)
        #only consider relevant part of Sym^2_0g:
        fprint(repstring+": Collecting relevant part of Sym20g...",f)
        RelevantSym20g=lie.poly_null(rankG)
        for j in range(0,lsg):
            LiEGisotype=LiESym20g.expon(j+1)
            branched=LiEGisotypesbranched[j]
            if any(lie(branched._name+'|'+x._name).sage()!=0 for x in RelevantSym20m):
                fprint("Adding: "+omega(LiEGisotype.sage()),f)
                RelevantSym20g+=lie.X(LiEGisotype)
        #only consider relevant part of repG*Sym^2_0g:
        #to be precise it would be Hom(repG,Sym^2_0g), but Casimir stays the same
        fprint(repstring+": Tensoring LiERepG with RelevantSym20g...",f)
        #this is without multiplicities
        TensorProd=lie.tensor(lie.X(LiERepG),RelevantSym20g,LiEG)
        fprint(repstring+": Branching terms of RepG*Sym20g...",f)
        tplength=TensorProd.length().sage()
        skipped=0
        RelevantTensorProd=[]
        for i in range(0,tplength):
            fprint(repstring+": Branching terms of RepG*Sym20g... ("+str(i+1)+" of "+str(tplength)+")",f)
            TensorProdTerm=TensorProd.expon(i+1)
            fprint("TensorProdTerm: "+omega(TensorProdTerm.sage()),f)
            if TensorProdTerm in branchdict:
                branched=branchdict[TensorProdTerm]
                skipped+=1
                fprint("Skipped.",f)
            else:
                fprint("Branching TensorProdTerm...",f)
                branched=lie.branch(TensorProdTerm,LiEH,LiErm,LiEG)
                branchdict[TensorProdTerm]=branched
            if lie(branched._name+'|'+LiEHtrivial._name).sage()!=0:
                RelevantTensorProd+=[TensorProdTerm]
        #MixedTerm:
        fprint(repstring+": Computing Casimirs...",f)
        CasMixed=[CasIrred(SageRepfromList(G,x.sage()),[1],[G]) for x in RelevantTensorProd]
        MixedTermMax=max(CasMixed)
        #Lower bound for LL:
        #formula: LL=3/2*CasG-1/2*MixedTerm+FibrewiseTerm-CasIsotropy
        #thus LL>=3/2*CasG-1/2*max(MixedTerm)+min(FibrewiseTerm)-CasIsotropy
        LLlower=3/2*CasG-1/2*MixedTermMax+FibrewiseTermMin-CasIsotropy
        goodbound=""
        if LLlower>2*E: goodbound=">"
        if LLlower==2*E: goodbound=">="
        fprint(repstring+": Checking for Killing tensors...",f)
        #check for Killing tensors
        killing=""
        dimHomSym3=0
        dimHomSym20=0
        for i in range(0,LiERepGbranched.length().sage()):
            summand=LiERepGbranched.expon(i+1)
            fprint("Term of LiERepGbranched: "+str(summand),f)
            coef=LiERepGbranched.coef(i+1).sage()
            dimHomSym20+=coef*lie(LiESym20m._name+'|'+summand._name).sage()
            dimHomSym3+=coef*lie(LiESym3m._name+'|'+summand._name).sage()
        if dimHomSym3<dimHomSym20:
            killing="!"
            fprint("Killing tensors found!",f)
        fprint(format_string2.format("Fourier mode","CasG","KT?","LL >=","LL >= 2E?"),f)
        fprint(format_string2.format(repstring,str(CasG),killing,str(LLlower),goodbound)\
              +" (tplength: "+str(tplength)+", skipped: "+str(skipped)+")",f)
    else: print((spacing+50)*" ",end="\r")
fprint("\nDone!",f)
f.close()
