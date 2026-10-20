# Monta um par de reads Sanger (F + R) em consenso unico, ponderado por qualidade.
# Uso: python Scripts/montar_par.py <F.ab1> <R.ab1> <nome> <saida.fasta> [seqPrimerF seqPrimerR] [--mott 0.05]
# seqPrimerF/R: sequências 5'->3' dos primers, para cortar os sítios das pontas ("-" = não cortar daquele lado)
# R = - : modo read unico (so apara e exporta).
# --mott: limite do aparo por qualidade (0.01 ~ Q20, rigoroso; 0.05 ~ Q13, recomendado p/ ITS longo)
# --trnlf: corta o consenso de trnL-F para comecar em GACTTAAATAATTTGAGC e terminar em ...AGAGGACTGAAAATC
# Apara por Mott (0.01), inverte o R, alinhamento de sobreposicao (pontas livres),
# e decide cada coluna pela qualidade. Discordancias com QV parecido viram IUPAC e sao listadas.
import struct, sys
IUPAC={frozenset('AG'):'R',frozenset('CT'):'Y',frozenset('GC'):'S',frozenset('AT'):'W',frozenset('GT'):'K',frozenset('AC'):'M'}
COMP=str.maketrans('ACGTRYSWKMN','TGCAYRSWMKN')
def read_abif(p):
    d=open(p,'rb').read()
    num=struct.unpack('>I',d[18:22])[0]; off=struct.unpack('>I',d[26:30])[0]; t={}
    for i in range(num):
        e=d[off+i*28:off+(i+1)*28]
        name,n,et,es,ne,ds,do,_=struct.unpack('>4sIHHIIII',e)
        t[(name.decode('latin1'),n)]=(e[20:24] if ds<=4 else d[do:do+ds])[:ds]
    return t['PBAS',2].decode('latin1').upper(), list(t['PCON',2])
def mott(q,cut=0.01):
    best=s=0;bs=be=start=0
    for i,x in enumerate(q):
        s+=cut-10**(-x/10)
        if s<0: s=0; start=i+1
        if s>best: best=s; bs=start; be=i+1
    return bs,be
def align(a,b,M=2,X=-3,G=-5):
    n,m=len(a),len(b)
    S=[[0]*(m+1) for _ in range(n+1)]; T=[[0]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1):
        ai=a[i-1]; Si=S[i]; Sp=S[i-1]; Ti=T[i]
        for j in range(1,m+1):
            d=Sp[j-1]+(M if ai==b[j-1] else X); u=Sp[j]+G; l=Si[j-1]+G
            if d>=u and d>=l: Si[j]=d; Ti[j]=0
            elif u>=l: Si[j]=u; Ti[j]=1
            else: Si[j]=l; Ti[j]=2
    # fim livre: melhor na ultima linha ou coluna
    best=max([(S[n][j],n,j) for j in range(m+1)]+[(S[i][m],i,m) for i in range(n+1)])
    _,i,j=best; A=[];B=[]
    A+=[(k) for k in range(n-1,i-1,-1)]; B+=[None]*(n-i)
    B2=[(k) for k in range(m-1,j-1,-1)]; A2=[None]*(m-j)
    A=A+A2; B=B+B2
    while i>0 and j>0:
        t=T[i][j]
        if t==0: i-=1;j-=1;A.append(i);B.append(j)
        elif t==1: i-=1;A.append(i);B.append(None)
        else: j-=1;A.append(None);B.append(j)
    while i>0: i-=1;A.append(i);B.append(None)
    while j>0: j-=1;A.append(None);B.append(j)
    return A[::-1],B[::-1],best[0]
CORTE=0.01
TRNLF='--trnlf' in sys.argv
if TRNLF: sys.argv.remove('--trnlf')
if '--mott' in sys.argv:
    i=sys.argv.index('--mott'); CORTE=float(sys.argv[i+1]); del sys.argv[i:i+2]
f,r,nome,out=sys.argv[1:5]
fs,fq=read_abif(f)
if r=='-':
    # read unico: so apara; Q<10 vira N; bases Q10-19 listadas
    a,b=mott(fq,CORTE); fs,fq=fs[a:b],fq[a:b]
    seq=''.join(c if q>=10 else 'N' for c,q in zip(fs,fq)).strip('N')
    baixas=[i+1 for i,q in enumerate(fq) if 10<=q<20]
    open(out,'w').write(f">{nome}\n"+'\n'.join(seq[i:i+70] for i in range(0,len(seq),70))+'\n')
    print(f"{nome}\tREAD_UNICO\taparado={len(fs)}\tconsenso={len(seq)}\tN/IUPAC={sum(c not in 'ACGT' for c in seq)}\tbases_Q10-19={len(baixas)}")
    if baixas: print("    posicoes Q10-19 (so 1 read, conferir):", baixas)
    sys.exit()
rs,rq=read_abif(r)
a,b=mott(fq,CORTE); fs,fq=fs[a:b],fq[a:b]
a,b=mott(rq,CORTE); rs,rq=rs[a:b],rq[a:b]
rs=rs.translate(COMP)[::-1]; rq=rq[::-1]
A,B,score=align(fs,rs)
cons=[];notas=[];ov=disc=0
for ia,ib in zip(A,B):
    ca=fs[ia] if ia is not None else None; cb=rs[ib] if ib is not None else None
    qa=fq[ia] if ia is not None else 0; qb=rq[ib] if ib is not None else 0
    if ca and cb:
        ov+=1
        if ca==cb: cons.append(ca); continue
        disc+=1; pos=len(cons)+1
        EXP={'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'GC','W':'AT','K':'GT','M':'AC','N':'ACGT'}
        inter=set(EXP.get(ca,'ACGT'))&set(EXP.get(cb,'ACGT'))
        if len(inter)==1 and (ca not in 'ACGT' or cb not in 'ACGT'):
            win=inter.pop(); cons.append(win); notas.append(f"pos {pos}: F={ca}(Q{qa}) R={cb}(Q{qb}) -> {win} (ambiguo compativel)"); continue
        if abs(qa-qb)>=10:
            win=ca if qa>qb else cb; cons.append(win); notas.append(f"pos {pos}: F={ca}(Q{qa}) R={cb}(Q{qb}) -> {win}")
        else:
            k=IUPAC.get(frozenset(ca+cb),'N'); cons.append(k); notas.append(f"pos {pos}: F={ca}(Q{qa}) R={cb}(Q{qb}) -> {k} REVISAR")
    elif ca or cb:
        c,q=(ca,qa) if ca else (cb,qb)
        if ov==0 or True:
            if q>=20: cons.append(c)
            elif q>=10: cons.append(c); notas.append(f"pos {len(cons)}: so 1 read, {c}(Q{q}) baixa")
            else: cons.append('N')
seq=''.join(cons).strip('N').replace('-','')
# remove sitios de primer opcionais (argumentos 5 e 6: primer F e primer R, 5'->3')
def fuzzy(s,p,maxmm=5):
    best=None
    for i in range(0,len(s)-len(p)+1):
        mm=sum(1 for x,y in zip(s[i:i+len(p)],p) if x!=y)
        if mm<=maxmm and (best is None or mm<best[1]): best=(i,mm)
    return best
corte=''
if len(sys.argv)>6:   # "-" = sem sequência para esse primer
    pf=sys.argv[5].upper(); pr=sys.argv[6].upper().translate(COMP)[::-1]
    mmF=max(1,len(pf)//5); mmR=max(1,len(pr)//5)   # até ~20% de diferenças (bases de baixa qualidade nas pontas)
    h=fuzzy(seq[:60],pf,mmF) if pf!='-' else None
    if h: seq=seq[h[0]+len(pf):]; corte+=f' primerF@{h[0]}({h[1]}mm)'
    elif pf!='-': corte+=' primerF_nao_achado'
    t=fuzzy(seq[-60:],pr,mmR) if pr!='-' else None
    if t: seq=seq[:len(seq)-60+t[0]]; corte+=f' primerR({t[1]}mm)'
    elif pr!='-': corte+=' primerR_nao_achado'
if TRNLF:
    # trnL-F: o consenso começa em GACTTAAATAATTTGAGC e termina em ...AGAGGACTGAAAATC (fora ficam os sítios dos primers c/f)
    ini=seq.find('GACTTAAATAATTTGAGC'); fim=seq.rfind('AGAGGACTGAAAATC')
    if 0<=ini<80: seq=seq[ini:]; corte+=' inicio_trnLF_ok'
    else: notas.append('AVISO trnL-F: inicio GACTTAAATAATTTGAGC nao encontrado nas primeiras 80 bases - conferir a ponta 5\' no MEGA')
    fim=seq.rfind('AGAGGACTGAAAATC')
    if fim>=0 and len(seq)-(fim+15)<80: seq=seq[:fim+15]; corte+=' fim_trnLF_ok'
    else: notas.append('AVISO trnL-F: fim ...AGAGGACTGAAAATC nao encontrado - read curto ou ponta 3\' a conferir')
open(out,'w').write(f">{nome}\n"+'\n'.join(seq[i:i+70] for i in range(0,len(seq),70))+'\n')
print(f"{nome}\tF_aparado={len(fs)}\tR_aparado={len(rs)}\tsobreposicao={ov}\tdiscord={disc}\tconsenso={len(seq)}\tN/IUPAC={sum(c not in 'ACGT' for c in seq)}{corte}")
for x in notas: print("   ",x)
