# Junta duas sequencias consenso que se sobrepoem (ex.: metade ITS1 [ITS1F+5.8R] e metade ITS2 [5.8F+ITS2R]).
# Uso: python3 Scripts/juntar_metades.py <metade5.fasta> <metade3.fasta> <nome> <saida.fasta> [--lacuna N]
# --lacuna N: as metades NAO se encostam (ex.: read faltando); insere N bases 'N' entre elas.
#             Estime N mapeando as metades numa sequencia proxima completa (blastn -subject).
# A ordem e o sentido das metades sao conferidos sozinhos: o script testa as combinacoes, escolhe a que
# sobrepoe e deixa o ITS no sentido 18S->26S (5.8S direto). Se nao houver sobreposicao, para com aviso.
# Cada arquivo deve ter UMA sequencia (o consenso da metade); gaps (-) e ? do MEGA sao removidos.
# Na sobreposicao: bases iguais ficam; N de um lado e preenchido pelo outro;
# discordancias viram IUPAC e sao listadas para revisar no cromatograma.
import sys
IUPAC={frozenset('AG'):'R',frozenset('CT'):'Y',frozenset('GC'):'S',frozenset('AT'):'W',frozenset('GT'):'K',frozenset('AC'):'M'}
COMP=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
def rc(s): return s.translate(COMP)[::-1]
def ler(p):
    nomes=[l[1:].strip() for l in open(p,encoding='utf-8',errors='replace') if l.startswith('>')]
    if len(nomes)>1:
        sys.exit(f"ERRO: {p} tem {len(nomes)} sequencias ({', '.join(nomes[:4])}). Ele deve ter so o CONSENSO da metade.\n"
                 "  Se for o alinhamento dos reads exportado do MEGA, gere o consenso antes:\n"
                 "  python3 Scripts/consenso_mega.py <alinhamento.fasta> <NOME> <consenso.fasta>")
    seq=''.join(l.strip().upper() for l in open(p,encoding='utf-8',errors='replace') if not l.startswith('>'))
    seq=seq.replace('-','').replace('?','').replace('.','').replace(' ','')
    if not seq: sys.exit(f'ERRO: {p} esta vazio')
    return seq
def sobrepor(a,b,M=2,X=-3,G=-5):
    # sufixo de a alinhado com prefixo de b (inicio de a livre, fim de b livre)
    n,m=len(a),len(b)
    S=[[0]*(m+1) for _ in range(n+1)]; T=[[0]*(m+1) for _ in range(n+1)]
    for j in range(1,m+1): S[0][j]=-10**9
    for i in range(1,n+1):
        for j in range(1,m+1):
            d=S[i-1][j-1]+(M if a[i-1]==b[j-1] else X); u=S[i-1][j]+G; l=S[i][j-1]+G
            S[i][j],T[i][j]=max((d,0),(u,1),(l,2))
    j=max(range(m+1),key=lambda k:S[n][k]); i=n; A=[];B=[]
    jfim=j
    while i>0 and j>0:
        t=T[i][j]
        if t==0: i-=1;j-=1;A.append(a[i]);B.append(b[j])
        elif t==1: i-=1;A.append(a[i]);B.append('-')
        else: j-=1;A.append('-');B.append(b[j])
    return i,A[::-1],B[::-1],jfim,S[n][jfim]
LACUNA=None
if '--lacuna' in sys.argv:
    k=sys.argv.index('--lacuna'); LACUNA=int(sys.argv[k+1]); del sys.argv[k:k+2]
a5,a3,nome,out=sys.argv[1:5]
a,b=ler(a5),ler(a3)
if LACUNA is not None:
    seq=a+'N'*LACUNA+b
    open(out,'w').write(f">{nome}\n"+'\n'.join(seq[i:i+70] for i in range(0,len(seq),70))+'\n')
    print(f"{nome}\tmetade5={len(a)}\tmetade3={len(b)}\tlacuna={LACUNA} N\tconsenso={len(seq)}\tN/IUPAC={sum(c not in 'ACGT' for c in seq)}")
    sys.exit()
# testa a ordem e o sentido das metades; fica com a combinacao de melhor sobreposicao
def ident(A,B): return sum(x==y for x,y in zip(A,B))/max(1,len(A))
combos=[(a,b,''),(b,a,"ordem trocada (a 2a metade era a do ITS1)"),(a,rc(b),'2a metade estava invertida'),(rc(b),a,'ordem trocada e 2a metade invertida')]
melhor=None
for x,y,obs in combos:
    r=sobrepor(x,y)
    if len(r[1])>=15 and ident(r[1],r[2])>=0.85 and (melhor is None or r[4]>melhor[1][4]): melhor=((x,y,obs),r)
if melhor is None:
    sys.exit("ERRO: as duas metades nao se sobrepoem (nenhuma combinacao de ordem/sentido teve >=15 bases parecidas).\n"
             "  Confira: (1) os dois arquivos sao da MESMA amostra e sao as metades ITS1F-5.8R e 5.8F-ITS2R;\n"
             "  (2) as pontas no 5.8S nao foram aparadas demais no MEGA (a sobreposicao fica no 5.8S, ~20-30 pb);\n"
             "  (3) se faltar um pedaco de verdade, use --lacuna N (manual, secao 6.4);\n"
             "  (4) se voce conferiu no MEGA que as metades se encostam sem sobrepor, use --lacuna 0.")
(a,b,obs),(i0,A,B,jfim,score)=melhor
if len(A)>0.8*min(len(a),len(b)):
    sys.exit(f"ERRO: as duas sequencias se sobrepoem em {len(A)} pb - parecem ser a MESMA metade (ou o ITS inteiro).\n"
             "  Use a metade ITS1F-5.8R e a metade 5.8F-ITS2R da mesma amostra.")
if obs: print(f"    obs: {obs} - corrigido automaticamente")
cons=list(a[:i0]); notas=[]; disc=0
for pa,(x,y) in enumerate(zip(A,B),1):  # pa = posicao dentro da sobreposicao
    if x==y: cons.append(x)
    elif x=='-': cons.append(y); notas.append(f"pos {len(cons)}: base so na metade 3' ({y}) REVISAR")
    elif y=='-': cons.append(x); notas.append(f"pos {len(cons)}: base so na metade 5' ({x}) REVISAR")
    elif x=='N': cons.append(y)
    elif y=='N': cons.append(x)
    else:
        disc+=1
        if pa>len(A)-5: cons.append(y); notas.append(f"pos {len(cons)}: 5'={x} 3'={y} -> {y} (ponta da metade 5', vale a 3')")
        elif pa<=5: cons.append(x); notas.append(f"pos {len(cons)}: 5'={x} 3'={y} -> {x} (ponta da metade 3', vale a 5')")
        else:
            k=IUPAC.get(frozenset(x+y),'N'); cons.append(k); notas.append(f"pos {len(cons)}: 5'={x} 3'={y} -> {k} REVISAR")
cons+=list(b[jfim:])
seq=''.join(cons)
ov=len(A)
if 'GCTGCGTTCTTC' in seq and 'GAAGAACGCAGC' not in seq:   # 5.8S ao contrario -> ITS inteiro invertido
    seq=rc(seq); print("    obs: ITS estava no sentido reverso (5.8S invertido) - gravado no sentido 18S->26S")
elif 'GAAGAACGCAGC' not in seq: print("    AVISO: trecho do 5.8S (GAAGAACGCAGC) nao encontrado - conferir o sentido no MEGA")
open(out,'w').write(f">{nome}\n"+'\n'.join(seq[i:i+70] for i in range(0,len(seq),70))+'\n')
print(f"{nome}\tmetade5={len(a)}\tmetade3={len(b)}\tsobreposicao={ov}\tdiscord={disc}\tconsenso={len(seq)}\tN/IUPAC={sum(c not in 'ACGT' for c in seq)}")
fim_sobrep=len(seq)-(len(b)-jfim)
if ov<30: print(f"    obs: sobreposicao curta (sitio do primer no 5.8S) - conferir a juncao: ...{seq[fim_sobrep-ov-10:fim_sobrep+10]}...")
for x in notas: print("   ",x)
