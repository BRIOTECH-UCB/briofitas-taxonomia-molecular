# QC de leituras Sanger (.ab1): comprimento, QV medio, bases >=Q20 e comprimento apos aparo (Mott, limite 0.01).
# Uso: python Scripts/qc_ab1.py <pasta_do_projeto> > qc.csv   (varre todas as pastas Organizado_por_amostra)
import struct, sys, os, glob, csv
def read_abif(p):
    d=open(p,'rb').read()
    if d[:4]!=b'ABIF': return {}
    num=struct.unpack('>I',d[18:22])[0]; off=struct.unpack('>I',d[26:30])[0]
    tags={}
    for i in range(num):
        e=d[off+i*28:off+(i+1)*28]
        name,n,et,es,ne,ds,do,_=struct.unpack('>4sIHHIIII',e)
        data=e[20:24] if ds<=4 else d[do:do+ds]
        tags[(name.decode('latin1'),n)]=data[:ds]
    return tags
def mott(q,cut=0.01):
    best=0;s=0;bs=be=0;start=0
    for i,x in enumerate(q):
        s+= cut-10**(-x/10)
        if s<0: s=0; start=i+1
        if s>best: best=s; bs=start; be=i+1
    return bs,be
rows=[]
root=sys.argv[1]
for p in sorted(glob.glob(root+'/**/Organizado_por_amostra/**/*.ab1',recursive=True)):
    t=read_abif(p)
    seq=t.get(('PBAS',2),b'').decode('latin1'); q=list(t.get(('PCON',2),b''))
    if not q: continue
    a,b=mott(q)
    hq=sum(1 for x in q if x>=20)
    rows.append([os.path.relpath(p,root).replace(os.sep,'/'),len(seq),round(sum(q)/len(q),1),hq,b-a, round(sum(q[a:b])/max(1,b-a),1)])
w=csv.writer(sys.stdout,lineterminator='\n'); w.writerow(['arquivo','len','QVmedio','bases_Q20','len_aparado','QV_aparado']); w.writerows(rows)
