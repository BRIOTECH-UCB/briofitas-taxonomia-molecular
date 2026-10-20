# Matriz concatenada (supermatriz) de dois ou mais marcadores, com partições para o IQ-TREE.
#
# Passo 1 - preparar: gera um FASTA de entrada por marcador, com os MESMOS rótulos em todos.
#   python3 Scripts/concatenado.py preparar <lista.tsv> <pasta_saida> ITS=Referencias/BLAST_ITS/bryophyta_ITS.fasta trnLF=Referencias/BLAST_trnLF/bryophyta_trnLF.fasta
#   lista.tsv: rotulo<TAB>marcador1<TAB>marcador2 (acesso, caminho de FASTA local, ou -), na ordem dos bancos dados.
#
# Passo 2 - (alinhar e aparar cada marcador separadamente: MAFFT + trimAl)
#
# Passo 3 - juntar: concatena os alinhamentos aparados; espécime sem um marcador recebe '?' nesse bloco.
#   python3 Scripts/concatenado.py juntar <saida_prefixo> ITS=ITS_aparado.fasta trnLF=trnLF_aparado.fasta
#   Gera <prefixo>.fasta (supermatriz), <prefixo>_particoes.nex e <prefixo>_cobertura.tsv.
import sys, os

def ler_fasta(p):
    d = {}; o = []; k = None
    for l in open(p, encoding='utf-8', errors='replace'):
        l = l.strip()
        if not l: continue
        if l.startswith('>'):
            k = l[1:].split()[0]
            if k.startswith('_R_'): k = k[3:]  # MAFFT --adjustdirection
            d[k] = []; o.append(k)
        elif k: d[k].append(l.upper())
    return {k: ''.join(v) for k, v in d.items()}, o

def preparar(lista, pasta, bancos):
    nomes = [b.split('=', 1)[0] for b in bancos]
    caminhos = [b.split('=', 1)[1] for b in bancos]
    linhas = []
    for l in open(lista, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#'): continue
        c = l.split('\t')
        linhas.append((c[0], c[1:1 + len(nomes)]))
    os.makedirs(pasta, exist_ok=True)
    for i, (m, banco) in enumerate(zip(nomes, caminhos)):
        pedidos = {fontes[i]: rot for rot, fontes in linhas if fontes[i] != '-' and not os.path.exists(fontes[i])}
        achados = {}
        if pedidos:
            cur = None
            for l in open(banco, encoding='utf-8', errors='replace'):
                if l.startswith('>'):
                    a = l[1:].split()[0]; cur = a if a in pedidos else None
                    if cur: achados[cur] = []
                elif cur: achados[cur].append(l.strip())
        faltam = [a for a in pedidos if a not in achados]
        if faltam: sys.exit(f'{m}: não encontrados no banco: ' + ', '.join(faltam))
        n = 0
        with open(os.path.join(pasta, f'{m}_entrada.fasta'), 'w') as out:
            for rot, fontes in linhas:
                f = fontes[i]
                if f == '-': continue
                if os.path.exists(f):
                    d, o = ler_fasta(f); s = d[o[0]]
                else:
                    s = ''.join(achados[f])
                out.write(f'>{rot}\n' + '\n'.join(s[j:j + 70] for j in range(0, len(s), 70)) + '\n'); n += 1
        print(f'{m}: {n} sequências -> {pasta}/{m}_entrada.fasta')
    print(f'Espécimes na lista: {len(linhas)}')

def juntar(pref, alns):
    blocos = []
    for a in alns:
        m, p = a.split('=', 1)
        d, o = ler_fasta(p)
        L = {len(s) for s in d.values()}
        if len(L) != 1: sys.exit(f'{p}: sequências com comprimentos diferentes (não é alinhamento?)')
        blocos.append((m, d, o, L.pop()))
    todos = []
    for _, d, o, _ in blocos:
        for k in o:
            if k not in todos: todos.append(k)
    with open(pref + '.fasta', 'w') as out:
        for k in todos:
            s = ''.join(d.get(k, '?' * L) for _, d, _, L in blocos)
            out.write(f'>{k}\n' + '\n'.join(s[j:j + 70] for j in range(0, len(s), 70)) + '\n')
    ini = 1
    with open(pref + '_particoes.nex', 'w') as nx:
        nx.write('#nexus\nbegin sets;\n')
        for m, _, _, L in blocos:
            nx.write(f'  charset {m} = {ini}-{ini + L - 1};\n'); ini += L
        nx.write('end;\n')
    with open(pref + '_cobertura.tsv', 'w') as cb:
        cb.write('rotulo\t' + '\t'.join(m for m, *_ in blocos) + '\n')
        for k in todos:
            cb.write(k + '\t' + '\t'.join('sim' if k in d else '-' for _, d, _, _ in blocos) + '\n')
    print(f'{len(todos)} espécimes, {ini - 1} colunas -> {pref}.fasta')
    for m, d, _, L in blocos:
        print(f'  {m}: {L} colunas, {len(d)} espécimes com dado, {len(todos) - len(d)} faltando (?)')
    print(f'Partições: {pref}_particoes.nex | cobertura: {pref}_cobertura.tsv')

if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('preparar', 'juntar'):
        sys.exit(__doc__ if __doc__ else 'uso: preparar | juntar')
    if sys.argv[1] == 'preparar': preparar(sys.argv[2], sys.argv[3], sys.argv[4:])
    else: juntar(sys.argv[2], sys.argv[3:])
