# Prepara uma árvore do IQ-TREE para virar figura no FigTree (sem precisar mexer em Python).
# - enraíza no grupo externo (um ou vários táxons, pelo início do nome)
# - ordena os clados (menores em cima)
# - mantém nos nós só os suportes que interessam
#
# Na 1a vez cria <saida>_nomes.tsv com nomes limpos sugeridos para a figura; edite e rode de novo.
# Uso:  python3 Scripts/preparar_figura.py <arvore.treefile> <inicio_do_nome_do_outgroup>[,outro,...] <saida.tree> [--mostrar fortes|qualquer|todos]
#   fortes   = só nós com SH-aLRT >= 80 E UFBoot >= 95
#   qualquer = nós em que pelo menos um dos dois passa do limite (padrão)
#   todos    = todos os valores, menos os de ramos de comprimento zero (SH-aLRT = 0)
# Exemplo: python3 Scripts/preparar_figura.py Arvores/Bryum_concat/Bryum_concat.treefile OUT_ Arvores/Bryum_concat/Bryum_concat_figura.tree
import sys, re, os
from collections import defaultdict

LIM_SH, LIM_UF = 80.0, 95.0

def ler_newick(s):
    s = s.strip().rstrip(';')
    adj = defaultdict(list); nomes = {}; cont = [0]; pos = [0]
    def novo():
        cont[0] += 1; return cont[0]
    def no():
        n = novo()
        if s[pos[0]] == '(':
            pos[0] += 1
            while True:
                f, lab, bl = no()
                adj[n].append([f, lab, bl]); adj[f].append([n, lab, bl])
                if s[pos[0]] == ',': pos[0] += 1; continue
                if s[pos[0]] == ')': pos[0] += 1; break
        m = re.match(r"('[^']*'|[^:,();]*)(?::([0-9.eE+-]+))?", s[pos[0]:])
        pos[0] += m.end()
        rot = m.group(1).strip("'"); bl = float(m.group(2)) if m.group(2) else 0.0
        if len(adj[n]) == 0: nomes[n] = rot; rot = ''
        return n, rot, bl
    no()
    return adj, nomes

def folhas_do_lado(adj, nomes, u, pai):
    pilha = [(u, pai)]; out = set()
    while pilha:
        a, p = pilha.pop()
        if a in nomes: out.add(a)
        for v, _, _ in adj[a]:
            if v != p: pilha.append((v, a))
    return out

def sugerir_nome(r):
    # 'P_pseudotriquetrum_Spse65_MK270162_Antarctica' -> 'P. pseudotriquetrum Spse65 MK270162 Antarctica'
    t = [x for x in r.split('_') if x]
    if t and t[0] == 'OUT': t = t[1:]
    if not t: return r
    if t[0] == 'AMOSTRA': return 'AMOSTRA ' + ' '.join(t[1:])
    if len(t) >= 2 and t[1][:1].islower():
        gen = t[0] + '.' if len(t[0]) == 1 else t[0]
        return ' '.join([gen, t[1]] + t[2:])
    return ' '.join(t)

def filtrar(rot, modo):
    if not rot: return ''
    m = re.match(r'^([0-9.]+)/([0-9.]+)$', rot)
    if not m: return rot
    sh, uf = float(m.group(1)), float(m.group(2))
    fmt = lambda x: f'{x:g}'
    if modo == 'fortes': ok = sh >= LIM_SH and uf >= LIM_UF
    elif modo == 'qualquer': ok = sh >= LIM_SH or uf >= LIM_UF
    else: ok = sh > 0
    return f'{fmt(sh)}/{fmt(uf)}' if ok else ''

def main():
    args = sys.argv[1:]
    modo = 'qualquer'
    if '--mostrar' in args:
        i = args.index('--mostrar'); modo = args[i + 1]; del args[i:i + 2]
    if modo not in ('fortes', 'qualquer', 'todos'): sys.exit('--mostrar deve ser fortes, qualquer ou todos')
    entrada, og, saida = args
    adj, nomes = ler_newick(open(entrada).read())
    prefixos = [p for p in og.split(',') if p]
    alvo = {n for n, r in nomes.items() if any(r.startswith(p) for p in prefixos)}
    if not alvo: sys.exit(f'Nenhum táxon começa com: {og}')
    todas = set(nomes)
    # procura o ramo que separa exatamente o grupo externo do resto
    raiz_ramo = None
    for u in adj:
        for v, lab, bl in adj[u]:
            lado = folhas_do_lado(adj, nomes, v, u)
            if lado == alvo or (todas - lado) == alvo:
                raiz_ramo = (u, v, lab, bl); break
        if raiz_ramo: break
    if not raiz_ramo:
        print('AVISO: o grupo externo não forma um grupo único na árvore; enraizei no primeiro táxon dele.')
        f = sorted(alvo)[0]; v, lab, bl = adj[f][0]; raiz_ramo = (f, v, lab, bl)
    u, v, lab, bl = raiz_ramo
    contagem = {}
    def nfolhas(a, p):
        if (a, p) in contagem: return contagem[(a, p)]
        n = 1 if a in nomes else sum(nfolhas(b, a) for b, _, _ in adj[a] if b != p)
        contagem[(a, p)] = n; return n
    # nomes da figura: tabela <saida>_nomes.tsv (criada com sugestões na 1ª vez; edite e rode de novo)
    arq_nomes = os.path.splitext(saida)[0] + '_nomes.tsv'
    if os.path.exists(arq_nomes):
        troca = {}
        for l in open(arq_nomes, encoding='utf-8'):
            c = l.rstrip('\n').split('\t')
            if len(c) >= 2 and not l.startswith('#') and c[1].strip(): troca[c[0]] = c[1].strip()
        print(f'Nomes da figura lidos de {arq_nomes}')
    else:
        troca = {r: sugerir_nome(r) for r in sorted(nomes.values())}
        with open(arq_nomes, 'w', encoding='utf-8') as f:
            f.write('# nome_na_arvore\tnome_na_figura   (edite a 2a coluna no Notepad++ e rode o comando de novo)\n')
            for r, s in troca.items(): f.write(f'{r}\t{s}\n')
        print(f'Criei {arq_nomes} com nomes sugeridos para a figura.')
    def escrever(a, p, rotulo, comp):
        if a in nomes:
            nm = troca.get(nomes[a], nomes[a])
            txt = "'" + nm.replace("'", "''") + "'"
        else:
            filhos = sorted([(b, l, c) for b, l, c in adj[a] if b != p], key=lambda x: nfolhas(x[0], a))
            txt = '(' + ','.join(escrever(b, a, l, c) for b, l, c in filhos) + ')' + filtrar(rotulo, modo)
        return f'{txt}:{comp:.6g}'
    # o lado do grupo externo fica sem suporte; o suporte do ramo-raiz vai para o grupo interno
    lado_u_e_out = folhas_do_lado(adj, nomes, u, v) == alvo
    ext, intg = (u, v) if lado_u_e_out else (v, u)
    s_ext = escrever(ext, intg, '', bl / 2)
    s_int = escrever(intg, ext, lab, bl / 2)
    arvore = f'({s_int},{s_ext});'
    open(saida, 'w').write(arvore + '\n')
    mostrados = len(re.findall(r"\)[0-9]", arvore))
    print(f'Árvore enraizada em {len(alvo)} táxon(s) do grupo externo; {len(nomes)} terminais; {mostrados} nós com suporte mostrado (modo: {modo}).')
    print(f'Árvore pronta: {saida}  -> desenhe com: python3 Scripts/desenhar_arvore.py {saida} <saida.svg> "<título>"')

if __name__ == '__main__':
    main()
