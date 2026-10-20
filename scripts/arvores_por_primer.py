# Faz uma árvore para CADA PAR DE PRIMERS (ITS1-ITS4, 18F-25R, FM-CM, rpl16Fb-rpl16R ...), do começo ao fim:
# monta a entrada, alinha (MAFFT), apara (trimAl), roda o IQ-TREE 3 e desenha a figura.
# Uso:  python3 Scripts/arvores_por_primer.py --genero Bryum
#       python3 Scripts/arvores_por_primer.py --genero Bryum --pares ITS1-ITS4,FM-CM     (só alguns pares)
#       python3 Scripts/arvores_por_primer.py --genero Bryum --so-preparar              (só monta as entradas)
#
# --genero = o nome da análise (a pasta Referencias/<genero>, criada pelo selecionar_referencias.py).
# Para cada par:
#   - referências do GenBank: as da lista_concatenado.tsv, do marcador do par (ITS para ITS1-ITS4 e 18F-25R, trnLF para FM-CM ...)
#   - amostras: o consenso DAQUELE par, Consensos/por_primer/<amostra>_<par>.fasta (montar_todos.py ou corrigido no MEGA).
#     Amostra sem esse par, ou com reads sem sobreposição, fica de fora (é avisado).
# Resultado em Arvores/<genero>/primer_<par>/ : <par>.treefile, <par>.iqtree e a figura <genero>_<par>.svg
# Pares usados: os "completo" do Scripts/primers.tsv (ou os de --pares).
import sys, os, re
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from ler_primers import ler_primers, achar_par
from concatenado import preparar
import rodar_arvore as ra

COLUNAS = ['ITS', 'trnLF', 'rpl16']   # ordem das colunas de marcador na lista_concatenado.tsv

def situacao(a, par):
    """(caminho, texto) do consenso da amostra para o par; caminho None = não usar."""
    f = f'Consensos/por_primer/{a}_{par}.fasta'
    if not os.path.exists(f): return None, 'sem consenso deste par'
    log = f.replace('.fasta', '_log.txt')
    if not os.path.exists(log): return f, 'consenso sem relatório (conferir)'
    t = open(log, encoding='utf-8', errors='replace').read()
    g = lambda k: int((re.search(k + r'=(\d+)', t) or [0, '0'])[1])
    if 'READ_UNICO' in t: return f, f'read único, {g("consenso")} pb (conferir no MEGA)'
    if g('sobreposicao') < 30: return None, 'reads sem sobreposição'
    return f, f'F+R, {g("consenso")} pb'

def main():
    a = sys.argv[1:]
    op = {'--so-preparar': '--so-preparar' in a}
    a = [x for x in a if x != '--so-preparar']
    while a:
        k = a.pop(0)
        if k in ('--genero', '--pares') and a: op[k] = a.pop(0)
        else: sys.exit('Uso: python3 Scripts/arvores_por_primer.py --genero GENERO [--pares ITS1-ITS4,FM-CM] [--so-preparar]')
    if '--genero' not in op: sys.exit('Falta --genero (ex.: --genero Bryum)')
    G = op['--genero']
    todos = ler_primers()
    if '--pares' in op:
        pares = []
        for n in op['--pares'].split(','):
            p = achar_par(n, todos)
            if not p: sys.exit(f'Par {n} não está no Scripts/primers.tsv')
            pares.append(p)
    else: pares = [p for p in todos if p['papel'] == 'completo']
    lista = f'Referencias/{G}/lista_concatenado.tsv'
    if not os.path.exists(lista):
        sys.exit(f'Falta {lista}. Rode antes o selecionar_referencias.py com --genero {G} (manual, passo 9).')
    if not op['--so-preparar']: ra.conferir_programas()
    ra.registrar(f'Arvores/{G}/registro_arvores_por_primer.txt')
    linhas = []
    for l in open(lista, encoding='utf-8'):
        if l.startswith('#') or not l.strip(): continue
        c = l.rstrip('\n').split('\t')
        linhas.append((c[0], c[1:]))
    feitos = []
    for p in pares:
        par, m = p['par'], p['marcador']
        print(f'\n== {par} (marcador {m})')
        if m not in COLUNAS: print(f'   marcador {m} desconhecido; pulei'); continue
        i = COLUNAS.index(m)
        banco = f'Referencias/BLAST_{m}/bryophyta_{m}.fasta'
        saida, n_ref, n_am, vistos = [], 0, 0, set()
        for rot, fontes in linhas:
            if rot in vistos: continue
            fonte = fontes[i] if i < len(fontes) else '-'
            mm = re.match(r'(?:OUT_)?AMOSTRA_(.+)$', rot)
            if mm:   # amostra: usa o consenso deste par
                am = mm.group(1)
                f, txt = situacao(am, par)
                print(f'   amostra {am:10s} {"entra" if f else "FORA "} ({txt})')
                if not f: continue
                saida.append((rot, f)); n_am += 1
            else:
                if fonte == '-': continue
                saida.append((rot, fonte)); n_ref += 1
            vistos.add(rot)
        if n_am == 0: print('   nenhuma amostra tem este par - árvore não feita'); continue
        if n_ref == 0: print(f'   a lista não tem referências de {m} - árvore não feita'); continue
        if not any(r.startswith('OUT_') for r, _ in saida): print('   AVISO: nenhum grupo externo (OUT_) com este marcador')
        if not os.path.exists(banco) and any(not os.path.exists(f) for _, f in saida):
            sys.exit(f'Falta o banco {banco}: rode python3 Scripts/blast_local.py --marcador {m}')
        pasta = f'Arvores/{G}/primer_{par}'
        os.makedirs(pasta, exist_ok=True)
        lp = f'{pasta}/lista_{par}.tsv'
        with open(lp, 'w', encoding='utf-8') as o:
            o.write(f'# Gerada por arvores_por_primer.py a partir de {lista}: referências de {m} + amostras com o par {par}\n')
            for rot, f in saida: o.write(f'{rot}\t{f}\n')
        print(f'   {n_ref} referências + {n_am} amostra(s)')
        preparar(lp, pasta, [f'{par}={banco}'])
        if op['--so-preparar']: feitos.append(par); continue
        base = f'{pasta}/{par}'
        print('   alinhando (MAFFT) e aparando (trimAl)...')
        ap, inv = ra.alinhar_aparar(f'{base}_entrada.fasta', base)
        if inv: print(f'   {inv} sequência(s) estavam invertidas e foram viradas pelo MAFFT (normal)')
        print('   rodando o IQ-TREE 3 (alguns minutos)...')
        tree = ra.iqtree(ap, base)
        svg = ra.figura(tree, f'{pasta}/{G}_{par}.svg', f'{G} — {par} (máxima verossimilhança, IQ-TREE 3)')
        print(f'   PRONTO: árvore {tree}' + (f' | figura {svg}' if svg else ''))
        feitos.append(par)
    if not feitos: sys.exit('\nNenhuma árvore feita.')
    print(f'\nÁrvores por primer de {G}: ' + ', '.join(feitos) + f'  (pasta Arvores/{G}/)')
    if not op['--so-preparar']: print('Abra os .svg no navegador do Windows (dois cliques).')

if __name__ == '__main__':
    main()
