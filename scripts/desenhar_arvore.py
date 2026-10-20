# Desenha uma árvore pronta (saída do preparar_figura.py) como figura SVG, sem precisar do FigTree.
# O SVG abre no navegador e no Inkscape (para ajustes finais e exportar PNG/PDF).
#
# Uso:  python3 Scripts/desenhar_arvore.py <arvore_figura.tree> <saida.svg> "<título>" [--nota "texto de rodapé"]
# Exemplo: python3 Scripts/desenhar_arvore.py Arvores/Sanionia_concat/Sanionia_concat_figura.tree Arvores/Sanionia_concat/Sanionia_concat.svg "Sanionia — ITS + trnL-F"
#
# Convenções da figura:
#   - nome da espécie em itálico; voucher, acesso e país em cinza
#   - amostra deste projeto (nome começando com AMOSTRA) em vermelho e negrito
#   - material antártico/subantártico em azul
#   - grupo externo em cinza (lido do <arvore>_nomes.tsv: linhas cujo nome original começa com OUT_)
#   - suportes SH-aLRT/UFBoot nos nós; em negrito quando fortes (SH-aLRT >= 80 e UFBoot >= 95)
import sys, re, html

ANTARTICO = re.compile(r'Antarc|Shetland|KingGeorge|King George|Marion|SouthGeorgia|South Georgia|Orkney|Sandwich|Victoria|'
                       r'Melbourne|Rittmann|Kerguelen|Falkland|Malvinas|Fuego|Macquarie|Heard|Crozet|Signy|Livingston|Deception', re.I)

def ler_newick(s):
    s = s.strip().rstrip(';')
    i = 0
    def nome():
        nonlocal i
        if i < len(s) and s[i] == "'":
            j = i + 1; out = ''
            while True:
                if s[j] == "'" and j + 1 < len(s) and s[j + 1] == "'": out += "'"; j += 2; continue
                if s[j] == "'": break
                out += s[j]; j += 1
            i = j + 1; return out
        m = re.match(r"[^:,();]*", s[i:]); i += m.end(); return m.group(0)
    def comp():
        nonlocal i
        if i < len(s) and s[i] == ':':
            m = re.match(r":([-0-9.eE+]+)", s[i:]); i += m.end(); return float(m.group(1))
        return 0.0
    def no():
        nonlocal i
        n = {'filhos': [], 'nome': '', 'comp': 0.0}
        if s[i] == '(':
            i += 1
            while True:
                n['filhos'].append(no())
                if s[i] == ',': i += 1; continue
                i += 1; break
        n['nome'] = nome(); n['comp'] = comp()
        return n
    return no()

def folhas(n):
    return [n] if not n['filhos'] else [f for c in n['filhos'] for f in folhas(c)]

def main():
    a = sys.argv[1:]
    if len(a) < 3: sys.exit('Uso: desenhar_arvore.py arvore_figura.tree saida.svg "titulo" [--nota "texto"]')
    entrada, saida, titulo = a[0], a[1], a[2]
    nota = a[a.index('--nota') + 1] if '--nota' in a else ''
    import os
    og = set()
    arq_nomes = os.path.splitext(entrada)[0] + '_nomes.tsv'
    if os.path.exists(arq_nomes):
        for l in open(arq_nomes, encoding='utf-8'):
            c = l.rstrip('\n').split('\t')
            if len(c) == 2 and c[0].startswith('OUT_'): og.add(c[1])
    raiz = ler_newick(open(entrada, encoding='utf-8').read())

    # ramo muito longo saindo da raiz (em geral o grupo externo) achata a figura: encurta e marca com //
    def prof(n): return max([prof(c) + c['comp'] for c in n['filhos']] or [0.0])
    resto = max(prof(c) for c in raiz['filhos'])
    for c in raiz['filhos']:
        if c['comp'] > 0.5 * resto:
            c['cortado'] = c['comp']; c['comp'] = 0.12 * resto
    # posições: x = distância da raiz; y = ordem das folhas
    def px(n, x0):
        n['x'] = x0 + n['comp'] if n is not raiz else 0.0
        for c in n['filhos']: px(c, n['x'])
    px(raiz, 0.0)
    fl = folhas(raiz)
    for k, f in enumerate(fl): f['y'] = k
    def py(n):
        if n['filhos']:
            for c in n['filhos']: py(c)
            n['y'] = (n['filhos'][0]['y'] + n['filhos'][-1]['y']) / 2
    py(raiz)
    xmax = max(f['x'] for f in fl) or 1.0

    LH, ML, MT, TW, FS = 16, 30, 70, 680, 11
    def X(v): return ML + v / xmax * TW
    def Y(v): return MT + v * LH
    largura_nome = max(len(f['nome']) for f in fl) * 6.4 + 40
    W = ML + TW + largura_nome + 20
    H = MT + len(fl) * LH + 80
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}" '
           f'font-family="Arial, Helvetica, sans-serif">',
           f'<rect width="100%" height="100%" fill="white"/>',
           f'<text x="{ML}" y="28" font-size="16" font-weight="bold">{html.escape(titulo)}</text>']

    def forte(lab):
        m = re.match(r'([0-9.]+)/([0-9.]+)', lab or '')
        return bool(m) and float(m.group(1)) >= 80 and float(m.group(2)) >= 95

    suportes = []  # desenhados por último, por cima dos ramos
    def ramos(n):
        for c in n['filhos']:
            out.append(f'<line x1="{X(n["x"]):.1f}" y1="{Y(c["y"]):.1f}" x2="{X(c["x"]):.1f}" y2="{Y(c["y"]):.1f}" stroke="black" stroke-width="1.2"/>')
            if c.get('cortado'):
                xm = (X(n['x']) + X(c['x'])) / 2; yc = Y(c['y'])
                for dx in (-2, 2):
                    out.append(f'<line x1="{xm + dx - 3:.1f}" y1="{yc + 5:.1f}" x2="{xm + dx + 3:.1f}" y2="{yc - 5:.1f}" stroke="black" stroke-width="1.2"/>')
                suportes.append(f'<text x="{xm:.1f}" y="{yc + 13:.1f}" font-size="8" text-anchor="middle" fill="#666">ramo encurtado ({c["cortado"]:.2g})</text>')
            ramos(c)
        if n['filhos']:
            out.append(f'<line x1="{X(n["x"]):.1f}" y1="{Y(n["filhos"][0]["y"]):.1f}" x2="{X(n["x"]):.1f}" y2="{Y(n["filhos"][-1]["y"]):.1f}" stroke="black" stroke-width="1.2"/>')
            lab = n['nome'].strip()
            if lab and n is not raiz and re.match(r'[0-9.]+/[0-9.]+', lab):
                peso = 'bold' if forte(lab) else 'normal'
                cor = 'black' if forte(lab) else '#666'
                suportes.append(f'<text x="{X(n["x"]) - 3:.1f}" y="{Y(n["y"]) - 3:.1f}" font-size="8.5" text-anchor="end" stroke="white" stroke-width="3" paint-order="stroke" '
                           f'font-weight="{peso}" fill="{cor}">{html.escape(lab)}</text>')
    out.append(f'<line x1="{X(0) - 8:.1f}" y1="{Y(raiz["y"]):.1f}" x2="{X(0):.1f}" y2="{Y(raiz["y"]):.1f}" stroke="black" stroke-width="1.2"/>')
    ramos(raiz)
    out.extend(suportes)

    for f in fl:
        nm = f['nome']; x = X(f['x']) + 5; y = Y(f['y']) + 4
        if nm.upper().startswith('AMOSTRA'):
            txt = nm.split(None, 1)[1] if ' ' in nm else nm
            out.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{FS + 1}" font-weight="bold" fill="#C0392B">'
                       f'&#9668; {html.escape(txt)} (este estudo)</text>')
            continue
        p = nm.split()
        n_tax = 3 if len(p) > 2 and p[1] in ('sp.', 'cf.', 'aff.', 'x') else 2
        tax, resto = ' '.join(p[:n_tax]), ' '.join(p[n_tax:])
        cor = '#888' if nm in og else ('#1F5FAD' if ANTARTICO.search(nm) else 'black')
        cor_resto = cor if cor != 'black' else '#777'
        out.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{FS}" fill="{cor}"><tspan font-style="italic">{html.escape(tax)}</tspan>'
                   f'<tspan fill="{cor_resto}" font-size="{FS - 1.5}"> {html.escape(resto)}</tspan></text>')

    # barra de escala
    alvo = xmax / 5
    pot = 10 ** __import__('math').floor(__import__('math').log10(alvo))
    esc = max(v * pot for v in (1, 2, 5) if v * pot <= alvo)
    yb = Y(len(fl)) + 15
    out.append(f'<line x1="{ML}" y1="{yb}" x2="{X(esc):.1f}" y2="{yb}" stroke="black" stroke-width="1.5"/>')
    out.append(f'<text x="{(ML + X(esc)) / 2:.1f}" y="{yb + 14}" font-size="10" text-anchor="middle">{esc:g} subst./sítio</text>')
    # legenda
    yl = yb + 38
    leg = [('#C0392B', 'amostra deste estudo'), ('#1F5FAD', 'material antártico/subantártico')]
    if og: leg.append(('#888', 'grupo externo'))
    xl = ML
    for cor, t in leg:
        out.append(f'<rect x="{xl}" y="{yl - 9}" width="10" height="10" fill="{cor}"/><text x="{xl + 14}" y="{yl}" font-size="10">{t}</text>')
        xl += 30 + len(t) * 5.6
    out.append(f'<text x="{xl}" y="{yl}" font-size="10">nós: SH-aLRT/UFBoot (negrito = forte, ≥ 80 e ≥ 95)</text>')
    if nota:
        out.append(f'<text x="{ML}" y="{yl + 18}" font-size="10" fill="#444">{html.escape(nota)}</text>')
    out.append('</svg>')
    open(saida, 'w', encoding='utf-8').write('\n'.join(out))
    print(f'Figura salva em {saida} ({len(fl)} terminais). Abra no navegador ou no Inkscape.')

if __name__ == '__main__':
    main()
