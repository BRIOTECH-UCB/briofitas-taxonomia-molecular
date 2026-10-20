# Escolhe as referências do GenBank para a árvore de um gênero, por regras fixas, e gera a lista e a Tabela S1.
# Uso:
#   python3 Scripts/selecionar_referencias.py --genero Sanionia --amostras 179 --outgroup Hamatocaulis
# Opções:
#   --genero G[,G2...]     gênero(s) do grupo interno, como no GenBank (ex.: Warnstorfia,Sarmentypnum,Straminergon)
#   --amostras A[,B...]    nossas amostras desse gênero (códigos das pastas)
#   --outgroup G           gênero do grupo externo (qualquer um, até fora de Bryophyta, ex.: Marchantia;
#                          o que não estiver no banco local é baixado do GenBank)
#   --amostra-outgroup A   (opcional) nossa amostra usada como grupo externo (ex.: 010)
#   --proximas N           quantas referências mais próximas da amostra (padrão 8)
#   --por-especie N        máximo por espécie não dominante (padrão 3)
#   --atualizar            baixa de novo os metadados do GenBank
# Precisa dos bancos locais de BLAST (Referencias/BLAST_ITS e BLAST_trnLF) e, para as "mais próximas",
# do BLAST das amostras (Consensos/<marcador>/BLAST_<marcador>_completo.tsv).
# Grava em Referencias/<genero>/: lista_concatenado.tsv (a anterior vira .bak), tabela_S1.tsv e metadados.
import sys, os, re, csv, time, collections, shutil, urllib.parse, urllib.request

MARC = ['ITS', 'trnLF']   # + 'rpl16' quando o banco de rpl16 existir (main)
BANCO = {'ITS': 'Referencias/BLAST_ITS/bryophyta_ITS.fasta', 'trnLF': 'Referencias/BLAST_trnLF/bryophyta_trnLF.fasta',
         'rpl16': 'Referencias/BLAST_rpl16/bryophyta_rpl16.fasta'}
TITULO_OK = {'ITS': re.compile(r'internal transcribed|ITS', re.I), 'trnLF': re.compile(r'trnL|trnF|tRNA-Leu|tRNA-Phe', re.I),
             'rpl16': re.compile(r'rpl16|ribosomal protein L16', re.I)}
ANT = re.compile(r'Antarc|Shetland|King Georg|Orkney|South Georgia|Bouvet|Crozet|Kerguelen|Marion|Heard|Prince Edward|'
                 r'Macquarie|Victoria Land|Signy|Livingston|Deception', re.I)
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

def args():
    # aceita listas separadas por vírgula OU por espaço: --amostras 010,046  ou  --amostras 010 046
    a = sys.argv[1:]; o = {'proximas': 8, 'por-especie': 3}
    validas = ('genero', 'amostras', 'outgroup', 'amostra-outgroup', 'proximas', 'por-especie', 'atualizar')
    k = None
    for x in a:
        if x.startswith('--'):
            k = x[2:]
            if k not in validas: sys.exit(f'Opção desconhecida: {x}. Opções: ' + ', '.join('--' + v for v in validas))
            if k == 'atualizar': o[k] = True; k = None
            else: o[k] = ''
        elif k is None: sys.exit(f'"{x}" está solto no comando: coloque depois de uma opção (ex.: --amostras 010 046)')
        else: o[k] = (o[k] + ',' + x) if o[k] else x
    for k in ('genero', 'amostras', 'outgroup'):
        if not o.get(k): sys.exit(f'Falta --{k}. Exemplo: python3 Scripts/selecionar_referencias.py --genero Sanionia --amostras 179 --outgroup Hamatocaulis')
    o['genero'] = [g.strip() for g in o['genero'].split(',')]
    o['amostras'] = [x.strip() for x in o['amostras'].split(',') if x.strip()]
    o['proximas'] = int(o['proximas']); o['por-especie'] = int(o['por-especie'])
    return o

def acessos(banco, generos):
    out = {}; cur = None
    for l in open(banco, encoding='utf-8', errors='replace'):
        if l.startswith('>'):
            p = l[1:].split()
            cur = p[0] if len(p) > 1 and p[1] in generos else None
            if cur: out[cur] = [l[1:].strip(), 0]
        elif cur: out[cur][1] += len(l.strip())
    return out  # acesso -> [titulo, comprimento]

def baixar_gb(accs, destino, atualizar):
    if os.path.exists(destino) and not atualizar: return
    accs = sorted(accs); tmp = destino + '.parcial'
    with open(tmp, 'w', encoding='utf-8') as f:
        for i in range(0, len(accs), 200):
            dados = urllib.parse.urlencode({'db': 'nuccore', 'rettype': 'gb', 'retmode': 'text', 'id': ','.join(accs[i:i + 200])}).encode()
            for t in range(5):
                try:
                    f.write(urllib.request.urlopen(EUTILS, dados, timeout=120).read().decode('utf-8', 'replace')); break
                except Exception as e:
                    print('   falha no download, tentando de novo...', e); time.sleep(3 * (t + 1))
            else: sys.exit('Não consegui baixar os metadados do GenBank. Verifique a internet e rode de novo.')
            time.sleep(0.5)
    os.replace(tmp, destino)

CONSULTA = {'ITS': '{g}[Organism] AND ("internal transcribed spacer"[All Fields] OR ITS1[All Fields] OR ITS2[All Fields]) AND 200:5000[SLEN]',
            'trnLF': '{g}[Organism] AND (trnL[All Fields] OR trnF[All Fields]) AND 150:5000[SLEN]',
            'rpl16': '{g}[Organism] AND rpl16[All Fields] AND 200:5000[SLEN]'}

def garantir_no_banco(m, genero):
    """Gênero que não está no banco local (ex.: hepáticas, fora de Bryophyta): baixa as sequências dele do
    GenBank e acrescenta ao FASTA do banco, para o selecionar_referencias e o concatenado.py acharem."""
    if acessos(BANCO[m], {genero}): return
    print(f'{genero} {m}: não está no banco local (que tem só musgos); baixando do GenBank...')
    q = urllib.parse.urlencode({'db': 'nuccore', 'retmax': 300, 'term': CONSULTA[m].format(g=genero)}).encode()
    for t in range(5):
        try:
            ids = re.findall(r'<Id>(\d+)</Id>', urllib.request.urlopen(EUTILS.replace('efetch', 'esearch'), q, timeout=120).read().decode()); break
        except Exception as e:
            print('   falha no download, tentando de novo...', e); time.sleep(3 * (t + 1))
    else: sys.exit('Não consegui consultar o GenBank. Verifique a internet e rode de novo.')
    if not ids: print(f'   nenhuma sequência de {m} de {genero} no GenBank'); return
    dados = urllib.parse.urlencode({'db': 'nuccore', 'rettype': 'fasta', 'retmode': 'text', 'id': ','.join(ids)}).encode()
    for t in range(5):
        try: fa = urllib.request.urlopen(EUTILS, dados, timeout=120).read().decode('utf-8', 'replace'); break
        except Exception as e: print('   falha no download, tentando de novo...', e); time.sleep(3 * (t + 1))
    else: sys.exit('Não consegui baixar as sequências do GenBank. Verifique a internet e rode de novo.')
    with open(BANCO[m], 'a', encoding='utf-8') as f: f.write('\n' + fa.strip() + '\n')
    print(f'   {fa.count(">")} sequências acrescentadas a {BANCO[m]}')

def metadados(gb):
    m = {}
    for r in open(gb, encoding='utf-8', errors='replace').read().split('\n//'):
        v = re.search(r'^VERSION\s+(\S+)', r, re.M)
        if not v: continue
        q = lambda k: re.sub(r'\s+', ' ', (re.search(r'/' + k + r'="([^"]*)"', r) or [None, ''])[1])
        m[v.group(1)] = dict(org=q('organism'), iso=q('isolate'), vou=q('specimen_voucher'), pais=q('country') or q('geo_loc_name'))
    return m

def tokens(md):
    s = (md['vou'] + ' ' + md['iso']).lower()
    return {t for t in re.split(r'[\s,;:()/\[\]<>]+', s) if len(t) >= 3 and re.search(r'\d', t)}

def ident_amostras(amostras):
    ident = collections.defaultdict(float)
    for m in MARC:
        f = f'Consensos/{m}/BLAST_{m}_completo.tsv'
        if not os.path.exists(f): continue
        for l in open(f, encoding='utf-8', errors='replace'):
            c = l.split('\t')
            if len(c) > 4 and c[0] in {f'{a}_{m}' for a in amostras} and int(c[4]) >= 80:
                ident[(m, c[1])] = max(ident[(m, c[1])], float(c[2]))
    return ident

def rotulo(md, acc, pre=''):
    g, *e = (md['org'] or 'sp').split()
    ep = e[1] if len(e) > 1 and e[0] in ('var.', 'subsp.') else (e[0] if e else 'sp')
    v = re.sub(r'[^A-Za-z0-9]', '', (md['vou'] or md['iso']).split()[-1] if (md['vou'] or md['iso']) else '')[:14]
    partes = md['pais'].split(':')
    lugar = partes[1].split(',')[0] if len(partes) > 1 and ANT.search(md['pais']) else partes[0]
    return pre + '_'.join(x for x in [g[0], ep, v, acc.split('.')[0], re.sub(r'[^A-Za-z]', '', lugar)[:16]] if x)

def main():
    o = args()
    G = o['genero'][0]; pasta = f'Referencias/{G}'; os.makedirs(pasta, exist_ok=True)
    R16 = os.path.exists(BANCO['rpl16'])
    if R16 and 'rpl16' not in MARC: MARC.append('rpl16')
    if R16: print('Banco de rpl16 encontrado: a lista terá também a coluna de rpl16.')
    for m in MARC:
        if not os.path.exists(BANCO[m]): sys.exit(f'Falta o banco local de {m}: rode python3 Scripts/blast_local.py --marcador {m}')
    for g in o['genero'] + [o['outgroup']]:   # qualquer gênero, até fora de Bryophyta: o que faltar é baixado do GenBank
        for m in MARC: garantir_no_banco(m, g)
    regs, excl = {}, collections.Counter()
    for grupo, gens in (('in', o['genero']), ('out', [o['outgroup']])):
        for m in MARC:
            acc = acessos(BANCO[m], set(gens))
            if grupo == 'out' and len(acc) > 150:   # grupo externo: bastam alguns; evita baixar milhares de registros
                acc = dict(sorted(acc.items())[:150])
            gb =f'{pasta}/{"_".join(gens)}_{m}.gb'
            print(f'{",".join(gens)} {m}: {len(acc)} sequências no banco local; metadados em {gb}')
            baixar_gb(acc, gb, o.get('atualizar'))
            md = metadados(gb)
            for a, (tit, n) in acc.items():
                d = md.get(a) or dict(org=' '.join(tit.split()[1:3]), iso='', vou='', pais='')
                if ' x ' in d['org'] or 'UNVERIFIED' in tit: excl['híbrido/UNVERIFIED'] += 1; continue
                if not TITULO_OK[m].search(tit): excl[f'{m}: título de outra região'] += 1; continue
                if m == 'ITS' and n < 300: excl['ITS parcial < 300 pb'] += 1; continue
                regs[(m, a)] = dict(d, acc=a, m=m, grupo=grupo, n=n)
    # pares ITS-trnLF pelo mesmo voucher (mesma espécie + mesmo número de voucher/isolado)
    idx = collections.defaultdict(list)
    for (m, a), r in regs.items():
        if m == 'trnLF':
            for t in tokens(r): idx[(r['org'], t)].append(r)
    pares, usados = [], set()
    for (m, a), r in regs.items():
        if m != 'ITS': continue
        for t in tokens(r):
            cand = [x for x in idx.get((r['org'], t), []) if x['acc'] not in usados]
            if cand: usados.add(cand[0]['acc']); pares.append((r, cand[0])); break
    ident = ident_amostras(o['amostras'])
    esp = lambda p: p[0]['org']
    ID = lambda p: max(ident.get(('ITS', p[0]['acc']), 0), ident.get(('trnLF', p[1]['acc']), 0) if p[1] else 0)
    pin = [p for p in pares if p[0]['grupo'] == 'in']; pout = [p for p in pares if p[0]['grupo'] == 'out']
    print(f'Pares ITS + trnL-F do mesmo voucher: {len(pin)} no grupo interno, {len(pout)} no externo')
    escolhidos = collections.OrderedDict()
    def add(p, regra):
        k = (p[0]['acc'] if p[0] else '-', p[1]['acc'] if p[1] else '-')
        if k not in escolhidos: escolhidos[k] = (p, regra)
    poucos_pares = len(pin) < 10
    if poucos_pares:
        print('Poucos pares no GenBank: o ITS e o trnL-F entram separados (árvores separadas, como em Pohlia).')
        singles = [((r, None) if r['m'] == 'ITS' else (None, r)) for r in regs.values() if r['grupo'] == 'in' and r['m'] in ('ITS', 'trnLF')]
        cont = collections.Counter()
        for p in sorted(singles, key=lambda p: -ident.get(((p[0] or p[1])['m'], (p[0] or p[1])['acc']), 0)):
            r = p[0] or p[1]
            if ANT.search(r['pais']): add(p, 'antártico/subantártico'); continue
            if cont[(r['m'], r['org'])] < o['por-especie'] + 1: cont[(r['m'], r['org'])] += 1; add(p, f'até {o["por-especie"] + 1} por espécie')
    else:
        dom = collections.Counter(esp(p) for p in pin).most_common(1)[0][0]
        for p in pin:
            if ANT.search(p[0]['pais']) or ANT.search(p[1]['pais']): add(p, 'antártico/subantártico')
        for p in sorted(pin, key=lambda p: -ID(p))[:o['proximas']]: add(p, 'mais próxima da amostra (BLAST)')
        cont = collections.Counter()
        for p in sorted(pin, key=lambda p: -ID(p)):
            if esp(p) != dom and cont[esp(p)] < o['por-especie']: cont[esp(p)] += 1; add(p, 'espécie irmã')
        paises = set()
        for p in sorted([p for p in pin if esp(p) == dom], key=lambda p: -ID(p)):
            c = p[0]['pais'].split(':')[0]
            if c and c not in paises: paises.add(c); add(p, f'{dom}: 1 por país')
    for p in pout[:2]: add(p, 'grupo externo')
    if not pout:   # nenhum espécime do outgroup com os dois marcadores do mesmo voucher: um de cada marcador, separados
        for r in [r for r in regs.values() if r['grupo'] == 'out' and r['m'] == 'ITS'][:2]: add((r, None), 'grupo externo (só ITS)')
        for r in [r for r in regs.values() if r['grupo'] == 'out' and r['m'] == 'trnLF'][:2]: add((None, r), 'grupo externo (só trnL-F)')
    # rpl16: o do mesmo voucher dos espécimes escolhidos e, à parte, referências só de rpl16 pelas mesmas regras
    r16_de, so_r16 = {}, []
    if R16:
        idx16 = collections.defaultdict(list)
        for (m, a), r in regs.items():
            if m == 'rpl16':
                for t in tokens(r): idx16[(r['org'], t)].append(r)
        usados16 = set()
        for k, (p, regra) in escolhidos.items():
            r = p[0] or p[1]
            for t in tokens(r):
                cand = [x for x in idx16.get((r['org'], t), []) if x['acc'] not in usados16]
                if cand: usados16.add(cand[0]['acc']); r16_de[k] = cand[0]['acc']; break
        livres = [r for (m, a), r in regs.items() if m == 'rpl16' and a not in usados16]
        def poe(r, regra):
            if r['acc'] not in usados16: usados16.add(r['acc']); so_r16.append((r, regra))
        ins = [r for r in livres if r['grupo'] == 'in']
        for r in ins:
            if ANT.search(r['pais']): poe(r, 'rpl16: antártico/subantártico')
        for r in sorted(ins, key=lambda r: -ident.get(('rpl16', r['acc']), 0))[:o['proximas']]: poe(r, 'rpl16: mais próxima da amostra (BLAST)')
        cont16 = collections.Counter()
        for r in sorted(ins, key=lambda r: -ident.get(('rpl16', r['acc']), 0)):
            if cont16[r['org']] < o['por-especie']: cont16[r['org']] += 1; poe(r, 'rpl16: até %d por espécie' % o['por-especie'])
        if not any(k in r16_de for k, (p, _) in escolhidos.items() if (p[0] or p[1])['grupo'] == 'out'):
            for r in [r for r in livres if r['grupo'] == 'out'][:2]: poe(r, 'rpl16: grupo externo')
        print(f'rpl16: {len(r16_de)} dos espécimes escolhidos têm rpl16 do mesmo voucher; {len(so_r16)} referências só de rpl16')
        tem_out16 = any(k in r16_de for k, (p, _) in escolhidos.items() if (p[0] or p[1])['grupo'] == 'out') or \
                    any(r['grupo'] == 'out' for r, _ in so_r16)
        if not tem_out16:
            print(f'ATENÇÃO: o GenBank não tem rpl16 de {o["outgroup"]}: a árvore de rpl16 ficará SEM grupo externo.\n'
                  '         Para ela, use outro --outgroup que tenha rpl16 (ou uma amostra nossa com --amostra-outgroup).')
    # escrever lista e Tabela S1
    lista = f'{pasta}/lista_concatenado.tsv'
    if os.path.exists(lista): shutil.copy2(lista, lista + '.bak')
    t1 = f'{pasta}/tabela_S1.tsv'
    regras = collections.Counter()
    with open(lista, 'w', encoding='utf-8') as f, open(t1, 'w', encoding='utf-8') as s:
        f.write(f'# Lista gerada por selecionar_referencias.py em {time.strftime("%d/%m/%Y")}. 1 linha = 1 espécime (mesmo voucher).\n')
        f.write(f'# Regras: todo material antártico/subantártico; {o["proximas"]} mais próximas da amostra; até {o["por-especie"]} por espécie irmã; '
                f'1 por país na espécie dominante; outgroup {o["outgroup"]}. Excluídos: híbridos, UNVERIFIED, ITS < 300 pb, títulos de outra região.\n')
        s.write('táxon\tvoucher/isolado\tlocalidade\tITS\ttrnL-F' + ('\trpl16' if R16 else '') + '\tregra de inclusão\n')
        f.write('# --- Amostras ---\n')
        def cols(a):
            fs = [f'Consensos/{m}/{a}_{m}.fasta' for m in MARC]
            return '\t'.join(x if os.path.exists(x) else '-' for x in fs)
        for a in o['amostras']: f.write(f'AMOSTRA_{a}\t{cols(a)}\n')
        if o.get('amostra-outgroup'): f.write(f'OUT_AMOSTRA_{o["amostra-outgroup"]}\t{cols(o["amostra-outgroup"])}\n')
        atual = None; vistos = set()
        for (ki, kt), (p, regra) in escolhidos.items():
            r = p[0] or p[1]; pre = 'OUT_' if r['grupo'] == 'out' else ''
            rot = rotulo(r, r['acc'], pre)
            if rot in vistos: continue
            vistos.add(rot); regras[regra] += 1
            if regra != atual: f.write(f'# --- {regra} ---\n'); atual = regra
            c16 = ('\t' + r16_de.get((ki, kt), '-')) if R16 else ''
            f.write(f'{rot}\t{ki}\t{kt}{c16}\n')
            s.write(f'{r["org"]}\t{r["vou"] or r["iso"]}\t{r["pais"]}\t{ki}\t{kt}{c16}\t{regra}\n')
        for r, regra in so_r16:
            rot = rotulo(r, r['acc'], 'OUT_' if r['grupo'] == 'out' else '')
            if rot in vistos: continue
            vistos.add(rot); regras[regra] += 1
            if regra != atual: f.write(f'# --- {regra} ---\n'); atual = regra
            f.write(f'{rot}\t-\t-\t{r["acc"]}\n')
            s.write(f'{r["org"]}\t{r["vou"] or r["iso"]}\t{r["pais"]}\t-\t-\t{r["acc"]}\t{regra}\n')
    if not any('grupo externo' in r for r in regras) and not o.get('amostra-outgroup'):
        print(f'\nATENÇÃO: nenhuma sequência de {o["outgroup"]} entrou como grupo externo (OUT_). A árvore ficaria sem raiz:'
              ' confira o nome do gênero do --outgroup.')
    n = sum(regras.values())
    print(f'\n{n} referências escolhidas -> {lista}  (anterior guardada como .bak)')
    print(f'Tabela S1 -> {t1}')
    print('\nCONFERÊNCIA:')
    for k, v in regras.items(): print(f'  {v:3d}  {k}')
    esp_c = collections.Counter((p[0] or p[1])['org'] for p, _ in escolhidos.values())
    print('  Espécies: ' + '; '.join(f'{k} ({v})' for k, v in esp_c.most_common()))
    if excl: print('  Excluídos dos bancos: ' + '; '.join(f'{k}: {v}' for k, v in excl.items()))
    print('Depois rode as árvores: arvores_por_primer.py e arvore_concatenada.py (manual, passos 10 e 11).')

if __name__ == '__main__':
    main()
