# Lê a tabela de pares de primers (Scripts/primers.tsv). Usado pelos outros scripts; não precisa rodar.
import os, sys, re

ARQUIVO = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'primers.tsv')

def ler_primers(arquivo=ARQUIVO):
    """Devolve a lista de pares, na ordem do arquivo, como dicionários."""
    if not os.path.exists(arquivo): sys.exit(f'Falta o arquivo de primers: {arquivo}')
    pares, nomes = [], set()
    for n, l in enumerate(open(arquivo, encoding='utf-8'), 1):
        if l.startswith('#') or not l.strip(): continue
        c = [x.strip() for x in l.rstrip('\n').split('\t')]
        if len(c) < 5: sys.exit(f'primers.tsv, linha {n}: precisa de pelo menos 5 colunas separadas por TAB (par, primer_F, primer_R, marcador, papel)')
        c += ['-'] * (9 - len(c))
        par, pf, pr, marcador, papel, opcoes, sf, sr, outros = c[:9]
        if papel not in ('completo', 'metade1', 'metade2'):
            sys.exit(f'primers.tsv, linha {n}: papel "{papel}" inválido (use completo, metade1 ou metade2)')
        for s in (sf, sr):
            if s != '-' and not re.fullmatch(r'[ACGTRYSWKMBDHVN]+', s.upper()):
                sys.exit(f'primers.tsv, linha {n}: sequência de primer inválida "{s}" (só letras de DNA, 5\'->3\')')
        if par.upper() in nomes: sys.exit(f'primers.tsv, linha {n}: par {par} repetido')
        nomes.add(par.upper())
        pares.append(dict(par=par, pf=pf, pr=pr, marcador=marcador, papel=papel,
                          opcoes=[] if opcoes == '-' else opcoes.split(),
                          seq_f=None if sf == '-' else sf.upper(), seq_r=None if sr == '-' else sr.upper(),
                          nomes=[par] + ([] if outros == '-' else [x.strip() for x in outros.split(',') if x.strip()])))
    return pares

def achar_par(nome, pares):
    """Acha o par pelo nome da planilha (ou um dos outros_nomes), sem diferenciar maiúsculas."""
    n = re.sub(r'\s+', '', nome).upper()
    return next((p for p in pares if n in (x.upper() for x in p['nomes'])), None)

def nomes_primers(pares):
    return {x.upper() for p in pares for x in (p['pf'], p['pr'])}
