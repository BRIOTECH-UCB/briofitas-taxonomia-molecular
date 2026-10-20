# BLAST local de consensos contra um banco de musgos (Bryophyta) do GenBank, por marcador.
#
# 1. Na primeira execucao de cada marcador baixa do NCBI (so download, nada e enviado)
#    as sequencias do marcador e cria o banco com makeblastdb em Referencias/BLAST_<marcador>/.
# 2. Roda blastn localmente com Consensos/<marcador>/<marcador>_todas_amostras.fasta.
# 3. Gera a tabela completa e um resumo dos melhores hits por amostra
#    (ordenado por: cobertura >= 90%, depois % identidade).
#
# Uso (na pasta Sequenciamentos):
#   python3 Scripts/blast_local.py --marcador ITS
#   python3 Scripts/blast_local.py --marcador trnLF
#   python3 Scripts/blast_local.py --marcador ITS --atualizar-banco   (baixa o banco de novo)
#   python3 Scripts/blast_local.py --marcador ITS --query outro.fasta (outra consulta)
import argparse, os, subprocess, sys, time, urllib.parse, urllib.request, re

MARCADORES = {
    'trnLF': 'Bryophyta[Organism] AND (trnL[All Fields] OR trnF[All Fields]) AND 150:5000[SLEN]',
    'ITS': 'Bryophyta[Organism] AND ("internal transcribed spacer"[All Fields] OR ITS1[All Fields] OR ITS2[All Fields]) AND 200:5000[SLEN]',
    'rpl16': 'Bryophyta[Organism] AND rpl16[All Fields] AND 200:5000[SLEN]',
}
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
LOTE = 500

def configurar(m):
    global TERMO, PASTA_DB, FASTA_DB, NOME_DB
    TERMO = MARCADORES[m]
    PASTA_DB = os.path.join('Referencias', 'BLAST_' + m)
    FASTA_DB = os.path.join(PASTA_DB, 'bryophyta_' + m + '.fasta')
    NOME_DB = os.path.join(PASTA_DB, 'bryophyta_' + m)

def pedir(url, dados, tentativas=5):
    corpo = urllib.parse.urlencode(dados).encode()
    for t in range(tentativas):
        try:
            with urllib.request.urlopen(url, corpo, timeout=120) as r:
                return r.read().decode('utf-8', 'replace')
        except Exception as e:
            print(f'   falha ({e}), tentando de novo...')
            time.sleep(3 * (t + 1))
    sys.exit('Nao consegui acessar o NCBI. Verifique a internet e rode de novo.')

def baixar_banco():
    os.makedirs(PASTA_DB, exist_ok=True)
    print('Buscando no GenBank:', TERMO)
    xml = pedir(EUTILS + 'esearch.fcgi', {'db': 'nuccore', 'term': TERMO, 'usehistory': 'y', 'retmax': 0, 'tool': 'briotech_blast'})
    total = int(re.search(r'<Count>(\d+)', xml).group(1))
    webenv = re.search(r'<WebEnv>(\S+?)</WebEnv>', xml).group(1)
    qkey = re.search(r'<QueryKey>(\d+)', xml).group(1)
    print(f'{total} sequencias. Baixando em lotes de {LOTE}...')
    tmp = FASTA_DB + '.parcial'
    with open(tmp, 'w', encoding='utf-8') as out:
        for ini in range(0, total, LOTE):
            fa = pedir(EUTILS + 'efetch.fcgi', {'db': 'nuccore', 'query_key': qkey, 'WebEnv': webenv,
                                                'rettype': 'fasta', 'retmode': 'text', 'retstart': ini,
                                                'retmax': LOTE, 'tool': 'briotech_blast'})
            if not fa.lstrip().startswith('>'):
                sys.exit(f'Resposta inesperada do NCBI no lote {ini}. Rode de novo.')
            out.write(fa.strip() + '\n')
            print(f'   {min(ini + LOTE, total)}/{total}')
            time.sleep(0.4)  # limite do NCBI: 3 requisicoes/s sem chave
    os.replace(tmp, FASTA_DB)
    with open(os.path.join(PASTA_DB, 'LEIAME.txt'), 'w', encoding='utf-8') as f:
        f.write(f'Banco baixado em {time.strftime("%Y-%m-%d")}\nTermo: {TERMO}\nSequencias: {total}\n')
    print('Criando banco BLAST...')
    subprocess.run(['makeblastdb', '-in', FASTA_DB, '-dbtype', 'nucl', '-out', NOME_DB,
                    '-title', 'Bryophyta ' + os.path.basename(PASTA_DB) + ' GenBank'], check=True)

def especie(titulo):
    p = titulo.split()
    return ' '.join(p[1:3]) if len(p) >= 3 else titulo

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--marcador', required=True, choices=sorted(MARCADORES))
    ap.add_argument('--query')
    ap.add_argument('--atualizar-banco', action='store_true')
    ap.add_argument('--saida')
    a = ap.parse_args()
    configurar(a.marcador)
    a.query = a.query or os.path.join('Consensos', a.marcador, a.marcador + '_todas_amostras.fasta')
    a.saida = a.saida or os.path.join('Consensos', a.marcador, 'BLAST_' + a.marcador)

    if a.atualizar_banco or not os.path.exists(NOME_DB + '.nsq') and not os.path.exists(NOME_DB + '.00.nsq') and not os.path.exists(NOME_DB + '.nal'):
        baixar_banco()

    if not os.path.exists(a.query) or os.path.getsize(a.query) == 0:
        sys.exit(f'Nenhuma amostra de {a.marcador} para o BLAST ({a.query} não existe ou está vazio).\n'
                 '  Monte os consensos (passo 5) e junte as amostras com: python3 Scripts/juntar_amostras.py (passo 6).\n'
                 '  O banco de BLAST já ficou pronto para quando houver amostras.')
    tabela = a.saida + '_completo.tsv'
    print('Rodando blastn local...')
    subprocess.run(['blastn', '-task', 'blastn', '-query', a.query, '-db', NOME_DB,
                    '-max_target_seqs', '5000', '-evalue', '1e-20',
                    '-outfmt', '6 qseqid sacc pident length qcovs mismatch gapopen evalue bitscore stitle',
                    '-out', tabela], check=True)

    hits = {}
    for l in open(tabela, encoding='utf-8', errors='replace'):
        c = l.rstrip('\n').split('\t')
        q, acc, pid, alen, qcov, bits, tit = c[0], c[1], float(c[2]), int(c[3]), int(c[4]), float(c[8]), c[9]
        lst = hits.setdefault(q, {})
        sp = especie(tit)
        # melhor hit por especie: cobertura >= 90% primeiro, depois identidade, depois bitscore
        chave = (qcov >= 90, pid, bits)
        if sp not in lst or chave > (lst[sp][3] >= 90, lst[sp][1], lst[sp][4]):
            lst[sp] = (acc, pid, alen, qcov, bits, tit)

    resumo = a.saida + '_resumo.tsv'
    with open(resumo, 'w', encoding='utf-8') as out:
        out.write('amostra\trank\tespecie\tacesso\tidentidade_%\talinhamento_pb\tcobertura_%\tbitscore\tdescricao\n')
        for q, d in hits.items():
            top = sorted(d.items(), key=lambda x: (x[1][3] < 90, -x[1][1], -x[1][4]))[:10]
            print(f'\n== {q}')
            for i, (sp, (acc, pid, alen, qcov, bits, tit)) in enumerate(top, 1):
                print(f'  {i:>2}. {sp:<35} {acc:<12} {pid:6.2f}%  {alen:>4} pb  cob {qcov:>3}%')
                out.write(f'{q}\t{i}\t{sp}\t{acc}\t{pid}\t{alen}\t{qcov}\t{bits}\t{tit}\n')
    print(f'\nResumo: {resumo}\nCompleto: {tabela}')

if __name__ == '__main__':
    main()
