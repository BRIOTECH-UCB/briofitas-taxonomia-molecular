# Funções comuns às árvores (usado pelo arvores_por_primer.py e pelo arvore_concatenada.py; não precisa rodar).
# Alinha (MAFFT L-INS-i), apara as pontas (trimAl), roda o IQ-TREE 3 e desenha a figura.
import os, re, sys, shutil, subprocess, time

AQUI = os.path.dirname(os.path.abspath(__file__))

def conferir_programas():
    falta = [p for p in ('mafft', 'trimal', 'iqtree3') if not shutil.which(p)]
    if falta:
        sys.exit('Programa não encontrado: ' + ', '.join(falta) + '. Feche e abra o terminal (ou rode: source ~/.bashrc); '
                 'se continuar, refaça a instalação (manual, parte 1).')

def rodar(cmd, saida=None):
    """Roda um programa; se der erro, mostra as últimas linhas e para."""
    with open(saida, 'w') if saida else open(os.devnull, 'w') as out:
        r = subprocess.run(cmd, stdout=out if saida else subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        msg = [l for l in (r.stderr or '').splitlines() if 'libarcher' not in l][-8:]
        sys.exit(f'ERRO ao rodar {cmd[0]}:\n   ' + '\n   '.join(msg))

def alinhar_aparar(entrada, prefixo):
    """<prefixo>_aln.fasta (MAFFT) e <prefixo>_aparado.fasta (trimAl). Tira o _R_ que o MAFFT põe nas invertidas."""
    aln = prefixo + '_aln.fasta'
    rodar(['mafft', '--quiet', '--adjustdirection', '--localpair', '--maxiterate', '1000', entrada], aln)
    txt = open(aln).read()
    invertidas = len(re.findall(r'^>_R_', txt, re.M))
    open(aln, 'w').write(re.sub(r'^>_R_', '>', txt, flags=re.M))
    ap = prefixo + '_aparado.fasta'
    rodar(['trimal', '-in', aln, '-out', ap, '-terminalonly', '-gt', '0.5', '-htmlout', prefixo + '_trimal.html'])
    return ap, invertidas

def outgroups(fasta):
    return [l[1:].split()[0] for l in open(fasta) if l.startswith('>OUT_')]

def iqtree(alinhamento, prefixo, particoes=None):
    out = outgroups(alinhamento)
    cmd = ['iqtree3', '-s', alinhamento, '-m', 'MFP+MERGE' if particoes else 'MFP', '-B', '5000', '--alrt', '1000',
           '--bnni', '-T', 'AUTO', '--prefix', prefixo, '--redo', '-quiet']
    if particoes: cmd[3:3] = ['-p', particoes]
    if out: cmd += ['-o', ','.join(out)]
    else: print('   AVISO: nenhum grupo externo (OUT_) neste alinhamento - árvore sem raiz definida')
    rodar(cmd)
    rel = open(prefixo + '.iqtree', encoding='utf-8', errors='replace').read()
    for k in ('Input data', 'Number of parsimony informative sites', 'Best-fit model according to BIC'):
        m = re.search(re.escape(k) + r'.*', rel)
        if m: print('   ' + m.group(0))
    return prefixo + '.treefile'

def figura(treefile, svg, titulo):
    fig = treefile.replace('.treefile', '_figura.tree')
    for cmd in ([sys.executable, os.path.join(AQUI, 'preparar_figura.py'), treefile, 'OUT_', fig],
                [sys.executable, os.path.join(AQUI, 'desenhar_arvore.py'), fig, svg, titulo]):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0: print('   AVISO: a figura não foi gerada: ' + (r.stderr.strip().splitlines() or ['?'])[-1]); return None
    return svg

class _Tee:
    def __init__(self, *fs): self.fs = fs
    def write(self, x):
        for f in self.fs: f.write(x); f.flush()
    def flush(self):
        for f in self.fs: f.flush()

def registrar(arquivo):
    """Grava tudo o que aparece na tela também em <arquivo> (acrescenta; cada rodada com data e hora)."""
    os.makedirs(os.path.dirname(arquivo), exist_ok=True)
    f = open(arquivo, 'a', encoding='utf-8')
    f.write(f"\n===== {time.strftime('%d/%m/%Y %H:%M')} | python3 {' '.join(sys.argv)}\n")
    sys.stdout = _Tee(sys.__stdout__, f)
    print(f'(registro desta rodada em {arquivo})')
