# Gera a sequência consenso a partir do alinhamento dos reads exportado do MEGA (FASTA).
# As decisões (aparar, corrigir discordâncias olhando o cromatograma) são feitas no MEGA;
# este script só junta as linhas e confere se sobrou alguma discordância.
#
# Uso: python3 Scripts/consenso_mega.py <reads_alinhados.fasta> <NOME> <saida.fasta>
#   Funciona com 2 reads (F + R) ou com 4 (as duas metades do ITS), já no mesmo sentido.
#   Coluna em que só um read tem base -> essa base.
#   Reads concordam -> a base.
#   Reads discordam (ou um tem gap no meio da sua leitura e outro tem base) -> código IUPAC + aviso REVISAR.
import sys

IUPAC = {frozenset('AG'): 'R', frozenset('CT'): 'Y', frozenset('GC'): 'S', frozenset('AT'): 'W',
         frozenset('GT'): 'K', frozenset('AC'): 'M', frozenset('CGT'): 'B', frozenset('AGT'): 'D',
         frozenset('ACT'): 'H', frozenset('ACG'): 'V'}
EXP = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
       'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG'}

def ler(p):
    d = {}; o = []; k = None
    for l in open(p, encoding='utf-8', errors='replace'):
        l = l.strip()
        if not l: continue
        if l.startswith('>'):
            k = l[1:].strip(); d[k] = []; o.append(k)
        elif k: d[k].append(l.upper().replace('?', '-').replace('.', '-'))
    return [(k, ''.join(d[k])) for k in o]

def main(entrada, nome, saida):
    reads = ler(entrada)
    if len(reads) < 2: sys.exit('O arquivo precisa ter pelo menos 2 reads alinhados.')
    L = {len(s) for _, s in reads}
    if len(L) != 1: sys.exit('Os reads têm comprimentos diferentes: exporte o ALINHAMENTO do MEGA (Data > Export Alignment > FASTA).')
    L = L.pop()
    faixa = []
    for _, s in reads:
        b = [i for i, c in enumerate(s) if c != '-']
        faixa.append((b[0], b[-1]) if b else (L, -1))
    cons = []; avisos = []
    for i in range(L):
        bases = []; gap_interno = []
        for (nm, s), (a, z) in zip(reads, faixa):
            if not (a <= i <= z): continue
            c = s[i]
            if c == '-': gap_interno.append(nm)
            elif c != 'N': bases.append((nm, c))
        if not bases and not gap_interno: continue
        if not bases:
            continue  # todos os reads que cobrem têm gap: coluna vazia
        conj = set()
        for _, c in bases: conj |= set(EXP.get(c, 'ACGT'))
        distintas = {c for _, c in bases}
        pos = len(cons) + 1
        if gap_interno:
            avisos.append(f'pos {pos} (coluna {i + 1}): base em {", ".join(n for n, _ in bases)} mas gap em {", ".join(gap_interno)} -> REVISAR (indel?)')
        if len(distintas) == 1:
            cons.append(distintas.pop())
        else:
            k = IUPAC.get(frozenset(conj), 'N') if len(conj) > 1 else conj.pop()
            cons.append(k)
            avisos.append(f'pos {pos} (coluna {i + 1}): ' + ', '.join(f'{n}={c}' for n, c in bases) + f' -> {k} REVISAR')
    seq = ''.join(cons).strip('N')
    with open(saida, 'w') as f:
        f.write(f'>{nome}\n' + '\n'.join(seq[j:j + 70] for j in range(0, len(seq), 70)) + '\n')
    print(f'{nome}: {len(reads)} reads, consenso {len(seq)} pb, ambíguos {sum(c not in "ACGT" for c in seq)}, avisos {len(avisos)} -> {saida}')
    for a in avisos: print('   ', a)
    if avisos:
        print('Há posições para REVISAR: volte ao MEGA, confira no cromatograma, corrija os reads e exporte de novo.')

if __name__ == '__main__':
    if len(sys.argv) != 4: sys.exit('Uso: python3 Scripts/consenso_mega.py <reads_alinhados.fasta> <NOME> <saida.fasta>')
    main(*sys.argv[1:4])
