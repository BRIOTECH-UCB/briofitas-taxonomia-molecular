# Monta automaticamente o consenso de TODOS os pares de primers de uma ou mais amostras.
# Uso:  python3 Scripts/montar_todos.py AMOSTRA [AMOSTRA2 ...]
#       python3 Scripts/montar_todos.py --todas          (todas as pastas de Organizado_por_amostra)
#       python3 Scripts/montar_todos.py AMOSTRA --pares rpl16Fb-rpl16R   (só esses pares; os outros arquivos ficam intactos)
#
# Os pares de primers vêm de Scripts/primers.tsv (para incluir um par novo, acrescente uma linha lá).
# Para cada par que tiver os .ab1 em Organizado_por_amostra/<amostra>/ (só o F ou só o R = read único):
#   monta Consensos/por_primer/<amostra>_<par>.fasta (+ _log.txt)  -> um por par (usado na planilha)
# Consensos finais, por marcador (Consensos/<marcador>/<amostra>_<marcador>.fasta):
#   o 1o par "completo" (na ordem do primers.tsv) com F+R que deu certo; se nenhum, as metades
#   (metade1 + metade2) juntas; se o marcador não tem metades, aceita read único do par completo.
#   Os consensos finais NÃO são sobrescritos se já existirem (ex.: montados no MEGA 7); use --sobrescrever.
import sys, os, re, subprocess, glob

AQUI = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
sys.path.insert(0, AQUI)
from ler_primers import ler_primers, achar_par
TODOS = ler_primers()
PARES = TODOS

def montar(a, p):
    par, pf, pr = p['par'], p['pf'], p['pr']
    opts = p['opcoes']
    if p['seq_f'] or p['seq_r']: opts = [p['seq_f'] or '-', p['seq_r'] or '-'] + opts
    d = os.path.join('Organizado_por_amostra', a)
    f, r = os.path.join(d, f'{a}_{pf}.ab1'), os.path.join(d, f'{a}_{pr}.ab1')
    tem_f, tem_r = os.path.exists(f), os.path.exists(r)
    if not (tem_f or tem_r): return None
    os.makedirs('Consensos/por_primer', exist_ok=True)
    saida = f'Consensos/por_primer/{a}_{par}.fasta'
    if tem_f and tem_r: args = [f, r]; tipo = 'F+R'
    else: args = [f if tem_f else r, '-']; tipo = 'read único (' + (pf if tem_f else pr) + ')'
    res = subprocess.run([PY, os.path.join(AQUI, 'montar_par.py')] + args + [f'{a}_{par}', saida] + opts,
                         capture_output=True, text=True, encoding='utf-8', errors='replace')
    out = res.stdout
    open(saida.replace('.fasta', '_log.txt'), 'w', encoding='utf-8').write(out + res.stderr)
    if res.returncode != 0 or not os.path.exists(saida):
        return dict(par=par, tipo=tipo, ok=False, txt='ERRO: ' + (res.stderr.strip().splitlines() or ['?'])[-1], saida=None)
    g = lambda k: (re.search(k + r'=(\d+)', out) or [None, '0'])[1]
    ov, cons = int(g('sobreposicao')), int(g('consenso'))
    revisar, avisos = out.count('REVISAR'), out.count('AVISO')
    ok = cons > 0 and (tipo != 'F+R' or ov >= 30)
    txt = f'{tipo}: consenso {cons} pb' + (f', sobreposição {ov}' if tipo == 'F+R' else '')
    if tipo == 'F+R' and ov < 30: txt += '  -> SEM SOBREPOSIÇÃO, não usar'
    if revisar: txt += f', {revisar} posição(ões) REVISAR'
    if avisos: txt += f', {avisos} AVISO'
    return dict(par=par, tipo=tipo, ok=ok, txt=txt, saida=saida)

def final(origem, destino, sobrescrever, nome):
    if os.path.exists(destino) and not sobrescrever: return f'{destino} já existe (mantido)'
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    seq = ''.join(l.strip() for l in open(origem) if not l.startswith('>'))
    open(destino, 'w').write(f'>{nome}\n' + '\n'.join(seq[i:i + 70] for i in range(0, len(seq), 70)) + '\n')
    return f'{destino} gravado ({len(seq)} pb, de {os.path.basename(origem)})'

def amostra(a, sobrescrever):
    print(f'== {a}')
    feitos = {}
    for p in PARES:
        r = montar(a, p)
        if r: feitos[p['par']] = r; print(f'   {p["par"]:15s} {r["txt"]}')
    conhecidos = {x.upper() for p in TODOS for x in (p['pf'], p['pr'])}
    for f in sorted(glob.glob(os.path.join('Organizado_por_amostra', a, f'{a}_*.ab1'))):
        primer = os.path.basename(f)[len(a) + 1:-4]
        if primer.upper() not in conhecidos:
            print(f'   AVISO: {os.path.basename(f)} - primer "{primer}" não está no Scripts/primers.tsv; não foi montado (acrescente o par lá)')
    if not feitos: print('   nenhum .ab1 encontrado em Organizado_por_amostra/' + a); return
    ok = lambda par: par in feitos and feitos[par]['ok']
    for m in dict.fromkeys(p['marcador'] for p in PARES):
        do_m = [p for p in PARES if p['marcador'] == m]
        if not any(p['par'] in feitos for p in do_m): continue
        destino = f'Consensos/{m}/{a}_{m}.fasta'
        completos = [p['par'] for p in do_m if p['papel'] == 'completo']
        m1 = next((p['par'] for p in do_m if p['papel'] == 'metade1'), None)
        m2 = next((p['par'] for p in do_m if p['papel'] == 'metade2'), None)
        bom = next((x for x in completos if ok(x) and feitos[x]['tipo'] == 'F+R'), None)
        if bom: print('   -> ' + final(feitos[bom]['saida'], destino, sobrescrever, f'{a}_{m}'))
        elif m1 and m2 and ok(m1) and ok(m2):
            if os.path.exists(destino) and not sobrescrever: print(f'   -> {destino} já existe (mantido)'); continue
            os.makedirs(os.path.dirname(destino), exist_ok=True)
            res = subprocess.run([PY, os.path.join(AQUI, 'juntar_metades.py'), feitos[m1]['saida'], feitos[m2]['saida'],
                                  f'{a}_{m}', destino], capture_output=True, text=True, encoding='utf-8', errors='replace')
            print(f'   -> {m} das duas metades: ' + (res.stdout.strip().splitlines() or [res.stderr.strip()])[0])
            print('      confira a junção: ' + next((l.strip() for l in res.stdout.splitlines() if 'juncao' in l or 'junção' in l), '-'))
        elif not (m1 or m2) and any(ok(x) for x in completos):
            x = next(x for x in completos if ok(x))
            print('   -> ' + final(feitos[x]['saida'], destino, sobrescrever, f'{a}_{m}'))
        else:
            print(f'   -> {m} final NÃO montado: nenhum par de {m} serviu. Veja no MEGA 7; se só um read de uma metade for bom,')
            print('      monte essa metade como read único e junte com --lacuna (manual, seções 6.3 e 6.4).')
    if any('REVISAR' in f['txt'] or 'AVISO' in f['txt'] or 'read único' in f['tipo'] for f in feitos.values()):
        print('   Olhe no MEGA 7 os pares com REVISAR/AVISO/read único (detalhes em Consensos/por_primer/<amostra>_<par>_log.txt).')

if __name__ == '__main__':
    args = sys.argv[1:]
    sobrescrever = '--sobrescrever' in args
    args = [x for x in args if x != '--sobrescrever']
    if '--pares' in args:
        i = args.index('--pares')
        if i + 1 >= len(args): sys.exit('Depois de --pares, diga os pares separados por vírgula (ex.: --pares rpl16Fb-rpl16R)')
        nomes = args[i + 1].split(','); del args[i:i + 2]
        PARES = []
        for n in nomes:
            p = achar_par(n, TODOS)
            if not p: sys.exit(f'Par {n} não está no Scripts/primers.tsv')
            PARES.append(p)
    if not args: sys.exit('Uso: python3 Scripts/montar_todos.py AMOSTRA [AMOSTRA2 ...]   ou   --todas')
    if args == ['--todas']: args = sorted(os.path.basename(p) for p in glob.glob('Organizado_por_amostra/*') if os.path.isdir(p))
    for a in args: amostra(a, sobrescrever)
