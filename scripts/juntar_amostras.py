# Junta os consensos finais de TODAS as amostras, um arquivo por marcador, para o BLAST (passo 6 do manual).
# Uso:  python3 Scripts/juntar_amostras.py
# Lê Consensos/<marcador>/<amostra>_<marcador>.fasta (ITS, trnLF, rpl16) e grava
#     Consensos/<marcador>/<marcador>_todas_amostras.fasta
# Marcador sem nenhuma amostra ainda (ex.: rpl16 antes de chegar) só é avisado.
import os, glob, sys

if not os.path.isdir('Consensos'): sys.exit('Rode na pasta do projeto: seq && python3 Scripts/juntar_amostras.py')
for m in ('ITS', 'trnLF', 'rpl16'):
    os.makedirs(f'Consensos/{m}', exist_ok=True)
    saida = f'Consensos/{m}/{m}_todas_amostras.fasta'
    arqs = sorted(f for f in glob.glob(f'Consensos/{m}/*_{m}.fasta') if not f.endswith('_todas_amostras.fasta'))
    if not arqs:
        print(f'{m:6s}: ainda nenhuma amostra')
        continue
    n = 0
    with open(saida, 'w', encoding='utf-8') as out:
        for f in arqs:
            txt = open(f, encoding='utf-8', errors='replace').read().strip()
            if txt:
                out.write(txt + '\n'); n += txt.count('>')
    print(f'{m:6s}: {n} amostras -> {saida}')
