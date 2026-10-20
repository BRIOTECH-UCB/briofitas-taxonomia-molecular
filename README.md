# Taxonomia molecular de briófitas — BRIOTECH

Método desenvolvido pelo **BRIOTECH**, projeto da **Universidade Católica de Brasília**, e utilizado na rotina de taxonomia molecular de briófitas do grupo. Os scripts, em Python 3, cobrem a montagem das sequências consenso (Sanger), a identificação por BLASTn local, a seleção das sequências de referência do GenBank, as análises de máxima verossimilhança e a preparação das figuras das árvores.

*Method developed by BRIOTECH, a project of the Catholic University of Brasília (Universidade Católica de Brasília, Brazil), and used in the group's routine molecular taxonomy of bryophytes. The Python 3 scripts cover Sanger consensus assembly, local BLASTn identification, selection of GenBank reference sequences, maximum likelihood analyses and tree figures. Code comments are in Portuguese.*

## Requisitos

- Python 3 (só a biblioteca padrão)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/), [MAFFT](https://mafft.cbrc.jp/alignment/software/) v7, [trimAl](https://github.com/inab/trimal) v1.5 e [IQ-TREE 3](https://github.com/iqtree/iqtree3) no PATH
- Testado no Ubuntu 24.04 (Windows Subsystem for Linux 2)

## Arquivos

| Arquivo | Etapa | O que faz |
| --- | --- | --- |
| `qc_ab1.py` | Qualidade das leituras | lê os valores de qualidade Phred dos `.ab1` e calcula o trecho de boa qualidade |
| `primers.tsv` | Pares de primers | pares de primers, região amplificada e opções de montagem |
| `ler_primers.py` | Pares de primers | lê o `primers.tsv` (usado pelos demais scripts) |
| `montar_par.py` | Consenso | apara as leituras (algoritmo de Mott modificado), inverte a reversa, alinha e decide cada base pela qualidade; discordâncias viram IUPAC; no trnL-F remove os sítios dos primers |
| `montar_todos.py` | Consenso | aplica o `montar_par.py` a todos os pares de primers de cada amostra |
| `juntar_metades.py` | Consenso | une as duas metades do ITS (ITS1F/5.8R e 5.8F/ITS2R) pela sobreposição no 5.8S |
| `consenso_mega.py` | Consenso | gera o consenso das leituras alinhadas e editadas no MEGA 7 |
| `juntar_amostras.py` | BLAST | reúne os consensos de todas as amostras, um arquivo por região |
| `blast_local.py` | BLAST | baixa do GenBank (NCBI E-utilities) as sequências de Bryophyta da região, cria o banco local e roda o BLASTn |
| `selecionar_referencias.py` | Referências | aplica os critérios de seleção das referências e do grupo externo; gera a lista de espécimes e a tabela suplementar |
| `concatenado.py` | Entradas | monta os FASTA de entrada a partir da lista de espécimes |
| `rodar_arvore.py` | Árvores | MAFFT (L-INS-i), trimAl (`-terminalonly -gt 0.5`) e IQ-TREE 3 (ModelFinder/BIC, 5.000 UFBoot2, 1.000 SH-aLRT) |
| `arvores_por_primer.py` | Árvores | uma análise por par de primers (referências + consensos das amostras) |
| `preparar_figura.py` | Figuras | enraíza no grupo externo, ordena os clados e mantém os suportes SH-aLRT/UFBoot |
| `desenhar_arvore.py` | Figuras | desenha a árvore em SVG |

## Ordem de execução

Na pasta do projeto (que contém a pasta `scripts`, os `.ab1` organizados em `Organizado_por_amostra/AMOSTRA/AMOSTRA_PRIMER.ab1`):

```
python3 scripts/qc_ab1.py . > qc.csv
python3 scripts/montar_todos.py AMOSTRA
python3 scripts/juntar_amostras.py
python3 scripts/blast_local.py --marcador ITS
python3 scripts/selecionar_referencias.py --genero GENERO --amostras AMOSTRA --outgroup Funaria
python3 scripts/arvores_por_primer.py --genero GENERO
```

Cada script descreve o próprio uso nas primeiras linhas.

## Como citar

BRIOTECH – Universidade Católica de Brasília. 2026. *Taxonomia molecular de briófitas — BRIOTECH* [Software]. Disponível em: https://github.com/BRIOTECH-UCB/briofitas-taxonomia-molecular.

## Licença

[MIT](LICENSE).
