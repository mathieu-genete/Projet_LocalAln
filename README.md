## Description
`localAln.py` est un programme d'alignement local de séquences. Il permet d'analyser des séquences au format FASTA et de calculer des alignements locaux en utilisant différents paramètres de configuration.

## Utilisation
Pour afficher l'aide :
```bash
python localAln.py -h
```

### Arguments
- `-i`, `--infasta` : Fichier FASTA contenant les séquences à analyser (obligatoire).
- `-c`, `--config` : Fichier de configuration des scores (défaut : `datas/score.txt`).
- `-o`, `--outstats` : Exporter les statistiques d'alignement.
- `-f`, `--outfasta` : Exporter l'alignement au format FASTA.
- `-l`, `--linesize` : Nombre de bases par ligne pour l'alignement (défaut : 50).
- `-n`, `--rndseq` : Nombre de séquences aléatoires à générer pour le calcul de la significativité (défaut : 100).
- `-g`, `--gapaffine` : Utiliser le système de score affine pour les indels.
- `-s`, `--significance` : Calcul de la significativité du score d'alignement.
- `-M`, `--showmatrix` : Afficher la matrice d'alignement.

### Exemple de commande
```bash
python localAln.py -i sequences.fasta -c config.txt -o stats.txt -f alignment.fasta -l 60 -n 200 -g -s -M
```


## Auteur
- **Mathieu Genete** - [mathieu.genete@univ-lille.fr](mailto:mathieu.genete@univ-lille.fr)

## Licence
Ce projet est sous licence copyleft.

## Version
1.0.0 (Janvier 2024)
