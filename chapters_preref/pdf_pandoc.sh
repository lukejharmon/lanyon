cat chapter1_introduction_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter1_introduction.pdf


cat chapter2_stats_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter2_stats.pdf

cat chapter3_bmintro_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter3_bmintro.pdf

cat chapter4_fitbm_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter4_fitbm.pdf

cat chapter5_mvbm_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter5_mvbm.pdf

cat chapter6_beyondbm_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter6_beyondbm.pdf

cat chapter7_introdiscrete_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter7_introdiscrete.pdf


"/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ../pdf/phylogeneticComparativeMethods.pdf ../pdf/chapter1_introduction.pdf ../pdf/chapter2_stats.pdf ../pdf/chapter3_bmintro.pdf ../pdf/chapter4_fitbm.pdf ../pdf/chapter5_mvbm.pdf ../pdf/chapter6_beyondbm.pdf ../pdf/chapter7_introdiscrete.pdf
