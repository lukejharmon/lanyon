cat chapter1_introduction_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter1_introduction.pdf


cat chapter2_stats_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter2_stats.pdf

cat chapter3_bmintro_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > temp.md
pandoc temp.md --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter3_bmintro.pdf


"/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ../pdf/phylogeneticComparativeMethods.pdf ../pdf/chapter1_introduction.pdf ../pdf/chapter2_stats.pdf ../pdf/chapter3_bmintro.pdf
