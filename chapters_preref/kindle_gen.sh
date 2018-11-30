#!/bin/bash
#
# requires : pandoc, pandoc-citeproc, pandoc-fignos 

cat chapter1_introduction_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t1.md
cat chapter2_stats_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t2.md
cat chapter3_bmintro_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t3.md
cat chapter4_fitbm_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t4.md
cat chapter5_mvbm_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t5.md
cat chapter6_beyondbm_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t6.md
cat chapter7_introdiscrete_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t7.md
cat chapter8_fitdiscrete_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t8.md
cat chapter9_beyondmk_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t9.md
cat chapter10_birthdeath_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t10.md
cat chapter11_fitbd_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t11.md
cat chapter12_beyondbd_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t12.md
cat chapter13_chardiv_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t13.md
cat chapter14_summary_pre.md | sed s@\{\{[[:space:]]site.baseurl[[:space:]]\}\}@..@ > t14.md

pandoc titlepage.md copyright.md acknowledgements.md toc.md t1.md t2.md t3.md t4.md t5.md t6.md t7.md t8.md t9.md t10.md t11.md t12.md t13.md t14.md --from=markdown-markdown_in_html_blocks-native_divs --pdf-engine=xelatex --filter pandoc-fignos --webtex https://latex.codecogs.com/svg.latex? --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -V geometry:margin=0.7in -V fontsize=12pt -V geometry:c5paper -o ../pdf/phylogeneticComparativeMethods.epub

rm t1.md t2.md t3.md t4.md t5.md t6.md t7.md t8.md t9.md t10.md t11.md t12.md t13.md t14.md
