
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

pandoc t1.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter1_introduction.pdf
pandoc t2.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter2_stats.pdf
pandoc t3.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter3_bmintro.pdf
pandoc t4.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 4." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter4_fitbm.pdf
pandoc t5.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 5." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter5_mvbm.pdf
pandoc t6.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 6." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter6_beyondbm.pdf
pandoc t7.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 7." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter7_introdiscrete.pdf
pandoc t8.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 8." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter8_fitdiscrete.pdf
pandoc t9.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 9." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter9_beyondmk.pdf
pandoc t10.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 10." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter10_birthdeath.pdf
pandoc t11.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 11." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter11_fitbd.pdf
pandoc t12.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 12." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter12_beyondbd.pdf
pandoc t13.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 13." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter13_chardiv.pdf
pandoc t14.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --filter pandoc-fignos -M fignos-caption-name="Figure 14." --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/chapter14_summary.pdf


pandoc titlepage.md copyright.md acknowledgements.md toc.md t1.md t2.md t3.md t4.md t5.md t6.md t7.md t8.md t9.md t10.md t11.md t12.md t13.md t14.md --from=markdown-markdown_in_html_blocks-native_divs --latex-engine=xelatex --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true -o ../pdf/phylogeneticComparativeMethods.pdf

rm t1.md t2.md t3.md t4.md t5.md t6.md t7.md t8.md t9.md t10.md t11.md t12.md t13.md t14.md
