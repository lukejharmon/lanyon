pandoc -s --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl chicago-fullnote-bibliography.csl chapter1_introduction_pre.md -o t1.html
cat metadata.yaml t1.html | sed s/%7B%7B%20/{{/ | sed s/%20%7D%7D/}}/ > ../chapter1_introduction.html
rm t1.html
