pandoc -s --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true chapter1_introduction_pre.md -o t1.html
cat metadata_ch1.yaml t1.html | sed s/%7B%7B%20/{{/ | sed s/%20%7D%7D/}}/ > ../chapter1_introduction.html
rm t1.html

pandoc -s --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true chapter2_stats_pre.md -o t1.html
cat metadata_ch2.yaml t1.html | sed s/%7B%7B%20/{{/ | sed s/%20%7D%7D/}}/ > ../chapter2_stats.html
rm t1.html

pandoc -s --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true chapter3_bmintro_pre.md -o t1.html
cat metadata_ch3.yaml t1.html | sed s/%7B%7B%20/{{/ | sed s/%20%7D%7D/}}/ > ../chapter3_bmintro.html
rm t1.html

pandoc -s --bibliography pcm_paperpile.bib --filter pandoc-citeproc --csl evolution.csl --metadata link-citations=true chapter4_fitbm_pre.md -o t1.html
cat metadata_ch4.yaml t1.html | sed s/%7B%7B%20/{{/ | sed s/%20%7D%7D/}}/ > ../chapter4_fitbm.html
rm t1.html
