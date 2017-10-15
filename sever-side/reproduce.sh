#Download all the files needed
cd annotation
wget http://merlot.lbl.gov/keliu/HeLab-AndrewProject.support.files/annotation/*.gtf
wait
cd ../

wget -r http://merlot.lbl.gov/keliu/HeLab-AndrewProject.support.files/fastq/
wait

#Reproduce the results
perl code/map_read.pl annotation/gao.chip-seq.exp.txt 
wait

perl code/call.peak.and.build.signal.pl annotation/gao.chip-seq.exp.txt 
wait

Rscript code/analyze.alignment.R
wait 
