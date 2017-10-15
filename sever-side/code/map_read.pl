use strict;
my $bowtie2_index="/srv/persistent/keliu/genomes/mm10/mm10.index";
my $bowtie2_option="--very-sensitive -p 1";

my @map_cmd_list;
my @picard_cmd_list;
my @filter_cmd_list;
my @feature_counts_cmd_list;
my @macs2_cmd_list;


my %chip_exp_name;
my $out_put = `ls fastq| grep fq |grep -v H3K27me3 |sed s/.fq.gz//g`;
grep{
    chomp $_;
    my @tmp = split /\./,$_;
    my $exp_name =$tmp[0].".".$tmp[1];
    if($_!~/input/) {$exp_name = $exp_name . "." .$tmp[2];}
    $chip_exp_name{$exp_name}=1;

} split /\n+/,$out_put;



my $out_put = `ls fastq| grep input |sed s/.fq.gz//g`;
grep{
    chomp $_;
    my @tmp = split /\./,$_;
    my $exp_name =$tmp[0].".".$tmp[1];
    if($_!~/input/) {$exp_name = $exp_name . "." .$tmp[2];}
    $chip_exp_name{$exp_name}=1;

} split /\n+/,$out_put;





foreach my $exp_name (keys %chip_exp_name){
    my ($map_cmd,$filter_cmd,$feature_counts_cmd,$macs2_cmd,$picard_cmd);
    my @tmp = split /\./,$exp_name;
    my $input_name = $tmp[0].".input.unique.bam";    
    if(-e "fastq/$exp_name.1.fq.gz"){ #paired-end sequencing
        $map_cmd    = "bowtie2 $bowtie2_option --no-discordant --no-mixed -x $bowtie2_index -1 fastq/$exp_name.1.fq.gz -2 fastq/$exp_name.2.fq.gz  -X 1000  | samtools view -bS -f 2  - |samtools sort -o  alignment/$exp_name.all.bam  -T sorted.all.$exp_name - \n";
    }else{
        $map_cmd    = "bowtie2 $bowtie2_option -x $bowtie2_index -U fastq/$exp_name.fq.gz                                                                   | samtools view -bS -F 4  - |samtools sort -o  alignment/$exp_name.all.bam  -T sorted.all.$exp_name - \n";
    }  
    $filter_cmd = "samtools view -b -q 3  alignment/$exp_name.all.bam | samtools sort -o alignment/$exp_name.unique.bam  -T sorted.$exp_name - \n";
    $picard_cmd = "java -jar /users/keliu/bin/picard/picard.jar MarkDuplicates INPUT=alignment/$exp_name.unique.bam OUTPUT=alignment/$exp_name.unique.rmdup.bam METRICS_FILE=tmp/$exp_name.unique.metrics REMOVE_DUPLICATES=true\n".
                  "java -jar /users/keliu/bin/picard/picard.jar MarkDuplicates INPUT=alignment/$exp_name.all.bam    OUTPUT=alignment/$exp_name.all.rmdup.bam    METRICS_FILE=tmp/$exp_name.all.metrics    REMOVE_DUPLICATES=true\n";

   push @map_cmd_list,$map_cmd       if not -e "alignment/$exp_name.all.bam";
   push @filter_cmd_list,$filter_cmd if not -e "alignment/$exp_name.unique.bam";
   push @picard_cmd_list,$picard_cmd if not -e "alignment/$exp_name.unique.rmdup.bam";

}

open(OUTFILE,"> tmp/map_cmd_list");
print OUTFILE join("",@map_cmd_list);
close OUTFILE;

open(OUTFILE,"> tmp/filter_cmd_list");
print OUTFILE join("",@filter_cmd_list);
close OUTFILE;

open(OUTFILE,"> tmp/picard_cmd_list");
print OUTFILE join("",@picard_cmd_list);
close OUTFILE;




`parallel --no-notice -j 10  < tmp/map_cmd_list`;
`parallel --no-notice -j 10  < tmp/filter_cmd_list`;
`parallel --no-notice -j 10 < tmp/picard_cmd_list`;










