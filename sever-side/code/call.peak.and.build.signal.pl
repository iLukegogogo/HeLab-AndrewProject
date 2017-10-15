use strict;
use Data::Dumper;
my  $chip_seq_exp_list = $ARGV[0];
my  %global_config_hash;
my  @chip_seq_exp_array;
my  $chip_seq_exp={};

open(INFILE,"<  $chip_seq_exp_list"  ) || die "Can not find  $chip_seq_exp_list\n";
while(my $line=<INFILE>){
    chomp $line;
    $line=~s/^\s+//;
    $line=~s/\s+$//;
    next if $line=~/^#/; # skip annotation line
    next if $line eq ""; # skip empty line
    my @tmp = split /\:/,$line;
    my $key = $tmp[0];
    my $value = $tmp[1];
    $value=~s/^\s+//;
    $value=~s/\s+$//;
    if($key =~/^GLOBAL/){
        $global_config_hash{$key} = $value;
    }else{
        $chip_seq_exp->{$key} = $value;    
        if($key eq "MODE"){
            push @chip_seq_exp_array,$chip_seq_exp;
            $chip_seq_exp = {};
        }
    }

}
close INFILE;

print(Dumper(\%global_config_hash));
my $global_bam_directory    = $global_config_hash{'GLOBAL_BAM_DIRECTORY'};
my $global_tmp_directory    = $global_config_hash{'GLOBAL_TMP_DIRECTORY'};
my $global_peak_directory   = $global_config_hash{'GLOBAL_PEAK_DIRECTORY'};
my $global_signal_directory = $global_config_hash{'GLOBAL_SIGNAL_DIRECTORY'};
my $global_species          = $global_config_hash{'GLOBAL_SPECIES'};
my $global_chr_size_file    = $global_config_hash{'GLOBAL_CHR_SIZE_FILE'};


my @macs2_call_peak_cmd;
my @macs2_bdgcmp_cmd;
my @sort_cmd;
my @bigwig_cmd;
my @mv_peak_cmd;
my @sort_peak_cmd;
my @bigbed_cmd;

foreach my $h (@chip_seq_exp_array) {
    print $h->{'MODE'},"\n";
    my $macs2_cmd = "macs2 callpeak --outdir $global_tmp_directory  -g $global_species ";
    if(exists $h->{'IP'}){
        my $value  = $h->{'IP'};
        my @tmp    = map{ $global_bam_directory."/".$_;}  split /\s+/,$value;
        $macs2_cmd = $macs2_cmd . " -t " . join(" ",@tmp) . " ";
    }
    if(exists $h->{'INPUT'}){
        my $value  = $h->{'INPUT'};
        my @tmp    = map{ $global_bam_directory."/".$_;}   split /\s+/,$value;
        $macs2_cmd = $macs2_cmd . " -c " . join(" ",@tmp) . " ";    
    }
    if(exists $h->{'MODE'}){
        my $value  = $h->{'MODE'};
        $macs2_cmd = $macs2_cmd . " " . $global_config_hash{$value} . " ";
    }
    if(exists $h->{'EXP_NAME'}){
        my $value  = $h->{'EXP_NAME'};
        $macs2_cmd = $macs2_cmd . " -n " . $value . " ";
    }
    $macs2_cmd = $macs2_cmd . "\n";
    
    my $treat_pileup_bdg   = $global_tmp_directory    . "/"  . $h->{'EXP_NAME'} . "_treat_pileup.bdg";
    my $control_lambda_bdg = $global_tmp_directory    . "/"  . $h->{'EXP_NAME'} . "_control_lambda.bdg";
    my $FE_bdg          = $global_tmp_directory    . "/"  . $h->{'EXP_NAME'} . "_FE.bdg";
    my $bigwig_file        = $global_signal_directory . "/"  . $h->{'EXP_NAME'} . ".FE.bigwig";
    my $bdgcmp_cmd         = " macs2  bdgcmp -t $treat_pileup_bdg -c $control_lambda_bdg -o $FE_bdg -m FE -p 0.01 \n";
    my $sort_cmd           = " sort -k1,1 -k2,2n -o $FE_bdg $FE_bdg \n";
    my $bigwig_cmd         = " bedGraphToBigWig $FE_bdg  $global_chr_size_file $bigwig_file \n";
    my $mv_peak_cmd        = " cat $global_tmp_directory/".$h->{'EXP_NAME'} ."_peaks.broadPeak | grep -v '#' |cut -f 1,2,3 > $global_peak_directory/".$h->{'EXP_NAME'} ."_peaks.bed\n";
    my $peak_file          = " $global_peak_directory/".$h->{'EXP_NAME'} ."_peaks.bed";
    my $sort_peak_cmd      = " sort -k1,1 -k2,2n -o $peak_file $peak_file\n";
    my $bigbed_file        = $global_signal_directory . "/"  . $h->{'EXP_NAME'} . ".peaks.bigbed";
    my $bigbed_cmd         = "bedToBigBed $peak_file $global_chr_size_file $bigbed_file\n";
    
    push @macs2_call_peak_cmd, $macs2_cmd;
    push @macs2_bdgcmp_cmd,    $bdgcmp_cmd;
    push @sort_cmd,            $sort_cmd;
    push @bigwig_cmd,          $bigwig_cmd;
    push @mv_peak_cmd,         $mv_peak_cmd;
    push @sort_peak_cmd,       $sort_peak_cmd;
    push @bigbed_cmd,          $bigbed_cmd;


}


open(OUTFILE, "> $global_tmp_directory/macs2_cmd_list");
print OUTFILE join("",@macs2_call_peak_cmd);
close OUTFILE;


open(OUTFILE,"> $global_tmp_directory/bdgcmp_cmd_list");
print OUTFILE join("",@macs2_bdgcmp_cmd);
close OUTFILE;

open(OUTFILE,"> $global_tmp_directory/sort_cmd_list");
print OUTFILE join("",@sort_cmd);
close OUTFILE;

open(OUTFILE,"> $global_tmp_directory/bigwig_cmd_list");
print OUTFILE join("",@bigwig_cmd);
close OUTFILE;


open(OUTFILE,"> $global_tmp_directory/mv_peak_cmd_list");
print OUTFILE join("",@mv_peak_cmd);
close OUTFILE;


open(OUTFILE,"> $global_tmp_directory/sort_peak_cmd_list");
print OUTFILE join("",@sort_peak_cmd);
close OUTFILE;


open(OUTFILE,"> $global_tmp_directory/bigbed_cmd_list");
print OUTFILE join("",@bigbed_cmd);
close OUTFILE;



`parallel --no-notice -j 20 < $global_tmp_directory/macs2_cmd_list`;
`parallel --no-notice -j 20 < $global_tmp_directory/bdgcmp_cmd_list`;
`parallel --no-notice -j 20 < $global_tmp_directory/sort_cmd_list`;
`parallel --no-notice -j 20 < $global_tmp_directory/bigwig_cmd_list`;
`parallel --no-notice -j 20 < $global_tmp_directory/mv_peak_cmd_list`;
`parallel --no-notice -j 20 < $global_tmp_directory/sort_peak_cmd_list`;
`parallel --no-notice -j 20 < $global_tmp_directory/bigbed_cmd_list`;



 











