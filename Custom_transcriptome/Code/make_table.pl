#!/bin/perl


use strict;

my @Table;


my @files = <../BAM/*.summary.txt>;
my $i = 0;
foreach my $file (@files) 
  {
    $Table[$i]{'File'} = $file;
    $Table[$i]{'File'} =~ s/^...//g;
    $Table[$i]{'Stem'} = $Table[$i]{'File'};
    $Table[$i]{'Stem'} =~ s/_trimmed.hisat2.summary.txt//g;
    $Table[$i]{'Stem'} =~ s/BAM.//g;

    my $unique = `grep 'aligned exactly 1 time' $file`;
    chomp($unique);
    $unique =~ s/.*\(//g;
    $unique =~ s/\).*//g;
    $Table[$i]{'unique'} = $unique;

    my $overall = `grep 'overall alignment rate' $file`;
    chomp($overall);
    $overall =~ s/.*\(//g;
    $overall =~ s/ overall.*//g;
    $Table[$i]{'overall'} = $overall;

    $Table[$i]{'FeatureCounts'} = "[FeatureCounts](FeatureCounts/" . $Table[$i]{'Stem'} . 
                                  "_trimmed.bam_featureCounts_counts.txt)";
    $Table[$i]{'GTF'}           = "[GTF](GTF/" . $Table[$i]{'Stem'} . 
                                  "_trimmed.gtf)";
    $i++;
  }



printf "| Sample | Alignment (Unique) | Alignment (Overall) | FeatureCounts | GTF |\n";
print  "| ------ | ------------------ | ------------------- | ------------- | --- |\n";

for($i=0; $i<=$#Table; $i++)
  {
    printf "| %s | %s | %s | [%s] | [%s] |\n", 
           $Table[$i]{'Stem'}, $Table[$i]{'unique'}, $Table[$i]{'overall'},
           $Table[$i]{'FeatureCounts'}, $Table[$i]{'GTF'}; 
  }


#------------------------------------------------------------------------------
# END OF SCRIPT
#------------------------------------------------------------------------------
