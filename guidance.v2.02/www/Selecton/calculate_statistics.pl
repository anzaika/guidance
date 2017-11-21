#!/usr/bin/perl -w
use strict;

my $stat_file = "/bioseq/Selecton/selecton_models_statistics.log";
my $out_file = "/bioseq/Selecton/total_models_statistics.csv";
#my $stat_file = "D:\\My Documents\\TAU\\toRemove\\selecton\\selecton_models_statistics.log";


my %stat_hash = ();
open STAT, $stat_file;
while (<STAT>){
    chomp;
    if (/^(\d+) total runTime: (.+)/){
        if (exists $stat_hash{$1}){
            $stat_hash{$1}{RUN_TIME} = $2;
            $stat_hash{$1}{FINISH} = 1;
        }
    }
    elsif(/^(\d+) (\S+) (\d+) (\d+)/){
        $stat_hash{$1} = {MODEL=>$2, SEQ_NUM=>$3, LENGTH=>$4};
        $stat_hash{$1}{FINISH} = 0;
    }
}
close STAT;

my $int=0;
my @m8 = ();
my @m7 = ();
my @m5 = ();
my @m8a = ();
my @mec = ();

while ( my ($key, $value) = each(%stat_hash) ) {        
    if ($stat_hash{$key}{FINISH}){
        my $t = $stat_hash{$key}{RUN_TIME};
        $t =~ /(\d+):(\d+):(\d+)/;
        my $hour = $1;
        my $min = $2;
        my @arr = ($stat_hash{$key}{SEQ_NUM}, $stat_hash{$key}{LENGTH}, $hour, $min);
        if ($stat_hash{$key}{MODEL} eq 'M8'){            
            push @m8, [ @arr ];
        }
        elsif ($stat_hash{$key}{MODEL} eq 'MEC'){            
            push @mec, [ @arr ];
        }
        elsif ($stat_hash{$key}{MODEL} eq 'M7'){
            push @m7, [ @arr ];
        }
        elsif ($stat_hash{$key}{MODEL} eq 'M5'){
            push @m5, [ @arr ];
        }
        elsif ($stat_hash{$key}{MODEL} eq 'M8a'){
            push @m8a, [ @arr ];
        }
    }        
}
open OUT, ">$out_file";
print OUT "----m5----\n";
print_arrays(\@m5);
print OUT "----m7----\n";
print_arrays(\@m7);
print OUT "----m8----\n";
print_arrays(\@m8);
print OUT "----m8a----\n";
print_arrays(\@m8a);
print OUT "----mec----\n";
print_arrays(\@mec);
close OUT;


sub print_arrays{
    my $arr_ref = shift;
    my @arr_to_print = @$arr_ref;
    print OUT "NUM_OF_SEQ,LENGTH,HOUR,MINUTE\n";
    for (my $i = 0; $i < @arr_to_print; $i++) {
        for (my $j=0; $j<@{$arr_to_print[$i]}; $j++){
            print OUT $arr_to_print[$i][$j].",";
        }
        print OUT "\n";
    }    
}
#print $m8[0][0];


#print $m8[0]." ";
#print $m8[1]." ";
#print $m8[2];

#for my $hash ( @m8 ) {
#    while ( my ($key, $value) = each(%$hash) ) {
#        print "$key: $value,";
#    }
#    print "\n";
#}

#while ( my ($key, $value) = each(%stat_hash) ) {        
#        if ($stat_hash{$key}{FINISH}){
#            print "$key ; ";
#            while (my ($key2, $val2) = each (%$value)){
#                if ($key2 ne "FINISH"){
#                    print "$key2: $val2, ";
#                }
#            }
#            print "\n";
#        }        
#}