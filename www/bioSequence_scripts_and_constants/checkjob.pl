#!/usr/bin/perl -w

my $ans_out_file = "checkjob.out";
my $ans_err_file = "checkjob.err";
$qsub_job_no = 124222;
$error = 0;
$cmd = "checkjob $qsub_job_no >$ans_out_file 2>$ans_err_file";
`$cmd`;
#open ANS, $ans_res_file;
#while (<ANS>){
#    if (/ERROR/){
#        print "error was found!\n";
#        $error = 1;
#        last;
#    }
#}
#if ($error==1){
#    print "job $qsub_job_no does not exist\n";
#}
#else{
#    print "job $qsub_job_no exist\n";
#}