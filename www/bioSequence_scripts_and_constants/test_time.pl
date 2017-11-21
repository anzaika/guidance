
#!/usr/bin/perl -w
use BIOSEQUENCE_FUNCTIONS;
use strict;


my $ans = &BIOSEQUENCE_FUNCTIONS::subtract_time_from_now("11:53:00", "29-10-2007");
print $ans;