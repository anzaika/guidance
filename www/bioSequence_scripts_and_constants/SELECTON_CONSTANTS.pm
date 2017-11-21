#!/usr/bin/perl

package SELECTON_CONSTANTS;

use GENERAL_CONSTANTS;
use constant STATISTICS_FILE => "/bioseq/Selecton/selecton_models_statistics.log";
# FOR BLAST
use constant MAXIMUM_MODIFIED_PERCENT => 0.15;
use constant FRAGMENT_REDUNDANCY_RATE => 95;
use constant FRAGMENT_OVERLAP => 0.10;
use constant FRAGMENT_MINIMUM_LENGTH => 0.60;
use constant MINIMUM_FRAGMENTS_FOR_MSA => 10;
use constant LOW_NUM_FRAGMENTS_FOR_MSA => 10;

use constant GB_CDS_INDEX_FILE => "/bioseq/data/results/GB_CDS/All_GB_CDS.Without_CON.bp_index";
1;
