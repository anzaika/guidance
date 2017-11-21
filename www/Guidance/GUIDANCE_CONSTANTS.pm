#!/usr/bin/perl

package GUIDANCE_CONSTANTS; #don't forget: a package must end with a return value (1; in the end)!!!!!

use GENERAL_CONSTANTS;

# constants to use when sending e-mails using the server admin's email address.
use constant PDB_FILE_NAME =>           "pdb_file.ent";
use constant USER_MSA_FILE_NAME =>      "msa_file.msa";
use constant USER_TREE_FILE_NAME =>     "user_tree.txt";
use constant MAX_UPLOAD_FILE_SIZE =>    20;

use constant GUIDANCE_RESULTS_DIR => GENERAL_CONSTANTS::SERVERS_RESULTS_DIR."Guidance/";
use constant GUIDANCE_LOGS_DIR => GENERAL_CONSTANTS::SERVERS_LOGS_DIR."Guidance/";
use constant GUIDANCE_RESULTS_URL => GENERAL_CONSTANTS::GUIDANCE_URL."results/";

use constant GUIDANCE_COL_FP_TP_SP_CUTOFFS => "/bioseq/Guidance/exec/mergeBpCos.columnScores.perf.100";


use constant BOOTSTRAPS=> 100;

use constant MAXIMUM_MODIFIED_PERCENT => 0.15;
use constant FRAGMENT_REDUNDANCY_RATE => 95;
use constant FRAGMENT_OVERLAP => 0.10;
use constant FRAGMENT_MINIMUM_LENGTH => 0.60;
use constant MINIMUM_FRAGMENTS_FOR_MSA => 5;
use constant LOW_NUM_FRAGMENTS_FOR_MSA => 10;

use constant PDB_FRAGMENT_MINIMUM_LENGTH => 0.50;
use constant PDB_FRAGMENT_MINIMUM_IDENTITY => 0.35;

use constant BAYES_INTERVAL => 3;

use constant CONSURF_PDF_LIB => "/var/www/html/ConSurf/php/ConSurf_PDF/ConSurf_PDF_Class.php";

use constant PROC_NUM => 2; #num of processors to use

1;
