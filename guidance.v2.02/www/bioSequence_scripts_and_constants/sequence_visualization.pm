#!/usr/local/bin/perl
package sequence_visualization; 
use strict;

sub create_sequence_visualization_html_page
{
    #input is:
    #1. server: ConSeq/Selecton/Epitopia
    #2. input_file: full path to the input file which is:
        #for ConSeq: query.grades file
        #for Selecton: selection4Site.txt file
        #for Epitopia: protein_seq_immunogenic_residues.txt or protein_seq_immunogenic_regions.txt, depending on argument 5
    #3.out_dir: where to print the output to. output file is always out_dir/sequence_colored.html. ends with a "/"
    #4.run_name: the run id name
    #5.output file name: the name of the html output file to be created

    my $server = $_[0];
    my $input_file = $_[1];
    my $out_dir = $_[2];
    my $run_name = $_[3];
    my $output_file = $_[4];

    my (@seq_indices);
    my (%aa_index2info,%category2color);
    if($server eq "ConSeq"){
	&get_conseq_data($input_file,\%aa_index2info,\@seq_indices);
    }
    if($server eq "Selecton"){
	&get_selecton_data($input_file,\%aa_index2info,\@seq_indices);
    }
    if($server eq "Epitopia"){
	if(@_ == 5){
	    &get_epitopia_data($input_file,\%aa_index2info,\@seq_indices);
	}
    }
    &get_colors($server,\%category2color);
    &create_html($server,$out_dir.$output_file,$run_name,\%aa_index2info,\@seq_indices,\%category2color);
    return (1);
}

sub get_conseq_data
{
    my ($line,$aa_index);
    my (@info);
    my $file = $_[0];
    my %aa_index2info = %{$_[1]};
    my @seq_indices = @{$_[2]};
    open INPUT, "< $file";
    while($line=<INPUT>){
	chomp($line);
	if($line =~ m/^((\s+)?)(\d+)(\s+)(\w)(\s+)(\-?\d+\.\d+)(\s+)(\d)(\s+)(\w)(\s+)(\w?)(\s+)(\d+)(\/)(\d+)(.*)/){
	    $aa_index = $3;
	    push(@info,$5);
	    push(@info,$9);
	    push(@info,$11);
	    if(($15 < 6)or(($17 > 60)and($15/$17 < 0.1))){
	       push(@info,"i");
	    }
	    $aa_index2info{$aa_index} = [ @info ];
	    @info = ();
	    push(@seq_indices,$aa_index);
	}
    }
    close(INPUT);
    %{$_[1]} = %aa_index2info;
    @{$_[2]} = @seq_indices;
    return (1);
}

sub get_selecton_data
{
    my ($line,$aa_index);
    my (@info);
    my $file = $_[0];
    my %aa_index2info = %{$_[1]};
    my @seq_indices = @{$_[2]};
    open INPUT, "< $file";
    while($line=<INPUT>){
	chomp($line);
	next if($line =~ m/^$/);
	@info = split(/\s+/,$line);
	$aa_index = shift(@info);
	$aa_index2info{$aa_index} = [ @info ];
	@info = ();
	push(@seq_indices,$aa_index);
    }
    close(INPUT);
    %{$_[1]} = %aa_index2info;
    @{$_[2]} = @seq_indices;
    return (1);
}

sub get_epitopia_data
{
    my ($line,$aa_index,$counter,$index,$border,$out,$i,$j,$colour,$num,$total,$limit,$residue,$div,$aa_name,$acc_pred,$category,$color_category);
    my (@info,@array);
    my (%residuesMap);
    my $file = $_[0];
    my %aa_index2info = %{$_[1]};
    my @seq_indices = @{$_[2]};
    $counter = 0;
    open INPUT, "< $file";
    while($line=<INPUT>){
	chomp($line);
	next if(($line =~ m/^$/)or($line =~ m/^Residue/));
	$counter += 1;
	@array = split(/\s+/,$line);
	$aa_index = substr($array[0],3,length($array[0])-3);
	$aa_name = &string2char(substr($array[0],0,3));
	$acc_pred = $array[@array-1];
	$acc_pred =~ tr/A-Z/a-z/;
	$residuesMap{$aa_index} = [ ($aa_name,$counter,$acc_pred) ];
    }
    close(INPUT);

    $div = $counter/5;
    if($div =~ m/^(\d+)(\.)(\d+)/){
	$border = $1;
	$border += 1;
    }
    else{
	$border = $div;
    }
    for($aa_index = 1;$aa_index <  keys %residuesMap;++$aa_index){
	$counter = @{$residuesMap{$aa_index}}[1];
	$category = 5;
	for(my $i = $border;$i <= 5*$border;$i+=$border){
	    if($counter <= $i){
		$color_category = $category;
		last;
	    }
	    $category -= 1;
	}
	$aa_name = @{$residuesMap{$aa_index}}[0];
	$acc_pred = @{$residuesMap{$aa_index}}[2];
	@info = ($aa_name,$color_category,$acc_pred);
	$aa_index2info{$aa_index} = [ @info ];
	@info = ();
	push(@seq_indices,$aa_index);
    }
    %{$_[1]} = %aa_index2info;
    @{$_[2]} = @seq_indices;
    return (1);
}

sub get_colors
{
    my $server = $_[0];
    my %category2color = %{$_[1]};
    if($server eq "ConSeq"){
	$category2color{1} = "#10C8D1";
	$category2color{2} = "#8CFFFF";
	$category2color{3} = "#D7FFFF";
	$category2color{4} = "#EAFFFF";
	$category2color{5} = "#FFFFFF";
	$category2color{6} = "#FCEDF4";
	$category2color{7} = "#FAC9DE";
	$category2color{8} = "#F07DAB";
	$category2color{9} = "#A02560";
    }
    elsif($server eq "Selecton"){
	$category2color{1} = "#FFBD00";
	$category2color{2} = "#FFFF78";
	$category2color{3} = "#FFFFFF";
	$category2color{4} = "#FCEDF4";
	$category2color{5} = "#FAC9DE";
	$category2color{6} = "#F07DAB";
	$category2color{7} = "#A02560";
    }
    elsif($server eq "Epitopia"){
	$category2color{1} = "#10C8D1";
	$category2color{2} = "#D7FFFF";
	$category2color{3} = "#FFFFFF";
	$category2color{4} = "#FAC9DE";
	$category2color{5} = "#A02560";
    }
    %{$_[1]} = %category2color;
    return 1;
}


sub create_html
{
    my ($aa_index,$count,$count_num,$space_num,$spaces,$aa_name,$color_category,$aa_color,$aa_print,$acc_print,$func_print,$category,$title,$width);
    
    my $server = $_[0];
    my $out_file = $_[1];
    my $run_name = $_[2];
    my %aa_index2info = %{$_[3]};
    my @seq_indices = @{$_[4]};
    my %category2color = %{$_[5]};
    
    open OUTPUT, "> $out_file";
    print OUTPUT "<html>\n<title>$server Results: $run_name</title>\n"; 
    print OUTPUT "<body bgcolor='white'>\n";
    print OUTPUT "<H1 align=center><u>$server Results</u></H1>\n\n";

    print OUTPUT "\n<table border=0 width=750>\n";
    print OUTPUT "<tr><td>\n";

    $aa_print = "";
    $acc_print = "";
    $func_print = "";
    $count = 1;
    foreach $aa_index (@seq_indices){
	$aa_name = @{$aa_index2info{$aa_index}}[0];
	$color_category = @{$aa_index2info{$aa_index}}[1];
	$aa_color = $category2color{$color_category};
	# print the counter above the beginning of each 10 characters
	if ($count % 50 == 1){
	    $count_num = $count;
	    for ($count_num ; $count_num < $count+50 ; $count_num += 10){
		if ($count_num <= @seq_indices){
		    $space_num = 11 - length($count_num);
		    $spaces = "&nbsp;" x $space_num;
		    print OUTPUT "<font face='Courier New' color='black' size=+1>".$count_num.$spaces."</font>";
		}
	    }
	    print OUTPUT "<br>\n";
	}
	### print the colored letters and 'b/e' for the buried/exposed residues (ConSeq and Epitopia)
	# after 50 characters - print newline
	if ($count % 50 == 0 or $count == @seq_indices){
	    #print aa name and corresponding background color
	    if($server eq "ConSeq"){
		if(@{$aa_index2info{$aa_index}} == 4){
		    $aa_print .= "<b><font face='Courier New' color='yellow' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b><br>";
		}
		else{
		    if($color_category == (keys %category2color)){
			$aa_print .= "<b><font face='Courier New' color='white' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b><br>\n";
		    }
		    else{
			$aa_print .= "<b><font face='Courier New' color='black' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b><br>\n";
		    }
		}
	    }
	    else{
		if($color_category == (keys %category2color)){
		    $aa_print .= "<b><font face='Courier New' color='white' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b><br>\n";
		}
		else{
		    $aa_print .= "<b><font face='Courier New' color='black' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b><br>\n";
		}
	    }
	    #print exposed/buried and functional/structural
	    if(($server eq "ConSeq")or($server eq "Epitopia")){
		if(@{$aa_index2info{$aa_index}}[2] eq 'e'){#exposed
		    $acc_print .= "<b><font face='Courier New' color='orange' size=+1>e</font></b><br>\n";
		    if($server eq "ConSeq"){
			if((@{$aa_index2info{$aa_index}}[1] == (keys %category2color))or(@{$aa_index2info{$aa_index}}[1] == ((keys %category2color)-1))){#functional
			    $func_print .= "<b><font face='Courier New' color='red' size=+1>f</font></b>\n";
			}
			else{
			    $func_print .= "<font face='Courier New' size=+1>&nbsp;</font>\n";
			}
		    }
		}
		else{#buried
		    $acc_print .= "<b><font face='Courier New' color='#00cc00' size=+1>b</font></b><br>\n";
		    if($server eq "ConSeq"){
			if(@{$aa_index2info{$aa_index}}[1] == (keys %category2color)){#structural
			    $func_print .= "<b><font face='Courier New' color='#000099' size=+1>s</font></b>\n";
			}
			else{
			    $func_print .= "<font face='Courier New' size=+1>&nbsp;</font>\n";
			}
		    }
		}
	    }
	    print OUTPUT $aa_print;
	    if($server ne "Selecton"){
		print OUTPUT $acc_print;
		if($server eq "ConSeq"){
		    print OUTPUT $func_print;
		}
	    }
	    print OUTPUT "</td></tr>\n";
	    print OUTPUT "<tr><td>\n";
	    $aa_print = "";
	    $acc_print = "";
	    $func_print = "";
	}
	# after 10 characters - print a space ('&nbsp;')
	elsif ($count % 10 == 0) {
	    if($server eq "ConSeq"){
		if(@{$aa_index2info{$aa_index}} == 4){
		    $aa_print .= "<b><font face='Courier New' color='yellow' size=+1><span style='background: $aa_color;'>$aa_name</span> </font></b>\n";
		}
		else{
		    if($color_category == (keys %category2color)){
			$aa_print .= "<b><font face='Courier New' color='white' size=+1><span style='background: $aa_color;'>$aa_name</span> </font></b>\n";
		    }
		    else{
			$aa_print .= "<b><font face='Courier New' color='black' size=+1><span style='background: $aa_color;'>$aa_name</span> </font></b>\n";
		    }
		}
	    }
	    else{
		if($color_category == (keys %category2color)){
		    $aa_print .= "<b><font face='Courier New' color='white' size=+1><span style='background: $aa_color;'>$aa_name</span> </font></b>\n";
		}
		else{
		    $aa_print .= "<b><font face='Courier New' color='black' size=+1><span style='background: $aa_color;'>$aa_name</span> </font></b>\n";
		}
	    }
	    if(($server eq "ConSeq")or($server eq "Epitopia")){
		if(@{$aa_index2info{$aa_index}}[2] eq 'e'){
		    $acc_print .= "<b><font face='Courier New' color='orange' size=+1>e&nbsp;</font></b>";
		    if($server eq "ConSeq"){
			if((@{$aa_index2info{$aa_index}}[1] == (keys %category2color))or(@{$aa_index2info{$aa_index}}[1] == ((keys %category2color)-1))){
			    $func_print .= "<b><font face='Courier New' color='red' size=+1>f&nbsp;</font></b>";
			}
			else{
			    $func_print .= "<font face='Courier New' size=+1>&nbsp;&nbsp;</font>";
			}
		    }
		    
		}
		else{
		    $acc_print .= "<b><font face='Courier New' color='#00cc00' size=+1>b&nbsp;</font></b>";
		    if($server eq "ConSeq"){
			if(@{$aa_index2info{$aa_index}}[1] == (keys %category2color)){
			    $func_print .= "<b><font face='Courier New' color='#000099' size=+1>s&nbsp;</font></b>";
			}
			else{
			    $func_print .= "<font face='Courier New' size=+1>&nbsp;&nbsp;</font>";
			}
		    }
		}
	    }
	}
	else{
	    if($server eq "ConSeq"){
		if(@{$aa_index2info{$aa_index}} == 4){
		    $aa_print .= "<b><font face='Courier New' color='yellow' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b>";
		}
		else{
		    if($color_category == (keys %category2color)){
			$aa_print .= "<b><font face='Courier New' color='white' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b>";
		    }
		    else{
			$aa_print .= "<b><font face='Courier New' color='black' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b>";
		    }
		}
	    }
	    else{
		if($color_category == (keys %category2color)){
		    $aa_print .= "<b><font face='Courier New' color='white' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b>";
		}
		else{
		    $aa_print .= "<b><font face='Courier New' color='black' size=+1><span style='background: $aa_color;'>$aa_name</span></font></b>";
		}
	    }
	    if(($server eq "ConSeq")or($server eq "Epitopia")){
		if(@{$aa_index2info{$aa_index}}[2] eq 'e'){
		    $acc_print .= "<b><font face='Courier New' color='orange' size=+1>e</font></b>";
		    if($server eq "ConSeq"){
			if((@{$aa_index2info{$aa_index}}[1] == (keys %category2color))or(@{$aa_index2info{$aa_index}}[1] == ((keys %category2color)-1))){
			    $func_print .= "<b><font face='Courier New' color='red' size=+1>f</font></b>";
			}
			else{
			    $func_print .= "<font face='Courier New' size=+1>&nbsp;</font>";
			}
		    }
		}
		else{
		    $acc_print .= "<b><font face='Courier New' color='#00cc00' size=+1>b</font></b>";
		    if($server eq "ConSeq"){
			if(@{$aa_index2info{$aa_index}}[1] == (keys %category2color)){
			    $func_print .= "<b><font face='Courier New' color='#000099' size=+1>s</font></b>";
			}
			else{
			    $func_print .= "<font face='Courier New' size=+1>&nbsp;</font>";
			}
		    }
		}
	    }
	}
    	$count++;
    }
    print OUTPUT "</td></tr>\n</table><br>\n";
    if($server eq "ConSeq"){
	$title ="The conservation scale";
	$width = 310;
    }
    elsif($server eq "Selecton"){
	$title ="The selection scale";
	$width = 240;
    }
    elsif($server eq "Epitopia"){
	$title ="The immunogenicity scale";
	$width = 180;
    }
    print OUTPUT "\n<br><b><u>Legend:</u><br><br>\n$title:</b><br>\n<table border=0 cols=1 width=$width>\n<tr><td align=center>\n<font face='Courier New' color='black' size=+1><center>\n";
    for(my $i = 1;$i <= keys %category2color;++$i){
	if($i == (keys %category2color)){
	    print OUTPUT "<font color='white'><span style='background: $category2color{$i};'>&nbsp;$i&nbsp;</span></font>";
	}
	else{
	    print OUTPUT "<span style='background: $category2color{$i};'>&nbsp;$i&nbsp;</span>";
	}
    }
    if($server eq "ConSeq"){
	print OUTPUT "</font></center>\n<center><table border=0 cols=3 width=$width>\n<tr>\n<td align=left><b>Variable</b></td>\n<td align=center><b>Average</b></td>\n<td align=right><b>Conserved</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";
    }
    if($server eq "Selecton"){
	print OUTPUT "</font></center>\n<center><table border=0 cols=3 width=$width>\n<tr>\n<td align=left><b>Positive</b></td>\n<td align=right><b>Purifying</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";
    }
    if($server eq "Epitopia"){
	print OUTPUT "</font></center>\n<center><table border=0 cols=3 width=$width>\n<tr>\n<td align=left><b>Low</b></td>\n<td align=center><b>Average</b></td>\n<td align=right><b>High</b></td>\n</tr>\n</table></center>\n</td>\n</tr>\n</table>\n";
    }
    if(($server eq "ConSeq")or($server eq "Epitopia")){
	print OUTPUT "<b><font face='Courier New' color='orange' size=+1>e</font> - A predicted exposed residue.</b><br>\n";
	print OUTPUT "<b><font face='Courier New' color='#00cc00' size=+1>b</font> - A predicted buried residue.</b><br>\n";
    }
    if($server eq "ConSeq"){
	print OUTPUT "<b><font face='Courier New' color='red' size=+1>f</font> - A predicted functional residue (highly conserved and exposed).</b><br>\n";
	print OUTPUT "<b><font face='Courier New' color='#000099' size=+1>s</font> - A predicted structural residue (highly conserved and buried).</b><br>\n";
	print OUTPUT "<b><font face='Courier New' color='yellow' size=+1><span style='background: $category2color{1};'>X</span></font> - Insufficient data - the calculation for this site was performed on less than 10% of the sequences.</font></b><br>\n"; 
    }
    print OUTPUT "</body>\n</html>\n";
    close (OUTPUT);
    return (1);
}


sub string2char
{
    my $result;
    if($_[0] eq "ALA"){
	$result = "A";
    }
    elsif($_[0] eq "ARG"){
	$result = "R";
    }
    elsif($_[0] eq "ASN"){
	$result = "N";
    }
    elsif($_[0] eq "ASP"){
	$result = "D";
    }
    elsif($_[0] eq "CYS"){
	$result = "C";
    }
    elsif($_[0] eq "GLU"){
	$result = "E";
    }
    elsif($_[0] eq "GLN"){
	$result = "Q";
    }
    elsif($_[0] eq "GLY"){
	$result = "G";
    }
    elsif($_[0] eq "HIS"){
	$result = "H";
    }
    elsif($_[0] eq "ILE"){
	$result = "I";
    }
    elsif($_[0] eq "LEU"){
	$result = "L";
    }
    elsif($_[0] eq "LYS"){
	$result = "K";
    }
    elsif($_[0] eq "MET"){
	$result = "M";
    }
    elsif($_[0] eq "PHE"){
	$result = "F";
    }
    elsif($_[0] eq "PRO"){
	$result = "P";
    }
    elsif($_[0] eq "SER"){
	$result = "S";
    }
    elsif($_[0] eq "THR"){
	$result = "T";
    }
    elsif($_[0] eq "TRP"){
	$result = "W";
    }
    elsif($_[0] eq "TYR"){
	$result = "Y";
    }
    elsif($_[0] eq "VAL"){
	$result = "V";
    }
    return ($result);
}
1;
