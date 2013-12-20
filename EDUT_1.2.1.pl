#!/usr/bin/perl

	##########################################################################################	
			# include the number of individuals in the actual binary data for each site.
			# standard number is the number of sequences. The number counted for the non-gapped and non-missing data sites.
			#gapped and missing should be used to adjust the priority score
    ##########################################################################################	

use Getopt::Std;
%options=();
getopts("a:i:bd:vmsg:n:e:", \%options) or usage();
	# -a followed by the filename for the ms formatted simulation output
	# -i followed by the filename for the input fasta format - required.
	#    the fasta file should be formatted so that there are no line
	#    breaks within the sequence data. A line break should appear after the
	#    sequence name, after each sequence and after the last sequence. This
	#    version also requires all capital letters for the nucleotides. 
	# either -a or -i input is required.
	# -b indicates print binary data
	# -d followed by the delimiter for the binary data output.
	#    tab is the default
	# 	 space
	#    none
	# -v indicates print verbose output
    # -m indicates printing a moderate/reasonable amount of output
	# -s indicates that the output should be only a short summary
	# The following options are new to EDUT_1.1:
	# -g use 1 to exclude the gapped sites 
    #    use 0 include the gapped sites (the default)
    # -n use 1 to exclude the sites where sequence is missing indicated by 'N'
    #	 use 0 include the sites where sequence is missing (the default)
	
##########################################################################
	
	#print $options{i};
	#print "\n";
	#print $options{n};
	#print "\n";
	$missing_counter = 0;
	$gap_counter = 0;
	$specific_excluded = 0;
	$info_array = collect_input(\%options);
	my ($ref_input, $length, $seq_id_count, $missing_data, $gapped_seq, $seq_name) = @$info_array;
	if (!($options{a}) and !($options{i})) 
	{
	 print "\nAn input file is required.\nInput may be a fasta file with option -i\nor an ms formatted file from a simulated data set with option -a\n";
	 usage();
	 
	}
	if ($options{a})
{
	print "#Simulation data was entered\n#\n";
	print "#The length of the simulated sequences is $length.\n";
	print "#The number of sequences in the simulation is $seq_id_count.\n";
	for (my $count1=0; $count1< $seq_id_count; $count1++)
		{
		for (my $count2= 0 ; $count2 < $length; $count2++)
			{
			$segsites[$count2] = $count2+1; #remove when position data is used. For sim data every site is segregating.
			#print $ref_input -> [$count1][$count2];#printing the binary
			$new_bin[$count1][$count2+1] = $ref_input -> [$count1][$count2];
			}
		#print "\n";#print a new line between the haplotypes
		}
}
	if (($options{v}) or ($options{m}))
	{
	print "#length of sequence is $length and there are $seq_id_count sequences.\n";
	if ($missing_data) 
		{
		print "Some base calls are represented with an N for missing data.\n\n"
		}	
	if ($gapped_seq) 
		{
		 print "Some base calls are gapped representing historical indel events.\n\n"
		}	
	}
	
 ############################
	my $counter = 0;
	for (my $position=0; $position < $length; $position++)
	{
		@current_position = ();
      	for ($i= 1; $i <= $seq_id_count; $i++)
	    	{	
        	#create an array at each site and find the segregating sites
	    	push(@current_position,$ref_input ->[$i][$position]);
			}
		$counter ++;
		push (@current_position, \%options);
		if (($options{e}) == ($position +1))
		{
			$specific_excluded = true;
		}
		else 
		{
			$specific_excluded = 0;
		}
		$site_info = &each_position(@current_position, $specific_excluded);
		$majority_base = pop (@$site_info);
		#print "The majority base is here:", $majority_base, "\n";
		#print "\n\n"; 
		#print "current position array: @current_position \n\n";
		($missing_excluded, $gap_excluded, $singleton, $illegal_char, $informative_seg_site) = @$site_info;
		if ($missing_excluded)
		{
			$missing_counter++;
		}
		if ($gap_excluded)
		{
			$gap_counter++;
		}
		if ($specific_excluded)
		{
			print "\nSite $options{e} has been excluded\n";
		}
		if ($majority_base)
		{
			$junk++;
			$column_info = &preliminary_binary($majority_base,\@current_position, $site_info, $position);
			$site_position = pop @$column_info;
		    $binary_column = pop @$column_info;
			push (@segsites, $site_position);
			@{$binary[$junk]} = @$binary_column;#starts at 1
			$num_elements_in_binary = @binary;
		}
		
	}	
	if ($options{v})
	{
	print "#Segregating sites are at positions: @segsites \n";
	}
foreach $pos(@segsites) 
# the scalar variable $pos represents the actual position (starting at 0) of the 
# segregating site in the original dataset. 
	{
	#print "The position is $pos\n";
    $count++;
	$hash_for_segsites{$pos} = $count; #count represents the number of the segsite.
	$another_segsite_hash{$count} =  $pos;#count starts at 1
	}
	if ($options{a})
	{ $num_seg_sites = $length;}
	else
	{
	$num_seg_sites = $count;
	}
if (($options{n}) or ($options{g}))
{
	print "\n\n";
	print "# $missing_counter missing sites or sites with illegal characters have been excluded from the analysis.\n\n";
	print "# $gap_counter gapped sites have been excluded from the analysis. \n\n";
	print "#  **Warning**: Excluding all missing or gapped sites may exclude some sites with segregating variation.\n\n";
}
print "\n\n# There are $num_seg_sites sites considered.\n#\n";


my $base;
$col_num = 0;
foreach $site (@binary)
{
	$seq_id = 0;
	foreach $base(@$site)
		{
		#print "A: sequence id $seq_id and position number $col_num. Store in the array $base.\n";
		$new_bin[$seq_id][$col_num] = $base;
		#print "Check here for the base $base\n";
		$seq_id++;
		}
    #print "A: The segsite number is $col_num.\n";
	$col_num++;
 }
$informative_site_counts = &process_columns($seq_id_count, $num_seg_sites,\@new_bin,\@segsites,\%hash_for_segsites);
($informative_sites, $non_informative, $singleton_count, $only_gapped_count, $new_segsite_hash_ref) = @$informative_site_counts;


print "# $informative_sites are parsimony informative and $non_informative sites are not parsimony informative.\n";
print "# $singleton_count are singleton mutations and $only_gapped_count sites are in indel polymorphisms with no segregating bases.\n";
#check if ($seq_id == $seq_id_count)
$count_info = &triplets(\@new_bin, \@segsites, $seq_id,\%hash_for_segsites);
#print "This is the big array:", @$count_info; #sites and counts from the triplets subroutine.
#print "\n";
# Look for the minimum cell count for the first triplet site pair
$a_count = pop @$count_info; #from the last array element passed from the triplets sub routine.
print "#\n# Number of pattern a triplets is $a_count\n";
if ($informative_sites == 0)
{
	print "\nThere are no informative SNP sites in the specified input data.\n";
	print "Check the input file format.\n";
	exit;
}
if (($informative_sites == 1) or ($informative_sites == 2))
{
	print "\n  #  There are no SNP site triplets in the specified input data\n  #  because there are fewer than 3 informative SNPs.\n";
	print "  #  Check the input file format.\n\n\n";
	exit;
}
$proportion = $a_count/(($informative_sites * ($informative_sites - 1)* ($informative_sites - 2))/ 6 );
print "# The proportion of pattern a triplets is $proportion\n#\n";
#if ($proportion > 0.28)
#{
#	print "# This proportion of pattern a is unlikely.\n#\n";
#}

$triplets_each_site = (($informative_sites * ($informative_sites - 1))/2);
$summary = &summarize($a_count,\@new_bin, $seq_id, \%hash_for_segsites, $seq_name, \%options);
$middle_correction = &how_much_middle_site($new_segsite_hash_ref, $informative_sites, $summary);
if ($options -> {s})
	{
	$summary_array_ref = &print_summary($triplets_each_site, $summary, \@segsites, $seq_name, $seq_id, $proportion,$middle_correction);
	}
#send a reference to the summary hash to a subroutine and then loop through all of the segregating sites.
#sort the summary numbers and then report the sites in the order that they should be checked.

#################################################################################################
sub how_much_middle_site
{
($hash_reference, $num_informative_sites, $summary_ref) = @_;
foreach $site(keys(%$summary_ref)) 
	{
	#print "$site is the site\n";
	#print "\n";
	$x1 = $hash_reference -> {$site};
	#print "$x1 out of $num_informative_sites sites\n"; 
	$segsite_middle_correction{$site} = (($x1 - 1) * ($num_informative_sites - $x1)); #where x is the position of the segsite (start counting at 1) in relation to all of the informative segsites and S is the total number of segsites.
	#print "The correction for the extra possible flags due to a middle triplet site is  $segsite_middle_correction{$site} for site number $site\n"; 
	}
return \%segsite_middle_correction;	
}
##########################################################################

sub summarize
{
	$CUTOFF_FREQ = 1;
	($a_count, $bin_ref, $seq_id, $hash_reference, $seq_name, $options) = @_;
	for (my $a_num = 0; $a_num < $a_count; $a_num++)
	{
		%pair_hash= ((2 + ($a_num * 12)) => '00', # Create a string representing the rare class pair
			(3 + ($a_num * 12)) => '10',  # match this string with every cell number in the big
			(4 + ($a_num * 12)) => '01',  # array skipping two of every six entries. 
			(5 + ($a_num * 12)) => '11',  # 
			(8 + ($a_num * 12)) => '00', 
			(9 + ($a_num * 12)) => '10', 
			(10 + ($a_num * 12)) => '01',
			(11 + ($a_num * 12)) => '11');
		 
		# Look for the minimum cell count for the first triplet site pair
		@pair_code = (0,6);
		foreach $pair (@pair_code)
		{
			for (my $cell_number = (2 + $pair + ($a_num)*12); $cell_number < (6 + $pair + ($a_num)*12) ;$cell_number++)
			{
				if (($count_info -> [$cell_number]) <=  $CUTOFF_FREQ)
				{
				#print "pattern a number $a_num\n";
				#print "Cell number $cell_number \n";
				my $the_rare_string = $pair_hash{$cell_number};
				my $seg_site1 = $segsites[$count_info -> [int(($cell_number-0.5)/6)*6]];
				my $seg_site2 = $segsites[$count_info -> [int(($cell_number-0.5)/6)*6 + 1]];
				#print "The two sites to the flag sequences routine one: $seg_site1 two: $seg_site2\n";
				$seq_num = &flag_sequences($the_rare_string, $seg_site1, $seg_site2, $bin_ref, $seq_id,$hash_reference);
				@flagged = ();
				if (($options -> {v}) or ($options -> {m}))
					{
					print "User should check sequence $label_key";
					#print  @flagged ;
					if ($options{a})
					{
					$name = $seq_num;
					print "number $name";
					}
					else
					{
					$name = $seq_name -> {$seq_num};
					print "labeled '$name'";
					}
					
					print " at site numbers: ", $seg_site1, " ", $seg_site2, " \n\n";
					}
				($summary{$seg_site1}{$seq_num})++;
				($summary{$seg_site2}{$seq_num})++;
				#print "For segregating site number $seg_site2 and sequence $seq_name{$label_key} ";
				#print "hash of hash holds ", $summary{$seg_site2}{$seq_num};
				#print "\n";
				}
			}
		}	
	}
	return \%summary;
}
#####################################################################################################

# The following subroutine will be called when the $option{s} (-s) is turned on.

# The following subroutine will be called when the $option{s} (-s) is turned on.

sub print_summary
{
my ($max_triplets, $summary_ref, $segsite_ref, $seq_name_ref, $seq_id, $proportion_pattern_a, $middle_hash) = @_;

	
 #print "summarize subroutine\n";
 #print "segsites are @$segsite_ref";
 #print "\n\tsummary of scored sites and sequences:\n";
 print "#\traw priority score\tcorrected priority score\tsite number \t";
 if (!($options -> {a}))
	 {
	 print "sequence name\t";
	 }
 print "sequence number\n";
 print "#\t__________________\t________________________\t_____________\t";
 if (!($options -> {a}))
	 {
	 print "_____________\t";
	 }
 print "_______________\n";
 print "#\n";
foreach $key1(sort hashValueAscending (keys(%$summary_ref))) 
	{
	#print $key1;
	#print "\n";

	for my $key2(sort hashValueAscending (keys(%{$summary_ref->{$key1}})))
		{
		$middle_correction = $middle_hash -> {$key1};
		#print "The mysterious key2 is $key2 end of the line\n";
		my $seq_label;
		$seq_num = $key2+1;
		$seq_label = $seq_name_ref-> {$key2};
		if (length($seq_label) < 8)
		{
			#print "$key2 should start at 0\n\n";
			$seq_label = ($seq_name_ref -> {$key2}) . "       ";
		}
		else
		{ 
			$seq_label = substr($seq_name_ref -> {$key2}, 0, 13);
		}
		$raw_priority_score = $summary_ref -> {$key1} -> {$key2};
		#$correction_factor = 0.1; #The correction factor can be a constant estimate of the expected proportion of pattern a.
		$correction_factor = $proportion_pattern_a;
		$corrected_priority_score = ((int ($raw_priority_score / (($max_triplets + $middle_correction) * $correction_factor)*10))/10);#rounded to the tenths place.
   		print "\t\t$raw_priority_score\t\t\t$corrected_priority_score\t\t\t     $key1\t\t";
		if (!($options -> {a}))
			{
			print "$seq_label\t\t";
			}
		print "$seq_num\n";
   		}
	}
}


####################################################################################################
sub hashValueAscending {
   $hash{$b} <=> $hash{$a};
}

#####################################################################################################

#The following subroutine called each_position tallies the states found at each position.
#
#Next, translate the information at the segregating sites into binary data (subroutine binary)
#send the binary data to the triplets code. 
#the four-gamete test then counts the three site configurations found.
#send the count information back to figure out the rarest state
#print "sequence $seq_id and position $position is $seq_position[$seq_id][$position]\n\n";

sub each_position
{
	#print "Each Postion\n";
	my $informative_seg_site = 0;
	$excluded_site = pop (@_);
	# print "The site is exclude: $excluded_site\n";
	$option_ref = pop (@_);
	#$junk = ($option_ref -> {$i});
	#print $junk , "\n";
	#print "The missing should be counted ", $option_ref -> {n}, "\n";
    my @position_array = @_;
    #print "This array is passed to the function: @position_array\n";
    my %count; 
    my %repeat;
    my $current_base;
    my $how_many_states = 0;
    my $gap = 0;
    my $gap_detected = 0;
    my $singleton = 0;
    my $missing = 0;
    my $illegal_char = 0;
    for my $i ('A','C','G','T') #initialize the values for 4 keys.
    {
	$count{$i} = 0;
	$repeat{$i} = 0;
    }
    
    #j is the sequence number.
    for (my $j = 0; $j < @position_array; $j++)
    {
	#consider the array of bases one position at a time
	#if an element of the array is not A, C, G, or T we need to output a message.
	# unless A, C, G, T : output a message. 
	# Maybe give the chance to include or exclude gaps.
	$current_base = $position_array[$j];
	#print "The current base is $current_base\n";
	if ($count{$current_base} > 0)
	{
	    $repeat{$current_base} = true;
	    #print "$position_array[$j] is a repeat of a previous base at this position.\n";
	}
	else
	{
	    $how_many_states++;
	}
	if (($current_base eq 'A') or ($current_base eq 'C') or ($current_base eq 'G') or ($current_base eq 'T')) #only put a count in the hash if the current base is a legal base call.
	{
	    $count{$current_base}++;
	}
	# Check if it appeared previously in the array	
	# Is it a third state?
	# Is it a singleton?
	for my $i ('A','C','G','T') #initialize the values for 4 keys.
	{
	    if ($count{$i} == 1)
	    {
		$singleton = true;
	    }
	    
	}
    }
    #print "For the current position there are $how_many_states different states.\n"; #This part works.
    if ($how_many_states > 1)
    {
	for (my $j=0; $j < @position_array; $j++)
	{
	    
	    $current_base = $position_array[$j];
	    if (!($current_base eq 'A' or $current_base eq 'C' or $current_base eq 'G' or $current_base eq 'T'))
	    { 
		if ($current_base eq 'N')
		{
		    $missing = true;
		}
		else
		{
		    if ($current_base eq '-')
		    {
			$gap = true;
		    }
		    else
		    {
		    $illegal_char = true;					
			#print "\n\n  **Warning**: An illegal character has been detected in the sequence: $current_base\n";
			#print "  Please check the sequence.\n";
			#print "  Use -n option set to 1 to exclude this site and other missing data from analysis";
			#print "\n  Use -e ? to exclude only this site from the analysis.";
		    }
		}
		#loop back through everything and check for illegal characters.
		#either a monomorphic site or a position with three states.
	    }	
	}
	#Somehow here I need to report whether there is a gap or missing.
	if (((($option_ref -> {g})) and $gap)) 
	#not counting the sites with missing and the site has missing or illegal character data 
	#not counting the gaps and site has a gap.
	{
		#print "\nExcluding a gapped site\n";
		$gap_detected = 1;
	}
	elsif ((($option_ref -> {n})) and (($missing) or ($illegal_char)))
	{
		#print "Excluding sites with missing data\n";
		$missing = 1;
	}
	#if we know the site is excluded we can bypass the each_position sub-routine altogether.
	elsif ($excluded_site)
	{
		$informative_seg_site = 0;
	}
	else 
	{
		$informative_seg_site = (&majority_base(\@position_array,\%count));
	}
	#print "\n";
    }
    else
    {
	$informative_seg_site = 0;   
    }
	#print "return from the eachpositions subroutine:  ",$informative_seg_site, "\n";
    @site_information = ($missing, $gap_detected, $singleton, $illegal_char, $informative_seg_site);
    return \@site_information;
}

#########################################################################################################

sub majority_base #pass in only the count array for the 4 bases
    #return the majority base which at this time is the same as the biggest count.	
{
    #print"In the majority base subroutine.";
    my $majority = 0;
    my $temp = 0;
    my ($one_segsite_array_ref, $base_count_hash_ref) = @_;
    for (my $i=0; $i < @$one_segsite_array_ref; $i++) 
    {
	#print  $one_segsite_array_ref ->[$i];
	#print "\n\n";
    }
    #print "The majority base sub-routine:"; 
    for my $base (keys %$base_count_hash_ref)
    {
	#print "Consider each base: $base is considered.\n";
	my $count = $base_count_hash_ref->{$base};
	if ($count > $temp)
		{
	    $majority = $base;
	    $temp = $count; 
		}
    }
    #print "\nmajority base is :  $majority\n\n";
    return $majority;
}

#########################################################################
#
#  This subroutine is used to find all triplets in the set of 
#  polymorphic sites in binary format.  The binary data is passed as a 
#  reference to an array of arrays  
#  The code determines if each
#  triplet found meets the criteria referred to as "pattern a" in Padhukasahasram 
#  et al. It returns 0 if no triplets are found and returns the appropriate
#  counts for each of the eight states making up each triplet.
# 00, 01, 11, 10 for both incompatible position combinations AB and BC.
#########################################################################

sub triplets
{
    my ($binary_ref, $segsite_ref, $seq_id, $segsite_hash_ref) = @_;
    my @sites;
    my @sites_and_counts;
    my $count =0;
	my $hap;
	my @binary_array = @$binary_ref;	
    for (my $hap_num= 0; $hap_num < $seq_id; $hap_num++)
    {
		$count =0;
		#print "haplotype number $hap_num\n";
		$temp[$hap_num] =();
		#print "Triplets: The haplotype number is $hap_num\n";
		foreach $pos(@$segsite_ref)     
			{
			$count++;
			#print "\nTriplets: position number $count\n";
			#print "temp is:", $temp[$hap_num];
			#print "\n";
			#print "add this";
			#print $binary_array[$hap_num][$count];
			#print "\n";
			$temp[$hap_num] .= ($binary_array[$hap_num][$count]);
			#print "temp now is: ",$temp[$hap_num];
    	    #print "\n";
			}
	$hap = $temp[$hap_num];
	push (@haplotypes, $hap);
	#print "Haplotype to split is: $hap\n";
	if ($options{b})
	{
	    @temp = split(//,$hap);
	    for ($i = 0; $i < @temp;$i++)
	    {	
		if ($i > 0)  #after the first character put the first delimiter.
			{
				if ($options{d} eq "space")
					{
						print " ";
					}
				elsif (!($options{d} eq "none"))
					{
						print "\t"; #tab is the default delimiter for the binary file
					}
			}
	        print $temp[$i];
	    } #end the for loop.
		#print $hap;
	    print "\n";
	} #end the conditional block for a command line request for binary data.
	#print the binary data here.
}
    
    
	#print "#The binary data has been read\n";
    @sites = split(//,$haplotypes[0]);        #split the first haplotype into an array of sites 
    my $N_sites = @sites;                
    my $N_haplotypes = @haplotypes; #the length of the array of haplotypes gives the num of haplotypes.
    my $a_count = 0;  #initialize the variable for number of times 
    #pattern a is found.
    
    my $pos_A = 0;
    my $pos_B = 0;
    my $pos_C = 0;
    my $triplet_num =0;
    
	#Create a four-gamete test lookup table 
    for ($first_site = 0; $first_site < $N_sites - 1; $first_site++)
    {
		for ($second_site = $first_site + 1; $second_site < $N_sites; $second_site++)
		{
			for (my $i= 0; $i < $N_haplotypes; $i++)
				{
				my @sites = split(//,$haplotypes[$i]);
				push(@first, $sites[$first_site]);
				push(@second, $sites[$second_site]);
				}
			#print "The first site is: $first_site array is: @first   \n";
			#print "The second site is: $second_site array is @second\n\n";
			my $temp = &four_gamete($N_haplotypes,\@first,\@second);
			if ($temp -> [0] > 0)
				{
				$lookup[$first_site][$second_site] = [@$temp];
				}
			else
				{
				$lookup[$first_site][$second_site] = 0;
				}
	    
			my $test_array_ref = $lookup[$first_site][$second_site];
			#print  "This is the test result for $first_site $second_site\n";
			#print "\n";
			@first = ();
			@second = ();
		}
    }
    
    #print "number of sites is $N_sites\n";
    for ($pos_A=0; $pos_A < ($N_sites -2); $pos_A++)  #position A must be at least 2 						
	#positions from the end of the 
	#binary data.
	{
		#print "First for loop";
		#print  $pos_A;
		$initial_B = $pos_A + 1;
		#position B must be at least 1 							
	    #character from the end of the 							
	    #binary data.
		
		for ($pos_B= $initial_B; $pos_B < ($N_sites-1); $pos_B++)
		{	
			#print $pos_B;
			$initial_C = $pos_B + 1;
			#position C can be the last position in the file.
			for ($pos_C= $initial_C; $pos_C < $N_sites; $pos_C++)
			{
				#print "Here are the three positions A, B, C:  ";
				#print  $pos_A," ", $pos_B," ",$pos_C;
				$triplet_num++;
				#print "  triplet number: ", $triplet_num, "\n";
				#print "\n\n";
				#print "\nOutside the for loop.\n";    
				#print "\nThis is the array holding this first triplet site\n", @triplet_first, "\n";
				#print "\nThis is the array holding this second triplet site\n", @triplet_second, "\n";
				#print "\nThis is the array holding this third triplet site\n", @triplet_third, "\n";
		
		# for each haplotype I make an array called sites.
		# save the contents of $sites[$pos_A] $sites[$pos_B] $sites[$pos_C]
		# each in its own array the three arrays are
		# called triplet (first,second and third).
		# determine the pairs in the triplet
		#- compare first and third
		#print "four gamete for $pos_A and $pos_C.\n";
				
				if (!$lookup[$pos_A][$pos_C])
					{
						#print "First and third are compatible.\n";
						#print "Compare four gamete 1 and 2 ";
						if ($lookup[$pos_A][$pos_B])
						{
							#print "First and second are not compatible.\n";
							#print "Compare four gamete 2 and 3 ";
							if ($lookup[$pos_B][$pos_C])
							{
								#print "Second and third are not compatible.\n";
								if ($options{v})#verbose option
								{
								print "Detected 'pattern a' for the triplet at the following positions.\n";
								print  $pos_A + 1," ", $pos_B + 1," ",$pos_C+1,"\n";#add one to start at 1 instead of 0 for convenience.
								}
								$a_count++;
			    
			    # every time you find a triplet create the array of sites for each position in the triplet. 
			    # Then find counts for all eight configurations of the four types for each site pair.			    
			    #print "The first site in this triplet is: @first.";
			    #print "Position A is $pos_A\n"; #These are the segregating site positions starting at 0.
			    #print "Position B is $pos_B\n";
			    #print "Position C is $pos_C\n";
								push (@sites_and_counts, $pos_A);
								push (@sites_and_counts, $pos_B);
								#print "There are $N_haplotypes haplotypes.";
			    
								if ($options{v})#verbose option
								{					
									print "\nThe four counts for the first pair:" ,"\n";
				
									print "00 -> ", $lookup[$pos_A][$pos_B]-> [0], "\n";
									print "10 -> ", $lookup[$pos_A][$pos_B] -> [1], "\n";
									print "01 -> ", $lookup[$pos_A][$pos_B] -> [2], "\n";
									print "11 -> ", $lookup[$pos_A][$pos_B] -> [3], "\n";
									print "\n";
								}
								#print "\n";
			    # Find all of the counts that are <= some constant CUTOFF.
			    # Find which individual(s) are resposnsible for the smallest counts.
			    #@temp_1 = @$count_ref_1;
								for (my $i = 0 ; $i < 4; $i++)
								{	
								push (@sites_and_counts, ($lookup[$pos_A][$pos_B]-> [$i]));
								#print "\n",@sites_and_counts;
								}
								push (@sites_and_counts, $pos_B);
								#print "Add site B:\n",@sites_and_counts;
								push (@sites_and_counts, $pos_C);
								#print "Add site C:\n",@sites_and_counts;	
								#Add if (@$count_ref1 != 4) die now.
								#$count_ref_2 = &four_gamete($N_haplotypes, \@second, \@third);
								#print "The positions for the second pair are $pos_B and $pos_C\n"; 
								if ($options{v})#verbose option
								{
									print "The four counts for the second pair:\n";
				
									print "00 -> ", $lookup[$pos_B][$pos_C] -> [0], "\n";
									print "10 -> ", $lookup[$pos_B][$pos_C] -> [1], "\n";
									print "01 -> ", $lookup[$pos_B][$pos_C] -> [2], "\n";
									print "11 -> ", $lookup[$pos_B][$pos_C] -> [3], "\n";
									print "\n";
								}
								#print "\n";
								$reference_to_counts = $lookup[$pos_B][$pos_C];
								for (my $i = 0 ; $i < 4; $i++)
								{
									#my $test = @$count_ref_2;
									#print "the number of elements in the array of counts is $test ";
									push (@sites_and_counts, ($lookup[$pos_B][$pos_C]-> [$i]));
									#print "\n",@sites_and_counts;
								}
#print "pattern a number: $a_count \n The arrays of 4 counts: \n", @temp_1, "\n", @temp_2 , "\n\n";
#These temp arrays are empty. The verbose output includes the four counts for each of the pairs occurring in triplets.
#output currently done in the triplets subroutine.
							}
						}
					}
			}		   
		}
	}
    push (@sites_and_counts, $a_count);
    return \@sites_and_counts;
}

#############################################################################################
# uses two site references to determine whether all four possibilities 
# exist for the pair.

sub four_gamete 
    
{
    #print "Four gamete test\n\n";
    my ($num, $site1_ref, $site2_ref);
    ($num, $site1_ref,$site2_ref) = @_;
	#print @$site1_ref;
	#print "\n";
	#print @$site2_ref;
	#print "\n";
    my $sum_two = 0;
    my $sum_zero = 0;
    my $both_one = 0;
    my $arrangement1_sum1 = 0;
    my $arrangement2_sum1 = 0;
    my $four_gamete = 0;
    my @sum = ();
    for (my $j = 0; $j < $num; $j++)
    {
	$sum[$j] = ($site1_ref->[$j]) + ($site2_ref->[$j]);
        if ($sum[$j] == 2)
	{
	    $sum_two++;
		#print "increment the 11 count.\n $sum_two\n";
	}
        if ($sum[$j] == 1)
	{
	    #print "\nThe sum is one.\n";
	    if ($site1_ref->[$j]) #if it is 1
	    {
		
		$arrangement1_sum1++;
		#print "increment 10 count.\n$arrangement1_sum1\n";
		#print "\narrangement 1 is present\n.";
	    }
	    if ($site2_ref->[$j])
	    { 
		$arrangement2_sum1++;
		#print "increment 01 count.\n$arrangement2_sum1\n";
		#print "\narrangement 2 is present\n.";
	    }
	    if (($arrangement1_sum1) and ($arrangement2_sum1))
	    {
		$both_one = "true";
		#print "\nBoth ones occur: $both_one\n";
	    }
	}
	else
	{
	    if ($sum[$j] == 0)
	    {
		$sum_zero++;
#print "increment 00 count.\n $sum_zero \n";
	    }
	}
	
    }
    
    if (($both_one) and ($sum_two) and ($sum_zero))
    {
	#print "I found all four gametes.";
	#print $both_one
	@counts = ($sum_zero, $arrangement1_sum1, $arrangement2_sum1, $sum_two);
	#print "Check this:\n";
	#print @counts;
	$four_gamete = \@counts; # make this an array reference for the array describing the counts of 00, 10, 01 and 11.
    }                               else 
    {
	$four_gamete = 0;
    }
}

########################################

sub flag_sequences
{
    #print "fed to flag_sequences sub_routine:  @_ \n";
    ($string_match, $pos1, $pos2, $binary_ref, $num, $hash_ref) = @_;
	 my @sequence_array;
    	#print "flag: The string to match is: $string_match It should be the pair with freq = 1\n";
    	#print "flag here: Positions $pos1 $pos2\n";
	my $site1 = $hash_ref -> {$pos1};
	my $site2 = $hash_ref -> {$pos2};
	#print "flag: Positions  ", $hash_ref -> {$pos1}," ", $hash_ref -> {$pos2} ,"\n";
	for (my $j= 0 ; $j < $num ; $j++)
    {
	#print "flag: These are the two bases: ", ($binary_ref -> [$j][$site1]), " and ",$binary_ref -> [$j][$site2], "\n";
	$test = ($binary_ref -> [$j][$site1]) . ($binary_ref -> [$j][$site2]);
	#print "\n flag: Test string: $test\n";
	if ($test eq $string_match)
	{
	   #print "flag: String matching @sequence_array \n";
	   #print "flag: match!\n";
	   #print "flag: $j is the sequence number for $string_match \n";
	    return $j;
	}
		#else
	  	#{
		#	return 0; #It should always find the rare pair.
		#	print "Something unexpected has happened.\n";
	  	#}
    }
}

################################################################################
sub usage
    {
    print "The expected usage of this program:\n";
    print "\n";
    print "usage: $0 -b [-d delimiter] -v -m -s -i input_filename\n\n";
    
	print " -i followed by the filename for the input fasta format - input is required.\n";
	print "    the fasta file should be formatted so that there are no line\n";
	print "    breaks within the sequence data. A line break should appear after the\n";
	print "    sequence name, after each sequence and after the last sequence. This\n";
	print "    version also requires all capital letters for the nucleotides. \n\n";
	print " -a followed by the filename for the ms formatted simulation output\n";
    print "    the simulation input may substitute for the -i fasta file input.\n\n";
	print " -b indicates print binary data\n";
	print " -d followed by the delimiter for the binary data output.\n";
	print "    tab (the default)\n";
	print "    space\n";
	print "    none\n";
	print " -v indicates print verbose output\n";
    print " -m indicates printing a moderate/reasonable amount of output\n";
	print " -s indicates that the output should be only a short summary. \n";
	print "    This is the recommended output option.\n\n";
	
 print "example: $0 -i a_fasta_file.txt -b -d none -s \n\n";

	print "    The following options are new to EDUT_1.1:\n";
	print " -g use 1 to exclude the gapped sites\n";
	print "    use 0 include the gapped sites (the default)\n";
	print " -n use 1 to exclude the sites where sequence is missing indicated by 'N'\n";
	print "    use 0 include the sites where sequence is missing (the default)\n\n";


        exit;
    }

################################################################################

sub collect_input
{
	my $illegal;
	my $seq_id = 0;
	my @info;
	my @array;
	my $prev_length = 0;
	#my $options_ref= @_; 
	#print $options_ref, "\n";
	#$filename = $options_ref -> {i};
	#print $filename, "\n";
	my($options_ref) = @_;
	#foreach (sort keys %$options_ref) {
	#print "$_ => $$options_ref{$_}\n";
	#}
	#print $options_ref -> {i};
	#print "\n";
	if ($options_ref -> {a})
		{
		$info_array = &process_sim_file(\%options);
		return $info_array;
		}
	if ($options_ref -> {i})
	{
	open (FASTA_FILE, $options_ref -> {i}) || usage() ;
	while ($read_line = <FASTA_FILE>)
	{
		$prev_length = $length;
		chomp $read_line;
		$line = uc $read_line;
		if ($line =~ /^>/) #sequence tag 
			{
			#save the tag
				$line =~ s/>//;
				$seq_name{$seq_id} = $line;
				#print "name of sequence $seq_id is $line\n";
				#print "\n";
				@temp = ();
				$seq_id++;
			}
    
		elsif ($line =~ /^#/) #commented line
			{
				$seq_comment{$seq_id} = $line;
				#skip the comment
				#print "comment is $line.\n";
			}
    
		elsif (!($line=~ /^ *$/)) #not a blank line and not a tag.
    						  #check if it is nucleotide sequence
			{
				#print "Not a blank line: $line\n";
				if (!($line =~ /^[ACGTN-]/i) or !($line =~ /[ACGTN-]$/i) )
					{
					#if the line fails to begin or end with a valid nucleotide sequence symbol then print an error message.
					print "\nInput includes characters other than accepted bases:\n\n";
					print "$line\n\n";
					die "Not a properly formatted aligned fasta file. Fix input and try again.\n\n";
					}
				else 
				#save the two-dimensional array for later processing.
					{
					if ( $line =~/^[ACGT-]*N[ACGT-]*$/i)
						{
						$missing = true;
						}
					if ( $line =~/^[ACGTN]*-[ACGTN]*$/i)
						{
						$gapped = true;
						}
					@bases_in_sample = split(//,$line);
					$array[$seq_id]=[@bases_in_sample];
					$current_length = @bases_in_sample;	
					}	
					#print "The current length is $current_length\n";
				$length = $current_length;
			}
		if ($illegal)
		{
			die 
		}
		if (($seq_id > 1) and ($length != $prev_length) and ($length != 0)) #check if you are looking at the first sequence
			 #if not, check if the sequence length is the same as the length of the previously considered sequence.
			{
				print "\n\nThe length of the sequences before sequence $seq_id was $prev_length \n";
				print "Sequence $seq_id is of length $current_length \n";
				print "\nThe sequences are not all of the same length. \nSequence data may not be properly aligned.\n";
				die "Please fix the alignment and try again.\n\n";    		
		 	}
		elsif (($length == 0) and ($seq_id > 1) and ($prev_length))
			{
				die "Problem reading one or more sequences. Check file format.\n";
			}
	}
	}
  #end the while loop
   #end each line of sequence    	
push (@info,\@array); #a reference to the above array.
push (@info, $length); #length of sequence
push (@info, $seq_id); #number of sequences
push (@info, $missing);
push (@info, $gapped);
push (@info, \%seq_name);
return (\@info);
}

####################################################################################################
sub preliminary_binary
{
    ($majority_base, $current_position, $site_info, $position) = @_;
    ($missing, $gap, $singleton, $illegal_char, $informative_seg_site) = @$site_info;
    #print "Position $position is represented by the array: @$current_position\n";
    $site_position = $position + 1;
	#print "\nAdd the site position $site_position to the list of segsites.\n";#always adding one to the list of segsites.
	#keep an array containing all of the positions of the segregating sites.
    # This line was previously being called before the each_position subroutine was 
    # called and it didn't work properly.
    
    for ($i=0; $i< @$current_position; $i++)
    {
		#print "Consider this: ", @position_array[$i]; #An array holding two references. Odd.
		#print " $i\n";
		#print "Inside this loop position $position is represented by the array: @$current_position\n";
		@position_array = @$current_position;
		if (($majority_base))
		{ 
	  	  if (@position_array[$i] eq $majority_base) 
	   		{	
			#print "majority found\n";
			$binary_column[$i]= 0;
	    	}	
	    	else	
	    	{
				if ((@position_array[$i] eq 'A') or (@position_array[$i] eq 'C') or (@position_array[$i] eq 'G') or (@position_array[$i] eq 'T') or (@position_array[$i] eq 'a') or (@position_array[$i] eq 'c') or (@position_array[$i] eq 'g') or (@position_array[$i] eq 't'))
				{
		    	#print "another legal base.\n";	
		    	$binary_column[$i]= 1;
				}
				else
				{
		    		#print "something other than a base.\n";
		    		$binary_column[$i] = 'N';
		    		#do the checking here for missing bases or illegal characters.
				}
	    	}		   
	}
	#@binary_sequence[$i]= @binary_column; #was this accidentally commented out?
	#print " $i:  ";
	#print "@binary_sequence[$i]\n";
    } 
    #print "\nThe first element returned by the preliminary binary subroutine is: @binary_column\n"; #put a reference to the binary array into an array
	#print "The second thing returned is the site_position $site_position for this binary data.\n";
    #print "\n";
	push (@info,\@binary_column);#Each time the subroutine returns a reference it seems to replace the referenced array returned the last time.
	push (@info,$site_position);
	return \@info;
}
##############################################################################################

sub process_columns #I believe this subroutine is currently not being used. See comment for intended usage.
{
    ($num_sequences, $num_segsites, $binary_ref, $seg_sites_ref, $hash_ref ) = @_ ;
    my $counter = 0;
    my $counter2 = 0;
	$singleton_count= 0;
	$only_gapped_count=0;
    #print "This is the array for the binary data: @$binary_ref";
    #print "\n";
    #print "\nThe number of segsites is $num_segsites\n";
    
	
# Look at the binary data column by column and determine if the column is informative.
# This means that the code will check that it has a polymorphism other than the missing data.
# Also make sure it does not have a singleton mutation
# Use a constant. Instead of excluding only singletons allow exclusion of user defined low frequency
# mutations.     
    $count_informative_seg_site=0;
	$count_non_informative=0;
	foreach $column_num(@$seg_sites_ref)

    {
    	$site_num = $hash_ref -> {$column_num};
		#print "The segsites are @$seg_sites_ref\n";
		#print "$column_num is the position number \n";
		my $singleton = 0;
		my $only_gapped = 0;
		my @count;
		$count[0] = 0;
  		$count[1] = 0;
  		$only_gapped = 0;
  		$singleton = 0;
  		#print "In the process columns subroutine. $num_sequences is the number of sequences to loop\n";
		my @temp =[];
		for (my $seq_num= 0; $seq_num < $num_sequences; $seq_num++)    
		{
		 #print "The sequence number in the process columns loop: $seq_num\n";  
		#print "$counter Use this number: ", ($binary_ref -> [$counter]), "\n";
	    my $temp_variable = ($binary_ref -> [$seq_num][$site_num]);
	    #print "The temporary variable from the binary data is $temp_variable\n";
	    $counter++;
	    if (!($temp_variable eq 'N'))
	    	{
	    	#print "The temporary variable from the binary data is $temp_variable\n";
			$count[$temp_variable]++;
	    	}
		#print "the counts for column $column_num are $count[0] for 0 and $count[1] for the 1 state.\n"; 
		}
		
		if (($count[0]==1 ) or ($count[1] == 1))
			{
	    	$singleton = true;
			$singleton_count++;
			}
		if ($count[1] == 0)
			{
		    $only_gapped = true;
			$only_gapped_count++;
			#print "Only a gap site. Do not count!\n";
			}
		if ($count[0] == 0)
			{
		    #$only_gapped = true;
		    #print "This message should not appear.";
			}
		if (($singleton) and ($only_gapped))
		{
			$singleton = 0;
			$singleton_count--;
			#A site can not be a singleton if it is not segregating for 2 bases.
		}
		#print "Is this column a singleton: $singleton \n and is this site only gapped with no other polymorphism?:  $only_gapped\n";
		if (!($singleton) and !($only_gapped))
			{
		    #for (my $seq_num= 0; $seq_num < $num_sequences; $seq_num++) 
			#    {   
			#	$temp .= ($binary_ref -> [$counter2]);
			#	$counter2++;
				#print "\n";
				#print $temp;
				#print "\n";
			#    }
		    #push(@revised_binary_array, $temp);
		    #$temp = '';
			$count_informative_seg_site++; #This counts all of the sites that are informative
			$informative_segsites_hash{$column_num} = $count_informative_seg_site;
											#The number is used to calibrate the priority scores.
			}
		else
			{
		    #for (my $seq_num= 0; $seq_num < $num_sequences; $seq_num++) 
			#    {
			# 	$counter2++;
			#    }	  
			$count_non_informative++;  
			}
		#print "\n";
		#return \@revised_binary_array;
    }
	push (@informative_sites, $count_informative_seg_site);
	push (@informative_sites, $count_non_informative);
	push (@informative_sites, $singleton_count);
	push (@informative_sites, $only_gapped_count);
	push (@informative_sites, \%informative_segsites_hash);
	return (\@informative_sites);
}

sub process_sim_file

{
	my ($options_ref) = @_;
	my @bases_in_sample;
	my @array;
	my $current_length;
	my $prev_length;
	my $length;
	open (SIM_FILE, ($options_ref -> {a})) || usage();
	while ($line = <SIM_FILE>)
		{
		chomp $line;
		if (!($line =~ /^\/\/$/) and !($line =~ /^segsites/) and !($line =~ /^positions/) and !(/^\s+$/ ))
			{
				#print "Processing the simulation data";
				$prev_length = $current_length;
				@bases_in_sample = split(//,$line);
				$array[$seq_id]=[@bases_in_sample];
				$current_length = @bases_in_sample;
				if ($current_length)
				{
					$seq_id++; #all part of the messy fix for my simulation processing
				}	
				#print "Length of simulated sequence is $current_length.\n"; 
				
			}
			if ($current_length != 0)
			{
			$length = $current_length;
			}
			else
			{
			$length = $prev_length;
			#print "length = $length (problem is here)\n\n";#messy way to fix the problem caused by the trailing newline 
			}
			$seq_name{$seq_id} = "seq_". $seq_id; #not used. not working
		}
		push (@info,\@array); #a reference to the above array. It is actually already binary if from a sim.
		#print "Length to push $length.\n";
		push (@info, $length); #length of sequence
		push (@info, $seq_id); #number of sequences
		push (@info, 0); #represents no missing data (never in sims)
		push (@info, 0); #represents no gapped data  (never in sims)
		push (@info, \%seq_name);
		return (\@info);
}


