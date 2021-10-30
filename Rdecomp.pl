# Rdecomp.pl derived from clean_for_hybridtrees.pl
#
# Dec. 3, 2020
# DG This program calls clean_for_hybridtrees recursively to do the full cluster decomposition of the input matrix. The output is in file "decomponly" where each
# leaf of the recursion tree has the string written below it, along with the KWARG upper bound (for the cluster problem, not the hybridization problem)
# for that matrix. Then to find the hybridization number of the original input matrix, we find the hybridization number
# of each of these subproblems, and sum their hybridization numbers. The output file ``decomp" also has these matrices, but has much more information
# that can be used to check correctness.
#
# Currently, Dec. 11, 2020, this program is taylored to examine the grass data. Each grass datafile has the name "grassx.ymatrix", where
# x is the number of trees in the input, and y is the largest leaf index. In the grass data, the leaves
# are numbered starting from 0, so y is one less than the number of leaves (in each tree). The program is gerry-rigged so that the file "outx.y" is written to and then read, and then
# written to. This needs to be cleaned up someday, but it works for now.
#
# Call Rdecomp.pl as: perl Rdecomp.pl grassx.ymatrix > outx.y
#
# The stand-alone program clean_for_hybridtrees.pl has been incorporated as a subroutine, and information about the stand-alone program is given below.
#
#
# DG Nov. 20, modified clean_for_hybridtrees so that duplicate columns inside the same tree are removed. Such duplicate cols. will never be part of the input, since
# each column corresponds to a split, and each split in a tree is recorded exactly once, but after the other part of cleanup (removing some rows and cols),
# duplicate cols. can be created. Note that, unlike the case for building a DAG for representing clusters (in a cluster-network), we must not
# remove duplicate cols. that come from different trees, because that split must be part of each tree specification.
#
#
# DG This program should clean a binary matrix as used in computing the hybridization number in the version when a given set of
# trees must be represented in the resulting DAG. 
#
# The input is the output of parse_Newick.pl where the FIRST TWO lines produced
# by parse_Newick.pl (the borders line and the tree_string line) are NOT removed. The output has a line with
# only the new borders between the columns representing the input trees, followed by the cleaned matrix.
# 
# July 29, 2020 Currently the program assumes that the number of trees is at most nine.
#
# the input matrix should have no spaces between the bits in the binary matrix.
#
# call on a terminal command line as: perl clean_for_hybridtrees.pl input-file-name
#
# derived from cleanup.pl
# July 28, 2020
#

# July 30, 2020 These programs do not use or alter border information. They just pass borders, as cborders or tborders or newcborders back and forth
# to make sure that passing is working correctly, in case we do want to use borders later. If we never do need it, we will delete it then.


$infile = $ARGV[0];
if ($infile =~ /(\d+\.\d+)/) {
	$extension = $1;
	$out = 'out' . $extension;
	print "The accumlated out file is $out \n";
}

open (IN, $infile);  # Here we read in the 0/1 input matrix, which will be used later.  No spaces between bits.
$outfile = $infile . 'clean';
system ("rm $outfile");

open (CLEAN, ">>$outfile"); # the file where the clean matrices are written - because of the recursive call structure, the order of these reflect the depth-first structure
                            # of the recursive calls.
open (DEBUG, '>>debug');
open (DECOMP, '>>decomp');
open (DECOMPONLY, '>>decomponly');

$borders = <IN>; # the first line in the input file has the border information for the trees represented in the input matrix
chomp $borders;
@tborders = $borders =~ /(\d+)/g;

$tree_string = <IN>; # the second line in the input file has the tree labels, 1 through the number of trees, written in each column.
chomp $tree_string;

#
print DECOMP "The full input from $infile is: \n";
print DECOMP "$borders\n";
print DECOMP "$tree_string\n";

print "The full input from $infile is: \n";
print "$borders\n";
print "$tree_string\n";

@inarray = <IN>;

foreach $line (@inarray) {
	chomp $line;
	print "$line\n";
	print DECOMP "$line\n";
}
print "\n";
$bits = @inarray;   # This is the number of lines in inarray.

clean_for_hybridization ($borders, $tree_string, @inarray);

 close (CLEAN);
 close (DEBUG);
 close (DECOMP);
 close (DECOMPONLY);
 close (IN);

open (OUTSCAN, $out);
@inscan = <OUTSCAN>;
$tswcount = 0;
$total_TSWtime = 0;
$total_unsat_value = 0;
$total_indeterm_value = 0;

print "\nExecution Summary\n";
foreach $line (@inscan) {
	if ($line =~ /took (\d+\.\d+) wall-time/) {
		$TSWtime = $1;
		$tswcount++;
		print "tsw execution $tswcount took: $TSWtime \n";
		$total_TSWtime += $TSWtime;
	}
	if ($line =~ /value (\d+) is UNSATISFIABLE/) {
		$unsat_value = $1 + 1;
		$total_unsat_value += $unsat_value;
		print "unsat values: $unsat_value, $total_unsat_value\n";
	}

	if ($line =~ /value (\d+) is INDETERMINATE/) {
		$indeterm_value = $1 + 1;
		$total_indeterm_value += $indeterm_value;
		print "indeterm_values: $indeterm_value, $total_indeterm_value\n";
	}
}

$total_value = $total_unsat_value + $total_indeterm_value;
print "\nThe obtained value (which might not be optimal if there are any indeterm values) is $total_value\n";
print "The total time for all the TSW computations is $total_TSWtime\n";

#
#
#
sub clean_for_hybridization {

my ($borders, $tree_string, @inarray) = @_;

#print CLEAN "Entering clean_for_hybridization with the passed matrix with borders: $borders \n";
print DEBUG "Entering clean_for_hybridization with the passed matrix with borders: $borders \n";

print DEBUG "\nThe data passed to clean_for_hybridization is: \n";
print DEBUG "$borders";
print DEBUG "$tree_string\n";
foreach $blat (@inarray) {
	chomp $blat;
	print DEBUG "$blat\n";
}

$allones = "";  # the all-ones string records the original set of rows that the program is working with. Initially,
                # that set consists of all the rows, but some 1 entries change to 0 as rows are removed.  #DG why is this set at every entry? Maybe it should
		#                                                                                         be set outside and passed to clean_for_hybridization.
foreach $i (1 ... $bits) {
       $allones .=  '1';
  }

    @cleanarray = ();
    ($newbinary, $newborders, $newtree_string, @cleanedarray) = clarray ($allones, $borders, $tree_string, @inarray);  # pass both the binary vector detailing the indices of the subset of remaining rows, 
                                                                 # and a string giving the boundries of the trees in the input, and then the row subset.
                                                                 # newbinary is the binary vector detailing the rows of the cleaned array.
                                                                 # newboundries is a string giving the new tree boundry values.

     print "returned from clarray. New binary is $newbinary \n";								 
     print "returned from clarray. New borders is $newborders \n";
     print "returned from clarray. New treestring is $newtree_string \n";

   # print "The cleaned array after return from clarray \n"; 

    #    foreach $line (@cleanedarray) {
         #print OUT "$line\n";
	 #}

$colsnow = length ($cleanedarray[0]);
$ntslen = length ($newtree_string);
#print "After clarray, but before removing duplicate cols, the number of cols is $colsnow, and the new treestring is $newtree_string of length $ntslen  \n";

   #print OUT "\n The binary vector showing which original rows are contained in the cleaned array - although some columns have also been removed \n";
   #print OUT "$newbinary\n";
   #print OUT "\n\n";

   unshift (@cleanedarray, $newtree_string); # DG Nov. 22, 2020 unshift does a pre-pend to a list, unlike push which does an append to a list
   $atree = @cleanedarray;
   #print "Having added the string $newtree_string, the array now has $atree rows \n";
   #   ($num_rows, $num_cols, @final_array) =  remove_dup_columns(@cleanedarray); # DG Jan 3, 2021 It seems that this was left in
                                                                              # by mistake in previous executions. When the input
									      # consists of explicit trees (as in TSW type computations), 
									      # there will be no duplications
									      # inside any single tree, and any duplications that there
									      # are across two or more trees must be left in so that all
									      # of the splits in the different trees are captured.
									      # This is different from SW type computations, where the
									      # input just consists of splits, without any explicit
									      # trees.
   @final_array = @cleanedarray; # DG Jan 3, 2021 Because I have commented out the above call to remove_dup_columns, the following three
                                 # variables need to be set explicitly here. Why not just avoid using them? For continuity with
				 # the existing code, not wanting to have to check if any problems would be created by their removal.
				 # Some day, this should be cleaned-up, but for now, it is the simplest and safest way to respond to
				 # the removed call to remove_dup_columns.
   $num_rows = @final_array;
   $num_cols = $colsnow;

   #print "The matrix after remove_dups has $num_rows rows and $num_cols columns: \n";

   $newtree_string = $final_array[0];

   chomp $newtree_string;
   #print "Back from remove_dups, newtree_string:\n$newtree_string\n";
   #print "Now we remove the added first row, newtree_string \n";
   shift (@final_array);
   $num_rows = @final_array;
   $num_cols = length ($newtree_string);
   #print "The final number of rows and cols is $num_rows, $num_cols \n";

   if ($num_rows < 3) {
	   return;
   }

   $num_trees = 1;
   @tree_string = split (//, $newtree_string);
   $newborders = "1 ";
   foreach $i (1 .. $num_cols - 1) {
       if ($tree_string[$i] ne $tree_string[$i-1]) {
          $newborders .= "$i ";
          $num_trees++;
       }
    }
    $newborders .= "$num_cols";
            
    $n = @cleanedarray;  # the number of rows in the cleaned array
    print "The number of rows and columns in the cleaned matrix:  $n, $num_cols\n";
    print DEBUG "The number of rows and columns in the cleaned matrix:  $n, $num_cols\n";
    print CLEAN "$newborders\n";
    print DEBUG"$newborders\n";
    print "$newborders\n";
    print "$newtree_string\n";
    foreach $line (@final_array) {
         print CLEAN "$line\n";
         print DEBUG "$line\n";
	 print "$line\n";
   }

   #print "The output is in file: $ARGV[0]clean\n";


 maxcluster($num_trees, $borders, $newborders, $newtree_string, @final_array);

 } #end of clean_for_hybridization



############
sub clarray {  # cleanarray calls cleaniter repeatedly until no changes are made in an iteration

%rowhash = ();
my $row = 0;

@cinarray = @_; 
# print OUT "cinarray before shift: \n @cinarray \n";

$cbinary= shift @cinarray;
$cborders = shift @cinarray;
$ctree_string = shift @cinarray;

#print OUT "In clarray, the string cbinary is $cbinary\n";  # cbinary is a binary string of length = bits, detailing the subset of rows in cinarray.
#print OUT "In clarray, the string cborders is $cborders \n";
#print OUT "In clarray, the string ctree_string is $ctree_string\n";  # ctree_string with a single integer for each column, detailing in which tree the column originates.

#print OUT "cinarray after shifts: \n@cinarray\n";

foreach $line (@cinarray) {
       chomp $line;
#       #print OUT "In clarray: $line \n";
       $rowhash{$row} = $line;
       $row++;
}

my $line_count = $row;   # a true count of the number of rows. But the largest index for the rows is $row-1.
my $line_length = length $cinarray[0];
##print OUT "Here is where line_length is assigned: $line_length \n";

my $old_line_count = 100000; # start with a large number bigger than the actual number of rows.
my $old_line_length = 1000000;


while (($line_count < $old_line_count) or ($line_length < $old_line_length))  {  # while we continue to make changes in the array, either removing a row
                                                                                 # or removing a column.
    $old_line_count  = $line_count;  # update the count and length
    $old_line_length = $line_length;

    @ccleanarray = ();
    #print OUT "About to call cleaniter with binary strings $cbinary and $cborders and $ctree_string \n";

    ($newcbinary, $newcborders, $newctree_string, @ccleanarray) = cleaniter ($cbinary,$cborders,$ctree_string,@cinarray);

    #print OUT "Returned from cleaniter, cbinary is now newcbinary, with value $newcbinary\n";
    #print OUT "Returned from cleaniter, cborders is now newcborders, with value $newcborders \n";
    #print OUT "Returned from cleaniter, ctree_string is now newctree_string, with value $newctree_string \n";

    if (! $ccleanarray[0]) {  # if ccleanarray is empty, set the count and length to 0
        $line_count = $line_length = 0;
    }

    else {
         $line_count = @ccleanarray; # update line_count with the new number of rows
         $line_length = length($ccleanarray[0]);  # update the new line length
    }

    $cbinary = $newcbinary;
    $cborders = $newcborders;
    $ctree_string = $newctree_string;
    @cinarray = @ccleanarray;
    
    $cinarray_len = @cinarray;
     #print OUT "cbinary; old_line_count, line_count, old_line_length, line_length: $cbinary; $old_line_count, $line_count, $old_line_length, $line_length \n";

     #print OUT "at the end of the while statement in clarray - returning to the top to compare old values with new values \n";
}  # end of while statement


#    #print OUT "\n The cleaned array has $line_count lines, each with $line_length bits\n";
#    foreach $line (@ccleanarray) {
#        #print OUT "$line\n";  
#    }

return ($cbinary, $cborders, $ctree_string, @cinarray);
}
################

sub cleaniter {
# remove duplicate rows and uninformative columns. This assumes that the root is the all zero sequence.
##print OUT "entering cleaniter \n";

my ($cbinary, $cborders, $ctree_string, @inlines) = @_;
##print OUT "Inside cleaniter. cbinary is $cbinary \n";
@cbinarybits = split (//,$cbinary);
%insequence = ();
@lines = ();
@outlines = ();
@cstring = split (//,$ctree_string);
$tstring = "";

$j = 0;
foreach $line (@inlines) {  # this block will remove duplicate rows in the input sequences
   chomp $line;
   while ($cbinarybits[$j] eq 0) { # find the next 1 in cbinarybits. The way that j and inlines are coordinated, there always should be
                                   # a next 1 at this point.
#        print "$j \n";
        $j++;
   }
      
#   #print OUT "$j, $cbinarybits[$j]: $line \n";
   if (!defined $insequence{$line}) {
      $insequence{$line} = 1;
      push (@lines, $line);
   }
   else {
        $cbinarybits[$j] = 0; # if the line is a duplicate, don't add it to @lines, and change bit j to 0 in cbinarybits
#        #print OUT "cbinarybits when line $j is a duplicate @cbinarybits \n";
   }
   $j++;
}

$cbinary = join ('', @cbinarybits);
# #print OUT "The new cbinary in cleaniter, after removing duplicate rows, is $cbinary \n";
#print OUT "In cleaniter, after removing duplicate rows. The lineas array is: \n";
#foreach  $line (@lines) {
   #print OUT "$line\n";
   #}

# Having removed duplicate rows, we now remove uninformative columns. But it is easier to work with rows, so we first transpose the matrix.

@tlines = transpose(@lines);  # transpose the matrix @lines
@translines = ();
$num_goodlines = 0;
$colnum = 0;
foreach $line (@tlines) {   # examine each row (which was originally a col.) to see if it has
                            # more than one entry of value 1, and at least one 0. If so, then add that row to the growing @translines.

      @bits = split (//,$line);
      $line_length = @bits;
      $onecount = 0;
      foreach $bit (@bits) {
         if ($bit eq '1') {
            $onecount++;
         }
      #print OUT "In cleaniter. colnum, bit and onecount: $colnum, $bit, $onecount \n";
      }

       if ($onecount > 1) {
#      if (($onecount > 1) and (!defined $columnstransposed{"$line"})) lbr 
#      if (($onecount > 1) and ($onecount < $line_length)) lbr # don't accumulate a line (actually a column in the inlines array) if it only has a single 1 or
                                                            # it has all ones. These cases are not informative. Note the asymmetry because we are
                                                            # assuming that the root sequence is the all-zero sequence. So, we will accumulate a line
                                                            # that has only a single 0.
         $goodline = join ('',@bits);
         $num_goodlines++;
         push (@translines, $goodline);

          $tstring .= $cstring[$colnum];
          } 
$colnum++;
}

##print OUT "The number of goodlines is $num_goodlines \n";
##print OUT "In cleaniter translines:  \n";
#foreach $lIne (@translines) {
   # print  OUT "$lIne \n";
#}

if ($num_goodlines == 0) {   # create and return cbinary as the all-zero string.
  $len_cbinary = length $cbinary;
  $cbinary = "";
  foreach $i (0 ... $len_cbinary - 1) {
       $cbinary .=  '0';
  }
}
else {
    @outlines = transpose(@translines);
}

#print OUT "About to exit cleaniter. The returned cbinary is $cbinary \n";
##print OUT "The resulting outlines matrix is \n ";
#foreach $outlin (@outlines) {
#      #print OUT "$outlin \n";
#}

$num_lines_in_outlines = @outlines;
#print OUT "The number of lines in outlines is $num_lines_in_outlines \n";

return ($cbinary, $cborders, $tstring,  @outlines);
}


#####################
sub transpose {

@lines = @_;

@transpose = ();
$i = 0;
foreach $line (@lines) {
  chomp $line;
  @row = split(//,$line);
  foreach $j (0 ... length($line) - 1) {
     $transpose[$j][$i] = $row[$j];
#     print "$i, $j, $row[$j], $transpose[$j][$i]\n";
  }
# print "\n";
$i = $i + 1;
$linelength = length($line);
}

# print "The number of rows in the transpose is $linelength, and the number of cols. is $i \n";

 @trans = ();
 foreach $p (0 ... $linelength - 1) {
    $newrow = "";
    foreach $q (0 ... $i-1) {
#      print "$transpose[$p][$q]";
      $newrow = $newrow . $transpose[$p][$q];
    }
#    print "\n";
    push (@trans, $newrow);
 }

# foreach $line (sort @trans) {
#      print CRES "$line\n";
# }

return @trans;
}


################
sub remove_dup_columns { # This should not be used in Rdecomp.pl or GRdecomp.pl

# ($cborders, @rows) = @_;
 @rows = @_;
 #print "entering remove_dup_columns, the matrix passed is: \n";
 #foreach $liner (@rows) {
	 #	print "$liner\n";
	 #}
#print "\n";

@out_rows = ();
@out_cols = ();
%unique_cols = ();

@cols = transpose(@rows);

$col_count = 0;
foreach $col (@cols) {
  $col_count++;
  if (!defined $unique_cols{$col}) {
     push (@out_cols,  $col);
     $unique_cols{$col} = 1;
  }
  else {
  #    print "In remove_dup_columns, just removed col $col_count \n";
  }
}        

$num_cols = @out_cols;
@out_rows = transpose(@out_cols);

$num_rows = @out_rows;
#print "Inside remove_dup_columns. The number of rows and columns after removal of duplicate cols. is $num_rows, $num_cols \n";

return ($num_rows, $num_cols, @out_rows);

}

##########

sub maxcluster {

my ($num_trees, $borders, $newborders, $newtree_string, @final_array) = @_;
#print "The array passed to maxcluster: \n";
#foreach $liner (@final_array) {
	#	 print "$liner \n";
	# }
#print "\n";

my @one_count = (); # counts of the numbers of ones in the columns


my $num_rows = @final_array;
my $num_cols = length($final_array[0]);
print "\nIn maxcluster with $newborders. The number of rows and cols in the input matrix is: $num_rows, $num_cols \n";
print CLEAN "In maxcluster with $newborders. The number of rows and cols in the input matrix is: $num_rows, $num_cols \n";
print DEBUG "In maxcluster with $newborders. The number of rows and cols in the input matrix is: $num_rows, $num_cols \n";

my @tfarray = transpose(@final_array);

#print "The transposed array in maxcluster: \n";
#foreach $liner (@tfarray) {
#	 print "$liner\n";
#}
#print "\n";

my $tnum_rows = $num_cols;
my $tnum_cols = $num_rows;
my %cvector = (); # hash whose keys are the distinct column vectors in final_array, and whose values are the number of times that vector appears in final_array.
my @cvectorcount = (); # array which lists for column k, the number of prior occurrences (in columns with index < k) of the vector in column k.

#print "In maxcluster, the number of rows and cols in the transposed array is: $tnum_rows, $tnum_cols \n";

foreach $j (0, $tnum_rows) {
    $one_count[$j] = 0;
}

foreach $trow (0 .. $tnum_rows - 1) {  # count the number of 1s in each row of tfarray, which is the number of 1s in the corresponding column of final_array.
    @line = split(//, $tfarray[$trow]);
    foreach $j (0 .. $tnum_cols - 1) {
       if ($line[$j] == 1) {
          $one_count[$trow]++;
       }
     }

my    $string = $tfarray[$trow];
#     print "In maxcluster, trow $trow in tfarray has string $string \n";
     if (defined $cvector{$string}) {
	     #  print "in maxcluster, there is repeated vector of tfarray at trow $trow \n";
               $cvector{$string} = $cvector{$string} +  1;
	       #print "$cvectorcount[$trow] \n";  
               $cvectorcount[$trow] = $cvector{$string};  
	       #print "Now cvectorcount for row $trow of tfarray has value $cvectorcount[$trow] \n";
     }

     else {
               $cvector{$string} = 1;  # the vector in row trow of tfarray is seen for the first time.
               $cvectorcount[$trow] = 1;
     }
}

# Now find a column in final_array (a row in tfarray) that is common to all input trees, with a one_count that is strictly less than tnum_cols, but with the maximum
# count over all the columns that are common to all input trees and have a one_count that is strictly less than tnum_cols.

my $selected_one_count = 0;
foreach $k (0 .. $tnum_rows-1) {   # print out the counts of 1s in the columns of final_array (rows in tfarray), along with that column and the number of prior times the vector
                                 # in that column has been seen.
				 #print "$k: $one_count[$k] $tfarray[$k] $cvector{$tfarray[$k]} $cvectorcount[$k] $selected_one_count $selectedk \n";

     if (($one_count[$k] < $tnum_cols) && ($cvectorcount[$k] == $num_trees) && ($one_count[$k] >= $selected_one_count)) {
             $selected_one_count = $one_count[$k];
             $selectedk = $k;
	     #print "There is a change of selectedk when k = $k \n"; # k is the column in final_array (row in tfarray) that defines the cluster we will extract.
     }
} # end of foreach k

if ($selected_one_count == 0) { # the dfs search for the next cluster hit a leaf without finding a next cluster
	print "There is no subcluster below the roots, so the only cluster in the tree is the entire set of leaves \n";
	print DECOMP "\nThere is no subcluster below the roots, so the only cluster in the tree is the entire set of leaves \n";
#	print CLEAN "\nThere is no subcluster below the roots, so the only cluster in the tree is the entire set of leaves \n";
	print DEBUG "\nThere is no subcluster below the roots, so the only cluster in the tree is the entire set of leaves \n";
	open (TSWGRASS, '>for_tswgrass');
	print DECOMPONLY "$newborders\n";
	print TSWGRASS "$newborders\n";
	print TSWGRASS "$newtree_string\n";

	open (FORKWARG, '>forkwarg');
	foreach $line (@final_array) { # write out a leaf-matrix found in the decomposition, to DECOMPONLY, and also to FORKWARG.
		chomp $line;
        	print DECOMPONLY "$line\n";
		print FORKWARG "$line\n";
		print TSWGRASS "$line\n";
	}

        close (FORKWARG);
        close (TSWGRASS);
	open (KWARGTEMP, '>kwargtemp');
        system ("./kwarg -k forkwarg > kwargtemp");
        close (KWARGTEMP);

	open (KWARGTEMP, 'kwargtemp');
	$line = <KWARGTEMP>;
	if ($line =~ /(\d+)/) {
		$kwargbound = $1;
	}
	close(KWARGTEMP);

        print DECOMPONLY "kwargbound: $kwargbound\n\n";

        system ("perl TSWgrass.pl for_tswgrass $kwargbound 5000"); # now find the hybridization number of the subproblem in file for_tswgrass
	return;
}

#print "The selected one count and selected k: $selected_one_count, $selectedk \n";

####
# Now we extract out the subproblem defined by the clusters identified by the cols. of M that are identical to the column selectedk.
# Note that this part of the code operates on final_array rather than tfarray, and so there is no need to transpose.

print DEBUG "Extracting clusters in maxcluster, the value of num_cols is $num_cols \n";

my @new_array = ();
my @remaining_array = ();
my @new_leaf = (); # this is the row in the remainder subproblem that defines the new leaf that the solution to the cluster subproblem will attach to.

foreach $i (0 .. $num_rows - 1) {
   my   $fstring = $final_array[$i];
   my   @row = split (//, $fstring);


   if (length ($fstring) != $num_cols) {
	   print "YIKES - mismatch in lengths in maxcluster \n"; #DG Dec. 4, 2020
   }

     $k = 0;
     if ($row[$selectedk] == 1) {  # check if this row identifies a leaf that is part of the cluster subproblem 
    my    @new_row = ();

          foreach $j (0 .. $num_cols - 1) { # scan the columns in this identified row
		  # print DEBUG "In the loop in extracting the newleaf: $j, $num_cols \n";

             if ($one_count[$j] < $one_count[$selectedk]) { # check if the 1s in column j are contained in the ones in column selectedk, so 
		                                            # the mutation represented by column j is in the subtree below the mutation represented by column selectedk
							    # Or there are no 1s in column j contained in the ones in column selectedk, so the mutation represented by
							    # column j has no relation to the mutation represented by column selectedk (actually in the copy of selectedk
							    # that is in the same input tree as column j
                $new_row[$j] = $row[$j];
		$new_leaf[$j] = 0;
             }
             if ($one_count[$j] == $one_count[$selectedk]) { # check if column j is identical to column selectedk (but in k-1 cases j will be in a different input tree)
                $new_row[$j] = $row[$j];
		$new_leaf[$j] = $row[$j]; # This is just a slick way to say "if (row[j] == 1 then new_leaf[j] = 1) else new_leaf[j] = 0. DG Dec. 7, 2020
             }
             if ($one_count[$j] > $one_count[$selectedk]) { # check if column j represents a mutation whose 1s contain all of the 1s in column selectedk, in which case
		                                            # the mutation represented by column j is on the path from the root to the subtree, in the input tree containing
							    # column j, whose root is the copy of selectedk that is in the same input tree as column j.
                $new_row[$j] = 0;
		$new_leaf[$j] = $row[$j]; # This is just a slick way to say "if (row[j] == 1 then new_leaf[j] = 1) else new_leaf[j] = 0. DG Dec. 7, 2020

		#if ($row[$j] == 1) {
			#	$new_leaf[$j] = 1;
			#}
	      }
          } #end of foreach j

      my  $newfstring = join ('',@new_row);
          chomp $newfstring; 
          $newlen = length ($newfstring);
	  #print "The length of string newfstring is $newlen\n";
          push (@new_array, $newfstring);
     } # end of if (row[selectk] == 1)

     else {
	  push (@remaining_array, $fstring);
      } 
} # end of foreach i


print "\nThe extracted subproblem from $borders is represented by the following matrix of $num_cols columns:\n";
print CLEAN "The extracted subproblem from $borders is represented by the following matrix of $num_cols columns:\n";
print DEBUG "The extracted subproblem from $borders is represented by the following matrix of $num_cols columns:\n";
print "\n$newborders\n";
print DECOMP "\n$newborders\n";
print DEBUG "\n$newborders\n";
print "$newtree_string\n";
print DECOMP "$newtree_string\n";
print DEBUG "$newtree_string\n";
foreach $string (@new_array) {
   chomp $string;
   print "$string\n";
   print DECOMP "$string\n";
   print DEBUG "$string\n";
}


print "\nThe remaining subproblem from $borders is represented by the following matrix: \n";
print "$newborders\n";
print DECOMP "\n$newborders\n";
print DEBUG "\n$newborders\n";
print "$newtree_string\n";
print DECOMP "$newtree_string\n";
print DEBUG "$newtree_string\n";
foreach $string (@remaining_array) {
   chomp $string;
   print "$string\n";
   print DECOMP "$string\n";
   print DEBUG "$string\n";
}

# Now create the needed new row providing the leaf that the
# solution to the cluster-subproblem attaches.
$anotherlenleaf = @new_leaf;
$nlstring = join (//,@new_leaf);
$lenleaf = length ($nlstring);

#print "The newleaf string has length $lenleaf $anotherlenleaf: \n$nlstring\n";
print CLEAN "The newleaf string has length $lenleaf $anotherlenleaf: \n$nlstring\n";
print DEBUG "The newleaf string has length $lenleaf $anotherlenleaf: \n$nlstring\n";
print DECOMP "The newleaf string has length $lenleaf $anotherlenleaf: \n$nlstring\n";


push (@remaining_array, $nlstring);

clean_for_hybridization ($newborders, $newtree_string, @new_array);
clean_for_hybridization ($newborders, $newtree_string, @remaining_array);

} # end of sub maxcluster



