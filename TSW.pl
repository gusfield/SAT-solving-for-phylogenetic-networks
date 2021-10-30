# TSW.pl
# Adapted from TSWhybpipe1.pl
# DG
# October 31, 2020
# Pipeline for exploring the SAT approach to computing the softwired hybridization number of a binary matrix,
# given a set of trees as input, in Newick format. The DAG must contain each of the given trees, not just their splits as
# in SWhybpipe1.pl and SWhybpipe3.pl. The root sequence is assumed to be the all-zero sequence.
#
# The input must contain a first line which gives the borders between trees, followed by a line with tree_string that explicitly
# labels the trees from 1 ... number of trees (at most 9)
#
# Call as: perl TSW.pl datafile at the fnamemstmatrix-stage target max_time
#
# When cryptominisat5 is used as the SAT-solver, max_time must be CPU time, which (on a four-core maching) could be four time
# the intended wall-time. However, this is only approximate, because the true core use will only be determined at run-time. It could
# be anything between one (unlikely) and four (more likely on large problem sizes). Often it is around three.
#
# With a target of t, the number of interior
# nodes can be bounded by n-2 + 2t, although that is generally larger than necessary. So, one can choose a smaller interior_num if willing
# to risk having the CNF be UNSAT, and then having to restart with a larger target.
#
  $dm = $ARGV[0];
  #$tswtrace = 'tswtrace.' . $dm;
  open (TRACE, '>>tswtrace');
  print "\nStarting TSW with input matrix below\n";
  print TRACE "\nStarting TSW with input matrix below\n";
  use Time::HiRes qw(gettimeofday);
  $before = gettimeofday;
  # print "The before time is $before \n";

  $total_sat_time = 0;
  $target = $ARGV[1]; # starting target
  print "The target passed to TSW.pl is $target \n";
  $mtime = $ARGV[2];
  # print "The input max_time value is $mtime \n";
  #print "The file name (passed to TSW from SWhybrid.pl, or at the direct call) is $ARGV[0] \n";

  system ("perl clean_no_cluster.pl $dm"); # DG Dec. 10, 2020 call the version of clean_for_hybridization that does not do the find-cluster computation,
                                         # which has already been done when TSWgrass.pl is called from Rdecomp or GRdecomp. Also, when TSW.pl is
					 # called from SWhybpipe.pl, we don't want to do any subcluster finding.
                                         #
   #system ("perl clean_for_hybridtrees.pl $dm");
$cdata  = $dm . 'clean';
#print "The file after clean_for_hybridtrees is: $cdata \n";
$gtimediff = 0;

# $interior_num = int (($n - 2 + 2 * $target) * 0.65); # now compute interior_num in Thybridnum.pl


# print "\n Just before call to search, The input max_time value is $mtime \n";
  $sat = search($cdata, $target, $mtime); # DG update Nov. 7, 2020
  #  print "Back from search. Variables sat, timeout, target equal to: $sat, $timeout, $target \n";

$time = localtime();
#print "Computation ending at time $time \n";
#print "Finito \n";
$after = gettimeofday;
#print "The after time is $after \n";
$elapsed_time = $after - $before;
print "The total TSW computation for this decision took $elapsed_time wall-time seconds \n";
print TRACE "The total TSW computation for this decision took $elapsed_time wall-time seconds \n";
print "\nThe total time used in SAT instances is $total_sat_time \n\n";
print TRACE "\nThe total time used in SAT instances is $total_sat_time \n\n";
#
close(TRACE); # DG Jan 22, 2021

##########

sub search {

         my ($cdata, $target, $mtime) = @_;
	 # print "\nIn search, cdata, target, max_time  = $cdata, $target, $mtime \n";

# The structure here is a bit odd since the code for sub search was repurposed from code that didn't use a sub.
# I could have streamlined a bit, but chose not to - if it ain't broke, don't fix it.

$timeout = 0;
$sat = 1;

 while ($sat)  { 

  print "\n Trying TSW Softwired (for tree network) target of $target, with max allowed CPU-time = $mtime seconds, \n";
  $walltime = $mtime/4;
  print "which could be as small as $walltime seconds, depending on the number of cores used on my macbook Pro (up to 4).\n";

  print TRACE "\n Trying TSW Softwired (for tree network) target of $target, with max allowed CPU-time = $mtime seconds \n";
  print TRACE "When cryptominisat5 is used as the SAT-solver, the max allowed time is CPU time, so the actual wall-time \
  will be between a fourth to a half of that time \n";
  system ("perl Thybridnum.pl $cdata $target");
  $time = localtime();
  $cd = $cdata . '.SAT';
  #  print "file $cd created.\n"; 

  #print "Starting SAT solver at time $time \n";
  $timeout = 0;

  $gtime0 =  gettimeofday;
  print TRACE "Starting SAT solver at epoch time: $gtime0 \n";

  $SIG{'ALRM'} = handler;
    alarm($max_time);
  $cdsol = $cdata . '.sol';
  #    system ("time ./cadical $cd > $cdsol");
  #    system ("./glucose-syrup -model $cdata.SAT > $cdsol");
  #system ("./plingeling $cd > $cdsol");
  system ("time cryptominisat5 -t 4 --maxtime $mtime $cd > $cdsol");
    alarm(0);

  $time1 = localtime();
  $gtime1 = gettimeofday;
  #print TRACE "SAT solvers for Softwired case finished at time $time1 \n";
  #print "\n Inside TSW: SAT solvers for Softwired case finished at time $time1 \n";

  $gtimediff = $gtime1 - $gtime0;

  print TRACE "SAT solvers for Softwired case finished at epoch time $gtime1 \n";
  print "SAT solvers for Softwired case finished at epoch time $gtime1 \n";

  print TRACE "The elapsed SAT-solver time for the Softwired case is $gtimediff \n";
  print "The elapsed SAT-solver time for the Softwired case is $gtimediff \n";

  $totalsearchtime += $gtimediff;


  print TRACE "The total time used by the solver is $totalsearchtime \n";
  print "\n Inside TSW: The total time used by the solver is $totalsearchtime \n";

  open (IN, "$cdsol");
  @lines = <IN>;
  close (IN);

  foreach $line (@lines) { 

    if ($line =~ /(\d+.\d+) seconds, /) {
        print "Execution time: $1 seconds \n";
    }

    if ($line =~ /real time : (\d+.\d+)/) {
        print "Execution time: $1 seconds \n";
    }

    if ($line =~ /UNKNOWN/) {
       print "The SAT-solver did not terminate in the allowed time \n"; # plingeling and cadical print UNKNOWN in this case.
       print TRACE "UNKNOWN: The SAT-solver did not terminate in the allowed time \n"; # plingeling and cadical print UNKNOWN in this case.
       $sat = 0;
       $timeout = 1;
#       return (0, 1, $target1, $interior_num1);
       return (0);
#       last;
    }

    if ($line =~ /INTERRUPTED/) {
       print "The SAT-solver did not terminate in the allowed time \n"; # Glucose prints INTERRUPTED in this case.
       print TRACE "INTERRUPTED: The SAT-solver did not terminate in the allowed time \n"; # Glucose prints INTERRUPTED in this case.
       $sat = 0;
       $timeout = 1;
       #return (0, 1, $target1, $interior_num1);
       return (0);
#       last;
    }

    if ($line =~ /INDETERMINATE/) {
       print "target value $target is INDETERMINATE with cryptominisat5. Time for this TIMEOUT execution \n";
       print "case is $gtimediff seconds \n";

       print TRACE "target value $target is INDETERMINATE with cryptominisat5. Time for this TIMEOUT execution \n";
       print TRACE "case is $gtimediff seconds \n";
       #       print "TIMEOUT in cryptominisat5 solver \n";
              print "cryptominisat5 did not terminate in the allowed time \n"; # cryptominisat5 prints INDETERMINATE in this case.
       $sat = 0;
       $timeout = 1;
       return (0);
#       last;
    }

    if ($line =~ /UNSATISFIABLE/) {
       #print "SAT solver finished \n";
       print "target value $target is UNSATISFIABLE \n";
       print TRACE "target value $target is UNSATISFIABLE \n";
       $value_obtained = $target + 1;
       $total_UNSATvalue += $value_obtained;

       $sat = 0;
       return (0);
#       last;
    }
    if (($line =~ / SATISFIABLE/) || ($line =~ /^SATISFIABLE/)) {
       #print "SAT solver finished \n";
       print "target value $target is SATISFIABLE \n";
       print TRACE "target value $target is SATISFIABLE \n";
       $sat = 1;
       $total_sat_time = $totalsearchtime;
       last;
     }
  }

  if (($sat == 1) and ($timeout == 0)) { # this confirms that the last execution ended with SATISFIABLE
	  #print "Before call to SAT_SOLconversion, the value of cdata is $cdata \n";

     system ("perl SAT_SOLconversion.pl $cdata $cdata");
     system ("perl tcheck-trees.pl $cdata.TRANSOL $cdata");
     #print "\n";
  
  
  $target = $target - 1;
  $interior_num = $interior_num - 1;
  }

    } # end while

    #  if ($timeout == 1) {
	  #       $SIG{'ALRM'} = handler;
	  #  alarm($max_time);
#    system("python3 pipeline.py $cdata -b $target -n $interior_num");
#      $cadfile = $cdata . '.SAT';
      #      system("./cadical --unsat -t 24 $cadfile");
      #      alarm(0);
     #  }


} # end sub search

########

sub handler {
  $sat = 0;
  $timeout = 1;
  print "TIMEOUT from handler in TSW\n";
  kill 2, -$$;
  return(0);
};

########

sub cleanup { # Note, in TSW.pl cleanup is not called, and cleanup is the only place where remove_dup_columns is called.
	      # so in TSW.pl (and also in TSWgrass.pl, remove_dup_columns is not called. So why have it here? Dunno.

@inarray = @_;

# DG This should clean a binary matrix as used in computing the history bound for an n by m matrix. 
# the input matrix should have no spaces between the bits.
#
#

open (OUT, ">$ARGV[0]clean"); # the file where the clean matrix is written

#@inarray = <IN>;

#print OUT "The input array is: \n";
#print OUT "@inarray \n";
#$bits = @inarray;   # This is the number of lines in inarray.

$allones = "";
foreach $i (0 ... $bits - 1) {
       $allones .=  '1';
  }

    @cleanarray = ();
    ($newresult, @cleanedarray) = clarray ($allones, @inarray);  # pass both the binary vector detailing the indexes of the subset, and then the subset itself.
                                                               # newresult is the binary vector detailing the rows of the cleaned array.
#    print OUT "The cleaned array after return from clarray \n"; 

#    foreach $line (@cleanedarray) {
#         print OUT "$line\n";
#    }

   print TRACE "\n The binary vector showing which original rows are contained in the cleaned array - although some columns have also been removed \n";
   print TRACE "$newresult \n";
   
   print TRACE "\n\n";
   ($num_rows, $num_cols, @final_array) =  remove_dup_columns(@cleanedarray);
   print "The final matrix after cleanup has $num_rows rows and $num_cols columns: \n";
   print TRACE "The final matrix after cleanup has $num_rows rows and $num_cols columns: \n";
   foreach $line (@final_array) {
         print TRACE "$line\n";
         print OUT "$line\n";
   }

   # print "Output is in file: $ARGV[0]clean\n";
 close (OUT);
 close (TRACE);

return ($num_rows, $num_cols, @final_array);
} #end of cleanup.pl

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
#print OUT "XX In clarray, the string ctree_string is $ctree_string\n";  # ctree_string with a single integer for each column, detailing in which tree the column originates.

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
foreach  $line (@lines) {
   #print OUT "$line\n";
}

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
#         print "XXX @bits \n";
         $goodline = join ('',@bits);
#         #print OUT "goodline: $goodline \n";
         $num_goodlines++;
         push (@translines, $goodline);

#         $columnstransposed{"$goodline"} = 1;
          $tstring .= $cstring[$colnum];
          #print OUT "colnum tstring $colnum $tstring\n";
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

# print "XXXXXX\n";
# foreach $line (sort @trans) {
#      print CRES "$line\n";
# }

return @trans;
}


################
sub remove_dup_columns {

($cborders, @rows) = @_;
# @rows = @_;
@out_rows = ();
@out_cols = ();
%unique_cols = ();

@cols = transpose(@rows);

foreach $col (@cols) {
  if (!defined $unique_cols{$col}) {
     push (@out_cols,  $col);
     $unique_cols{$col} = 1;
  }
  else {
    #print OUT "In remove_dups, just removed a cols \n";
  }
}        

$num_cols = @out_cols;
@out_rows = transpose(@out_cols);

$num_rows = @out_rows;
# print "The number of rows and columns after removal of duplicate cols. is $num_rows, $num_cols \n";

return ($num_rows, $num_cols, $cborders, @out_rows);

}
