# SWgrass.pl
# Nov. 10, 2020 addapted for use with tree data that is already in Newick form, so msTclean.pl is not needed, and
# there shouldn't be ``mst" in any file name.
#
# call as: SWgrass.pl file-name-of-Newickdata maxtime
#
# SWhybpipe.pl  # this is the correct one
# Nov. 9, 2020
#
# DG October 31, 2020
#
# Sept. 13 modify to add msTclean.pl at the beginning. msTclean.pl takes in the output from MS -T
# and extracts up to ten distinct trees in Newick format.
#
# DG July 16 This version is in summer19/glucose-syrup-4.1/parallel. Another version is in summer19/hybridization-number-sat-main in order to be able to 
# call python3 pipeline.py when looking for UNSAT.
# For some reason I can't get pipeling.py to run here 
# If major changes are made to either one, change the other as well.


# SWhybpipe1.pl
# DG July 3, 2020
# Pipeline for exploring the SAT approach to computing the softwired hybridization number of a binary matrix. 
# With a target of t, the number of interior
# nodes can be bounded by n-2 + 2t, although that is generally larger than necessary. So, one can choose a smaller interior_num if willing
# to risk having the CNF be UNSAT, and then having to restart with a larger one.

$fname = $ARGV[0]; # $ARGV[0] is the input file of trees in Newick format.
$maxtime = $ARGV[1];

use Time::HiRes qw(gettimeofday);

$before = gettimeofday;

#
system ("rm trace");
open (TRACE, '>trace');

# system ("date +%s"); # number of seconds since the start of the epoch

$fname = $ARGV[0];

  parse_Newick ($fname); # the output goes to fnamematrix 
  $fname = $fname . 'matrix';

    print "The input wall-maxtime value is $maxtime seconds\n";

open(IN, "$fname");
open (SHRUB, '>tempshrub');
open (KWARG, '>tempkwarg');

$line1 = <IN>; # input matrices from parse_Newick.pl have two lines at the start for the case that the input comes from trees
               # comment out these two read IN lines in the case that the input does not come from parse_Newick. Also, Bhybridnum.pl does not
               # use these two lines.
$line2 = <IN>;

@data = <IN>;
chomp @data;
close(IN);

($n, $m, @data) = pipeclean(@data);
print TRACE "After pipeclean, the cleaned array is of shape $n, $m \n";
print "After pipeclean, the cleaned array is of shape $n, $m \n";

$fname = $fname . 'clean';
@data_shrub = ();
@data_kwarg = ();
push (@data_shrub, "$fname");
$first_line = '0' x $m;
push (@data_shrub, $first_line);
push (@data_kwarg, $first_line);

foreach  $line (@data) {
     push (@data_shrub, $line);
     push (@data_kwarg, $line);
}

foreach $line (@data_shrub) {
	   print SHRUB "$line \n";
	}
close(SHRUB);

foreach $line (@data_kwarg) {
    print KWARG "$line \n";
}
close(KWARG);


print "The SHRUB fast upper bound and full upper bound are ";
system ("./shrub-mac -r tempshrub > outshrub");
system ("./kwarg -k tempkwarg > outkwarg");
close(SHRUB);
close(KWARG);

open (OUTSHRUB, 'outshrub');

@shrubfile = <OUTSHRUB>;
  foreach $line (@shrubfile) {
	       if ($line =~ /full upper bound is (\d+)/) {
		   $targetSHRUB = $1;
		   print "The best upper bound from SHRUB is $targetSHRUB \n";
		}
	}
close(OUTSHRUB);

open (OUTKWARG, 'outkwarg');
   
@kwargfile = <OUTKWARG>;
  foreach $line (@kwargfile) {
    if ($line =~ /Found history with (\d+)/) {
       $targetKWARG = $1;
       print "The best upper bound from KWARG is $targetKWARG \n";
    }
  }
close(OUTKWARG);

if ($targetKWARG < $targetSHRUB) {
	 $target = $targetKWARG;
	  }
  else {
	 	 $target = $targetSHRUB;
	 }
   
 print "The initial target is $target \n";

$interior_num = int (($n - 2 + 2 * $target) * 0.75);

  $cdata = $fname; 
  #  print "XXXX Assigning fname to cdata: $fname, $cdata \n";
#  ($sat, $timeout, $target, $interior_num) = search($target, $interior_num, $maxtime);
  #print "\n Just before call to Search, The input maxtime value is $maxtime \n";
  $sat = search($cdata);
  #print "Back from Search. Variables sat, timeout, target and interior_num equal to: $sat, $timeout, $target, $interior_num \n";

  if ($timeout == 1) {
  $after = gettimeofday;
  $elapsed = $after - $before;
  print "Search timed out. The time since the start of the program is $elapsed seconds. \n";
  #  print "Now we call cadical with --unsat option.  \n";
#    system("python3 pipeline.py $cdata -b $target -n $interior_num");
#    $cadfile = $cdata . '.SAT';
#    $unsatmaxtime = 4 * $maxtime;
#     system("./cadical --unsat -t $unsatmaxtime $cadfile");
    #        system("cryptominisat5 -t 4 --maxtime $unsatmaxtime $cadfile");
  }

  #  $date1 = date +%s;
  # print "date1 is $date1 \n";
$time1 = localtime();
print TRACE "Computation ending at time $time1 \n";
# print "Finito \n";

$after = gettimeofday;
$elapsed = $after - $before;
#print "Before was $before. After was $after. The total time for the SW problem is $elapsed seconds \n";
print TRACE "Before was $before. After was $after. The total wall-time for the SW problem is $elapsed seconds \n";
#
#
# print "\n Now we recompute using the same input, but solve the TSW variant of the problem \n";
$TSWtarget =  2 * $target;
$TSWmaxtime = 24 * $maxtime; # The true wall-time TSWmaxtime is intended to be six times the SW maxtime, but cryptominisat interprets
                             # maxtimes as CPU times, which are typically from three to four times as large as wall-times.

print TRACE "The TSW target and maxtime are $TSWtarget, $TSWmaxtime \n";
#print "The TSW target and maxtime are $TSWtarget, $TSWmaxtime \n";

$TSWfile = $ARGV[0] . 'matrix';
close(TRACE);
system ("perl TSWgrass.pl $TSWfile $TSWtarget $TSWmaxtime > outT");

open(TRACE, '>>trace');
$after = gettimeofday;
$elapsed = $after - $before;
print "The total wall-time for the SW+TSW problem instances is $elapsed seconds \n";
print TRACE "Before was $before. After was $after. The total wall-time for the SW+TSW problem instances is $elapsed seconds \n";
system ("cat out outT > outTT");
close (TRACE);
#
##########

sub search {

$cdata = $_[0];
#print "\nIn Search: file = $cdata; target and interior_num are $target, $interior_num \n";

# The structure here is a bit odd since the code for sub search was repurposed from code that didn't use a sub.
# I could have streamlined a bit, but chose not to - if it ain't broke, don't fix it.

$timeout = 0;
$sat = 1;
$totalsolvetime = 0;
while ($sat)  { 

  print TRACE "\n TRYING Softwired (cluster network) target of $target, with $interior_num interior nodes, and max allowed wall-time = $maxtime seconds \n";
  print "\n Trying Softwired (cluster network) target of $target, with $interior_num interior nodes, and max allowed wall-time = $maxtime seconds \n";
  system ("perl Bhybridnum.pl $cdata $m $target $interior_num");
  #  print TRACE "file $cdata.SAT created.\n"; 

  $time0 = localtime();
  $gtime0 =  gettimeofday;
  print TRACE "Starting SAT solver at time $time0 \n";
  #print "Starting SAT solver at time $time0 \n";
  $timeout = 0;

    $cdatag = $cdata . 'g';
    $cdatap = $cdata . 'p';
    $cdataCAD = $cdata . 'CAD';
    $cdataCRYP = $cdata . 'CRYP';


  #  print "\n SOLVING with Glucose-Syrup \n";
  #  system ("time ./glucose-syrup -model $cdata.SAT > $cdatag.sol");
    #    alarm(0);

    # $SIG{'ALRM'} = handler;
  # alarm($maxtime);

  #  print "\n SOLVING with Plingeling\n";
  #  system ("time ./plingeling $cdata.SAT > $cdatap.sol");
    #    alarm(0);


    # $SIG{'ALRM'} = handler;
  # alarm($maxtime);
  #print "\n SOLVING with Cadical \n";
  #system ("time ./cadical $cdata.SAT > $cdataCAD.sol");
  #  alarm(0);

  $SIG{'ALRM'} = handler;
  alarm($maxtime);
  # print "\n SOLVING with cryptominisat5\n";
  print TRACE "\n SOLVING with cryptominisat5\n";
  system ("time cryptominisat5 -t 4 $cdata.SAT > $cdata.sol");
  alarm(0);

  $time1 = localtime();
  $gtime1 = gettimeofday;
  print TRACE "SAT solvers for Softwired case finished at time $time1 \n";
  #print "SAT solvers for Softwired case finished at time $time1 \n";

  $gtimediff = $gtime1 - $gtime0;

  print TRACE "The elapsed time for the Softwired case is $gtimediff \n";
  print "The elapsed time for the Softwired case is $gtimediff \n";

  $totalsearchtime += $gtimediff;

  print TRACE "The total time used by the solver is $totalsearchtime \n";
  print "Inside SWhybpipe: The total time used by the solver is $totalsearchtime \n";

  open (IN, "$cdata.sol");
  @lines = <IN>;
  close (IN);

  foreach $line (@lines) { 

    if ($line =~ /(\d+.\d+) seconds, /) {
        print TRACE "Execution time: $1 seconds \n";
        print "Execution time: $1 seconds \n";
    }

    if ($line =~ /real time : (\d+.\d+)/) {
        print TRACE "Execution time: $1 seconds \n";
        print "Execution time: $1 seconds \n";
    }

    if ($line =~ /UNKNOWN/) {
       print "The SAT-solver did not terminate in the allowed time \n"; # plingeling and cadical print UNKNOWN in this case.
       $sat = 0;
       $timeout = 1;
#       return (0, 1, $target1, $interior_num1);
       return (0);
#       break;
    }

    if ($line =~ /INTERRUPTED/) {
       print "The SAT-solver did not terminate in the allowed time \n"; # Glucose prints INTERRUPTED in this case.
       $sat = 0;
       $timeout = 1;
       #return (0, 1, $target1, $interior_num1);
       return (0);
#       break;
    }

   if ($line =~ /INDETERMINATE/) {
       print "TIMEOUT in cryptominisat5 solver \n";
       print "The SAT-solver did not terminate in the allowed time \n"; # cryptominisat5 prints INDETERMINATE in this case.
       $sat = 0;
       $timeout = 1;
       return (0);
#       break;
    }

    if ($line =~ /UNSATISFIABLE/) {
       #print "SAT solver finished \n";
       print "target value $target is UNSATISFIABLE \n";
       $sat = 0;
       #return (0, 0, $target1, $interior_num1);
       return (0);
#       break;
    }
    if ($line =~ / SATISFIABLE/) {
       #print "SAT solver finished \n";
       print "target value $target is SATISFIABLE \n";
       $sat = 1;
       break;
     }
  }

  if (($sat == 1) and ($timeout == 0)) { # this confirms that the last execution ended with SATISFIABLE
       #system ("perl SAT_SOLconversion.pl $cdata $cdata");
     SAT_SOLconversion ($cdata);
     #     print TRACE "Inside Search. cdata is $cdata \n";
     check_trees ($cdata);
     #     print "\n";

  
  $target = $target - 1;
  $interior_num = $interior_num - 1;
  }

   } # end while
} # end sub search

###########

sub handler {
  $sat = 0;
  $timeout = 1;
  print "TIMEOUT \n";
  kill 2, -$$;
  return;
};

########

sub pipeclean {  # This is basically the program cleanup.pl modified so that it only outputs the cleaned matrix, while the program cleanup.pl
	         # produces other information before the cleaned matrix. So when hand-walking the programs using cleanup.pl instead of pipeclean
		 # one has to extracted the cleaned matrix from the output of cleanup.pl before inputting the matrix to Bhybridnum.pl

@inarray = @_;
$bits = @inarray; # This is the number of entries in inarray.

# DG This should clean a binary matrix as used in computing the history bound for an n by m matrix. 
# the input matrix should have no spaces between the bits.
#
#

my $fname = $fname . 'clean';
open (OUT, ">$fname"); # the file where the clean matrix is written

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
   #print TRACE "The final matrix after pipeclean has $num_rows rows and $num_cols columns: \n";
   #print "The final matrix after pipeclean has $num_rows rows and $num_cols columns: \n";
   foreach $line (@final_array) {
         print TRACE "$line\n";
         print OUT "$line\n";
	 print "$line\n"; # DG Nov 9
   }

 print TRACE "Output from Pipeclean is in file: $fname\n";
 close (OUT);

return ($num_rows, $num_cols, @final_array);
}

############
sub clarray {  # cleanarray calls cleaniter repeatedly until no changes are made in an iteration

%rowhash = ();
my $row = 0;

@cinarray = @_; 
# print "cinarray before shift: \n @cinarray \n";

$cbinary= shift @cinarray;
#print OUT "In clarray, the string cbinary is $cbinary\n";  # cbinary is a binary string of length = bits, detailing the subset of rows in cinarray.
#print "In clarray, the string cbinary is $cbinary\n";  # cbinary is a binary string of length = bits, detailing the subset of rows in cinarray.
#print "cinarray after shift: @cinarray \n";

foreach $line (@cinarray) {
       chomp $line;
#       print OUT "In clarray: $line \n";
       $rowhash{$row} = $line;
       $row++;
}

my $line_count = $row;   # a true count of the number of rows. But the largest index for the rows is $row-1.
my $line_length = length $cinarray[0];
#print OUT "Here is where line_length is assigned: $line_length \n";

my $old_line_count = 100000; # start with a large number bigger than the actual number of rows.
my $old_line_length = 1000000;


while (($line_count < $old_line_count) or ($line_length < $old_line_length))  {  # while we continue to make changes in the array, either removing a row
                                                                                 # or removing a column.
    $old_line_count  = $line_count;  # update the count and length
    $old_line_length = $line_length;

    @ccleanarray = ();
#    print OUT "About to call cleaniter with binary string $cbinary\n";

    ($newcbinary, @ccleanarray) = cleaniter ($cbinary,@cinarray);

    if (! $ccleanarray[0]) {  # if ccleanarray is empty, set the count and length to 0
        $line_count = $line_length = 0;
    }

    else {
         $line_count = @ccleanarray; # update line_count with the new number of rows
         $line_length = length($ccleanarray[0]);  # update the new line length
    }

    $cbinary = $newcbinary;
    @cinarray = @ccleanarray;
    $cinarray_len = @cinarray;
#     print OUT "cbinary; old_line_count, line_count, old_line_length, line_length: $cbinary; $old_line_count, $line_count, $old_line_length, $line_length \n";

#     print "at the end of the while statement - returning to the top \n";
#     print OUT "at the end of the while statement in clarray - returning to the top to compare old values with new values \n";
}  # end of while statement


#    print OUT "\n The cleaned array has $line_count lines, each with $line_length bits\n";
#    foreach $line (@ccleanarray) {
#        print OUT "$line\n";  
#    }

return ($cbinary, @cinarray);
}
################

sub cleaniter {
# remove duplicate rows and uninformative columns. This assumes that the root is the all zero sequence.
#print OUT "entering cleaniter \n";

my ($cbinary, @inlines) = @_;
#print OUT "Inside cleaniter. cbinary is $cbinary \n";
@cbinarybits = split (//,$cbinary);
%insequence = ();
@lines = ();
@outlines = ();
# %columnstranposed = ();

$j = 0;
foreach $line (@inlines) {  # this block will remove duplicate rows in the input sequences
   chomp $line;
   while ($cbinarybits[$j] eq 0) { # find the next 1 in cbinarybits. The way that j and inlines are coordinated, there always should be
                                   # a next 1 at this point.
#        print "$j \n";
        $j++;
   }
      
#   print OUT "$j, $cbinarybits[$j]: $line \n";
   if (!defined $insequence{$line}) {
      $insequence{$line} = 1;
      push (@lines, $line);
   }
   else {
        $cbinarybits[$j] = 0; # if the line is a duplicate, don't add it to @lines, and change bit j to 0 in cbinarybits
#        print OUT "cbinarybits when line $j is a duplicate @cbinarybits \n";
   }
   $j++;
}

$cbinary = join ('', @cbinarybits);
# print OUT "The new cbinary in cleaniter, after removing duplicate rows, is $cbinary \n";

# Having removed duplicate rows, we now remove uninformative columns. But it is easier to work with rows, so we first transpose the matrix.

@tlines = transpose(@lines);  # transpose the matrix @lines
@translines = ();
$num_goodlines = 0;

foreach $line (@tlines) {   # examine each row (which was originally a col.) to see if it has
                            # more than one entry of value 1, and at least one 0. If so, then add that row to the growing @translines.

      @bits = split (//,$line);
      $line_length = @bits;
      $onecount = 0;
      foreach $bit (@bits) {
         if ($bit == 1) {
            $onecount++;
         }
#      print OUT "bit and onecount: $bit $onecount \n";
      }

      if ($onecount > 1) {
#      if (($onecount > 1) and (!defined $columnstransposed{"$line"})) 
#      if (($onecount > 1) and ($onecount < $line_length))  # don't accumulate a line (actually a column in the inlines array) if it only has a single 1 or
                                                            # it has all ones. These cases are not informative. Note the asymmetry because we are
                                                            # assuming that the root sequence is the all-zero sequence. So, we will accumulate a line
                                                            # that has only a single 0.
#         print "XXX @bits \n";
         $goodline = join ('',@bits);
#         print OUT "goodline: $goodline \n";
         $num_goodlines++;
         push (@translines, $goodline);

#         $columnstransposed{"$goodline"} = 1;
      } 
#print "\n";
} # end of foreach line

#print OUT "The number of goodlines is $num_goodlines \n";
#print OUT "In cleaniter translines:  \n";
foreach $lIne (@translines) {
   # print  OUT "$lIne \n";
}

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
#print OUT "The resulting outlines matrix is \n ";
foreach $outlin (@outlines) {
#      print OUT "$outlin \n";
}

$num_lines_in_outlines = @outlines;
#print OUT "The number of lines in outlines is $num_lines_in_outlines \n";

return ($cbinary, @outlines);
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

@rows = @_;
@out_rows = ();
@out_cols = ();
%unique_cols = ();

@cols = transpose(@rows);

foreach $col (@cols) {
  if (!defined $unique_cols{$col}) {
     push (@out_cols,  $col);
     $unique_cols{$col} = 1;
  }
}        

$num_cols = @out_cols;
@out_rows = transpose(@out_cols);

$num_rows = @out_rows;
# print "The number of rows and columns after removal of duplicate cols. is $num_rows, $num_cols \n";

return ($num_rows, $num_cols, @out_rows);

}


#########
#
sub SAT_SOLconversion {

	# open (TRACE, '>>trace');
# translate the Dimacs formatted Solution into a solution using named variables
my $cdata = $_[0];
#print TRACE "In SAT_SOLconversion, cdata is $cdata \n";

open (KEY, "$cdata.KEY"); # this is the file that shows the correspondence between variable names and the integers given in the dimacs format. It is a .KEY file. Enter only the file name without the .KEY extension

open (INSOL, "$cdata.SOL"); # enter a file name of a .SOL file to process the CNF formulas. Enter only the file name without the  .SAT extension

# Typically, I run plingeling as ``plingeling test.SAT > test.SOL", where ``test.SAT" is the name of the file containing the
# CNF formula. The output goes to test.SOL. The lines that actually contain the solution all begin with a "v", so the program
# below can identify and extract them from test.SOL.
# For Glucose-syrup, one must use the option -model to get a list of variables set true.

#print "Starting SAT_SOLconversion \n";
open (OUT, ">$cdata.TRANSOL");  #DG June 28 changed 0 to 1
#print "The output should go into file $cdata[0].TRANSOL \n";
#print OUT "First is a list of CNF integers and their corresponding variable names \n";
#print OUT "The values of the variables in the SAT solution follow \n";

%chr = ();

@keylines = ();

@keylines = <KEY>;
$numkeylines = @keylines;
%conversion = ();

$cols = 1;  # this initialized the value of cols, which will be updated below 
$rows = 1;

foreach $j (0  .. $numkeylines - 1)  {
   if ($keylines[$j] =~ /(\w.+) (\d+)/) {
#      print "$j \n";
      $variable = $1;
      $key = $2;
#      $key = $1;
#      $variable = $2;
#      print "Accumulating the conversion hash with key and variable: $key $variable \n";
      $conversion{"$key"} = "$variable";
      $conversion{"-$key"} = "-$variable";
      if ($variable =~ /(\d+),(\d+)/) {
           $first = $1;
           $second = $2;
           if ($first > $rows) {
             $rows = $first;
            }
           if ($second > $cols) {
             $cols = $second;
            }
       }
   }
}

# print "Finished Accumulating Keys \n";

#foreach $key (sort {$a<=>$b} keys %conversion) {
#    print OUT "$key $conversion{$key} \n";
#}
#print OUT "End of variables list \n\n";
#print OUT "Now starts the list of the variable values in the solution \n";

@satlines = <INSOL>; # read in the .SAT or .SOL file containing the CNF clauses or the truth values assigned to the variables in a SAT solution.
$numlines = @satlines;
# print "The number of lines in the .SAT or .SOL file = $numlines \n";

foreach $i (0 .. $numlines -1) {

   if ($satlines[$i] =~ /^v [ ,-0123456789]+$/) {    # use this to convert the literals set true in a satisfiable solution, so -v set true means v set false

#   if ($satlines[$i] =~ /^v [ ,0123456789]+$/) (lb)  # DG this is wrong (March 28). use this to convert the variables set true in a satisfiable solution, so only v set true is reported

#       print OUT "$satlines[$i]";
       chomp $satlines[$i];
       @words = split(/ +/,$satlines[$i]);
       foreach $word (@words) {
	          chomp $word;
		 #	          print "the word in the split is: $word \n";
		 #	          print OUT "$word \n";

          if ($word eq '0') {
                print OUT " 0";
          }
          else {
            if ((not $word =~ /-/) and ($word ne 'v')) {
             print OUT " $conversion{$word} \n";
            } # end of if not -
          } # end of else
       } # end of foreach word
       #      print OUT "\n";
    } # end of if satlines[i] =~ /[ ,-...]/
} # end of foreach i


# print "The output is in file $cdata.TRANSOL \n";  # DG June 28

close(INSOL);
close(KEY);
close(OUT);
#close(TRACE);
}

###########
#
# check-trees.pl 
# June 10, 2020 DG
# check the solution produced by ILP or SAT to an instance of the problem of computing the Softwired hybridization number
# of a binary matrix, where we are only representing clusters. That is, the DAG is a cluster-network.
# 
# call as: perl check-trees.pl file.TRANSOL datafile
use warnings;
use diagnostics;

sub check_trees {

my	$fname = $_[0];
$INTfile = $fname . ".TRANSOL";
open (INT, "$INTfile"); # this is the TRANSOL file that has the Dimacs variables that are set true, 
                     # translated to their mneumonic variable names
open (CHARDATA, "$fname"); # this is the original binary data for the problem instance
# this is the original binary data for the problem instance.
# In the context of
# these experiments, it is a "...matrixclean" file.

open (OUT, ">DD");
#open (TRACE, '>>trace');
#print TRACE "YYY In check_trees. fname is $fname \n";
open (DAG, ">$fname.gml");
#print "Inside check_trees. The value of fname is $fname \n";

@clines = <CHARDATA>;
chomp @clines;
$num_leaves = @clines;
print TRACE "The number of leaves is $num_leaves \n";

@tlines = transpose(@clines);  # transpose the matrix @clines

# foreach $line (@tlines) {
	#      print TRACE "$line\n";
	#}


@lines = <INT>;

@D = ();
@R = ();

%RCT = ();
%chars = ();
$charnum = 0;
$maxtail = 0;
$maxhead = 0; # this is the largest node number which has an edge out of it in the solution. 

print DAG "graph [directed 1 \n";

foreach $line (@lines) {
  @words = split (/ /,$line);
   foreach $word (@words) {
     if ($word =~/X\((\d+),(\d+),(\d+)\)/) {
        $char = $1;
        $head = $2;
        $tail = $3;

        if (not defined $chars{$char}) {
           $chars{$char} = 1;
           $charnum++;
        }
        if ($tail > $maxtail) {
           $maxtail = $tail;
        }
        if ($head > $maxhead) {
             $maxhead = $head;
        }
        
        print DAG qq(edge [ label " "
        source $head target $tail 
        ]\n);

       if (! defined $nodehash{$head}) {
	      $nodehash{$head} = 1;
              print DAG qq(node[ id $head label "$head"
                  graphics[ Image[ Type "File" Location "white-ball.gif"] ]
                  vgj[ labelPosition "below" ]
             ]\n);
       }

       if (! defined $nodehash{$tail}) {
	      $nodehash{$tail} = 1;
              print DAG qq(node[ id $tail label "$tail"
                  graphics[ Image[ Type "File" Location "white-ball.gif"] ]
                  vgj[ labelPosition "below" ]
             ]\n);
       }

    }
      else {
        if ($word =~ /RCT\((\d+),(\d+)\)/) { # put RCT information into the RCT hash.
           $char = $1;
           $node = $2;
           $RCT{$char} = $node;
#           $RCT{"$char"} = $node;
        } # end of if word == RCT
      } # end of else

       if ($word =~ /T\((\d+)\)/) { # DG Aug.12 - I think that T takes two parameters: char and num, so this is never
	                            # executed. This does not cause a problem, because nodes are picked up as heads or tails
				    # from the X variables.
           $node = $1;
           print DAG qq(node[ id $node label "$node"
                  graphics[ Image[ Type "File" Location "white-ball.gif"] ]
                  vgj[ labelPosition "below" ]
             ]\n);
        } # end of if word == T
   }
} # end of foreach line

 print TRACE "The number of characters is $charnum \n";
 # print TRACE "The largest head index is $maxhead, and the largest tail index is $maxtail \n";
 $first_leaf = $maxtail - $num_leaves + 1;

foreach $char (1 ... $charnum) {
   foreach $head (0 ... $maxtail) {  # maxtail is always larger than maxhead because of the n leaves which have no out edges
                                     # we need to make the upper limit equal to maxtail in order to make a square array
      foreach $tail (0 ... $maxtail) {
        $R[$char][$head][$tail] = $D[$char][$head][$tail] = 0;
      }
    }
}

#$maxtailp1 = $maxtail + 1;
# print "The number of nodes = $maxtailp1 \n";

foreach $line (@lines) {
   @words = split (/ /,$line);
   foreach $word (@words) {
     if ($word =~/X\((\d+),(\d+),(\d+)\)/) {

	     #      print TRACE "$line";
      $char = $1;
      $head = $2;
      $tail = $3;
      $R[$char][$head][$tail] = $D[$char][$head][$tail] = 1;
      #      print TRACE "$char, $head, $tail $R[$char][$head][$tail] \n";
     }
   }
}

foreach $char (1 ... $charnum) {
   print OUT "INPUT TABLE $char \n";

   foreach $head (0 ... $maxtail) {
      foreach $tail (0 ... $maxtail) {
        print OUT  "$R[$char][$head][$tail] ";
      }
      print OUT "\n";
    }
print OUT "\n\n";
}

foreach $char (1 ... $charnum) {
    foreach $k (0 .. $maxtail) {
             $R[$char][$k][$k] = 1;
           }
}

foreach $char (1 ... $charnum) {
   print OUT "INTERMEDIAT TABLE $char \n";

   foreach $head (0 ... $maxtail) {
      foreach $tail (0 ... $maxtail) {
        print OUT  "$R[$char][$head][$tail] ";
      }
      print OUT "\n";
    }
print OUT "\n\n";
}

foreach $char (1 ... $charnum) { # compute the directed transitive closure
    foreach $k (0 .. $maxtail) {
	    #        print TRACE "\n Iteration $k for character $char \n";
        foreach $i (0 .. $maxtail) {
          foreach $j (0 .. $maxtail) {
             $E = (($R[$char][$i][$j]) || ($R[$char][$i][$k] && $R[$char][$k][$j]));
             $R[$char][$i][$j] = $E;
	     #             print TRACE "XX $char, $i, $j: $R[$char][$i][$j],  $E \n";
           }
        }
      }
} # end of foreach char


foreach $char (1 ... $charnum) {
   print OUT "REACHABILITY TABLE $char \n";

   foreach $head (0 ... $maxtail) {
      foreach $tail (0 ... $maxtail) {
        print OUT  "$R[$char][$head][$tail] ";
      }
      print OUT "\n";
    }
print OUT "\n\n";
}

# Now check that the reachability matrices are correct

foreach $char (1 ... $charnum) {
  $spanning = 1;
  foreach $leaf ($first_leaf ... $maxtail) {
     if ($R[$char][0][$leaf] != 1) {
       print TRACE "Incomplete spanning tree for character $char and leaf $leaf \n";
       $spanning = 0;
     }
  }

  if ($spanning == 1) {
	  #       print TRACE "The spanning tree for character $char reaches all of the leaves \n";
  }
}




$Dis = 1;
$char = '1';
foreach $line (@tlines) {
   chomp $line;
   @bits = split (//,$line); 

   $rootnode = $RCT{$char};
   #      print TRACE "The rootnode of the character cluster-tree $char is node $rootnode \n";

       foreach $j (0 .. $num_leaves - 1) {
           $leaf = $first_leaf + $j;
           if ($R[$char][$rootnode][$leaf] != $bits[$j]) {
               print TRACE "Disagreement \n";
               print TRACE "char = $char; j = $j; leaf = $leaf; bits = $bits[$j] \n";
               $Dis = 0;
           }
	   #           else {
		   #    print TRACE "Agreement \n";
		   #}
       }
   $char++;
}

if (($Dis) && ($spanning)) {
   print TRACE "The solution works \n";
   print "The solution works \n\n";

   print DAG "]\n";
  # print "The code for the graphics display of the DAG is in file $fname.gml \n";
   close (DAG);
}
else {
   print "The solution does NOT work \n";
   print "The Dis value is $Dis, and the spanning value is $spanning \n";

   print TRACE "The solution does NOT work \n";
   print TRACE "The Dis value is $Dis, and the spanning value is $spanning \n";

   print DAG "]\n";
   print TRACE "The code for the graphics display of the DAG is in file $fname.gml \n";
   close (DAG);
}

# close(TRACE);
close(INT);
close(CHARDATA);
}

########
sub parse_Newick {

# parse_Newick.pl
# DG, July 16, 2020
# 
# This program reads in (directed) trees (or just a single one) in Newick format and creates matrices (one for each tree) that specify the splits (clusters) in the trees.

my $file = $_[0];
#my $file = @_;
open (IN, "$file"); 
$mfile = $file . 'matrix';
open (OUT, ">$mfile");

$k = 0;  
@C = ();

@start_list = (1);
$tree_list = "";

$treenum = 0;
while ($tree = <IN>) { 
  chomp $tree;
  $treenum++;
  #  print OUT "Tree $treenum \n";
  # print OUT "$tree\n";

  $len = length $tree;
  #  print OUT "$len \n";

@t1 = ();
@t1 = splitwell ($tree);  # split the characters, but group consecutive digits into a single number
# print OUT "@t1";
# print  OUT "\n";
$len = @t1;
# print OUT "$len \n";


$Rtotal = 0;
$pcount = 0;
$maxleaf = 0;
$minleaf = 1;
foreach $i (0..$len-1) {
   if ($t1[$i] eq '(' ) {
      $pcount++;
      $t1[$i] = "L$pcount;";
      next;
   }
     
  if ($t1[$i] eq ')' ) {
      $t1[$i] = "R$pcount;";
      $pcount--;
      $Rtotal++;
      next;
  }

  if ($t1[$i] =~ /(\d+)/) {
      $t1[$i] = "$1";
      if ($1 > $maxleaf) {
          $maxleaf = $1;
      }
      if ($1 < $minleaf) {
          $minleaf = $1;
      }
  }

#  if ($t1[$i] eq ',') {
#      $t1[$i] = ' ';
#  }
}

#foreach $i (0..$len-1) {
	#     print OUT "$t1[$i]";
	#}
#print OUT "\n\n";

#print OUT "The number of columns = $Rtotal \n";
#print OUT "The number of rows = $maxleaf \n";

#if ($treenum > 1) {
#   $k++;
#}

foreach $j ($k .. $k - 1 + $Rtotal) { # DG July 21
    $tree_list .= $treenum;
    foreach $i (0 .. $maxleaf) { 
      $C[$i][$j] = 0;
   }
}

foreach $i (0 .. $maxleaf) {
   foreach $j (0 .. $k - 1 + $Rtotal) {
	   #print "$C[$i][$j]";
   }
   #print "\n";
}

 
# scan the t1 list to find matching R and L pairs of parentheses
foreach $i (0 .. $len-1) {
    if ($t1[$i] =~ /R(\d+)/) {
       $left = "L$1;";
       #print OUT "$t1[$i] $left $k \n";
       $j = $i-1;
       $k++; 

       # scan left to collect the leaf numbers between paired R and L parens
       until ($t1[$j] eq $left) {
          if ($t1[$j] =~ /^(\d+)/) {
             #print OUT "$k:$1 ";
             $C[$1][$k-1] = 1; #DG July 21
          }
          $j--;
       }
       #print OUT "\n";
       #print "$j $t1[$j] \n";
     }
}

push (@start_list, $k);

} # end of while tree = <IN>

print OUT "@start_list\n";
print OUT "$tree_list\n";

#foreach $i (1 .. $maxleaf) {  
foreach $i (0 .. $maxleaf) {  # change Dec. 18, 2020 to see if this solves the problem of not writing the first line of the matrix.
   foreach $j (0 .. $k-1) { # DG July 21
      print OUT "$C[$i][$j]";
   }
   print OUT "\n";
}
close(IN);
close(OUT);
} # end of parse_Newick


################
sub splitwell{

@NA = ();
$S = $_[0];
@SA = split (//,$S);
#print "Inside splitwell @SA \n";
$len = @SA;
# print "Inside splitwell length $len \n";

$num = "";

$working = 0;
$j = 0;
foreach $i (0 .. $len -1) {
#  print "Testing for  j num, SA: $j, $num, $SA[$i] \n";
  if ($SA[$i]  =~ /[(),]/) {
     if ( $working == 1){
       $NA[$j] = int $num;
       $num = "";
       $j++;
       $working = 0;
     }
     $NA[$j] = $SA[$i];
     $j++;
  }
  else {
    $working = 1;
    $num .= $SA[$i];
  }
}

$len = @NA;
return @NA;
}
