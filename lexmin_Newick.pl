# lexmin_Newick.pl
# May 3, 2021 Modify to reduce (if needed) the leaf numbers so that they start from 0. Still assume that the
# leaf numbers are consecutive numbers.
#
# DG April 17, 2021
# This takes in leaf-labeled tree described in Newick notation and for each input tree, 
# produces the unique Newick notation for the tree that is
# lexicographically minimum among all Newick notations for the tree. This then makes it easy to
# determine if two Newick notations describe the same tree, since they do if and only if their lexmin_Newick notations
# are identical. So, we have reduced the isomorphism problem to a string comparison problem. Properly implemented,
# creating the lexmin_Newick notation can be done in linear time.
#
# DG corrected Nov. 20, 2020 to include all rows - the prior version omitted the first row in the output.
#
# DG, July 16, 2020
# 
# This program reads in (directed) trees (or just a single one) in Newick format and creates matrices (one for each tree) that specify the splits (clusters) in the trees.

open (IN, "$ARGV[0]"); 
open (OUT, ">$ARGV[0]trace");
open (FINAL, ">$ARGV[0]lexmin");

$k = 0;  # the number of columns created so far
@C = (); # leaf-char binary matrix output

@start_list = (1); # left end
$tree_list = "";

$treenum = 0;
while ($tree = <IN>) { 
  chomp $tree;
  $treenum++; # count the number of trees input
    print OUT "Tree $treenum \n";
   print OUT "$tree\n";
  %lp = (); # hash for left position
  $rp = (); # list for right position

  $len = length $tree;  # of chars in Newick format of input tree
    print OUT "Number of symbols in input tree: $len \n";

@t1 = ();
@t1 = splitwell ($tree);  # split the characters, but group consecutive digits into a single number
 print OUT "@t1"; # individual tokens in Newick input
 print  OUT "\n";
$len = @t1;  # number of tokens in @t1
 print OUT "The number of tokens in Newick input: $len \n";

$Rtotal = 0;
$pcount = 0;
$maxleaf = 0; # largest leaf number encountered
$minleaf = 1000; # smallest leaf number encountered

foreach $i (0..$len-1) {  # scan through the Newick input
   if ($t1[$i] eq '(' ) {
      $pcount++; # counter for number of unmatched left parens 
      $t1[$i] = "L$pcount"; # replace the '(' with L$pcount
      print OUT "@t1\n";
      next;
   }
     
  if ($t1[$i] eq ')' ) {
      $t1[$i] = "R$pcount"; # replace the ')' with R$pcount
      print OUT "@t1\n";
      $pcount--; # decrement pcount
      $Rtotal++; # increase the number of right parrens seen - why increment here?
      next;
  }

  if ($t1[$i] =~ /(\d+)/) { # if entry i is a number
      $t1[$i] = "$1";  # replace the number with the string representation for that number
      print OUT "@t1\n";
      if ($1 < $minleaf) {
          $minleaf = $1; # update maxleaf if appropriate
      }
      if ($1 > $maxleaf) {
          $maxleaf = $1; # update maxleaf if appropriate
      }
  }
} # end of foreach

  print OUT "The minleaf is $minleaf \n";
  if ($minleaf > 0) {  # if the smallest leaf number is larger than 0, then subtract minleaf from each
                       # leafnumber so that they start with 0, which is assumed in the rest of the program.
                       # It is assumed, also, that the leafnumbers are consecutive.
           foreach $i (0..$len-1) {  # scan through the Newick input
                 if ($t1[$i] =~ /^(\d+)/) {
                    $t1[$i] = $t1[$i] - $minleaf;
                 }
            } # end of foreach i
            $maxleaf--;
  } # end of if minleaf

print OUT "The (possibly reduced) string \n";
foreach $i (0..$len-1) {
	     print "$t1[$i] ";
	     print OUT "$t1[$i] ";
	}
print "\n\n";

print OUT "The number of columns = $Rtotal \n"; # Rtotal gives the number of right parens. Each right
                                                # paren specifies a subtree, and induces one column specifying the
                                                # leaves in that subtree. So Rtotal is one less than the number of
                                                # leaves.
print "The number of columns = $Rtotal \n";
#print OUT "The number of rows = $maxleaf \n";
$maxleafp1 = $maxleaf + 1;
print OUT "The number of rows = $maxleafp1 \n";
print "The number of rows = $maxleafp1 \n";

#if ($treenum > 1) {
#   $k++;
#}

print OUT "About to add more columns to the matrix, initialized to 0\n";
print OUT "The value of k is $k, and the value of Rtotal is $Rtotal\n";
foreach $j ($k .. $k + $Rtotal - 1) { # DG July 21 j is a column number. Rtotal - 1 counts the
                                      # number of columns created in each tree. k seems to be the largest
                                      # column index created so far. where is k updated? 
    $tree_list .= $treenum;
    foreach $i (0 .. $maxleaf) { # is this correct for data where the leaf numbers start at 1?
      $C[$i][$j] = 0; # create the matrix with all zeros.
   }
}

foreach $i (0 .. $maxleaf) {
   foreach $j (0 .. $k + $Rtotal - 1) {
      print "$C[$i][$j]";
      print OUT "$C[$i][$j]";
   }
   print "\n";
   print OUT "\n";
}

 
# scan the t1 list to find matching R and L pairs of parentheses
foreach $i (0 .. $len-1) {
    if ($t1[$i] =~ /R(\d+)/) { # found a right parren, so must search left to find the matching left parren
       $right = "R$1";
       $left = "L$1";
       print OUT "XX Found an R: i:$i, $t1[$i], $left, k:$k \n";
       print "XX Found an R: $i $t1[$i] $left $k \n";
       $caseXY = $caseX = $caseG = 0; # flag to check if a real case has been found

##################### added section April 28, 2021 to find the min between right and left before we
##################### do anything else.

       $q = $i;
       $scan_min = $maxleaf;
       # scan left to collect the leaf numbers between paired R and L parens
       until ($t1[$q] eq $left) {
          if ($t1[$q] =~ /^(\d+)/) {
             print OUT "k, leaf: $k:$1 and update a single cell entry in C to 1\n";
             $C[$1][$k-1] = 1; #DG July 21
	     if ($t1[$q] < $scan_min) {
		     $scan_min = $t1[$q];
	     }
          }
          $q--;
        } # end of until. At this point, left is in position q.
        print OUT "QQQ The min between $left (at position $q)  and $right (at position $i), is $scan_min \n";
        $Lend= $q; # this is the position in t1 of the L paren in this matched pair
	

###############

       $j = $i-1; # j is one position to the left of the position i where the R parren was found
       print OUT "KK j = $j, i = $i, k = $k. incrementing k: $k \n";
       $k++; # DG counting the number of right parens found so far.

       print OUT "The token to the left of $right is $t1[$j] \n";

       # This next takes care of the special case that the right parren matches a left parren
	                          # with only two nodes between them, i.e., the smallest possible pair of subtrees

       if (($t1[$j] =~ /\d+/) && ($t1[$j-1] == ',') && ($t1[$j-2] =~ /^\d+/)  &&  ($t1[$j-3] =~ /L\d+/))  {  
          print OUT "Yes, we are in case (x,y)\n";
          $caseXY = 1;

       if ($t1[$j] =~ /(\d+)/) {
	       $digitR = $1;
               print OUT "Yes, we found digitR, it is $digitR\n"; 				  
               print "Yes, we found digitR, it is $digitR\n"; 				  
               }
 
               if (($t1[$j-1] == ',') && ($t1[$j-2] =~ /^(\d+)/)) {  
		       $digitL = $1;
                       print "Yes, we found digitL, it is $digitL\n";
                       print OUT "Yes, we found digitL, it is $digitL\n";
               }

#                    if ($1 != $2) {
#                        print OUT "Error ($1 $2) in identifying the case of right and left parens matching with only
#                                   two nodes between them \n";
#                    }


		       if ($digitR < $digitL) {
			       $t1[$j] = $digitL;
			       $t1[$j-2] = $digitR;
		       }
		       print OUT "ZZ @t1\n";
       } # end of if in case (x,y) 

       else { # *
            if (($t1[$j-1] eq ',') && ($t1[$j-2] =~ /^R/ )) { #** this is the case of ( ...),x)
                  print OUT "Yes, we are in the case of (...),x) scan_min is $scan_min\n";
                  $caseX = 1;
                  if ($t1[$j] == $scan_min) { #*** exchange the number to the left of right with the
                                              # string from left to the comma to the left of it.
                        print OUT "HH j $j, t1[j] $t1[$j], scan_min $scan_min q $q \n";
                             for ($jj = $j-2; $jj >= $q+1; $jj--) {
                                   $t1[$jj+2] = $t1[$jj]; # move the string two places to the right
  	                           print OUT "ZZZ @t1\n";
                             } #end of for
                             $t1[$q+1] = $scan_min;
                             $t1[$q+2] = ',';
		             print OUT "ZZZZ @t1\n";
                    } # match to ***
              } # match to **
        } # match to *


             # general case
            if ($t1[$j] =~ /^R(\d+)/) { # the token to the left of the right parren $right, should also be an R.
                    print OUT "Yes, we are in the general case ((  ),(  )) or (x,(   ))\n";
                    $caseG = 1;
                    print OUT "GG Token to the left of the right paren $right is also an R \n";
		    $SR = "R$1"; # record what that token is.
		    $SL = "L$1";  # set SL to the left parren that matches that R.
                    print OUT "The matched SL and RL parens are $SR and $SL \n";
              $q = $j;
              $SLscan_min = $maxleaf;
              # scan left to find the min leaf number between paired SR and SL parens

              until ($t1[$q] eq $SL) {
                  if (($t1[$q] =~ /^(\d+)/) && ($t1[$q] < $SLscan_min)) {
         	       $SLscan_min = $t1[$q];
                       print OUT "SLscan_min changed to $SLscan_min. New min is in position $q \n";
	          }
                  $q--;
              } # end of until

              if ($scan_min == $SLscan_min) {
                  print OUT "The min between $SL and $SR is in the second term, and so the two terms should be interchanged.\n"; 


                        @sav = ();
                        print OUT "HHH j $j, t1[j] $t1[$j], SLscan_min $SLscan_min q $q \n";
                             for ($jj = $q; $jj <= $j; $jj++) {
                                   $sav[$jj] = $t1[$jj]; # move the second term to sav
                                   print OUT "A jj $jj, t1[jj] $t1[$jj] \n";
                             } #end of for
                             
                             for ($jj = $q - 2; $jj >= $Lend+ 1; $jj--) {
                                   $newplace = $jj + ($j - $q + 2);
                                   $t1[$newplace] = $t1[$jj];
                                   print OUT "B jj $jj, newplace $newplace, t1[jj] $t1[$jj] \n";
                             }
		             print OUT "W @t1\n";
                             # Now we have to insert what we saved into the empty space at the left.

                             for ($offset = 0; $offset <= $j - $q; $offset++) {
                                  $t1[$Lend+ 1 + $offset] = $sav[$q + $offset];
                                   print OUT "C q $q, offset $offset, Lend $Lend, sav[q + offset] $sav[$q + $offset]   \n";
                             }
                             $t1[$Lend + $offset + 1] = ',';
		             print OUT "WWW @t1\n";

              } # end of if (scan_min
	    } # end of if for general case

	    if ($caseXY + $caseX + $caseG == 0) {
		    print "YIKES unmatched parens, \n"; # as a safety - make sure the above logic
		    print OUT "YIKES unmatched parens, \n"; # as a safety - make sure the above logic
		                                                            # matches reality
	    }

       } # end of found a right paren 

     #  print OUT "PP The min node between $right and $left is $scan_min \n";

     } # end of foreach

push (@start_list, $k); # add the current column index to the start_list, which shows the starting positions
                        # in the matrix (and hence starting from 1)

print OUT "$tree\n"; # tree is the original input string for the tree in Newick format.
print OUT "@t1\n";

foreach $element (@t1) {
	if ($element =~ /L/) {
		print OUT "(";
		print FINAL "(";
	}
	else {
		if ($element =~ /R/) {
	          print OUT ")";
	          print FINAL ")";
                }
        	else {
                   print OUT "$element";
                   print FINAL "$element";
	        }
	}
}
print OUT "\n";
print FINAL "\n";

} # end of while tree = <IN>

print OUT "@start_list\n";
print OUT "$tree_list\n";

foreach $ii (0 .. $maxleaf) { # 
   foreach $jj (0 .. $k-1) { # DG July 21
      print OUT "$C[$ii][$jj]";
   }
   print OUT "\n";
}

close (FINAL);
close (OUT);

################
sub splitwell{

@NA = ();
$S = $_[0];
@SA = split (//,$S);
 print "Inside splitwell @SA \n";
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
