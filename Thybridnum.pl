# Thybridnum.pl 
# July 29 modified for input from clean_for_hybridtrees. The input has a FIRST line 
# that gives the borders between the trees output from clean_for_hybridtrees (after cleaning the input to clean_for_hybridtrees).
#
# Then the ILP and CNF produced is for the problem where the DAG created must display each of the original trees, and not just
# their splits, as in the case of Bhybridnum.pl
#
#
# DG Aug. 13, 2020 the comments are left over from Bhybridnum.pl and so need to be updated.
# Call as perl Thybridnum.pl filename  reticulation_target 
#
# July 18, 2020
# DG
# June 6 modify to set the number of interior nodes as an input parameter. This is so we can try out reducing the size of the complete graph to 
# see if it effects efficiency and correctness.
#
# The underlying super-graph consists of a root node 0, plus a complete graph consisting of $interior_num nodes, and then $n leaves. Node 0 has
# a directed edge to each of the other nodes, and each interior node has a directed edge to each of the leaves.
#
#
#May 30, 2020 Add in the CNF to determine if there is a feasible solution where the number of reticulation nodes is at most equal to the
# input target value.
#
# DG I think the program will have a bug if the input contains a column with only a single value of 1
#
# DG  May 30, 2020 The logic in this program is somewhat involved and there is no doubt that there is redundancy, both in the logic and in the programming.
# Sometime later, perhaps, this should be streamlined. But since it works, we will leave it for now
#
# Bhybridnum.pl May 20, 2020 convert to remove the F variables, so that all the variables are binary.
#
# hybridnum.pl May 16, 2020 convert to compute the hybridization number - just add the constraint that the
# indegree of any node is at most two
#
# historyilp3.pl This seems to be the right one to use - April 24, 2020
# use strict
#
# Feb. 25, 2017
# converted from historyilp2.pl
# change the hypercube of dimension m, to a complete graph with m + n nodes and another n leaves.
#

%SATvariable = (); # This is a hash indexed by ILP variables, with value equal to the integer used to encode that variable in Dimacs format.


@lines = ();

open (IN, $ARGV[0]);  # the file with a first line showing the borders of the trees, and then the binary input detailing the splits.1G
print "On entry to Thybridnum.pl, the input file is $ARGV[0] and the target is $ARGV[1]\n";
open (OUT, ">$ARGV[0].lp");
open (SAT_KEY, ">$ARGV[0].KEY");
open (SAT_FINAL, ">$ARGV[0].SAT");


$target = $ARGV[1]; # target is the maximum number of allowed reticulation nodes.

$borders = <IN>; # borders lists the indices separating the columns in the different trees
print "The input border line in Thybridnum:  $borders";
@tborders = ($borders =~ /(\d+)/g);

#print "@tborders\n";
$numt = @tborders;
$numtrees = $numt - 1;
print "There are $numtrees trees \n";
%treenum = ();

foreach $i (0 .. $numtrees-1) {
      if ($i > 0) {
      $left = $tborders[$i] + 1;
      }
      else {
         $left = 1; 
      }
      $right = $tborders[$i + 1];

      foreach $character  ($left .. $right) {
         $treenum{$character} = $i+1;
      }
$lastright= $right;
}

$m = $lastright;
#print "The last right value is $lastright, which should equal to $m \n";

foreach $character (sort keys %treenum) {
#    print "C character $character is in tree $treenum{$character} \n";
}

while ($line = <IN>) {

  push(@lines, $line);
 # print "$line";
  chomp $line;
  } # end while line = <IN>


%inline = ();
$n = @lines;
# print "$n, $m\n";  # the number of rows and colums in the input matrix

$maxinterior_num = $n - 2 + 2 * $target; # this is the maximum number of interior nodes (in the complete graph) 
                                         # ever needed, given  the target.

#$interior_num = int ($maxinterior_num * 0.75); # start with a smaller number of interior nodes to speed up
                                                # start with a smaller number of interior nodes to speed up
                                               # computation. However, if you get UNSAT, then recompute
					       # using maxinterior_num.
$interior_num = int ($maxinterior_num * 1); # April 4, 2021 - change to 1 because we were running into cases
                                            # where we were getting UNSAT with too few interior nodes.

print "The maxinterior value is $maxinterior_num; the interior number we will try is $interior_num, and the target is $target \n";
$base = $interior_num + 1;
foreach $offset (0 .. $n-1) {
  $index = $base + $offset; 
  $inline{$index} = $lines[$offset];
}

# print (keys %inline);

$binaries = "";
$clausecount = 0;


$SATinteger = 1;
print OUT "Minimize \n"; 

 foreach $num (1 .. $n + $interior_num) {   # generate the R variables for the objective function. R stands for Reticulation, i.e., recording whether that
                              # node is a reticulation node or not. DG we might be able to make this 1 .. $n+$m, but this is safer for now. This
                              # will count a leaf as a reticulation node if it has two or more edges into it. If we change to 1 .. $n+$m, we must
                              # also forbid a leaf from being a reticulation node. The way it is here will allow a leaf to be a reticulation node,
                              # and count it as one.
     $binaries .= "R($num)\n";
     print OUT "+ R($num) ";  # this is printing to the objective function.
     $SATvariable{"R($num)"} = $SATinteger;
     $SATinteger++;

}

 foreach $num (0 .. $interior_num) { # generate the RCT variables for the roots of the cluster trees. No leaf is allowed to be a root of a cluster tree.
                              # Also use this in the larger block to create the X variables.
   foreach $character (1 .. $m) {
     $binaries .= "RCT($character,$num)\n";
     $SATvariable{"RCT($character,$num)"} = $SATinteger;
     $SATinteger++;
   }  # end of foreach character

    foreach $tree (1 .. $numtrees) {
     foreach $fnum ($num + 1 .. $n + $interior_num) { # create the X variables for edges from 0 through the internal nodes, to higher index nodes through the leaves.
        $binaries .= "X($tree,$num,$fnum)\n";
        $SATvariable{"X($tree,$num,$fnum)"} = $SATinteger;
        $SATinteger++;
     } # end of foreach fnum
   } # end of foreach tree 
} # end of foreach num
      

foreach $num (0 .. $n + $interior_num) { # generate the T and CT variables for the character-spanning trees  and the cluster trees
   foreach $tree (1 .. $numtrees) {
     $binaries .= "T($tree,$num)\n";
     $SATvariable{"T($tree,$num)"} = $SATinteger;
     $SATinteger++;

     # generate variables to record whether there is an edge into node num in the character-spanning tree
     $SATvariable{"Ti($tree,$num)"} = $SATinteger;
     $SATinteger++;

      # generate variables to record whether there is an edge out of node num in the character-spanning tree
     $SATvariable{"To($tree,$num)"} = $SATinteger;
     $SATinteger++;

     }  # end of foreach tree

     foreach $character (1 .. $m) {
        $tree = $treenum{$character};

        $binaries .= "CT($character,$num)\n"; 
        $SATvariable{"CT($character,$num)"} = $SATinteger;
        $SATinteger++;
     }
} # end of foreach num


foreach $head (0 .. $interior_num) {  # generate the Y and Z variables
   foreach $tail ($head + 1 .. $n + $interior_num) {
      $binaries .= "Y($head,$tail)\n";
      $SATvariable{"Y($head,$tail)"} = $SATinteger;
      $SATinteger++;
     
      foreach $character (1 .. $m) {
           $binaries .= "Z($character,$head,$tail)\n";
           $SATvariable{"Z($character,$head,$tail)"} = $SATinteger;
           $SATinteger++;
      } # end of foreach character
   }  # end of foreach tail
} # end of foreach head
# print "After generating the Z variables and before generating the I variables, SATinteger = $SATinteger \n";

$SATvariable{"RT(1,0)"} = $SATinteger; # initial value for RT
$SATinteger++;

foreach $num (2 .. $n + $interior_num + 1) { # create the CNF variables to count the number of reticulation nodes. # DG June 2, 2020
   foreach $sum (0 .. $num) {
        if ($sum < $target + 1) {
           $bt = $sum;
        }
        else {
           $bt = $target + 1;
        }

        $SATvariable{"RT($num,$bt)"} = $SATinteger; # the meaning of RT(num,bt) is that the total number of reticulation nodes counted before
                                                    # examining node $num, is at most $bt. The max before examining node num cannot be larger than
                                                    # num - 2, and we are only interested in whether the 
                                                    # number of reticulation nodes reaches target+1.
        $SATinteger++;
   }
} # end of foreach num

$CNFstring = "";
$after_leaves = $n + $interior_num + 1;
$targetp1 = $target + 1;

$v = $SATvariable{"RT($after_leaves,$targetp1)"};
$CNFstring .= "-$v 0 \n"; # forbid the number of reticulation nodes to exceed the value $target.
$clausecount++;

foreach $num (1 .. $n + $interior_num) {    # generate the I variables.
   $SATvariable{"I(0,$num)"} = $SATinteger;
   $CNFstring .= "-$SATinteger 0 \n";  # there are no edges coming into node 0
   $clausecount++;
   $SATinteger++;

   foreach $baltnum (1 .. $num - 1) {
      if ($baltnum <= $interior_num + 1)  { # don't include edges between leaves, but allow 
                                      # I(baltnum,num) and RI(baltnum,num) to have baltnum
                                      # equal to the first leaf (this is a kluge just for technical
                                      # programming reasons, not logical reasons).
          $SATvariable{"I($baltnum,$num)"} = $SATinteger;
          $SATinteger++;

          $SATvariable{"RI($baltnum,$num)"} = $SATinteger;
          $SATinteger++;
      }
    }
} # end of foreach num
      

%charedgehash = ();
%edgehash = ();
print OUT "\n\nsubject to \n";

foreach $tree (1 .. $numtrees) {
   print OUT "T($tree,0) = 1 \n"; # Node 0, the root of the super-network, is in the spanning tree for each character.
   $binaries .= "T($tree,0) \n";

   $Tvar = $SATvariable{"T($tree,0)"};
   $CNFstring .= "$Tvar 0 \n";
   $clausecount++;

 foreach $num (1 .. $interior_num) { # all the internal nodes.  #DG this is more complex than it should be, because num is used as a second parameter in the
                               # first part of this block, and as a first parameter in the second part of this block. This should be cleaned up someday.
                               # change 0 to 1 in the foreach statement June 2, 2020
      $Tvar = $SATvariable{"T($tree,$num)"};
      $Tivar = $SATvariable{"Ti($tree,$num)"};
      $Tovar = $SATvariable{"To($tree,$num)"};

      $Ti_string = "";
      $To_string = "";
      $Xinto = "";
      $Xoutof = "";

     foreach $baltnum (0 .. $num - 1) { # generate the edge-into node $num inequalities. These go from 0 and the internal nodes to higher index internal nodes.
                                        # note that when num = 0, this sub-block is skipped, which is what we want.
 
            $key = "$tree,$baltnum,$num"; # there is redundancy in setting up charedgehash, with the code in the
                                               # next loop, but it works.
            $Xinto .= "+ X($key) "; # for the X into constraint in the spanning-character tree
            $v1 = $SATvariable{"X($key)"};
            $CNFstring .= "$Tvar -$v1 0 \n"; # if there is an edge into a node, then the node is in the spanning-character tree
            $clausecount++;
            $Ti_string .= "$v1 "; # accumulate all of the variable integers for the edges incoming to node num.
             
            $CNFstring .= "$Tivar -$v1 0 \n"; # if there is an edge into a node, then set the Ti variable for the node for the spanning-character tree
            $clausecount++;

            $Xkey = "$baltnum,$num";
            $charedgehash{$key} = 1;
            $edgehash{$Xkey} = 1;


     } # end of foreach baltnum

     $CNFstring .= "$Ti_string -$Tivar 0 \n"; # Tivar for char,num can only be set to 1 if there is some X(char,baltnum,num) set to 1
     $clausecount++;


      foreach $faltnum ($num + 1 .. $n + $interior_num) { # generate the edge-out-of inequalities, from node $num to node $faltnum.
                                                   # these edges go from node 0 and the internal nodes to higher index
                                                   # internal nodes and to leaves.

             $key = "$tree,$num,$faltnum";
             $Xoutof .= "- X($key) "; # edge for the spanning $character tree. 


             $Xkey = "$num,$faltnum";
             $charedgehash{$key} = 1;
             $edgehash{$Xkey} = 1;

             $binaries .= "X($key)\n";
             print OUT "X($key) - T($tree,$num) <= 0 \n"; # DG May 21 if there is an edge out of node num, then node num must be in the spanning-character tree.
            
             $v1 = $SATvariable{"X($key)"};
             $CNFstring .= "$Tvar -$v1 0 \n";
             $clausecount++;
             $To_string .= "$v1 "; # accumulate all of the variable integers for the edges out of node num.

             $Tovar = $SATvariable{"To($tree,$num)"};
             $CNFstring .= "$Tovar -$v1 0 \n"; # if there is an edge out of a node, then set the To variable for the node for the spanning-character tree
             $clausecount++;
      } # end of foreach faltnum

     $CNFstring .= "$To_string -$Tovar 0 \n"; # Tovar for char,num can only be set to 1 if there is some X(char,num,faltnum) set to 1
     $clausecount++;

 if ($num > 0) {
     print OUT "$Xinto <= 1 \n";
     print OUT "$Xinto - T($tree,$num) = 0\n"; # node num (for the internal nodes) is in the spanning-character-tree iff there is an edge into node num

     print OUT "$Xinto  $Xoutof <= 0\n\n"; # If there is an edge into an internal node, there must be an edge out of it.
                                          # DG we also need a constraint to say that there can be edge out of an internal node
                                          # only if there is an edge into it. Actually, this is implied by the clauses establishing that a) if
                                          # a node is in the tree, then there must be an edge into it - unless is is the root; b) if there is
                                          # an edge out of a node, then it is in the tree.


     $Tivar = $SATvariable{"Ti($tree,$num)"}; # for any internal node and any character, there is an edge into the node in the character-spanning tree
                                                # iff there is an edge out of the node in the character-spanning tree
     $Tovar = $SATvariable{"To($tree,$num)"};
     $CNFstring .= "$Tivar -$Tovar 0 \n";
     $CNFstring .= "$Tovar -$Tivar 0 \n";
     $clausecount = $clausecount + 2;
  } # end of if num > 0 
} # end of foreach num 

# foreach $tree (1 .. $numtrees) (left squiggle)  # DG June 3, 2020 foreach character, prohibit more than one incoming edge into  
                               # node num in the spanning tree for that  character 

  foreach $num (1 .. $n + $interior_num) {  
     foreach $baltnum (0 .. $num - 1) {
         foreach $altnum ($baltnum + 1 .. $num - 1) {
            if ($altnum <= $interior_num) {
                $v1 = $SATvariable{"X($tree,$baltnum,$num)"};
                $v2 = $SATvariable{"X($tree,$altnum,$num)"};
                $CNFstring .= "-$v1 -$v2 0 \n";
                $clausecount++;
            }
         } # end of foreach altnum
     } # end of foreach baltnum
   } # end of foreach num

 foreach $num ($interior_num + 1 .. $n + $interior_num) { # for each of the leaves
     print OUT "T($tree,$num) = 1 \n"; #  every leaf must be in each character-spanning tree
     $Tvar = $SATvariable{"T($tree,$num)"};
     $CNFstring .= "$Tvar 0 \n";
     $clausecount++;


     $into = "";
     $Xinto = "";
     foreach $baltnum (0 .. $interior_num) { # generate the into leaf $num inequalities. #changed the initial value from 1 to 0 Feb. 28
#            print ("Into $num from $baltnum \n");
             $Xinto .= "+ X($tree,$baltnum,$num) "; # for the X constraint for the spanning $character tree

             $v1 = $SATvariable{"X($tree,$baltnum,$num)"}; # DG this might be redundant, but I can try to take it out later
             $CNFstring .= "$v1 ";
      }

 print OUT "$Xinto = 1\n\n";
 $CNFstring .= "0 \n";
 $clausecount++;

 } # end of "foreach leaf"


                   # This implements the CNF that if a node is in the spanning tree for some character, then the node (other than node 0) must
                   # have an edge into it in that spanning tree.
# foreach $tree (1 .. $numtrees) (left squiggle) # DG July 21 
     foreach $num (1 .. $n + $interior_num) {
        $Xstring = "";
        foreach $baltnum (0 .. $num - 1) {
              if ($baltnum < $interior_num + 1) {
                  $v = $SATvariable{"X($tree,$baltnum,$num)"};
                  $Xstring .= "$v ";
              }
         }
         $v1 = $SATvariable{"T($tree,$num)"};
         $CNFstring .= "$Xstring -$v1 0 \n";
         $clausecount++;
     }
} # end of "foreach tree" 

 foreach $num ($interior_num + 1 .. $n + $interior_num) { # for each of the leaves # DG May 21

     $v = $SATvariable{"Y(0,$num)"};
     $Xstring = "";
     $XstringY = "";

     foreach $tree (1 .. $numtrees) {
         $Xstring .= " + X($tree,0,$num)";  # This is a kluge, because the 0-leaf case for X didn't work before adding it.

## May 25      print OUT " + X($tree,0,$num)";  # This is a kluge, because the 0-leaf case for X didn't work before adding it.
         $v1 = $SATvariable{"X($tree,0,$num)"};
         $XstringY .= "$v1 ";
         $CNFstring .= "$v -$v1 0 \n"; # if for the given char, there is an edge from root 0 to leaf num, then set Y(0,num) to 1.
         $clausecount++;
     }

     print OUT "$Xstring - Y(0,$num) >= 0 \n\n";  # for each leaf, num, set Y(0,num) to 1 only if there is some char such that
                                          # (0,num) is an edge in the spanning tree for char, i.e. only if X(char,0,num) is 1.

     print OUT "$Xstring - $m Y(0,$num) <= 0 \n\n"; # for each leaf num, if there is a char such that edge (0,num) is an edge in the
                                            # spanning tree for that char, then set Y(0,num) to 1.


         $CNFstring .= "$XstringY -$v 0 \n"; # for each leaf, num, set Y(0,num) to 1 only if there is some char such that
                                              # (0,num) is an edge in the spanning tree for char, i.e. only if X(char,0,num) is 1.

         $clausecount++;


 } # end of foreach num, for the leaves,


foreach $edge (sort keys %edgehash) {    # create the inequalities for the Y(i,j) variables.
        $Xstring = "";
        $XstringY = "";
        $v = $SATvariable{"Y($edge)"};
        foreach $tree (1 .. $numtrees) {
                  $key = $tree . ',' . $edge;
                  $Xstring .= "+ X($key) ";  # sum the X variables for $edge over the different $characters
                  $v1 = $SATvariable{"X($key)"};
                  $XstringY .= "$v1 ";
                  $CNFstring .= "$v -$v1 0 \n"; # if for the given char, there is an specified edge, then set Y(edge) to 1.
                  $clausecount++;

            }

         print OUT "$Xstring - Y($edge) >= 0\n"; # DG May 17 if edge $edge is not used for any $character, set Y to 0. 
         print OUT "$Xstring - $m Y($edge) <= 0 \n"; # if edge $edge is used for some $character, set Y to 1.


         $CNFstring .= "$XstringY -$v 0 \n"; # for each edge, set Y(edge) to 1 only if there is some char such that
                                          # (0,num) is an edge in the spanning tree for char, i.e. only if X(char,edge) is 1.
         $clausecount++;
} # end of foreach edge


foreach $num (1 .. $n + $interior_num) {    # create the inequalities that set the I and R variables, based on the Y variables.
   $Ystring = ""; #DG May 18
   foreach $baltnum (0 .. $num - 1) {
     if ($baltnum < $interior_num + 1)  { # don't include edges between leaves.
        $Ystring .= "+ Y($baltnum,$num) "; # DG May 18

                          # Y(baltnum,num) => I(baltnum+1,num)
        $baltnump1 = $baltnum + 1;

            $v = $SATvariable{"Y($baltnum,$num)"};
            $v2 = $SATvariable{"I($baltnum,$num)"}; 

        if ($baltnump1 < $num) { 
            $v1 = $SATvariable{"I($baltnump1,$num)"};
            $CNFstring .= "$v1 -$v 0 \n";
            $clausecount++;
       

                         # I(baltnum,num) => I(baltnum+1,num)
            $CNFstring .= "$v1 -$v2 0 \n";
            $clausecount++;
        }

                        # Y(baltnum,num) AND I(baltnum,num) => R(num)
         if ($num > 1) {
            $v3 = $SATvariable{"R($num)"};
            $CNFstring .= "$v3 -$v2 -$v 0 \n";
            $clausecount++;
         }

                        # Y(baltnum,num) AND I(baltnum,num) => RI(batlnum +1,num)
         if (($baltnump1 > 1) and ($baltnum < $interior_num + 1) and ($baltnump1 < $num)) {  # DG May 27
              $v4 = $SATvariable{"RI($baltnump1,$num)"};  
              $CNFstring .= "$v4 -$v2 -$v 0 \n";
               $clausecount++;
          }

          $v4 = $SATvariable{"RI($baltnump1,$num)"};  
          $v5 = $SATvariable{"RI($baltnum,$num)"};  
 
                        # RI(baltnum,num) => RI(baltnum+1,num) but baltnum must be larger than 0, since (0,num) can't be the first edge into num in the
                        # edge ordering
         if (($baltnum > 0) and ($baltnum < $interior_num) and ($baltnump1 < $num)) { # DG May 27,28
              $CNFstring .= "$v4 -$v5 0 \n";
              $clausecount++;

         }
         
         if ($baltnum > 0) { # DG May 27

                      # RI(batlnum,num) => NOT Y(baltnum,num)
               $CNFstring .= "-$v -$v5 0 \n";
               $clausecount++;
          }


     } # end of if baltnum < interior_num + 1 
   } # end of foreach baltnum

   print OUT "$Ystring <= 2 \n";  # DG May 18 # DG May 16, change to compute the hybridization number. This is actually redundant give the next constraint.
   print OUT "$Ystring - R($num) <= 1 \n"; # DG May 18 # DG May 16, 2020 change to compute the hybrization number # June 21 change from 2 R(num)

} # end of foreach num

# for each character, pick one root of the cluster-tree for that character.

foreach $character (1 .. $m) {    # a node can be the root of the cluster-tree for $character only if it is in the spanning-tree for
                                  # that $character.
   $RCTstring = "";
   $RCTstringSAT = "";
   $tree = $treenum{$character}; # DG July 18

   foreach $num (0 .. $interior_num) {

     print OUT "RCT($character,$num) - T($tree,$num) <= 0\n"; # if node num is made the root of the cluster tree for a character,
                                                                   # then it must be in the spanning-tree for that character. 

     $v1 = $SATvariable{"RCT($character,$num)"}; # if a node is the root of a cluster tree for a character, then it must be in the spanning-tree for
                                                 # that character, and it must be in the cluster-tree for that character
     $v2 = $SATvariable{"T($tree,$num)"};
     $v3 = $SATvariable{"CT($character,$num)"};
     $CNFstring .= "$v2 -$v1 0 \n";
     $CNFstring .= "$v3 -$v1 0 \n";
     $clausecount = $clausecount + 2;


     print OUT "RCT($character,$num) - CT($character,$num) <= 0\n"; # if a node num is made the root of the cluster tree for a character
                                                                    # then it must be in the cluster tree for that character.

     $RCTstring .= "+ RCT($character,$num) ";
     $RCTstringSAT .= "$v1 ";
   }

  print OUT "$RCTstring = 1\n";  # for each character, pick one root of the cluster-tree for that character.
  $CNFstring .= "$RCTstringSAT 0 \n";
  $clausecount++;


  print OUT "CT($character,0) - RCT($character,0) = 0 \n"; # DG June 10. the 0 node is the root of a cluster-tree for a character if and only if it is in
                                                            # the cluster tree for that character.
  $v1 = $SATvariable{"CT($character,0)"};
  $v2 = $SATvariable{"RCT($character,0)"};
  $CNFstring .= "$v2 -$v1 0 \n"; # if the 0 node is in the cluster tree for a character, then it is the root of that tree. The clause for the other direction
                                     # was created above.
  $clausecount++;

} # end of foreach character

#DG June 5 for each char, generate all pairs of possible roots of the cluster tree for that character, and then
# a clause that says that at least one of those possible roots is not chosen.

foreach $character (1 .. $m) { # DG June 5, 2020 foreach character, prohibit more than one choice of root for the cluster tree for that char.

  foreach $num (1 .. $interior_num) {  # DG June 5 only allow interior nodes, plus 0, to be the root of a cluster tree. DG later, we should disallow 0 from
                                 # being the root. We can do that with an explicit clause implementing -RCT(c,0)

     foreach $baltnum (0 .. $num - 1) { # DG June 5 if we disallow 0 from being a root,  I think we can change 0 to 1 here.
                $v1 = $SATvariable{"RCT($character,$num)"};
                $v2 = $SATvariable{"RCT($character,$baltnum)"};
                $CNFstring .= "-$v1 -$v2 0 \n";
                $clausecount++;
     }
   }
} # end foreach character


foreach $character (1 .. $m) {
    $tree = $treenum{$character} ;
    foreach $baltnum (0 .. $interior_num) {
        foreach $num ($baltnum + 1 .. $n + $interior_num) { # generate the inequalities to forward propogate the CT relations 
            print OUT "CT($character,$baltnum) + X($tree,$baltnum,$num) - CT($character,$num) <= 1\n"; # forward propogation

            $v1 = $SATvariable{"CT($character,$baltnum)"};
            $v2 = $SATvariable{"X($tree,$baltnum,$num)"};
            $v3 = $SATvariable{"CT($character,$num)"};

            $CNFstring .= "$v3 -$v1 -$v2 0 \n"; # forward propogation
            $clausecount = $clausecount + 1;
        }        
     }


     foreach $num ($interior_num + 1 .. $n + $interior_num) {  #  look through the leaves  

                  $v1 = $SATvariable{"CT($character,$num)"};

                  $line = $inline{$num};
                  @chars = split (//,$line);
                  if ($chars[$character - 1] eq '1') {
#                     print "HEY $character $num \n";
                     print OUT "CT($character,$num) = 1 \n";  # generate the inequality that says a leaf with the
                                                            # given character must be reached in the cluster-tree
                                                            # for that character.
                     $CNFstring .= "$v1 0 \n";
                     $clausecount++;
                  }
                  else {
                     print OUT "CT($character,$num) = 0 \n";  # generate the inequality that says a leaf without  the
                                                              # given character must not be reached in the cluster-tree
                                                              # for that character.
                     $CNFstring .= "-$v1 0 \n";
                     $clausecount++;
                  }
     }

} # end of foreach $character


foreach $character (1 .. $m) {  # create the inequalities for the Z relations 
   $tree = $treenum {$character};
   foreach $num (1 .. $interior_num) {  # only for the interior nodes # DG May 31, change 1 to 0 to include the zero node. June 4 changed back to 1
       $CZsum = "";  # string to sum the Z variables in the only-if constraints for CT
       $CZstring = "";

       foreach $baltnum (0 .. $num - 1) {
          $key = "$character,$baltnum,$num"; 
          $CZsum .= "- Z($key) ";
          print OUT  "2 Z($key) - X($tree,$baltnum,$num) - CT($character,$baltnum) <= 0 \n"; # set Z to 1 only if CT and X are both 1.

          $v1 = $SATvariable{"Z($key)"}; # DG June 4 This is actually forward implication. Z is set to 1 if
                                                             # CT is 1 for a node with an edge into node num. 
          $CZstring .= "$v1 ";

          $v2 = $SATvariable{"CT($character,$baltnum)"};
          $v3 = $SATvariable{"X($tree,$baltnum,$num)"};
          $CNFstring .= "$v1 -$v2 -$v3 0 \n";
          $clausecount++;

          $CNFstring .= "-$v1 $v2 0 \n";  # DG June 4.  this is the only if direction. Z(c,balt,num)  is set to 1 only if CT(c,balt) and X(c,balt,num)
          $CNFstring .= "-$v1 $v3 0 \n";
          $clausecount = $clausecount + 2;

 
 
          $v5 = $SATvariable{"CT($character,$num)"};
          $CNFstring .= "-$v1 $v5 0 \n"; # June 4, for symmetry, and because it can't hurt, add the clause that Z(c,balt,num) implies CT(c,num).
                                         # the clause for Z(c,balt,num) implies CT(c,balt) is generated above.
          $clausecount++;
         
        } # end of foreach baltnum

    print OUT "CT($character,$num)  $CZsum - RCT($character,$num) <= 0 \n"; # notice that each entry in CZnum is negated.

    $v5 = $SATvariable{"CT($character,$num)"};
    $v4 = $SATvariable{"RCT($character,$num)"};
    $CNFstring .= "$CZstring $v4 -$v5 0 \n";  # DG June 4 when node num is in the cluster tree for char, then either node num has an edge into it from
                                              # a node baltnum that is in the cluster tree for char, or node num is the root of that cluster-tree.
    $clausecount++;
  } # end of foreach num  for interior nodes


 foreach $num ($interior_num + 1 .. $n + $interior_num) { # Here we do the only if logic for the leaves
    $CZsum = "";
    $CZstring = "";
 

   foreach $baltnum (0 .. $interior_num) {
     $key = "$character,$baltnum,$num"; 
     print OUT  "2 Z($key) - X($tree,$baltnum,$num) - CT($character,$baltnum) <= 0 \n"; # set Z to 1 only if CT and X are both set to 1.
     $CZsum .= "- Z($key) "; # accumulate the Z values. Note the negative sign here

     $v1 = $SATvariable{"Z($character,$baltnum,$num)"};

     $CZstring .= "$v1 ";

     $v2 = $SATvariable{"CT($character,$baltnum)"};
     $v3 = $SATvariable{"X($tree,$baltnum,$num)"};
     $CNFstring .= "$v1 -$v2 -$v3 0 \n";
     $clausecount++;


     $CNFstring .= "-$v1 $v2 0 \n"; # Z implies CT #DG June 1
     $CNFstring .= "-$v1 $v3 0 \n"; # Z implies X
     $clausecount = $clausecount + 2;
     

     $v5 = $SATvariable{"CT($character,$num)"};
     $CNFstring .= "-$v1 $v5 0 \n";
     $clausecount++;


   }

  print OUT "CT($character,$num)  $CZsum <= 0\n"; # if leaf-node num is placed in the cluster-tree for a character, then
                                                                         # node num must be the tail of an edge
                                                                         # in the spanning tree for that character, where the head of that edge is
                                                                         # in the cluster-tree for that character. Recall that all the terms that
                                                                         # contribute to CZsum are negative.
    $v5 = $SATvariable{"CT($character,$num)"};
    $CNFstring .= "$CZstring -$v5 0 \n";  # DG June 4 when leaf num is in the cluster tree for char, then since it can't be the root of that cluster
                                              # tree, there must be an edge into it from
                                              # a node baltnum that is in the cluster tree for char. 
    $clausecount++;

 } # end foreach num that labels a leaf
} # end foreach character 

print OUT "binary\n";
print OUT "$binaries";
print OUT "end\n";

print SAT_KEY "variable - key correspondence \n"; #DG this block is placed after all the variable and clause generating code.
%RevSATvariable = ();
foreach $key (keys %SATvariable) {
        $SATint = $SATvariable{$key};
        $RevSATvariable{"$SATint"} = $key;
}


foreach $SATint (sort {$a <=> $b} keys %RevSATvariable) {
        print SAT_KEY "$RevSATvariable{$SATint} $SATint \n";
}

$v = $SATvariable{"RT(1,0)"}; # initialize the RT count. Node 2 is the first node that could be a reticulation node, so the total number of reticulation
                              # nodes possible before node 2 is examined is zero.

$v1 = $SATvariable{"RT(2,0)"}; # initialize the RT count. Node 2 is the first node that could be a reticulation node, so the total number of reticulation
                              # nodes possible before node 2 is examined is zero.
$CNFstring .= "$v 0 \n";
$CNFstring .= "$v1 0 \n";
$clausecount = $clausecount + 2;

foreach $num (2 .. $n + $interior_num) { # create the CNF variables to count the number of reticulation nodes.
     if ($num - 2 < $target) {
        $bt = $num - 2;
     }
     else {
         $bt = $target + 1;   # DG June 3, 2020 if num is large enough that the number of reticulations observed before 
                              # node num is examined
                              # is as larger as target, then if a reticulation is observed at node num, the total number of
                              # reticulations will go to target+1, and that information will be propogated to the dummy node
                              # after the leaves, where target+1 reticulations is prohibited, so that will make the CNF formula
                              # unsatisfiable. This is how we restrict the total number of reticulations. 
     }
        
     foreach $sum (0 .. $bt) {
        $nump1 = $num + 1;      # RT(num,sum) => RT(num+1,sum)
        $v1 = $SATvariable{"RT($num,$sum)"}; 
        $v2 = $SATvariable{"RT($nump1,$sum)"}; 
        $CNFstring .= "$v2 -$v1 0 \n";
        $clausecount++;

        if ($sum <= $target) {  # DG June 3, 2020
             $sump1 = $sum + 1;
             $v3 = $SATvariable{"R($num)"}; 
             $v4 = $SATvariable{"RT($nump1,$sump1)"}; 
         
             $CNFstring .= "$v4 -$v1 -$v3 0 \n";
             $clausecount++;
        }
     } # end foreach sum
} # end foreach num


$varcount = $SATinteger - 1;
#print "The total variable and clause counts are $varcount, $clausecount \n";

print SAT_FINAL "p cnf $varcount $clausecount \n";
print SAT_FINAL "$CNFstring";



################
sub splitwell{

@NA = ();
$S = $_[0];
@SA = split (//,$S);
 print "Inside splitwell @SA \n";
$len = @SA;
 print "Inside splitwell length $len \n";

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
"About to exit splitwell, the retuned array is @NA \n";
return @NA;
}
