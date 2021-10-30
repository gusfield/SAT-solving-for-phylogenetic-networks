# parse_Newick.pl
# DG corrected Nov. 20, 2020 to include all rows - the prior version omitted the first row in the output.
#
# DG, July 16, 2020
# 
# This program reads in (directed) trees (or just a single one) in Newick format and creates matrices (one for each tree) that specify the splits (clusters) in the trees.

# call as parse_Newick.pl file_with_Newick_formated_tree

open (IN, "$ARGV[0]"); 
open (OUT, ">$ARGV[0]matrix");

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
  %lp = ();
  $rp = ();

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
$minleaf = 0;  # Nov. 20, 2020 we seem never to use this - what is it for?
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
      print "$C[$i][$j]";
   }
   print "\n";
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
       print "$j $t1[$j] \n";
     }
}

push (@start_list, $k);

} # end of while tree = <IN>

print OUT "@start_list\n";
print OUT "$tree_list\n";

#foreach $i (1 .. $maxleaf) {
foreach $i (0 .. $maxleaf) { # DG Nov. 20, 2020
   foreach $j (0 .. $k-1) { # DG July 21
      print OUT "$C[$i][$j]";
   }
   print OUT "\n";
}


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
