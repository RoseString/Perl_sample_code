#!/usr/bin/perl
#by: Dan Sun

use strict;
use warnings;

#declare and initialize some variables
my @P = ();
my %protein;
my $USAGE = "USAGE:  $0  file_name  protein_1  protein_2\n";

unless ($ARGV[0] && $ARGV[1] && $ARGV[2]){
    
    die "Wrong!\n$USAGE\n"; #if arguments are less than 3, report error and exit
}

#sub Create_Adjacency_Matrix
my @A = &Create_Adjacency_Matrix();

#sub All_Pairs_Shortest_Paths
&All_Pairs_Shortest_Paths(\@A, \@P);

#sub Get_Shortest_Path
my @shortest_path = &Get_Shortest_Path(@P);

#print the all pairs shortest paths matrix
print "\nThe all pairs shortest paths matrix is:\n\n";

for (my $i = 0; $i < scalar @A; ++$i){
    for (my $j = 0; $j < scalar @A; ++$j){
        print $A[$i][$j], "\t";
    }
    print "\n";
}

#print the shortest path between two chosen proteins
print "\nThe shortest path passing from $ARGV[1] to $ARGV[2] is: @shortest_path\n\n";

exit;

#############################################################################################

sub Create_Adjacency_Matrix {
    
   
    my $i = 0;
    my @A = ();
    
    unless (open(INPUT, $ARGV[0])){
        
        die "Cannot open the input file!\n";
    }

    while(my $line = <INPUT>){
        
        chomp $line;
        my @cols = split(/\s+/, $line);
        #map the protein ids into indexes
        if (! exists $protein{$cols[0]}){
            $protein{$cols[0]} = $i;
            ++$i;
        }
        if (! exists $protein{$cols[2]}){
            $protein{$cols[2]} = $i;
            ++$i;
        }
        #if there is an edge between two proteins, the value of two dimensional array is 1
        $A[$protein{$cols[0]}][$protein{$cols[2]}] = 1;
        $A[$protein{$cols[2]}][$protein{$cols[0]}] = 1;
    }
    
    for (my $m = 0; $m < values %protein; ++$m){
        $A[$m][$m] = 0;#if the same protein, give it a value 0
        for (my $n = 0; $n < values %protein; ++$n){
            $A[$m][$n] = 100 if not defined $A[$m][$n];#if there is no edge between two proteins, the value of two dimensional array is 100(very big)
        }
    }
    
    close INPUT;
    
    return @A;
    
}

#############################################################################################

sub All_Pairs_Shortest_Paths{
    
    my ($A, $P) = @_;
    #initialize every element of array P
    for (my $m = 0; $m < scalar @$A; ++$m){
        for (my $n = 0; $n < scalar @$A; ++$n){
            if($$A[$m][$n] == 0 || $$A[$m][$n] == 100){
                $$P[$m][$n] = "N";#"N" means no predecessor
            }else{
                $$P[$m][$n] = $m;
            }
        }
    }
    for (my $k = 0; $k < scalar @$A; ++$k){
        for (my $i = 0; $i < scalar @$A; ++$i){
            for (my $j = 0; $j < scalar @$A; ++$j){
                if ($$A[$i][$k] + $$A[$k][$j] < $$A[$i][$j]){
                    $$P[$i][$j] = $$P[$k][$j];#update the value of array P if it satisfies the condition
                    $$A[$i][$j] = $$A[$i][$k] + $$A[$k][$j];#update the value of array A if it satisfies the condition
                }
            }
        }
    }
}

#############################################################################################

sub Get_Shortest_Path{

    my (@P) = @_;
    my @shortest_path_index = ();
    my $i = $protein{$ARGV[1]};
    my $j = $protein{$ARGV[2]};
    die "The two proteins are the same!\n" if $i == $j;#if the input protein names are the same, report error and exit
    
    #trace the shortest path and unshift every node to the front of an array
    unshift(@shortest_path_index, $j);
    do{
        $j = $P[$i][$j];
        unshift(@shortest_path_index, $j);
    }until($j == $i);
    
    my %rev_protein = reverse %protein;#reverse the keys and values of hash %protein
    
    foreach(@shortest_path_index){
        push(@shortest_path, $rev_protein{$_});#push every value(protein) of the reversed hash to an array
    }
    
    return @shortest_path;
}