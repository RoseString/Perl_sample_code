#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw[min max];

# touch an output file
system("touch output.txt");

# initialize variables & hashes
my $pct = $ARGV[2];
my $setpoint = 0;
my $join = 0;
my $fixed_chr = "";
my $output = "output.txt";
my $USAGE = "Usage: $0  InputFile1  InputFile2  Percentage% [-J if you want to join]\n";
my (%entries, %starts, %stops);

# if any argument is not given, report error and exit the program
unless ($ARGV[0] && $ARGV[1] && $ARGV[2]){
    print $USAGE;
    exit;
}

# if input files cannot be opened, report error and exit the program
unless (open(FILE1, $ARGV[0]) && open(FILE2, $ARGV[1])){
    print "Cannot open inputfiles! Please check!\n";
    exit;
}

# if the user choose to join two strands, give $join a value 1
if (defined $ARGV[3]){
    if ($ARGV[3] eq ("-J") || $ARGV[3] eq ("-j")){$join = 1};
}

# seperately put every line, start & stop into different hashes
while (my $line = <FILE2>){
    chomp $line;
    my @col = split (/\s+/, $line);
    my $chr = $col[0];
    # if the first file is not sorted, report error and exit
    if ($chr lt $fixed_chr){
        print "Input files need to be sorted!\n";
        exit;
    }else{
        $fixed_chr = $chr;
    }
    push (@{$entries{$chr}}, $line);
    push (@{$starts{$chr}}, $col[1]);
    push (@{$stops{$chr}}, $col[2]);
}

close FILE2;
print "Writing results to output.txt...\n";
open (OUTPUT, ">$output");

# this is a big loop for the first strand
$fixed_chr = "";
while (my $line = <FILE1>){
    chomp $line;
    my @col = split (/\s+/, $line);
    my $chr = $col[0];
    
    # give setpoint a value zero everytime chr is not the same; if the second file is not sorted, report error and exit
    if ($chr ne $fixed_chr){
        if ($chr lt $fixed_chr){
            print "Input files need to be sorted!\n";
            exit;
        }else{
        $setpoint = 0;
        $fixed_chr = $chr;
        }
    }
    # if stop of the first strand is smaller than start of the second strand, do the next line
    if ($col[2] < $starts{$chr}[$setpoint]) {next}
    # if TE is 100% contained in the intron, it is not necessary to calculate percentage which may waste some time
    if ($col[1] >= $starts{$chr}[$setpoint] && $col[2] <= $stops{$chr}[$setpoint]){
       if ($join == 0) {print OUTPUT "$line\n"}
       if ($join == 1) {print OUTPUT "$line\t$entries{$chr}[$setpoint]\n"}
       next;
    }
    
    # this is a loop for the second strand
    for (my $i = $setpoint; $i < scalar(@{$entries{$chr}})-1; $i++){
        # still, if TE is 100% contained in the intron, it is not necessary to calculate percentage which may waste some time
        if ($col[1] >= $starts{$chr}[$setpoint] && $col[2] <= $stops{$chr}[$setpoint]){
            if ($join == 0) {print OUTPUT "$line\n"}
            if ($join == 1) {print OUTPUT "$line\t$entries{$chr}[$setpoint]\n"}
            last;
        }
        
        # if there is overlapped region between two strands, calculate percentage
        if (min($stops{$chr}[$setpoint], $col[2]) - max($starts{$chr}[$setpoint], $col[1]) >= 0){
            my $length = $col[2] - $col[1];
            my $coverage = min($stops{$chr}[$setpoint], $col[2])- max($starts{$chr}[$setpoint], $col[1]);
            my $percent = int( $coverage / $length * 100);
            # if the percentage is greater than the argument, then print and exit the second loop
            if ($percent >= $pct && $join == 0){
                print OUTPUT "$line\n";
            }elsif(($percent >= $pct && $join == 1)){
                print OUTPUT "$line\t$entries{$chr}[$setpoint]\n";
            }
            last;
        # if there is no overlapped region, and stop of the first strand is smaller than start of the second strand, exit the second loop
        }elsif($col[2] < $starts{$chr}[$setpoint]){
            last;
        # if there is no overlapped region, but stop of the first strand is greater than start of the second strand, do the next line of Intron.bed and add 1 to setpoint
        }else{
            $setpoint++;
            next;
        }
    }
}
    
close FILE1;
close OUTPUT;

print "Finished!\nResults have been written to output.txt!\n";
exit;