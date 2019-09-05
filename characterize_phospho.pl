#!/usr/bin/perl

use strict;
use warnings;

#initialize some variables, arrays and hashes
my $USAGE = "USAGE:  $0  protein_name  pdb_inputfile\n";
my $num = 0;
my @aas = ();
my %phospho_seq;
my %phospho_spa;

unless ($ARGV[0] && $ARGV[1]){
    die "Wrong!\n$USAGE"; #if any argument is not given, exit the program and report error
}

unless (open(FILE, $ARGV[1])){
    die "Cannot open the input file!\n"; #if the input files cannot be opened, exit the program and report error
}

while (my $line = <FILE>){
    chomp $line;
    my @cols = split(/\s+/, $line);
    if ($cols[0] eq "SEQRES"){ #the "if" here could help us get protein sequences 
        for (my $i = 0; $i < 4; $i++){ #this for loop is to shift the first four coloums to get pure protein sequences
            shift @cols;
        }
        foreach (@cols){ #push amino acids to an array
            push (@aas, $_);
        }
    }
    
    if ($cols[0] eq "MODRES"){ #the "if" here could help us get the index of the first appeared phosphorylation record
        $num = $cols[4];
        last;
    }
}
close FILE;

%phospho_seq = &Find_Phospho_seq; #go to the corresponding subroutine and get sequential neighborhoods
foreach (sort {$a <=> $b}(keys (%phospho_seq))){ #print sequential neighborhoods
    
    print "These are the sequential neighborhoods for MODEL $_ of protein $ARGV[0]:@{$phospho_seq{$_}}\n";
}

foreach my $key (sort {$a <=> $b}(keys (%phospho_seq))){ 
    my %Counters_phospho_seq = &Frequencies(@{$phospho_seq{$key}}); #go to the corresponding subroutine and get frequencies of each amino acid of sequential neighborhoods
    
    foreach (keys(%Counters_phospho_seq)){ #print frequencies
        print "Model $key: The frequency of amino acid $_ of sequential neighborhoods for protein $ARGV[0] is $Counters_phospho_seq{$_}.\n";
    }
}

%phospho_spa = &Find_Phospho_spa;#go to the corresponding subroutine and get spatial neighborhoods
foreach (sort {$a <=> $b}(keys (%phospho_spa))){ #print spatial neighborhoods
    print "These are the spatial neighborhoods for MODEL $_ of protein $ARGV[0]:@{$phospho_spa{$_}}\n";
}

foreach my $key (sort {$a <=> $b}(keys (%phospho_seq))){
    my %Counters_phospho_spa = &Frequencies(@{$phospho_spa{$key}}); #go to the corresponding subroutine and get frequencies of each amino acid of spatial neighborhoods
    
    foreach (keys(%Counters_phospho_spa)){ #print frequencies
        print "Model $key: The frequency of amino acid $_ of spatial neighborhoods for protein $ARGV[0] is $Counters_phospho_spa{$_}.\n";
    }
}

my %Counters_complete_seq = &Frequencies(@aas); #go to the corresponding subroutine and get frequencies of each amino acid of the complete protein sequence
foreach (keys(%Counters_complete_seq)){ #print frequencies
    print "The frequency of amino acid $_ for complete protein $ARGV[0] sequence is $Counters_complete_seq{$_}.\n";
}

exit; #exit the program
####################################################################################

sub Find_Phospho_seq { #this subroutine is to find sequential neighborhoods 
    
    my $i = 1;
    
    open (FILE, $ARGV[1]);
    while (my $line = <FILE>){
        chomp $line;
        my @cols = split(/\s+/, $line);
        if ($cols[0] eq "MODEL"){ 
            $i = $cols[1];
        }
        if (defined $cols[5]){
            for (my $j = -7; $j <= 7; $j++){ #this for loop is to get 7 amino acids preceding and following the phosporylated amino acid(excluding itself), respectively
                if ( ($cols[0] eq "HETATM" || $cols[0] eq "ATOM") && $cols[2] eq "CA" && $cols[5] == $num + $j ){
                    next if $cols[5] == $num; #to exclude the phosporylated amino acid
                    push (@{$phospho_seq{$i}}, $cols[3]); #keys of the hash are different models
                    next;
                }else{
                    next;
                }
            }
        }
    }    
    return %phospho_seq;
    close FILE;
}

####################################################################################

sub Find_Phospho_spa { #this subroutine is to find spatial neighborhoods 
    
    my (@x0, @y0, @z0) = ();
    my ($x, $y, $z, $distance) = 0;
    my $i = 1;
    
    open (FILE, $ARGV[1]);
    while (my $line = <FILE>){
        chomp $line;
        my @cols = split(/\s+/, $line);
        if ($cols[0] eq "MODEL"){
            $i = $cols[1];
        }
        if (defined $cols[5]){
            if ($cols[0] eq "HETATM" && $cols[2] eq "CA" && $cols[5] == $num){ #get x0 y0 z0 of different models
                push (@x0, $cols[6]);
                push (@y0, $cols[7]);
                push (@z0, $cols[8]);
            }
        }
    }
    close FILE;
    
    open (FILE, $ARGV[1]);
    
    while (my $line = <FILE>){
        chomp $line;
        my @cols = split(/\s+/, $line);
        if ($cols[0] eq "MODEL"){
            $i = $cols[1];
        }
        if (defined $cols[5]){
            next if $cols[0] eq "HETATM" && $cols[5] == $num; #to exclude the phosporylated amino acid
            if (($cols[0] eq "HETATM" || $cols[0] eq "ATOM") && $cols[2] eq "CA"){ #calculate distances between each amino acid and the phosporylated one
                $x = $cols[6];
                $y = $cols[7];
                $z = $cols[8];
                $distance = (($x - $x0[$i-1]) ** 2 + ($y - $y0[$i-1]) ** 2 + ($z - $z0[$i-1]) ** 2) ** 0.5;
                
                push (@{$phospho_spa{$i}},$cols[3]) if $distance <=10; #if the distance is not greater than 10, push the result to the hash; keys of the hash are different models
            }
        }
        
    }
    close FILE;
    return %phospho_spa;
}

####################################################################################

sub Frequencies{ #this subroutine is to get frequencies of each amino acid
    my @aas = @_;
    my @AAs = ("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"); #store 20 amino acids names in an array
    my %Counters;
    
    foreach (@AAs){
        $Counters{$_} = 0; #intialize counter for each amino acid
    }
    
    foreach my $aa (@aas){ #this foreach loop is to replace modified amino acids with their original names
        $aa =~ s/SEP/SER/g if $aa eq "SEP";
        $aa =~ s/TPO/THR/g if $aa eq "TPO";
        $aa =~ s/MSE/MET/g if $aa eq "MSE";
        $aa =~ s/CSP/CYS/g if $aa eq "CSP";
        
        foreach my $AA (@AAs){            
            if ($aa eq $AA){ #count every amino acid
                $Counters{$AA}++;
                last;
            }
            
        }
    }
    return %Counters;
}