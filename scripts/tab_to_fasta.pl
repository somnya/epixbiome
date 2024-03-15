# inspired in "https://github.com/gungorbudak/ped2fasta/blob/master/ped2fasta.pl"
# use: perl tab_to_fa path_to_tsv

use strict;
my $file = $ARGV[0]; # Path to file
my @columns;
my $sequence;
open(my $inf, "<", $file) or die $!;
open(my $ouf, ">", $file.".fa") or die $!;
while (my $row = <$inf>) {
	# Each row belongs to one individual
	chomp $row;
	@columns = split(" ", $row); # Splits columns of each line into an array
	print $ouf ">", join(" ", @columns[0]), "\n"; # Prints first column of each line ah ID (starts at 0)
	$sequence = join("", @columns[1..$#columns]); # Joins all DMRs
	print $ouf $sequence, "\n"; # Prints the sequence
}
close $inf;
close $ouf;