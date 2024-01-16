#!/usr/bin/env perl
use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::SeqIO;
## adding simple whitespace trimming just in case
sub trim {my $s = shift; $s =~ s/^\s+|\s+$//g; return $s};

my $out = Bio::SeqIO->new( -format => "Fasta" ); 
# Create database from a single Fasta file or even a whole directory of files
# Creates a persistent index file on disk the first time it is run by parsing the sequence identifiers
# This index makes it possible to randomly access the sequences without being loaded into memory
# BioPython doesn't have anything comparable
my $db = Bio::DB::Fasta->new($ARGV[0]);
open my $buscos, $ARGV[1] || die "could not open busco file $!";
while (<$buscos>) {
  chomp;
  next if /^\s*\#|^\s*$/; # skip comments and blank lines
  my @p = split "\t";
  @p = map {trim($_)} @p; # trim trailing whitespace from the array
  #print join ',', @p;
  next unless $p[1] eq 'Complete';
  my $seq =  $db->seq($p[2],$p[3],$p[4],$p[5]); # check if the index is 0 or 1 based in the busco output!
  $p[9] =~ s/\s/_/g;
  my $id = 'SCR0023_' . $p[9] . '_' . $p[0];
  $seq = Bio::Seq->new(-seq => $seq, -display_id => $id, -description => $p[9]);
  $out->write_seq($seq); ## translation of unspliced genomic sequence is moot
}

__END__