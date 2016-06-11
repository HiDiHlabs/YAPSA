#!/usr/bin/env perl

###############################################################################
#                                                                 
#  Copyright 2015 Matthias Schlesner and Daniel Huebschmann.
#   
#  This file is part of YAPS.
#
#   YAPS is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   any later version.
#
#   YAPS is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this script.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

  use strict;
  use warnings;
  use Getopt::Std;
  $Getopt::Std::STANDARD_HELP_VERSION = 1;
  use v5.10;
  use List::Util qw(sum);
  
  my %opts;
  getopts( 'r:w:h', \%opts);
  if ($opts{h}) {
    print "\nUsage: ./optimizer.pl [Option...] [Value...]\n\n";
    print "\tOptions for input parameters:\n";
    print "\t-r\tfasta_file; file name; default: /ibios/co02/reference/Reference_1KG/hs37d5.fa\n";
    print "\t-w\tword_length; integer; default: 3\n";
    print "\t-h\tdisplay this message\n";
    print "\n";
    exit;
  }

  my $ref = $opts{r} || "/ibios/co02/reference/Reference_1KG/hs37d5.fa";
  my $word_length = $opts{w} || 3;

  my $seq;
  my %kmers;
  my $chrcount;
  open (REF, $ref);
  while(<REF>) {
    if (/^>/) {
      $chrcount++;
      count_kmers(\$seq) if ($chrcount > 1);  # \$seq is a hard reference on $seq
      $seq = '';
      next;
    }
    chomp; # by standard applied to $_; cut all \n
    $seq .= uc($_); # .= is append
  }
  
  
  foreach (sort keys %kmers) {
    say "$_\t$kmers{$_}";
  }
  
  #my $total = sum(values %kmers);
  #say $total;
  
  sub count_kmers {
    my $seq = $_[0]; # argument of the subroutine, call by reference
    my $kmer;
    #my $i;
    #my $length = length($$seq);
    while ($$seq =~ /(?=([^MNR]{$word_length}))/ig) {
      # ignore kmers with 'N' 
      # $$seq is dereferencing the hard reference given to the subroutine
      # (?=...) : equality is only successful, if ... comes afterwards
      # //g : global, i.e. return a list with all hits
      # //i : case insensitive
      $kmer = $1;
      # $1 contains the pattern in the first parenthesis
      #if (substr($kmer, 1,1) =~ /[GT]/) {
      #  $kmer =~ tr/ACGT/TGCA/;
      #  $kmer = reverse $kmer;
      #}
      $kmers{$kmer}++;
    }
  }

