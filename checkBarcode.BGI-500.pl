#!/usr/bin/env perl


=head1

        Author: Hao Yu (yuhao@genomics.cn)
        Date:   2017-04-25 (v1.0)

=cut


use strict;
use warnings;


die "USAGE: perl $0 <BARCODE_file1> <BARCODE_file2>\n" if @ARGV < 2;

print "#1.ID1\t#2.BARCODE1\t#3.ID2\t#4.BARCODE2\t#5.STATE\n";

my %index1;
open IN, '<', "$ARGV[0]" or die "$!\n";
while (<IN>)
{
        chomp;

        my @line = split m/\s+/;
        $index1{$line[0]}{'seq'} = $line[1];
}
close IN;

my %index2;
open IN, '<', "$ARGV[1]" or die "$!\n";
while (<IN>)
{
        chomp;

        my @line = split m/\s+/;
        $index2{$line[0]}{'seq'} = $line[1];
}
close IN;

foreach (sort {$a <=> $b} keys %index1)
{
        my $id1 = $_;

         my $BC   = $index1{$id1}{'seq'};
         my $rBC  = reverse $BC;
        (my $cBC  = $BC) =~ tr/ATCG/TAGC/;
         my $rcBC = reverse $cBC;

        foreach (sort {$a <=> $b} keys %index2)
        {
                my $id2 = $_;

                if    ($BC   eq $index2{$id2}{'seq'})
                {
                        printf "%d\t%s\t%d\t%s\tSTD\n",     $id1, $index1{$id1}{'seq'}, $id2, $index2{$id2}{'seq'};
                }
                elsif ($rBC  eq $index2{$id2}{'seq'})
                {
                        printf "%d\t%s\t%d\t%s\tREV\n",     $id1, $index1{$id1}{'seq'}, $id2, $index2{$id2}{'seq'};
                }
                elsif ($cBC  eq $index2{$id2}{'seq'})
                {
                        printf "%d\t%s\t%d\t%s\tCOM\n",     $id1, $index1{$id1}{'seq'}, $id2, $index2{$id2}{'seq'};
                }
                elsif ($rcBC eq $index2{$id2}{'seq'})
                {
                        printf "%d\t%s\t%d\t%s\tREV-COM\n", $id1, $index1{$id1}{'seq'}, $id2, $index2{$id2}{'seq'};
                }
        }
}
