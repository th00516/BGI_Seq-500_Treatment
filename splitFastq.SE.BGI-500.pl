#!/usr/bin/env perl


=head1

        Author: Hao Yu (yuhao@genomics.cn)
        Date:   2018-01-09 (v0.1)

=cut


use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/basename/;
use IO::File;


my ($idx, $lane, $phred, $outdir, $help);
GetOptions 'i:s' => \$idx, 'l:s' => \$lane, 'p:i' => \$phred, 'o:s' => \$outdir, 'h' => \$help;
die "ERROR: Incorrect barcode file or incorrect lane input\nUSAGE: $0 <-i id-barcode_file> -l <lane> [-p 33|64] [-o output_dir] [-h]\n"
if not defined $idx or not defined $lane;

if (defined $help) {print "USAGE: $0 <-i id-barcode_file> -l <lane> [-o output_dir] [-h]\n";}

$phred ||= 33; die "-p means Phred+\n" if $phred != 33 and $phred != 64;

$outdir ||= './';
my $prefix = basename $lane;
mkdir "$outdir/$prefix.cut";

my $bc_in;
$bc_in = IO::File -> new ("< $idx") or die "$!\n";

my (%box, @fq1_out, @sta_out);
while (<$bc_in>)
{
        my ($id, $bc) = (split /\s+/)[0,1];
        for (my $i = 0; $i <= 9; $i ++)
        {
                my @line = split //, $bc;

                $line[$i] = 'A'; my $bc1 = join '', @line;
                $line[$i] = 'T'; my $bc2 = join '', @line;
                $line[$i] = 'G'; my $bc3 = join '', @line;
                $line[$i] = 'C'; my $bc4 = join '', @line;
                $line[$i] = 'N'; my $bc5 = join '', @line;

                $box{$bc1} = $id if not defined $box{$bc1};
                $box{$bc2} = $id if not defined $box{$bc2};
                $box{$bc3} = $id if not defined $box{$bc3};
                $box{$bc4} = $id if not defined $box{$bc4};
                $box{$bc5} = $id if not defined $box{$bc5};
        }

        $fq1_out[$id - 1] = IO::File -> new ("| gzip -c > $outdir/$prefix.cut/$prefix\_$id.fq.gz") if not defined $fq1_out[$id - 1];
        $sta_out[$id - 1] = IO::File -> new ("          > $outdir/$prefix.cut/$prefix\_$id.stat")  if not defined $sta_out[$id - 1];
}

close $bc_in;

my $fq1_in;
if (-e "$lane.fq.gz")
{
        $fq1_in = IO::File -> new ("gzip -cd $lane.fq.gz |");
}
else {die "Incorrect Files\n";}

my (@buf1, @stat, %unsure_bc);
while (my $fq1_id = <$fq1_in>)
{
        chomp $fq1_id;

        my $id1 = $fq1_id;

        my $fq1_se = <$fq1_in>; chomp $fq1_se; <$fq1_in>;
        my $fq1_qs = <$fq1_in>; chomp $fq1_qs;

        my $bc_se; $fq1_se =~ m/(\S*)(\S{10})$/; $fq1_se = $1; $bc_se = $2;
        my $bc_qs; $fq1_qs =~ m/(\S*)(\S{10})$/; $fq1_qs = $1; $bc_qs = $2;

        my $bc_qu = 0; map {$bc_qu += ((ord $_) - $phred) / length $bc_se;} split m//, $bc_qs;

        next if $bc_qu < 10;

        if (defined $box{$bc_se})
        {
                push @{$buf1[$box{$bc_se} - 1]}, "$fq1_id $bc_se\n$fq1_se\n+\n$fq1_qs\n";

                $stat[$box{$bc_se} - 1]{$bc_se} ++;

                if (@{$buf1[$box{$bc_se} - 1]} == 400000)
                {
                        $fq1_out[$box{$bc_se} - 1] -> print (join '', @{$buf1[$box{$bc_se} - 1]}); undef $buf1[$box{$bc_se} - 1];
                }
        }
        else {$unsure_bc{$bc_se} ++;}
}

$fq1_in -> close;

map
{
        if (defined $buf1[$box{$_} - 1])
        {
                $fq1_out[$box{$_} - 1] -> print (join '', @{$buf1[$box{$_} - 1]}); undef $buf1[$box{$_} - 1];
        }

        $sta_out[$box{$_} - 1] -> printf ("=> %10s\t%8d\n", $_, $stat[$box{$_} - 1]{$_}) if defined $stat[$box{$_} - 1]{$_};
} keys %box;

map
{
        $fq1_out[$box{$_} - 1] -> close if defined $fq1_out[$box{$_} - 1];
        $sta_out[$box{$_} - 1] -> close if defined $sta_out[$box{$_} - 1];
} keys %box;

my $unsure_bc_out = IO::File -> new ("> $outdir/$prefix.cut/unsure_barcode.lst");
map {$unsure_bc_out -> print ("$_\t$unsure_bc{$_}\n");} sort {$unsure_bc{$b} <=> $unsure_bc{$a}} keys %unsure_bc;
$unsure_bc_out -> close;
