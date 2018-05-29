#!/usr/bin/env perl


=head1

        Author: Hao Yu (yuhao@genomics.cn)
        Date:   2017-05-05 (v0.1)
        Date:   2017-05-07 (v1.0)  // To create a file buffer for improving performance but getting more memory expanditure
        Date:   2017-05-09 (v2.0)  // To allow 1 bp mismatch
        Date:   2017-05-10 (v2.1)
        Date:   2017-05-19 (v2.2)  // To improve performance further
        Date:   2018-05-29 (v2.3)  // To modify the OUTPUT cache size from 200k to 10k
                                      To add a new function of check/debug data without OUTPUT
                                      To change the suffix of OUTPUT dirctory from '.cut' into '.splitted'

=cut


use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/basename/;
use IO::File;


my ($idx, $lane, $phred, $outdir, $check, $debug, $help);
GetOptions 'i:s' => \$idx, 'l:s' => \$lane, 'p:i' => \$phred, 'o:s' => \$outdir, 'c:s' => \$check, 'd:s' => \$debug, 'h' => \$help;
die "ERROR: Incorrect barcode file or incorrect lane input\nUSAGE: $0 -i <id-barcode_file> -l <lane> [-p 33|64] [-o output_dir] [-c|-d] [-h]\n"
if not defined $idx or not defined $lane;

if (defined $help) {print "USAGE: $0 -i <id-barcode_file> -l <lane> [-o output_dir] [-c|-d] [-h]\n";}

$phred ||= 33; die "-p means Phred+\n" if $phred != 33 and $phred != 64;

$outdir ||= './';
my $prefix = basename $lane;
mkdir "$outdir/$prefix.splitted";

my $bc_in;
$bc_in = IO::File -> new ("< $idx") or die "$!\n";

my (%box, @fq1_out, @fq2_out, @sta_out);
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

        unless (defined $check or defined $debug)
        {
                $fq1_out[$id - 1] = IO::File -> new ("| gzip -c > $outdir/$prefix.splitted/$prefix\_$id\_1.fq.gz") if not defined $fq1_out[$id - 1];
                $fq2_out[$id - 1] = IO::File -> new ("| gzip -c > $outdir/$prefix.splitted/$prefix\_$id\_2.fq.gz") if not defined $fq2_out[$id - 1];
        }
        $sta_out[$id - 1] = IO::File -> new ("> $outdir/$prefix.splitted/$prefix\_$id.stat") if not defined $sta_out[$id - 1];
}

close $bc_in;

my ($fq1_in, $fq2_in);
if (-e "$lane\_1.fq.gz" and -e "$lane\_2.fq.gz")
{
        $fq1_in = IO::File -> new ("gzip -cd $lane\_1.fq.gz |");
        $fq2_in = IO::File -> new ("gzip -cd $lane\_2.fq.gz |");
}
else {die "Incorrect Files\n";}

my (@buf1, @buf2, @stat, %unsure_bc);
while (my $fq1_id = <$fq1_in>, my $fq2_id = <$fq2_in>)
{
        chomp $fq1_id; chomp $fq2_id;

        (my $id1 = $fq1_id) =~ s/\/1//; (my $id2 = $fq2_id) =~ s/\/2//;
        next if $id1 ne $id2;

        my $fq1_se = <$fq1_in>; my $fq2_se = <$fq2_in>; chomp $fq1_se; chomp $fq2_se; <$fq1_in>; <$fq2_in>;
        my $fq1_qs = <$fq1_in>; my $fq2_qs = <$fq2_in>; chomp $fq1_qs; chomp $fq2_qs;

        my $bc_se; $fq2_se =~ m/(\S*)(\S{10})$/; $fq2_se = $1; $bc_se = $2;
        my $bc_qs; $fq2_qs =~ m/(\S*)(\S{10})$/; $fq2_qs = $1; $bc_qs = $2;

        my $bc_qu = 0; map {$bc_qu += ((ord $_) - $phred) / length $bc_se;} split m//, $bc_qs;

        next if $bc_qu < 10;

        if (defined $box{$bc_se})
        {
                push @{$buf1[$box{$bc_se} - 1]}, "$fq1_id $bc_se\n$fq1_se\n+\n$fq1_qs\n";
                push @{$buf2[$box{$bc_se} - 1]}, "$fq2_id $bc_se\n$fq2_se\n+\n$fq2_qs\n";

                $stat[$box{$bc_se} - 1]{$bc_se} ++;

                unless (defined $check or defined $debug)
                {
                        if (@{$buf1[$box{$bc_se} - 1]} == 10000 and @{$buf2[$box{$bc_se} - 1]} == 10000)
                        {
                                $fq1_out[$box{$bc_se} - 1] -> print (join '', @{$buf1[$box{$bc_se} - 1]}); undef $buf1[$box{$bc_se} - 1];
                                $fq2_out[$box{$bc_se} - 1] -> print (join '', @{$buf2[$box{$bc_se} - 1]}); undef $buf2[$box{$bc_se} - 1];
                        }
                }
        }
        else {$unsure_bc{$bc_se} ++;}
}

$fq1_in -> close;
$fq2_in -> close;

map
{
        unless (defined $check or defined $debug)
        {
                if (defined $buf1[$box{$_} - 1] and defined $buf2[$box{$_} - 1])
                {
                        $fq1_out[$box{$_} - 1] -> print (join '', @{$buf1[$box{$_} - 1]}); undef $buf1[$box{$_} - 1];
                        $fq2_out[$box{$_} - 1] -> print (join '', @{$buf2[$box{$_} - 1]}); undef $buf2[$box{$_} - 1];
                }
        }

        $sta_out[$box{$_} - 1] -> printf ("=> %10s\t%8d\n", $_, $stat[$box{$_} - 1]{$_}) if defined $stat[$box{$_} - 1]{$_};
} keys %box;

map
{
        unless (defined $check or defined $debug)
        {
                $fq1_out[$box{$_} - 1] -> close if defined $fq1_out[$box{$_} - 1];
                $fq2_out[$box{$_} - 1] -> close if defined $fq2_out[$box{$_} - 1];
        }
        $sta_out[$box{$_} - 1] -> close if defined $sta_out[$box{$_} - 1];
} keys %box;

my $unsure_bc_out = IO::File -> new ("> $outdir/$prefix.splitted/unsure_barcode.lst");
map {$unsure_bc_out -> print ("$_\t$unsure_bc{$_}\n");} sort {$unsure_bc{$b} <=> $unsure_bc{$a}} keys %unsure_bc;
$unsure_bc_out -> close;

