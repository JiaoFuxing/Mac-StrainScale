#!/usr/bin/perl
# get_site.pl — RUST-get-site/main.rs 的逐行等价 Perl 版（仅用于本地无 Rust 环境下测试）
# 用法: perl get_site.pl <input_file> <error_base_number> <min_depth>
# 输出: output/really-ref.fa (含'-'占位), output/ref-coords.txt, output/poc.txt, output/poc2.txt
#       filter-result-with-counts.txt（过滤后至少两种碱基的原始五列记录）
#       poc.txt:比对列号（比对列号，1 起始）
#       poc2.txt:输入文件第一列原始ID
use strict; use warnings;

my ($input, $ebn, $md) = @ARGV;
die "Usage: perl get_site.pl <input_file> <error_base_number> <min_depth>\n" unless defined $md;

mkdir "output" unless -d "output";
open my $in, '<', $input or die "cannot open $input: $!\n";
open my $ref, '>', "output/really-ref.fa" or die;
open my $poc, '>', "output/poc.txt" or die;
open my $poc2, '>', "output/poc2.txt" or die;
open my $crd, '>', "output/ref-coords.txt" or die;
open my $fil, '>', "filter-result-with-counts.txt" or die "cannot create filter-result-with-counts.txt: $!\n";

print $ref ">really\n";
my $seq = '';
my ($poc_cnt, $kept) = (0, 0);

while (my $line = <$in>) {
    chomp $line;
    next if $line eq '';
    my @cols = split /\t/, $line;
    next if @cols < 5;
    my ($col0_id, $bc, $err, $bases, $depth) = @cols[0..4]; # 取出第1列 $col0_id
    my $col = length($seq) + 1;   # 比对列号（1 起始）

    if (substr($bases, 0, 1) eq '*' || substr($bc, 0, 1) eq '*') {
        $seq .= '-';
        next;
    }

    my @filtered;
    if ($err <= $ebn) {
        for my $part (split /,/, $bc) {
            next unless length $part;
            if ($part =~ /^(.):(-?\d+)/) {
                push @filtered, $1 if $2 > $ebn;
            }
        }
    } else {
        @filtered = split //, substr($bases, 0, 16);
    }

    if ($depth > $md && @filtered == 1) {
        $seq .= $filtered[0];
        print $crd "$col\n";
        $kept++;
    } else {
        $seq .= '-';
    }

    if (@filtered >= 2) {
        print $poc "$col\n";
        print $poc2 "$col0_id\n";
        print $fil "$line\n";
        $poc_cnt++;
    }
}
close $in;
print $ref "$seq\n";
close $ref; close $poc; close $poc2; close $crd; close $fil;

printf STDERR "Done!\n  output/really-ref.fa: %d bases (alignment columns)\n  kept (non-'-'): %d  poc sites: %d\n",
    length($seq), $kept, $poc_cnt;