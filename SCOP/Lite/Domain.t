#!/usr/bin/perl
#===============================================================================
#
#         FILE:  Domain.t
#
#  DESCRIPTION:  
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  YOUR NAME (), 
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  17/12/09 15:23:08
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

use Test::More;                      # last test to print

use SCOP::Lite::Domain;

my $cla_file = 'dir.cla.scop.txt_1.75';
open my $CLA, '<', $cla_file or die "$!: $cla_file";
while (<$CLA>){
    next if /^#/;
    my $domain = SCOP::Lite::Domain->new($_);
    can_ok($domain,qw(domid pdbid defn sccs sunid class fold superfamily family));
    ok($domain->domid eq $domain->sid,'domid eq sid');
    while (<$CLA>){
        ok(eval { SCOP::Lite::Domain->new($_) },'read domain '.(split q{ },$_)[0]);
    }
}
done_testing();

