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

use Benchmark qw(:all);
use Test::More;                      # last test to print

use SCOP::Lite::Database;

my $cla_file = 'dir.cla.scop.txt_1.75';
my $t0 = Benchmark->new;
my $scop = SCOP::Lite::Database->new(file => $cla_file);
can_ok($scop, qw(sid release sunid domain_count));
ok($scop->domain_count == 110800,"domain count OK");
ok(keys %{$scop->sid} == $scop->domain_count,'sid hash contains correct no. of domains');
ok(keys %{$scop->sunid} == 0,'sunid hash empty');
ok(@{$scop->domains} == $scop->domain_count,'domains() returns correct no. of domains');
my $t1 = Benchmark->new;
print 'time: ',timestr(timediff($t1,$t0)),"\n";

my @domain_list = qw(d1fltx_ d2rhea_);
$scop = SCOP::Lite::Database->new(file => $cla_file, domain_list => \@domain_list);
ok($scop->domain_count == @domain_list,"domain count == ".@domain_list);
my $t2 =  Benchmark->new;
print 'time: ',timestr(timediff($t2,$t1)),"\n";

done_testing();

