package STAMP::Hit;

# $Id: Hit.pm,v 1.10 2007/05/10 15:55:56 www-stamp Exp $

use strict;
use Carp;
use STAMP::Domain;

our @ISA = qw(STAMP::Domain);
our $verbose = 0; 


# The regular expression that the 'data' line for a hit in a STAMP
# scan file matches.

our $hit_line_re = qr/^\#\sSc=\s{1,3}\d{1,2}\.\d{3}
  \s
  RMS=\s{1,3}(\d{1,3}\.\d{2,3})
  \s{1,2}
  len=\s+(\d+)
  \s
  nfit=\s+(\d+)
  \s
  seq_?id=\s+(\d+\.\d{2})
  \s
  sec_?id=\s+(\d+\.\d{2})
  \s
  q_len=\s+(\d+)
	\s
  d_len=\s+(\d+)
  \s
  n_sec=\s+(\d+)
  \s
  (?:n_equiv=?\s+(\d+)\s)?
  fit_pos=\s(\w?\s+-?\d+\s\w?)
  /ix;

=head1 CONSTRUCTOR

=head2 new(@arguments)

The constructor can be called in two ways:

new($score_string,$domain_string)

First argument is the "# Sc=..." line from the STAMP scan file.
Second argument is the domain descriptor (multiline string)

new($string)

Argument is the "# Sc=..." line and the domain descriptor
as a single multiline string.

=cut

sub new {
    my $type = shift;
    my ($string,$domstr);

    if (@_ == 2){
	($string,$domstr) = @_;
    }
    else {
	($string,$domstr) = split "\n",$_[0],2;
    } 
    chomp $string;
 
    # Check that the data line for the hit is kosher.
    if ($string !~ /$hit_line_re/){
	confess("STAMP::Hit::new(): String\n$string\ndoes not match expected pattern:\n$hit_line_re");
    }

    my $self = new STAMP::Domain($domstr);

    $self->{_to_string} = "$string\n$domstr";

    chomp $string;

    my @data;

    if ($string =~ /n_equiv/){
	@data = grep !/\#|=/, split q{ },$string,23;
    }
    else {
	# SORTTRANS output
	@data = grep !/\#|=/, split q{ },$string,21;
    }

    my @attr = grep /=/, split q{ },lc($string);

    @attr = map { chop $_; "_$_"} @attr;

    # This field is absent in files output by sorttrans so the attribute
    # has to be defined explicitly.
    $self->{_n_equiv} = 0 unless defined $self->{_n_equiv};

    @$self{@attr} = @data;

    $self->{_data} = \@data;

    bless $self,$type;
}

sub data {
    @{$_[0]->{_data}};
}

=head2 d_len

=cut

sub d_len {
    $_[0]->{_d_len};
}

=head2 fit_pos

=cut

sub fit_pos {
    $_[0]->{_fit_pos};
}

=head2 len

=cut

sub len {
    $_[0]->{_len};
}

=head2 n_equiv

=cut

sub n_equiv {
    $_[0]->{_n_equiv};
}

=head2 n_fit

=cut

sub n_fit {
    $_[0]->{_nfit};
}

=head2 n_sec

=cut

sub n_sec {
    $_[0]->{_n_sec};
}

sub property {
    my ($self,$name) = @_;
    if ($name eq "score"){
	$self->{_sc};
    }
    elsif ($name eq "n_fit"){
	$self->{_nfit};
    }
    else {
	my $val = $self->{"_$name"};
	if (defined $val){
	    return $val;
	}
	else {
	    confess "No such property $name in STAMP::Hit";
	}
    }
}

=head2 q_len

=cut

sub q_len {
    $_[0]->{_q_len};
}

=head2 rms

=cut

sub rms {
    $_[0]->{_rms};
}

=head2

=cut

sub score {
    $_[0]->{_sc};
}

=head2 sec_id

=cut

sub sec_id {
    $_[0]->{_sec_id};
}

=head2 seq_id

=cut

sub seq_id {
    $_[0]->{_seq_id};
}

1;

=head1 NAME

  STAMP::Hit

=head1 SYNOPSIS

Class representing a STAMP database hit.

=head1 INHERITANCE

  Subclass of STAMP::Domain.

=head1 DEPENDENCIES

STAMP::Domain

=head1 AUTHOR

Tom Walsh

=cut

=head1 PACKAGE VARIABLES

=head2 $hit_line_re

The regular expression that the 'data' line for each hit in a STAMP
output file matches.

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT & LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut
