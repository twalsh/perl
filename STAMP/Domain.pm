package STAMP::Domain;

# $Id: Domain.pm,v 1.9 2005/02/21 14:00:48 tom Exp $

use strict;
use Carp;

our $verbose = 0;

our $coord_definition_re = 'ALL|(CHAIN (\w))|((\w) (-?\d+) (\w) (?:TO|to) \w (-?\d+) (\w))';
# Use \s to match white space so that the pattern can be used with
# the /x modifier.
$coord_definition_re =~ s/ /\\s/g;

our $descriptor_re = qr/^(\S+) # Location of PDB file
  \s+
  (\S+) # Identifier
  \s+
  ({\s
    (?:(.*\n.*) # Coordinate description + transformation
     |          #  or
     (.*))      # Coordinate description only
    \s
   })/isox;

=head1 NAME

STAMP::Domain

=head1 SYNOPSIS

Implements a class to represent STAMP domain descriptors.

Provides functions for parsing STAMP descriptors.

=head1 AUTHOR

Tom Walsh

=head1 CONSTRUCTOR

=head2 new($string)

  Create a STAMP::Domain object from the specified STAMP domain
descriptor.

=cut

sub new {
    my ($type,$string) = @_;
    my $self = {_to_string => $string};

    print "Parsing string:\n$string\n" if $verbose;

    unless ($string =~ /$descriptor_re/){
	print STDERR "The string:\n---\n$string\n---\n";
	print STDERR "doesn't match the pattern expected for a STAMP domain definition:\n";
	print STDERR "$descriptor_re\n";
	confess;
    }

    @$self{qw(_pdb _id)} = ($1,$2);
    
    my ($coord,$trans);
    if ($4){
	($coord,$trans) = split /\n/,$4,2;
    }
    elsif ($5){
	$coord = $5;
    }
    else {
	print STDERR "No coordinate description in string $string";
	confess;
    }

    $self->{_coord} = $coord;
    $self->{_trans} = $trans;

    # Get list of segments from the coordinate description
    $self->{_segments} = [];
  
    unless ($coord =~ /($coord_definition_re)+/){
	confess "Coordinate descriptor:\n$coord\ndoes not match the expected format:\n$coord_definition_re\n---\n";
    }
    
    while ($coord =~ /($coord_definition_re)/ig){
	push @{$self->{_segments}},$1;
    }
    confess "No segments defined in STAMP domain descriptor $coord" unless @{$self->{_segments}};

    $self->{_matrix} = undef;
    $self->{_vector} = undef;

    if ($trans){
	# Extract tranformation matrix and vector
	$self->{_matrix} = [];
	$self->{_vector} = [];
	
	my @lines = split /\n/,$trans;
	for my $i (0..2){
	    (@{$self->{_matrix}[$i]}[0..2],$self->{_vector}[$i]) = split ' ',$lines[$i];
	}
	
	if ($verbose){
	    print "matrix:\n";
	    for my $i (0..2){
		print join ' ',@{$self->{_matrix}[$i]},"\n";
	    }
	    print "vector:\n";
	    print join ' ',@{$self->{_vector}},"\n";
	}
    }
    else {
	for my $i (0..2){
	    (@{$self->{_matrix}[$i]}[0..2],$self->{_vector}[$i]) = (0,0,0,0);
	    $self->{_matrix}[$i][$i] = 1.0;
	}
    }

    bless $self,$type;
}

=head1 METHODS

=cut

sub coord {
    $_[0]->{_coord};
}

=head2 id

Returns the domain identifier.

=cut

sub id {
    $_[0]->{_id};
}

=head2 matrix

Returns the transformation matrix for the domain. The returned value
is an array reference to a two-dimensional array

=cut

sub matrix {
    $_[0]->{_matrix};
}

=head2 pdb

The PDB file specified in the descriptor.

=cut

sub pdb {
    $_[0]->{_pdb};
}

=head2 segments

Returns a list of the segments in the domain, e.g. if the
descriptor is:

{ CHAIN A B 1 _ TO B 100 _ C 2 _ TO C 100 _ } 

the returned list will be:

["CHAIN A","B 1 _ TO B 100 _","C 2 _ TO C 100 _"]

=cut

sub segments {
    @{$_[0]->{_segments}};
}

=head2 to_string

Returns the domain descriptor as a string.

=cut

sub to_string {
    $_[0]->{_to_string};
}

=head2 transform

Returns the domain transform as a string.

=cut

sub trans {
    $_[0]->{_trans};
}

sub transform { $_[0]->{_trans} }

=head2 vector

Returns the a reference to the translation vector for the domain. 

=cut

sub vector {
    $_[0]->{_vector};
}

=head1 FUNCTIONS

=head2 match($string)

If $string matches the pattern expected for a STAMP domain description,
return a list comprising the matched component of the string and the
following elements of the domain definition:

=over

=item 

the PDB file

=item 

the domain identifier

=item 

the coordinate description if it includes a transformation

=item 

the coordinate description if it does not include a transformation.

=back

Only one of the last two elements on the list will be defined,
depending on whether or not there is a transformation defined.

If the string does not match the expected pattern, the function
returns undef.

=cut

sub match { 
    my $string = shift;

    if ($string =~ /$descriptor_re/){
	my $coor_str = $&;
	# Extract the coordinate description and check it's OK
	my @matches = ($&,$1,$2,$4,$5);
	my $coor_desc = $4 || $5;

        unless ($coor_desc){
	    confess "No coordinate description in string:\n$coor_str\n";
	}
	else {
	    if ($coor_desc =~ /$coord_definition_re/){
		return @matches;
	    }
	    else {
		return undef;
	    }
	}
    }
    else {
	return undef;
    }
}

1;

=head1 PACKAGE VARIABLES

=head2 $coord_definition_re

The regular expression that a correctly-formatted coordinate descriptor 
must match. The captured substrings are:

$1 - a descriptor matching CHAIN \w

$2 - the chain identifier from the descriptor matched by $1

$3 - a descriptor matching 

  ((\w) (-?\d+) (\w) (?:TO|to) \w (-?\d+) (\w)),

  e.g. B 1 _ TO B 100 P

$4 - the first chain identifier in $3

$5 - the first sequence number in $3

$6 - the first insertion code in $3

$7 - the second sequence number in $3

$8 - the second insertion code in $3

If the descriptor is 'ALL', the match string will be $&.

The pattern uses \s to match whitespace. Therefore, it is
safe to use with the /x modifier.

=head2 $descriptor_re

The regular expression that a valid STAMP domain descriptor matches.
The pattern uses \s to match whitespace and is safe to use with the /x
modifier.  The back references defined are:

$1 - Location of the PDB file.

$2 - Domain identifier

$3 - The complete coordinate description and (if defined) transformation.

$4 - The coordinate description and transformation. This is defined
only if the transformation is defined. It is equivalent to $3 with the
enclosing braces removed.

$5 - The coordinate description if the transformation is not defined.

=head2 $verbose

Set to 1 to turn on logging output.

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT & LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut

