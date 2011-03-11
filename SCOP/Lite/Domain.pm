package SCOP::Lite::Domain;

use strict;
use warnings;
use Carp;

my @attr = qw(sid pdbid defn sccs sunid class fold superfamily family protein species domain);

use Object::Tiny qw(sid pdbid defn sccs sunid class fold superfamily family protein species domain);

sub new {
    my ( $class, $line ) = @_;

    if (
        $line =~ /
        (\w[s\d]\w{3}[\w\.]\w)
        \s
        ([s\d]\w{3})
        \t
        ((?:\w:(?:-?\d+\w?-\d+\w?)?,?)+|-)
        \t
        (\w\.\d+\.\d+\.\d+)
        \t
        (\d+)
        \t
        cl=(\d+),cf=(\d+),sf=(\d+),fa=(\d+),dm=(\d+),sp=(\d+),px=(\d+)
        /ox
      )
    {
        my ( $sid, $pdbid, $defn, $sccs, $sunid, $cl, $fo, $sf, $fa ) =
          ( $1, $2, $3, $4, $5, $6, $7, $8, $9 );
        my $self = $class->SUPER::new();
        @{$self}{@attr} = ( $1, $2, $3, $4, $5, $6, $7, $8, $9, $10,$11,$12 );
        return $self;
    }
    else {
        croak("Bad domain data line: $line");
    }
}

sub domid {
    my $self = shift;
    return $self->{sid};
}

sub parse_coordinate_definition {
    my $self = shift;
    my $defn = $self->defn;

    my @defn;

    # Remove the PDB code if it is present. Note that non-PDB entries in
    # SCOP are assigned codes beginning with 's'.
    $defn =~ s/^[\ds]\w{3} //;

    if ( $defn ne q{-} ) {
        my @chains = split /,/, $defn;

        for my $chn (@chains) {
            my $cdef;
            if ( $chn =~ /^(\w):$/ ) {
                $cdef = [$1];
            }
            elsif ( $chn =~ /^(?:(\w):)?(-?\d+)([A-Z]?)-(\d+)([A-Z]?)$/o ) {
                my $c;
                if ($1) {
                    $c = $1;
                }
                else {
                    $c = q{ };
                }

                # Insertion codes.
                my $ins1 = $3 || q{ };
                my $ins2 = $5 || q{ };

                $cdef = [ $c, $2, $ins1, $4, $ins2 ];
            }
            else {
                confess "Cannot parse chain description $chn\n";
            }
            push @defn, $cdef;
        }
    }
    return @defn;
}

1;

=pod 

=head1 NAME

SCOP::Lite::Domain

=head1 SYNOPSIS

Lightweight object for representing SCOP domains. Created by reading the SCOP
classification (.cla) file.

=head1 CONSTRUCTOR

$domain = SCOP::Lite->new($line)

where $line is a line from the SCOP file. The constructor will croak if the 
line is wrongly formatted.

=head1 METHODS

Methods for retrieving properties of the domain. Read the SCOP documentation
for details of what each property is.

sunid is the SCOP unique identifer (the integer identifier assigned to
every object in SCOP).

=head2 class

sunid of the domain's SCOP class.

=head2 defn

The domain coordinate definition (e.g. 'A:', 'X:1-100').

=head2 family

sunid of the domain's SCOP family.

=head2 fold 

sunid of the domain's SCOP fold.

=head2 parse_coordinate_definition 

 Function: Parses the domain definition into a list of the segments that 
 comprise the domain.

  Return values:

If the domain definition is '-', i.e. the domain includes the entire
PDB entry, the value returned is an empty list.

Otherwise it returns a list of one or more elements, each of which is
itself a list, corresponding to a segment of the domain (most domains
contain only a single segment).

For each segment:

- if the segment is a chain (e.g. 'A:'), the corresponding list
will consist of a single element, the chain letter.

- if the segment is a section of a chain, the corresponding list will
contain the chain letter, the sequence number and insert code of the
first residue, and the sequence number and insert code of the last
residue in the segment. Blank chain letters and insert codes are
represented by ' ' [space].

 Examples:

Descriptor             Return value
-                      ()
A:                     (['A'])
A:4-153                (['A',4,' ',153,' '])
3-219                  ([' ',3,' ',219,' '])
H:114-223B             (['H',114,' ',223,'B'])
C:468-586,C:704-871    (['C',468,' ',586,' '],['C',704,' ',871,' '])

=head2 pdbid

PDB code of the domain's PDB entry.

=head2 sccs

SCOP classification string 

=head2 sid

The SCOP identifier (e.g. d2rhe__).

=head2 sunid

SCOP unique identifier 

=head2 superfamily

sunid of the domain's SCOP superfamily.

=head1 DEPENDENCIES

Object::Tiny, which is available from CPAN.

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT AND LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut
