package PDB::Atom;  

# $Id: Atom.pm,v 1.9 2005/04/07 10:22:42 tom Exp $

use strict;

our $record_format = "ATOM  %5d %4s%1s%3.3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s";

sub new {
    my $type = shift;

    my ($res,$rec,$self);
    if (@_ == 2){
	$res = shift;
    }
    $self = {};
    $self->{_res} = $res;
    
    if (@_) {
	my $rec = shift;
	$self->{_ser}     = int(substr($rec,6,5));
	my $atom = substr($rec,12,4);
	
	($self->{_name} = $atom) =~ s/ //g;
	# Store the atom name from the PDB file so that it can be used
	# in the to_record() method.
	$self->{_name_field} = $atom;
	($self->{_element} = substr($rec,12,2))=~ s/ //g;
	$self->{_altloc}  = substr($rec,16,1);
	# Add 0 to the coordinates so that the string representation
	# corresponds to the numeric value stored in each field,
	# rather than the field itself, which may contain leading
	# blanks.
	$self->{_coor}[0] = substr($rec,30,8) + 0;
	$self->{_coor}[1] = substr($rec,38,8) + 0;
	$self->{_coor}[2] = substr($rec,46,8) + 0;
	$self->{_occup}   = substr($rec,54,6) + 0;
	$self->{_bval}    = substr($rec,60,6) + 0;
    }
    return bless $self,$type;
}

sub make {
    my ($type,%arg) = @_;
    my $self = {};
    $self->{_name}       = $arg{-name} or croak("make(): No atom name\n");
    $self->{_name_field} = $self->{_name};
    $self->{_element}    = substr($self->{_name},0,1);
    $self->{_altloc}     = $arg{-altloc} || q{ };
    unless ($arg{-coor}){
	croak("make(): No coordinates defined\n");
    }
    @{$self->{_coor}}    = @{$arg{-coor}};
    $self->{_occup}      = $arg{-occ}  || 1.0;
    $self->{_bval}       = $arg{-bval} || 0.0;
    return bless $self,$type;
}

sub altloc {
    $_[0]->{_altloc};
}

sub bval {
    $_[0]->{_bval};
}

sub element {
    $_[0]->{_element};
}

sub coor {
    $_[0]->{_coor};
}   

sub distance {
    my ($a1,$a2) = @_;

    my @coor = ($a1->coor,$a2->coor);

    my $dist = ($coor[0][0]-$coor[1][0])**2 + ($coor[0][1]-$coor[1][1])**2 +
	($coor[0][2]-$coor[1][2])**2;

    return ($dist > 0) ? sqrt($dist) : 0;
}

sub name {
    $_[0]->{_name};
}

sub occup {
   $_[0]->{_occup};
}

sub res {
    $_[0]->{_res};
}

sub ser {
    $_[0]->{_ser};
}

sub to_record {
    my $self = shift;

    my $res = $self->res;

    sprintf 
	$record_format,$self->{_ser},$self->{_name_field},$self->{_altloc},
	  $res->type,$res->chain,$res->number,$res->ins,
	    @{$self->{_coor}}[0..2], $self->{_occup}, $self->{_bval},
	      $self->res->chain,$self->element,"";
}

sub to_string {
    my $self = shift;

    my $res = $self->res;

    sprintf 
	"%3.3s %1s%4s%1s %-4s%s %4d",
	$res->type,$res->chain,$res->resid,$res->ins,$self->name,$self->altloc,$self->ser;
}

1;	

__END__

=head1 NAME

PDB::Atom

=head1 AUTHOR

Tom Walsh

=head1 SYNOPSIS

Represents atom records from PDB files

=head1 METHODS

=head3 altloc()

Alternate location indicator

=head3 bval()

    Atom B-value

=head3 element()

Returns the element type of the atom.

=head3 coor

    Return atom coordinates (array reference)

=head2 Atom::distance($atom1,$atom2)
   
Function : Return the distance between two atoms

Arguments : PDB::Atom objects

Returns : Distance in Angstroms.

=head3 name()
    
   Atom name

=head3 occup()

    Atom occupancy

=head3 res

    Residue to which the Atom belongs

=head3 ser
    
    Return atom serial number

=head3 to_record

    Return Atom as a string in PDB ATOM record format.

=head3 to_string

    Return a string representation of an Atom. The fields are:

=over

=item Residue name

=item Chain identifier

=item Residue number and insert code

=item Atom name

=item Alternate location indicator 

=item Serial number

=back    

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT & LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut

