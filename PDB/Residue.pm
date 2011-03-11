package PDB::Residue;

# $Id: Residue.pm,v 1.12 2005/04/01 16:41:25 tom Exp $

=pod

=head1 NAME
    
Residue::PDB
    
=head1 DESCRIPTION

PDB residue class    

=head1 METHODS

=cut

use strict;
use Carp;
use PDB::Atom;

our @resname = qw (ALA CYS ASP ASX GLU GLX PHE GLY HIS ILE LYS LEU MET
	ASN PRO GLN ARG SER THR VAL TRP TYR CYS UNK);

our %throne = (
	ALA => 'A',
	ASX => 'B',
	CYS => 'C',
	ASP => 'D',
	GLU => 'E',
	PHE => 'F',
	GLY => 'G',
	HIS => 'H',
	ILE => 'I',
	LYS => 'K',
	LEU => 'L',
	MET => 'M',
	ASN => 'N',
	PRO => 'P',
	GLN => 'Q',
	ARG => 'R',
	SER => 'S',
	THR => 'T',
	UNK => ' ',
	VAL => 'V',
	TRP => 'W',
	TYR => 'Y',
	GLX => 'Z'
);

our $amino_codes = join '', values %throne;

our $Disulphide_bond_length = 3.0;

# Bond lengths and angle values for backbone bonds
our %Bond_geometry = (    #"C N CA" => {len => 1.46,the => 117.6},
	"N CA C" => { len => 1.515, the => 112.8 },
	"CA C O" => { len => 1.225, the => 122.6 },

	#"CA C N" => {len => 1.335,the => 113.9}
);

# Bond length and angle tolerances for check_geometry.  Warnings are
# issued if the bond length/angle is not within the specified range of
# the value in %Bond_geometry
our ( $Bond_tolerance, $Angle_tolerance ) = ( 0.05, 5 );

# Names of atoms for calculating chi angles
our @chi_names = qw (N CA CB CG CD CE CZ);

# Additional chi atom definitions for specific sidechains. If a
# key-value pair for one of the atoms in chi_names exists in
# chi_names2, the initial atom name is replaced with the corresponding
# residue-specific name from this list
our %chi_names2 = (
	ASP => { CD => 'OD1' },
	GLU => { CE => 'OE1' },
	PHE => { CD => 'CD1' },
	HIS => { CD => 'ND1' },
	ILE => { CG => 'CG1', CD => 'CD1' },
	LYS => { CZ => 'NZ' },
	LEU => { CD => 'CD1' },
	MET => { CD => 'SD' },
	ASN => { CD => 'ND2' },
	GLN => { CE => 'NE2' },
	ARG => { CE => 'NE' },
	SER => { CG => 'OG' },
	THR => { CG => 'OG1' },
	VAL => { CG => 'CG1' },
	TRP => { CD => 'CD1' },
	TYR => { CD => 'CD1' },
	CYS => { CG => 'SG' }
);

# Number of chi angles in sidechain type
our %chi_angles = (
	ALA => 0,
	CYS => 1,
	ASP => 2,
	GLU => 3,
	PHE => 2,
	GLY => 0,
	HIS => 2,
	ILE => 2,
	LYS => 4,
	LEU => 2,
	MET => 3,
	ASN => 2,
	PRO => 2,
	GLN => 3,
	ARG => 4,
	SER => 1,
	THR => 1,
	VAL => 1,
	TRP => 2,
	TYR => 2
);

################################################################################

sub add_atom {
	my ( $self, $atom ) = @_;
	push @{ $self->{_atoms} }, $atom;
}

################################################################################

=head2 atom()

    Usage     : atom(NUMBER) 
                 - Returns atom by its order in the residue. Atoms are numbered from 0.
                atom(NAME)
                 - Retrieves the atom using its name
                atom(-name => NAME, -altloc => ALTLOC)
                 - Retrieves the atom using its name and alternate location indicator.
                   The default alternate location indicator is ' '.

   Returns    : PDB::Atom object or undef if no atom is found.

=cut

sub atom {
	my ( $self, @args ) = @_;
	my $atom;
	if ( @args == 1 ) {
		if ( $args[0] =~ /^(\d+)$/ ) {    # It's an index
			my $ix = $args[0];
			if ( defined $self->{_atoms}[$ix] ) {
				$atom = $self->{_atoms}[$ix];
			}
		}
		else {
			my $name = $args[0];
			$atom = $self->{_atom_map}{ $name . '- ' };
		}
	}
	else {
		my $name;
		my $alt = ' ';

		my $argstr = join ' ', @args;

		while (@args) {
			my ( $k, $v ) = splice @args, 0, 2;
			if ( $k eq '-name' ) {
				$name = $v;
			}
			elsif ( $k eq '-altloc' ) {
				$alt = $v;
			}
			else {
				confess("unrecognised key $k in atom selection\n");
			}
		}

		confess("No atom name given\n") unless $name;

		$atom = $self->{_atom_map}{ $name . "-" . $alt };
	}
	return $atom;
}

################################################################################

=head2 atoms()

=over 4

    Usage    : atoms($arg1,$arg2)
    Function : Return atoms as a list reference
    Arguments: Two numeric arguments specify array bounds in the atom list.
               'heavy' - return heavy atoms only.
               'backbone' - return N,Ca,C and O atoms only.
               'sidechain' - return sidechain atoms only 

=back

=cut

sub atoms {
	my $self = shift;
	my @atoms;

	if ( @_ == 2 ) {
		if ( $_[0] =~ /^(\d+)$/ && $_[1] =~ /^(\d+)$/o ) {
			@atoms = @{ $self->{_atoms} }[ $_[0] .. $_[1] ];
			splice @_, 0, 2;
		}
	}
	else {
		@atoms = @{ $self->{_atoms} };
	}

	for (@_) {
		if (/backbone/i) {
			@atoms =
				grep { $_->name eq 'N' || $_->name eq 'CA' || $_->name eq 'C' || $_->name eq 'O' && $_->element ne 'H' } @atoms;
			last;
		}
		elsif (/sidechain/i) {
			@atoms =
				grep { $_->name ne 'N' && $_->name ne 'CA' && $_->name ne 'C' && $_->name ne 'O' && $_->element ne 'H' } @atoms;
		}

		@atoms = grep { $_->element ne 'H' } @atoms if /heavy/i;
	}
	\@atoms;
}

################################################################################

=head2 backbone_angles()

    Usage    : backbone_angles()
    Function : Return phi, psi and omega angles 

=cut

sub backbone_angles {
	my $self = shift;

	my @atom;

	my ( $omega, $phi );

	# Skip initial residues and beginnings of chains
	unless ( $self->{_nterm} ) {
		if ( $self->{_chain} eq $self->{_prev}->{_chain} ) {
			push @atom, $self->{_prev}->get_coor('CA');
			push @atom, $self->{_prev}->get_coor('C');
			push @atom, $self->get_coor('N');
			push @atom, $self->get_coor('CA');
			$omega = dihedral(@atom);
			splice @atom, 0;

			push @atom, $self->{_prev}->get_coor('C');
			push @atom, $self->get_coor('N');
			push @atom, $self->get_coor('CA');
			push @atom, $self->get_coor('C');
			$phi = dihedral(@atom);
			splice @atom, 0;
		}
	}
	else {
		( $omega, $phi ) = ( 999, 999 );
	}

	my $psi;
	unless ( $self->{_cterm} ) {
		if ( $self->{_chain} eq $self->{_next}->{_chain} ) {
			push @atom, $self->get_coor('N');
			push @atom, $self->get_coor('CA');
			push @atom, $self->get_coor('C');
			push @atom, $self->{_next}->get_coor('N');
			$psi = dihedral(@atom);
		}
	}
	else {
		$psi = 999;
	}

	return ( $phi, $psi, $omega );
}

################################################################################

sub chain {
	$_[0]->{_chain};
}

################################################################################

sub check_geometry {
	my $self = shift;

	for my $a ( keys %Bond_geometry ) {
		my @atom = split ' ', $a;

		my $len = Atom::distance( $self->get_atom( $atom[1] ), $self->get_atom( $atom[2] ) );

		if ( abs( $len - $Bond_geometry{$a}{len} ) > $Bond_tolerance ) {
			printf "check_geometry: bond length in %s exceeds tolerance\n", $self->to_string;
			printf
				"$atom[1] $atom[2] length: %8.3f expected: %8.3f tolerance: %8.3f\n\n",
				$len, $Bond_geometry{$a}{len}, $Bond_tolerance;
		}

		my $ang = angle( $self->get_coor( $atom[0] ), $self->get_coor( $atom[1] ), $self->get_coor( $atom[2] ) );

		if ( abs( $ang - $Bond_geometry{$a}{the} ) > $Angle_tolerance ) {
			printf "check_geometry: angle value in %s exceeds tolerance\n", $self->to_string;
			printf
				"$atom[0] $atom[1] $atom[2] angle: %8.3f expected: %8.3f tolerance: %8.3f\n\n",
				$ang, $Bond_geometry{$a}{the}, $Angle_tolerance;
		}
	}
}

################################################################################

=head2 disul()

    Usage     : disul($res)
    Function  : returns true if both residues are CYS/CSS and the 
                distance between the SG atoms is < PDB::Disulphide_bond_length.
    Arguments : $res - a Residue object

=cut

sub disul {
	my ( $self, $res ) = @_;

	unless ( $self->type =~ /CYS|CSS/ && $res->type =~ /CYS|CSS/ ) {
		return undef;
	}
	else {
		$self->atom('SG')->in_range( $res->atom('SG'), $Disulphide_bond_length );
	}
}

################################################################################

=head2 get_chi()
 
    Usage    : get_chi()
    Function : Get sidechain chi angles

=cut

sub get_chi {
	my $self = shift;

	my @chi;
	for my $i ( 1 .. $chi_angles{ $self->{_type} } ) {
		my @names = @chi_names[ $i - 1 .. $i + 2 ];

		# Check if we need to replace atom names with ones from the
		# residue-specific list
		if ( exists $chi_names2{ $self->{_type} } ) {

			# If the name has a replacement, change it
			@names = map { $_ = ${ $chi_names2{ $self->{_type} } }{$_} || $_ } @names;
		}

		my @atom = map { $self->get_atom("NAME=>$_")->coor } @names;

		push @chi, dihedral(@atom);
	}

	return @chi;
}

sub index { return $_[0]->{_index} }
sub set_index { $_[0]->{_index} = $_[1] }

################################################################################

=head2 ins()

 Returns the insert code

=cut

sub ins {
	$_[0]->{_ins};
}

################################################################################

=head2 is_amino

    Usage    : is_amino()
    Function : Returns true if the residue is an amino acid

=cut

sub is_amino {
	my $self = shift;

	return index( $amino_codes, $self->seq ) != -1 ? 1 : 0;
}

################################################################################

=head2 is_het
 
    Usage    : is_het()
    Function : Returns true if the residue is a HETATM residue

=cut

sub is_het {
	$_[0]->{_is_het};
}

################################################################################

=head2 is_term()
 
    Usage    : is_term()
    Function : returns 1 if the residue is at either the N- or C-terminus

=cut

sub is_term {
	my $self = shift;
	return $self->{_nterm} || $self->{_cterm};
}

################################################################################

=head2 new()
 
    Usage    : new($record)
    Function : Create a new PDB Residue. If the optional argument 
               begins with 'ATOM' or 'HETATM' it is assumed to be
               a PDB record and the residue number and identity
               are read from it.
    Argument : A PDB atom record

=cut 

sub new {
	my $type = shift;
	my $self = {};

	my $rec = shift;
	if ( $rec =~ /^ATOM/ ) {
		$self->{_is_het} = 0;
	}
	elsif ( $rec =~ /^HETATM/ ) {
		$self->{_is_het} = 1;
	}
	else {
		confess("PDB::Residue::new(): not an ATOM record:\n$rec\n");
	}

	( $self->{_type} = substr( $rec, 17, 3 ) ) =~ s/ //g;
	$self->{_chain} = substr( $rec, 21, 1 );
	( $self->{_resnum} = substr( $_, 22, 4 ) ) =~ s/ //g;
	$self->{_ins} = substr( $rec, 26, 1 );
	( $self->{_resid} = $self->{_resnum} . $self->{_ins} ) =~ s/ //g;
	$self->{_seq}   = $throne{ $self->{_type} } || '-';
	$self->{_atoms} = [];
	$self->{_prev}  = undef;
	$self->{_next}  = undef;

	return bless $self, $type;
}

################################################################################

=head2 next()

Returns the next residue in the PDB file.

=cut

sub next {
	$_[0]->{_next};
}

################################################################################

=head2 number

 Usage : number()
  Function : Return the residue number

=cut

sub number {
	$_[0]->{_resnum};
}

################################################################################

=head2 prev()

Returns the preceding residue in the PDB file.

=cut

sub prev {
	$_[0]->{_prev};
}

################################################################################

=head2 resid

 Usage : resid()
 Function : Return the residue identifier (sequence number + insert code)

=cut

sub resid {
	$_[0]->{_resid};
}

################################################################################

=head2 seq

Usage    : seq
Function : Return the one letter code of the residue

=cut

sub seq {
	$_[0]->{_seq};
}

################################################################################

sub set_next {
	$_[0]->{_next} = $_[1];
}

sub set_prev {
	$_[0]->{_prev} = $_[1];
}

################################################################################

=head2 to_string()
 
    Usage    : to_string()
    Function : return Residue::PDB object as a formatted string

=cut

sub to_string {
	my $str;
	my $self = shift;

	$str = sprintf
		"%3s %1s %1s%5s",
		$self->{_type}, $self->{_seq}, $self->{_chain}, $self->{_resid};

	$str =~ /^\s*(.*?)\s*$/;    # Trim white space
	$&;                         # Return trimmed string
}

################################################################################

=head2 type()

 Returns the residue type.

=cut

sub type {
	$_[0]->{_type};
}

################################################################################

# DESTROY must be explicitly defined because there are circular
# references between adjacent residues, and between each
# residue and its atoms.
sub DESTROY {
	my $self = shift;

	delete $self->{_atoms};
	delete $self->{_atom_map};
	delete $self->{_prev};
	delete $self->{_next};
}

1;

=head1 DEPENDENCIES

PDB::Atom

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT & LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut
