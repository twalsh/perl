package PDB::Entry;

use strict;
use warnings;
use Carp;
use Log::Log4perl;
use PDB::Atom;
use PDB::Residue;

my $log = Log::Log4perl->get_logger('PDB::Entry');

# Minimum distance between adjacent Ca atoms for defining a chain break
our $chain_break_threshold = 5.0;


sub new {
    my ( $class, $pdbfile, %args ) = @_;
    my $self = {
        _expdta     => 'UNKNOWN',
        _nchains    => 1,
        _file       => $pdbfile,
        _resolution => undef,
    };

    # Attributes from the header records
    my @attributes =
      qw(_title _expdta _header _compnd _method _seqres _secstruc);
    @{$self}{@attributes} = (q{}) x @attributes;

    my $fh;
    if ( ref($pdbfile) eq 'GLOB' ) {
        $fh = $pdbfile;
    }
    else {
        open $fh, '<', $pdbfile
          or confess("PDB::Entry::new: cannot open $pdbfile\n");
    }

    my $prev_resnum  = undef;
    my $prev_resname = undef;
    local $/ = "\n";
    my @lines = <$fh>;
    my $chain;
    my ( $next, $prev );
    for (@lines) {

        # Stop reading if we've reached the end of the first model
        # in an NMR file.
        last if /^ENDMDL/;

        if (/^HEADER/) {
            $self->{_header} = substr( $_, 10, 56 );
        }
        if (/^TITLE/) {
            $self->{_title} .= substr( $_, 10, 60 );
        }
        elsif (/^EXPDTA/) {
            $self->{_expdta} .= substr( $_, 10, 60 );
        }
        elsif (/^COMPND/) {
            $self->{_compnd} .= substr( $_, 10, 60 ) . "\n";
        }
        elsif (/^REMARK 2\d{2}  EXPERIMENT TYPE                : (\S+)/) {
            $self->{_expdta} = $1;
        }
        elsif (/^REMARK   2 RESOLUTION. (\S+)\s+ANGSTROMS./) {
            $self->{_resolution} = $1;
        }
        elsif (/^(SHEET|TURN|HELIX)/) {
            $self->{_secstruc} .= $_;
        }
        elsif (/^(SEQRES)/) {
            $self->{_seqres} .= $_;
        }
        elsif (/^TER/) {

            # End of a chain
            $prev_resnum  = undef;
            $prev_resname = undef;
        }
        elsif (/^(ATOM|HETATM)/) {
            my $line = $_;

            # Residue sequence number
            my $seqn;
            ( $seqn = substr( $line, 22, 5 ) ) =~ s/ //g;

            # Residue sequence name
            my $resname = substr( $line, 17, 3 );

            if (
                !defined($prev_resnum)    # First residue
                || !defined($chain)       # Start of the first chain
                || $seqn    ne $prev_resnum           # Check residue numbers
                || $chain   ne substr( $_, 21, 1 )    # Check chain identifier
                || $resname ne $prev_resname          # Check residue name
              )
            {
                $prev = $next;
                $next = PDB::Residue->new($_);    # Make object for new residue
                $prev->set_next($next)
                  if $prev;    # Link the new residue to the previous one
                $next->set_prev($prev);    # ... and vice versa
                push @{ $self->{_residues} },
                  $next;                   # Add the new residue to residue list
                $next->set_index( $#{ $self->{_residues} } );

                # Add the residue to the list of residues for the current chain
                if ( !defined($chain) || $chain ne substr( $_, 21, 1 ) ) {
                    $chain = substr( $_, 21, 1 );
                    $self->{_nchains}++;
                }
                push @{ $self->{_chains}{$chain} }, $next;

                # Add the residue to the residue map
                unless (
                    exists $self->{_residue_map}{$chain}
                    { $next->number . '-' . $next->ins } )
                {
                    $self->{_residue_map}{$chain}
                      { $next->number . '-' . $next->ins } = $next;
                }
                else {
                    $log->warn( 'Duplicate residue number for '
                          . $next->to_string
                          . " in $pdbfile" );
                }

                $prev_resnum  = $seqn;
                $prev_resname = $resname;
            }

            # Create a new atom from the current record
            my $atom = PDB::Atom->new( $next, $_ );

# If the -trace argument is used and the atom isn't a CA in an amino acid, then ignore it.
            unless (
                $args{-trace}
                && ( $atom->name ne 'CA'
                    || ( $atom->name eq 'CA' && $next->type eq 'CA' ) )
              )
            {
                push @{ $next->{_atoms} }, $atom;
                unless (
                    exists $next->{_atom_map}
                    { $atom->name . '-' . $atom->altloc } )
                {
                    $next->{_atom_map}{ $atom->name . '-' . $atom->altloc } =
                      $atom;
                }
                else {
                    $log->warn( 'Ignoring duplicate atom location '
                          . $atom->to_string );
                }
            }
        }
    }
    $self->{_method} = $self->{_expdta} || q{};
    for my $attr (@attributes) {
        $self->{$attr} =~ s/\s+$//;
    }

    $self->{_nres} = scalar @{ $self->{_residues} };
    return bless $self, $class;
}

################################################################################


sub chain {
    my ( $self, $chain ) = @_;

    confess("No chain specified\n") unless defined $chain;

    if ( $self->{_chains}{$chain} ) {
        return @{ $self->{_chains}{$chain} };
    }
    return ();
}

################################################################################


sub chains {
    return $_[0]->{_chains};
}

################################################################################


sub compnd {
    return $_[0]->{_compnd};
}

################################################################################

sub expdta {
    return $_[0]->{_expdta};
}

################################################################################

sub file {
    return $_[0]->{_file};
}

################################################################################

sub find_chain_breaks {
    my $self = shift;
    $self->{_residues}[0]->set_nterm(1);
    $self->{_residues}[ $#{ $self->{_residues} } ]->set_cterm(1);

    for my $i ( 0 .. $#{ $self->{_residues} } - 1 ) {
        if (   $self->{_residues}[$i]->is_amino
            && $self->{_residues}[ $i + 1 ]->is_amino )
        {
            my $ca1 = $self->{_residues}[$i]->get_atom('CA');

            my $ca2 = $self->{_residues}[ $i + 1 ]->get_atom('CA');

            if ( !( $ca1->in_range( $ca1, $ca2, $chain_break_threshold ) ) ) {
                $self->{_residues}[$i]->set_cterm(1);
                $self->{_residues}[ $i + 1 ]->set_nterm(1);
            }
            else {
                $self->{_residues}[$i]->set_cterm(0);
                $self->{_residues}[ $i + 1 ]->set_nterm(0);
            }
        }
        else {
            $self->{_residues}[$i]->set_nterm(0);
            $self->{_residues}[$i]->set_cterm(0);
            $self->{_residues}[ $i + 1 ]->set_nterm(0);
            $self->{_residues}[ $i + 1 ]->set_cterm(0);
        }
    }
    return;
}

################################################################################

# Deprecated but kept for backwards compatibility

sub get_residue {
    return residue(@_);
}

################################################################################

sub has_amino_residues {
    my $self = shift;
    for my $res ( $self->residues ) {
        return 1 if $res->is_amino;
    }
    return 0;
}

################################################################################

sub header {
    my $self = shift;
    return $self->{_header};
}

################################################################################

sub length {
    my $self = shift;
    return scalar @{ $self->{_residues} };
}

################################################################################

sub method {
    return $_[0]->{_method};
}

################################################################################

sub nchains {
    return $_[0]->{_nchains};
}

################################################################################


sub residue {
    my ( $self, %args ) = @_;

    my $chain = $args{-chain};
    my $ins   = $args{-insert};
    my $num   = $args{-number};

    unless ( defined($num) ) {
        confess
"No sequence number specified in arguments to residue()\nArguments are: "
          . join q{ }, %args;
    }

    if ( !defined $chain ) {
        $chain = q{ };
    }
    if ( !defined $ins ) {
        $ins = q{ };
    }

    return $self->{_residue_map}{$chain}{"$num-$ins"}
      or croak("Cannot find residue $num$ins in chain $chain\n");
}

################################################################################

sub residues {
    my $self = shift;

    if (@_) {
        my ( $nt, $ct );
        if ( ref( $_[0] ) eq 'PDB::Residue' && ref( $_[1] ) eq 'PDB::Residue' )
        {
            $nt = $_[0]->index;
            $ct = $_[1]->index;
        }
        else {
            ( $nt, $ct ) = @_;
        }
        return @{ $self->{_residues} }[ $nt .. $ct ];
    }
    else {
        return @{ $self->{_residues} };
    }
}

################################################################################

sub resolution {
    return $_[0]->{_resolution};
}

################################################################################

sub title {
    return $_[0]->{_title};
}

################################################################################

sub write {
    my ( $self, $pdbfh ) = @_;

    $self->write_pdb_header($pdbfh);

    for my $res ( @{ $self->{_residues} } ) {
        for my $atom ( @{ $res->atoms } ) {
            print {$pdbfh} $atom->to_record, "\n";
        }
    }
    print {$pdbfh} "TER\n";
    print {$pdbfh} "END\n";
    return;
}

sub write_header {
    my ( $self, $pdbfh ) = @_;
    for my $title (qw(header title compnd expdta)) {
        if ( $self->{"_$title"} ) {
            my $header = $self->{"_$title"};
            for my $line ( split /\n/, $header ) {
                printf {$pdbfh} "%-6s    %s\n", uc($title), $line;
            }
        }
    }
    if ( $self->{_seqres} ) {
        print {$pdbfh} $self->{_seqres}, "\n";
    }
    if ( $self->{_secstruc} ) {
        print {$pdbfh} $self->{_secstruc}, "\n";
    }
    return;
}

sub write_chain {
    my ( $self, $pdbfh, $chain ) = @_;
    for my $res ( @{$chain} ) {
        for my $atom ( @{ $res->atoms } ) {
            print {$pdbfh} $atom->to_record, "\n";
        }
    }
    print {$pdbfh} "TER\n";
    return;
}

################################################################################

# This must be explicitly defined because there are circular references
# between adjacent residues.

sub DESTROY {
    my $self = shift;

    for my $r ( $self->residues ) {
        $r->DESTROY if defined($r);
    }

    delete $self->{_residues};
    return;
}

1;

=head1 NAME

  PDB::Entry

=head1 AUTHOR

  Tom Walsh

=head1 SYNOPSIS

 Class for representing a PDB entry

=head1 DESCRIPTION

 Creates objects representing PDB files.

 Note that in NMR structures, only the first model is read.

=head1 CONSTRUCTOR

=head2 new

Usage    : $pdb = new PDB::Entry($pdbfile,%options);
           $pdb = new PDB::Entry($pdb_filehandle,%options);

Function : Read a PDB entry from the specified file or filehandle

Options  : -trace
            Read only the backbone CA atoms.

Comment  : If the file contains multiple NMR models of a 
           single structure, only the first model is read.

           Confesses if the PDB file cannot be opened.

=head1 METHODS

=head2 chain($chn)

Returns : List of the  residues in the specified chain or an empty list
          if no such chain exists.

=head2 chains

Returns a reference to a hash in which the keys are chain
letters and the values are references to the lists of
residues in each chain.

=head2 compnd

The COMPND record from the PDB file.

=head2 expdta

 The EXPDTA record from the PDB file

The file from which the object was created

=head2 has_amino_residues

=over

=item Usage

$pdb->has_amino_residues()

=item Function 

Returns true if the structure contains standard amino acid residues

=back

=head2 header

Returns the PDB header record

=head2 nchains

The number of chains in the PDB file.

=head2 residue

=over

=item Usage

$residue = residue(number => $number,chain => $chain, insert => $insert);

=item Function 

Get the residue with the specified sequence number, chain identifier
and insert code.

=item Comment  

The chain identifier and insert code default to blank spaces if they
are not specified.

If there are two consecutive residues with the same sequence number
and insert code, this method will return the first one. Use the next()
method of the returned residue to retrieve the second residue.

=back

=head2 residues

=over

=item Usage

residues(nt,ct)

=item Function 

Return a list of the residues in the PDB object

=item  Argument 

Optional: beginning and end of the required segment

=back

=head2 resolution

Return the resolution of the PDB file. If resolution is not
applicable to this file, the value is undefined.

=head2 title

The TITLE record from the PDB file.

=head2 write

=over

=item Usage 

$pdb->write($pdbfh);

=item Function 

Write the PDB structure to the specified filehandle

=back

=head2 write_chain

=over

=item Usage 

$pdb->write($pdbfh,$chain);

=item Function 

Write a chain of residues to the specified filehandle.

$chain is an ARRAYREF of PDB::Residue objects.

=back

=head1 CLASS VARIABLE

=head2 $chain_break_threshold

The minimum distance between consecutive alpha carbons for a chain
break to be defined between the corresponding residues.

=head1 DEPENDENCIES

=over

=item Log::Log4perl

=item PDB::Atom

=item PDB::Residue

=back

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 LICENSE AND COPYRIGHT

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the LICENSE file included with this module.

=cut

