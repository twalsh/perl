package STAMP::Tools;

# $Id: Tools.pm,v 1.21 2009/01/23 10:40:49 tom Exp $

use strict;
use warnings;
use Carp;
use Exporter;
use FileHandle;
use Log::Log4perl;
use STAMP::Domain;

use base qw(Exporter);
our @EXPORT_OK = qw(descriptor make_descriptor);
our $verbose   = 0;
my $log = Log::Log4perl->get_logger("STAMP::Tools");
our $NULL_CHAIN_CHAR = 'A';

sub descriptor {
    my $dom = shift;
    return sprintf
      "%s %s %s",
      $dom->pdbid, substr( $dom->sid, 1 ), make_descriptor($dom);
}

sub make_descriptor {
    my ( $dom, $options ) = @_;
    my $keep_null_ids = $options->{-keep_null_chain_ids} || 0;
    my @defn          = $dom->parse_coordinate_definition;
    my $stamp         = q{};
    if ( @defn == 0 ) {
        $stamp = 'ALL';
    }
    else {
        for my $chn (@defn) {
            if ( @{$chn} == 1 ) {
                $stamp = join q{ }, $stamp, "CHAIN $chn->[0]";
            }
            else {
                my $c;
                if ( $chn->[0] ne q{ } ) {
                    $c = $chn->[0];
                }
                elsif ($keep_null_ids) {
                    $c = '_';
                }
                else {
                    $c = $NULL_CHAIN_CHAR;
                }

                # Insertion codes. Blank spaces are represented by '_'
                my $ins1 = $chn->[2] eq q{ } ? '_' : $chn->[2];
                my $ins2 = $chn->[4] eq q{ } ? '_' : $chn->[4];

                $stamp .= uc " $c $chn->[1] $ins1 to $c $chn->[3] $ins2";
            }
        }
    }

    $stamp =~ s/^\s//;

    return "{ $stamp }";
}

################################################################################

my $require_ss;

sub check_descriptor {
    my ( $domain, $pdb_entry, $dssp_entry, $options ) = @_;

    my $require_ss = $options->{-require_ss};

    my $keep_null_ids = $options->{-keep_null_chain_ids} || 0;
    $require_ss ||= 0;

    my $desc      = undef;
    my $ndesc     = q{};
    my $nresidues = 0;
    $desc = make_descriptor($domain);

    my $domid        = $domain->sid;
    my $stamp_domain = STAMP::Domain->new("unk unk $desc");
    my $chain        = uc substr( $domid, 5, 1 );

    if ( $chain eq '.' ) {
        $log->debug("multichain descriptor: $domid $desc");
        $desc = undef;
    }
    else {
        if ( $chain eq '_' ) {
            if ($keep_null_ids) {
                $chain = ' ';
            }
            else {
                $chain = $NULL_CHAIN_CHAR;
            }
        }

        for my $seg ( $stamp_domain->segments ) {
            my $records = [];
            if ($dssp_entry) {
                $records = ( $dssp_entry->get_domain($seg) )[0];
            }
            $nresidues = @{$records};
            if ( $seg !~ /(\w) (-?\d+) (\w?) TO (\w) (-?\d+) (\w?)/o ) {
                if ( !( grep { substr( $_, 11, 1 ) eq $chain } @{$records} ) ) {

                    # No records. Could be either a nonprotein chain
                    # or a CA-only structure. STAMP can handle the
                    # latter. In any case, the only way to sort this
                    # out is to read the PDB file.
                    $log->debug("No DSSP records found for $domid $seg");
                    $log->debug(
                        "Searching for $seg in PDB file " . $pdb_entry->file );
                    my $found = 0;

                    my @residues;
                    if ( $seg =~ /ALL/ ) {
                        @residues = $pdb_entry->residues;
                    }
                    elsif ( $seg =~ /CHAIN (\w)/ ) {
                        my $chn = $1;
                        @residues = $pdb_entry->chain($chn);
                    }

                    $nresidues = @residues;
                    unless (@residues) {
                        $log->debug( "Cannot find residues for $seg in "
                              . $pdb_entry->file );
                    }
                    else {
                        for my $r (@residues) {
                            if ( $r->is_amino && !$r->is_het ) {
                                $found = 1;
                                last;
                            }
                        }
                    }

                    unless ($found) {
                        $log->debug( "No amino acid residues found for $seg in "
                              . $pdb_entry->file );
                        $seg = undef;
                    }
                }
            }
            else {

                # Domain is a segment of a chain. Check that the termini
                # are OK
                $seg = check_termini( $seg, $pdb_entry, $records );

                # Estimate number of residues in segment
                my ( $first, $last ) = ( split q{ }, $seg )[ ( 1, 5 ) ];
                $nresidues = $last - $first;
            }
            $ndesc .= $seg if $seg;
        }
        if ($ndesc) {
            $ndesc = "{ $ndesc }";
            $ndesc =~ s/\s+/ /g;
        }
        else {
            $ndesc = undef;
        }
    }

    return ( $ndesc, $log, $nresidues );
}

################################################################################

sub check_termini {
    my ( $desc, $pdb_entry, $dssp_records ) = @_;

    my ( $bc, $bn, $bi, $ec, $en, $ei );
    my $desc_re = qr/(\w)\s+(-?\d+)\s+(\w?)\s+TO\s+(\w)\s+(-?\d+)\s+(\w?)/io;
    if ( $desc =~ $desc_re ) {
        ( $bc, $bn, $bi, $ec, $en, $ei ) = ( $1, $2, $3, $4, $5, $6 );
    }
    else {
        confess "Descriptor $desc doesn't match regex $desc_re";
    }

    for ( $bc, $bi, $ec, $ei ) {
        $_ = ' ' if $_ eq '_';
    }

    my @records = grep { substr( $_, 11, 1 ) =~ /$bc|$ec/ } @$dssp_records;

    my $term_ok;
    my $modified = 0;
    my ( $new_bn, $new_en );

    if ( @records == 0 ) {
        $log->debug( 'No DSSP records for termini found. Searching in PDB file '
              . $pdb_entry->file );

        # Try to find the termini in the PDB file
        ( $term_ok, $bc, $bn, $bi, $ec, $en, $ei, $modified ) =
          check_pdb_termini( $desc, $pdb_entry, $bc, $bn, $bi, $ec, $en, $ei );
        $log->debug('Termini not found in PDB file') unless $term_ok;
    }
    else {
        my $nt;
        my $found_nt;
        $new_bn = $bn;

      NTERM: while ( $new_bn < $en ) {

            for my $i ( 0 .. $#records ) {
                my $residue = $records[$i];
                my $seq = substr( $residue, 13, 1 );
                next if $seq =~ /X|!/;
                my $resid = make_resid($residue);
                if ( $resid eq "$bc$new_bn$bi" ) {
                    if ( $new_bn != $bn ) {
                        $log->debug("modified start residue $bc $new_bn$bi");
                        $modified = 1;
                    }
                    $nt       = $i;
                    $found_nt = 1;
                    last NTERM;
                }
            }
            $new_bn++;
        }

        my $ct;
        my $found_ct;
        $new_en = $en;

      CTERM: while ( $new_en > $new_bn ) {
            for my $i ( reverse( 0 .. $#records ) ) {
                my $residue = $records[$i];
                my $seq = substr( $residue, 13, 1 );
                next if $seq =~ /X|!/;
                if ( make_resid($residue) eq "$ec$new_en$ei" ) {
                    if ( $new_en != $en ) {
                        $log->debug("modified end residue $ec $new_en$ei");
                        $modified = 1;
                    }
                    $ct       = $i;
                    $found_ct = 1;
                    last CTERM;
                }
            }
            $new_en--;
        }

        if ( $found_nt && $found_ct ) {
            if ($require_ss) {
                my $has_ss = 0;
                for my $r ( $nt .. $ct ) {
                    my $ss = substr( $records[$r], 16, 1 );
                    if ( $ss ne ' ' ) {
                        $has_ss++;
                    }
                }

                if ( $has_ss == 0 ) {
                    $log->debug(
                        'Domain appears to have no secondary structure');
                    $log->debug('No start residue found') unless $found_nt;
                    $log->debug('No end residue found')   unless $found_ct;
                    $term_ok = 0;
                }
                else {
                    $term_ok = 1;
                }
            }
            else {
                $bn      = $new_bn;
                $en      = $new_en;
                $term_ok = 1;
            }
        }
        else {
            $log->debug('No start residue for domain') unless $found_nt;
            $log->debug('No end residue for domain')   unless $found_ct;
        }
    }

    if ($modified) {
        $log->debug("Modified descriptor: $desc");
    }

    if ($term_ok) {
        ( $bc, $bi, $ec, $ei ) =
          map { $_ eq q{ } ? '_' : $_ } ( $bc, $bi, $ec, $ei );

        return " $bc $bn $bi TO $ec $en $ei ";
    }
    else {
        return;
    }
}

###############################################################################

sub find_residue {
    my ( $pdbentry, $chn, $num, $ins ) = @_;

    my $res = $pdbentry->get_residue(
        -chain  => $chn,
        -number => $num,
        -insert => $ins
    );
    if ( $res && !$res->is_het && $res->atom('CA') ) {
        return $res;
    }
    else {
        return;
    }
}

###############################################################################

sub check_pdb_termini {
    my ( $domid, $pdb_entry, $bc, $bn, $bi, $ec, $en, $ei ) = @_;

    my $term_ok  = 0;
    my $new_bn   = $bn;
    my $found_nt = 0;
    $log->debug("searching for $bc $new_bn $bi");
    while ( $new_bn < $en ) {
        if ( $found_nt = find_residue( $pdb_entry, $bc, $new_bn, $bi ) ) {
            last;
        }
        $new_bn++;
    }
    if ($found_nt) {
        if ( $new_bn != $bn ) {
            $log->debug("modified initial residue $bc $bn$bi");
        }
    }
    else {
        $log->debug(
            "Cannot find start residue $bc $bn $bi in " . $pdb_entry->file );
    }

    my $found_ct;
    my $new_en = $en;

    while ( $new_en > $new_bn ) {
        last if ( $found_ct = find_residue( $pdb_entry, $ec, $new_en, $ei ) );
        $new_en--;
    }

    if ($found_ct) {
        if ( $new_en != $en ) {
            $log->debug("modified final residue $ec $new_en$ei");
        }
    }
    else {
        $log->debug( "Cannot find end residue $ec $en $ei in ",
            $pdb_entry->file );
    }

    my $modified = ( $new_en != $en ) || ( $new_bn != $bn );

    # print STDERR "$new_en $en $new_bn $bn $modified\n";

    return ( ( $term_ok = $found_nt && $found_ct ),
        $bc, $new_bn, $bi, $ec, $new_en, $ei, $modified );
}

###############################################################################

sub make_resid {
    my $residue = shift;
    my $resid = substr( $residue, 6, 5 );
    $resid =~ s/^\s*//g;
    my $chn = substr( $residue, 11, 1 );

    return "$chn$resid";
}

################################################################################

sub descriptor_ok {
    my ( $stamp_domain, $pdb_entry ) = @_;

    confess("No PDB file specified") unless $pdb_entry;

    my $ndesc;

    my $pdb = $stamp_domain->pdb;
    my $id  = $stamp_domain->id;

    for my $seg ( $stamp_domain->segments ) {
        if ( $seg !~ /(\w) (-?\d+) (\w?) TO (\w) (-?\d+) (\w?)/o ) {
            $log->debug( "Searching for $seg in PDB file " . $pdb_entry->file );
            my $found = 0;

            my @residues;
            if ( $seg =~ /ALL/ ) {
                @residues = $pdb_entry->residues;
            }
            elsif ( $seg =~ /CHAIN (\w)/ ) {
                my $chn = $1;
                @residues = $pdb_entry->chain($chn);
            }

            unless (@residues) {
                croak "% Cannot find residues for $seg in "
                  . $pdb_entry->file . "\n";
            }
            else {
                for my $r (@residues) {
                    if ( $r->is_amino && !$r->is_het ) {
                        $found = 1;
                        last;
                    }
                }
            }

            unless ($found) {
                croak "% No amino acid residues found for $seg in "
                  . $pdb_entry->file . "\n";
                $seg = undef;
            }
        }
        else {

            # Domain is a segment of a chain. Check that the termini
            # are OK
            $seg = check_termini( $seg, $pdb_entry );
        }
        $ndesc .= "$seg " if $seg;
    }

    if ($ndesc) {
        $ndesc = "{ $ndesc }";
        $ndesc =~ s/\s+/ /g;
        return "$pdb $id $ndesc";
    }
    else {
        return q{};
    }
}

sub pdb2stamp {
    my ( $pdb, $pdbid ) = @_;

    my ( $desc, $nres ) = ( [], [] );
    while ( my ( $chain, $residues ) = each %{ $pdb->chains } ) {
        my $is_protein = 0;
        for my $r (@$residues) {
            if ( $r->is_amino && !$r->is_het ) {
                $is_protein = 1;
                last;
            }
        }
        if ($is_protein) {
            my $coord;
            if ( $chain eq q{ } ) {
                $coord = "ALL";
            }
            else {
                $coord = "CHAIN $chain";
            }
            $chain = "_" if $chain eq q{ };
            $chain = lc $chain;
            push @$desc, "$pdbid $pdbid${chain}_ { $coord }";
            push @$nres, scalar(@$residues);
        }
    }
    return ( $desc, $nres );
}

sub stamp2pdb {
    my ( $domain, $pdbentry ) = @_;

    my @residues;
    for my $seg ( $domain->segments ) {
        if ( $seg =~ /(\w) (-?\d+) (\w?) TO (\w) (-?\d+) (\w?)/io ) {
            my ( $bc, $ec ) =
              map { $_ eq '_' ? $NULL_CHAIN_CHAR : $_ } ( $1, $4 );
            my ( $bn, $en ) = ( $2, $5 );
            my ( $bi, $ei ) =
              map { $_ eq '_' ? $NULL_CHAIN_CHAR : $_ } ( $3, $6 );
            my $start = $pdbentry->residue(
                -chain  => $bc,
                -number => $bn,
                -insert => $bi
            );
            if ( !$start ) {
                confess "Cannot find residue $bc $bn $bi";
            }
            my $end = $pdbentry->residue(
                -chain  => $ec,
                -number => $en,
                -insert => $ei
            );
            if ( !$end ) {
                confess "Cannot find residue $ec $en $ei";
            }
            push @residues, $pdbentry->residues( $start, $end );
        }
        elsif ( $seg =~ /ALL/ ) {
            @residues = $pdbentry->residues;
        }
        elsif ( $seg =~ /CHAIN (\w)/ ) {
            my $chn = $1;
            push @residues, $pdbentry->chain($chn);
        }
        else {
            confess "Cannot match segment $seg\n";
        }
    }
    @residues = grep { $_->is_amino && !$_->is_het } @residues;
    return \@residues;
}

1;

=head1 NAME

STAMP::Tools

=head1 SYNOPSIS

Utility functions for working with STAMP.

=head1 FUNCTIONS

=head2 descriptor($SCOP_Domain)

Returns a full STAMP descriptor for the specified SCOP domain.

=head2 make_descriptor($domain,$options)

=head3 Function 

Returns a STAMP coordinate descriptor for a SCOP::Domain object.

=head3 Argument 

SCOP::Domain object or SCOP domain descriptor string.

=head3 Options

-keep_null_chain_ids => 0

Translate null chain ids in SCOP identifers

=head3 Returns 

STAMP coordinate descriptor string.

=head2 check_descriptor($domain,$pdb_entry,$dssp_entry)

=head3 Function

Generates a descriptor for the specified SCOP::Domain object and
checks it using the DSSP and PDB files. It checks that the
coordinate range specified in the descriptor exists in the PDB
structure. If start and end residues are specified, it checks that
these residues have CA atoms that are not in HETATM records. This
ensures that the coordinates can be loaded by STAMP. The DSSP file is
used first; if the residues cannot be found in this file, the PDB file
is opened.

=head3 Arguments

$domain - SCOP::Domain object

$pdb_entry - PDB:::Entry object representing the PDB file containing the domain

$dssp_entry - DSSP::File object representing the corresponding DSSP file.

$options

=head4 -keep_null_chain_ids => 1

Null chain ids (' ') are not converted to 'A' in the returned descriptor.

=head4 -require_ss => 1

Require that the input structure must have secondary structure.

Returns

A list consisting of the descriptor and a list reference for the
logging output generated during checking. If the segment specified by
the descriptor could not be found, the returned descriptor will be
undef.

=head2 inspect_descriptor($desc,$pdbfile)

=head3 Function 

Inspects a STAMP descriptor and cross-checks it with the PDB file to
ensure that the specified domain can be found and is a protein chain.

Returns

The input descriptor, possibly modified to take account of changes in
the start and end residues of a segment, owing to CA atoms missing
from the PDB files.  

=head2 pdb2stamp($pdbfile,$pdbid)

=head3 Function

Create STAMP descriptors for the protein chains in a PDB file.

=head3 Arguments

$pdbfile - PDB file

$pdbid   - PDB id code

=head3 Returns

List of STAMP descriptors

List of the sizes of the corresponding chains.

=head2 stamp2pdb($stamp_domain,$pdb_entry)

=head3 Function

Extract a domain specified by a STAMP domain descriptor
from a PDB file.

=head3 Arguments

$stamp_domain - STAMP::Domain object

$pdb_entry  - PDB::Entry object

=head3 Returns

Return an ARRAYREF of the PDB::Residue objects from the
PDB::Entry object that correspond to the
domain defined by the STAMP descriptor $stamp_domain.

=head1 DEPENDENCIES

=over

=item Log::Log4perl

=item STAMP::Domain

=back

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT AND LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut
