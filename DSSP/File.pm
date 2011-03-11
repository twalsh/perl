package DSSP::File;

# $Id: File.pm,v 1.6 2006/12/20 15:50:40 twalsh Exp $
# $Log: File.pm,v $
# Revision 1.6  2006/12/20 15:50:40  twalsh
# *** empty log message ***
#
#
# Revision 1.1  2003/12/03 17:03:27  tom
# First checkin.
#

=head1 NAME

DSSP::File

=head1 SYNOPSIS

DSSP file object

=cut

use strict;
use Carp;

=head1 CONSTRUCTOR

=head2 new($dsspfile)

Create a new DSSP object by reading a DSSP file

=cut

sub new {
    my ($type,$dsspfile) = @_;

    my $self = {};

    open DSSP,$dsspfile or confess "Cannot open DSSP file: $dsspfile";

    @{$self->{_records}} = <DSSP>;
    chomp @{$self->{_records}};

    while (@{$self->{_records}} && $self->{_records}[0] !~ /^  \#  RESIDUE AA/) {	
	push @{$self->{_header_lines}}, shift @{$self->{_records}};
    }

    push @{$self->{_header_lines}}, shift @{$self->{_records}};

    $self->{_header} = (split q{ },$self->{_header_lines}[2],2)[1];
    $self->{_compnd} = (split q{ },$self->{_header_lines}[3],2)[1];
    $self->{_source} = (split q{ },$self->{_header_lines}[4],2)[1];
    $self->{_file} = $dsspfile;
    _make_seq($self);
  
    bless $self, $type;
}

sub _make_seq {
    my $self = shift;
    $self->{_serial} = {};
    my $offset = 0;
    for my $i (0 .. @{$self->{_records}}) {
        my $record = $self->{_records}->[$i];
        my $oseq = substr($record,0,5) + 0;
        my $insert = 0;
        if (substr($record,13,1) eq "!") {
            if (record_chain($self,$self->{_records}->[$i-1]) eq record_chain($self,$self->{_records}->[$i+1])){
        #        my $brk1 = substr($self->{_records}->[$i-1],6,5);
        #        my $brk2 = substr($self->{_records}->[$i+1],6,5);
        #        $insert = $brk2 - $brk1 - 2;
        #        $offset += $insert;   
                $offset--;
            }
            else {
                $offset--;
            }
        }
        $self->{_serial}{$oseq} = $oseq + $offset; 
    }
}

sub record_chain {
    my ($self,$record) = @_;
    return substr($record,11,1);
}

sub record_ss {
    my ($self,$record) = @_;
    return substr($record,16,1); 
}

sub serial {
    my ($self,$record) = @_;
    my $num = substr($record,0,5) + 0;
    return $self->{_serial}{$num};
}

=head1 METHODS

=cut

################################################################################

sub compnd {
    my $self = shift;
    $self->{_compnd};
}

################################################################################

=head2 file()

  Returns the name of the DSSP file.

=cut

sub file {
    $_[0]->{_file};
}

################################################################################

=head2 get_chain

 Usage: $dssp->get_chain($chain)

 Returns: reference to a list of records for the specified chain

=cut

sub get_chain {
    # Return reference to a list of records for the specified chain
    my ($self,$chain) = @_;

    # Find start of chain
    my $offset = 0;
    for my $i (0 .. $#{$self->{_records}}){
	if (substr($self->{_records}[$i],11,1) eq $chain){
	    $offset = $i;
	    last;
	}
    }
    return ([grep { substr($_,11,1) eq $chain } @{$self->{_records}}],$offset);
}

################################################################################

=head2 get_domain

 Usage: $dssp->get_domain($domain_descriptor)

 Returns: reference to a list of records for a domain specified using a
          STAMP-style descriptor.

=cut

sub get_domain {
    # Return reference to a list of records for a domain specified using a
    # STAMP-style descriptor.
    my ($self,$stamp_desc) = @_;

    if ($stamp_desc =~ /CHAIN (\w)/) {
	# Return chain
	return $self->get_chain($1);
    }
    elsif ($stamp_desc =~ /ALL/) {
	# Return all records
	return ($self->{_records},0);
    }
    elsif ($stamp_desc =~ /(\w) (-?\d+) (\w) TO \w (-?\d+) (\w)/) {
	# Get specific range of residues
	my $chain = $1;
	my ($ins1,$ins2) = ($3,$5);
	my ($start,$end) = ($2,$4);

	# Make DSSP-style residue ids.
	$start = sprintf "%6s", "$start$ins1$chain";
	$end = sprintf "%6s", "$end$ins2$chain";

	$start =~ s/_/ /g;
	$end =~ s/_/ /g;

	# Get indices for start and end of the domain.
	my $first = $self->get_index($start);
	my $last = $self->get_index($end);
	
	if (defined($first) && defined($last)){
	    return ([@{$self->{_records}}[$first-1 .. $last-1]],$first-1);
	}
	else {
	    return ([]);
	}
    }
    else {
	confess "Not a valid STAMP domain descriptor: $stamp_desc\n";
    }
}

################################################################################

=head2 get_index($res_spec)

Returns the index of the record for the specified residue. The
specification must have the same format as in the DSSP record i.e. the
sequence number with the insert code and the chain identifier
appended. Records are indexed from 1.

=cut

sub get_index {
    my ($self,$res_spec) = @_;
    $res_spec = sprintf "%6s",$res_spec;
    for my $i (1 .. @{$self->{_records}}) {
	my $record = $self->{_records}[$i-1];
	my $residue = substr($record,6,6);
	if ($residue eq $res_spec) {
	    return $i;
	}
    }
    #carp "Cannot find residue <$res_spec> in DSSP file " . $self->{_file};
    return undef;
}

################################################################################

=head2 get_residue($res_spec)

Returns the record for the specified residue. The specification must
have the same format as in the DSSP record.i.e. the sequence number
with the insert code and the chain identifier appended.

=cut

sub get_residue {
    my ($self,$res_spec) = @_;

    $res_spec = sprintf "%6s",$res_spec;
    for my $record (@{$self->{_records}}) {
	my $residue = substr($record,6,6);
	if ($residue eq $res_spec) {
	    return $record;
	}
    }
    carp "Cannot find residue <$res_spec> in DSSP file " . $self->{_file};
    return undef;
}

################################################################################

=head2 header()

Returns the header section from the DSSP file.

=cut

sub header {
    my $self = shift;
    $self->{_header};
}

################################################################################

=head2 num_records()

Returns the number of records in the DSSP file.

=cut

sub num_records {
    my $self = shift;
    scalar @{$self->{_records}};
}

################################################################################

sub source {
    $_[0]->{_source};
}

################################################################################

1;

=head1 DEPENDENCIES

None

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT & LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut


