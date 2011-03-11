package STAMP::FileFinder;
use strict;
use warnings;
use Carp;
use Log::Log4perl;
use Object::Tiny qw(patterns dirfile);

my $log = Log::Log4perl->get_logger("STAMP::FileFinder");

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(@_, patterns => []);
    open my $FILE, '<', $self->dirfile or croak 'Cannot read directory file ',$self->dirfile;
    while (<$FILE>){
        chomp;
        next if /^%/;
        next if /^\s*$/; 
        my ($dir, $prefix, $suffix) = split;
        $prefix = q{} if $prefix eq '_';
        $suffix = q{} if $suffix eq '_';
        push @{$self->patterns}, { dir => $dir, prefix => $prefix, suffix => $suffix}; 
    }
    return $self;
}

sub find {
    my ($self,$code) = @_;
    for my $p (@{$self->patterns}){
        my $filename = $p->{dir}.'/'.$p->{prefix}.$code.$p->{suffix};
        $log->trace("searching for $filename");
        if (-e $filename){
            $log->debug("found $filename");
            return $filename;
        }
    }
    return;
}

1;

=head1 NAME

STAMP::FileFinder

=head1 SYNOPSIS

    use STAMP::FileFinder;
    my $finder = STAMP::FileFinder->new("$ENV{STAMPDIR}/pdb.directories");
    my $pdbcode = "1flt";
    my $file = $finder->find($pdbcode) or die "cannot find PDB file for $pdbcode\n";

=head1 DESCRIPTION

Finds PDB/DSSP files using the directory files used by STAMP. See the STAMP
manual for an explanation of the format of directory files.

=head1 DEPENDENCIES

    Log::Log4perl
    Object::Tiny

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT & LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut

