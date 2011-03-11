package STAMP::Scan;

# $Id: Scan.pm,v 1.11 2009/01/24 16:28:41 tom Exp $

use strict;
use Carp;
use FileHandle;
use STAMP::Hit;

our $debug = 0;
our $VERBOSE = 0;

sub new {
    my ($type,$arg) = @_;
    
    unless (ref($arg) eq 'HASH'){
        confess("Argument to new() must be an option hash");
    }

    my $file    = $arg->{file} || q{};
    my $fh      = $arg->{fh};
    my $getline = $arg->{getline};

    my $next_line;

    if ($getline){
	$next_line = $getline;
    }
    elsif ($fh) {
	$next_line = sub { $fh->getline() };
    }
    else {
	open my $nfh, $file or confess("Cannot open STAMP scan file: $file");
	$next_line = sub { return <$nfh> };
    }
    
    print STDERR "STAMP::Scan::new: File $file opened for reading\n" if $debug;
    
    my $self = {_file => $file,
		_domains => [],
		_hits => [], 
		_query => undef,
		_slide => undef};
    
    my $query = 0;
    my $cregex = $STAMP::Domain::coord_definition_re;
    
    while ($_ = $next_line->()){
	if (/^%/){
	    if (/at every\s+(\d+)\s+residue of the database sequences/){
		# Slide parameter
		$self->{_slide} = $1;
	    }
	    elsif (/Transformations were output for Sc=\s+(\S+)/){
		# Score cutoff
		$self->{_cutoff}{Sc} = $1;
	    }
	} elsif (/^\# Sc=/){
	    # Start reading a hit
	    my $score = $_;
	    
	    # Get the domain descriptor. The while loop is required to skip
	    # over empty lines that SORTTRANS adds to files that are
	    # themselves products of SORTTRANS. 
	    my $domain = q{};
	    my $i = 0;
	    while (!$domain && ++$i < 3){
		$domain = $next_line->();
		chomp $domain if $domain;
	    }
	    
	    confess "Error reading file $file: no data\n" unless $domain;
	    
	    $domain .= "\n";
	    
	    if ($domain !~ /\}\n$/){
		# The descriptor contains a transformation, so read the
		# next 3 lines
		for (1..3){
		    $domain .= $next_line->();
		}
	    }
	    
	    my $hit = new STAMP::Hit($score,$domain);
	    
	    # First hit is actually the query domain. 
	    if ($self->{_query}){
		push @{$self->{_hits}},$hit;
	    }
	    else {
		$self->{_query} = $hit;
	    }
	    push @{$self->{_domains}},$hit;
	}
	elsif (/^\S+ \S+ { $cregex/){
	    my $domain = $_;
	    
	    if ($domain !~ /\}\n$/){
		# The descriptor contains a transformation, so
		# read the next 3 lines
		for (1..3){
		    $domain .= $next_line->();
		}
	    }
	    
	    my $dom = STAMP::Domain->new($domain);
	    unless ($self->{_query}){
		$self->{_query} = $dom;
	    }
	    
	    push @{$self->{_domains}},$dom;
	}
	else {
	    # It's not a comment and it can't be part of a descriptor, so
	    # something is wrong with the file.
	    print STDERR "Error reading $file\n";
	    print STDERR "Line: $_\n";
	    print STDERR "This line isn't a comment or part of a descriptor\n";
	    confess;
	}
    }
    
    $self->{_nhits} = @{ $self->{_hits} };

    bless $self,$type;
}

################################################################################

=head2 domains 

=cut

sub domains { 
    @{$_[0]->{_domains}};
}

################################################################################

=head2 file

Returns the name of the scan file.

=cut

sub file {
    $_[0]->{_file};
}

################################################################################

=head2 hits

  Return the hits from STAMP file as an arrayref of STAMP::Hit objects.

=cut

sub hits {
    return [ @{ $_[0]->{_hits} } ];
}

=head2 nhits

Returns the number of hits in the STAMP output.

=cut

sub nhits {
    $_[0]->{_nhits};
}

################################################################################

=head2 query()

  Returns a STAMP::Domain object that represents the query
domain in the scan file.

=cut

sub query {
    $_[0]->{_query};
}

################################################################################

=head2 score_cutoff

Returns the score cutoff used. This is taken from the header of the
scan file.

=cut

sub score_cutoff {
    $_[0]->{_cutoff}{Sc};
}

################################################################################

=head2 slide

Returns the slide parameter value

=cut

sub slide {
    $_[0]->{_slide};
}

################################################################################

=head1 CLASS VARIABLE

$debug

Set to 1 to turn on debugging output.

=cut

################################################################################

1;

=pod

=head1 NAME

  STAMP::Scan

=head1 SYNOPSIS

  Object-oriented interface to STAMP domain files

=head1 METHODS 

=head2 new({ file => $file} )

  Create STAMP::Scan object from a STAMP domain file. 

=head2 new({ fh => $FileHandle} )

  Create STAMP::Scan object by reading from the given FileHandle object.

=head2 new({ getline => $GetNextLineFunce} )

  Create STAMP::Scan object by calling a subroutine which returns the
next line in a STAMP scan file at each iteration, and undef when EOF
is reached.

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT AND LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut

