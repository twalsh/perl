package SCOP::Lite::Database;

use strict;
use warnings;

use Carp;
use Log::Log4perl;
use Object::Tiny qw(sid release sunid domain_count);
use SCOP::Lite::Domain;

my $log = Log::Log4perl->get_logger("SCOP::Lite::Database");

sub new {
    my ( $class, %arg ) = @_;

    my $key = $arg{key} || 'sid';
    my $file = $arg{file} or croak "No file option give to constructor\n";
    my $domain_list = $arg{domain_list};
    my %get_domain;
    if ($domain_list) {
        for my $domain (@{$domain_list}){
            if ($domain !~ /^d\d\w{5}$/){
                $log->error_die("domain list contains illegally formatted domain id: $domain.\n");
            }
        }
        %get_domain = map { $_ => 1 } @{$domain_list};
    }
    my $self = $class->SUPER::new( sid => {}, sunid => {}, domain_count => 0 );
    open my $IN, '<', $file or croak "$!: $file";
    while (<$IN>) {
        if (/^#/) {
            if (/^# SCOP release (\S+)/o) {
                $self->{release} = $1;
            }
            next;
        }
        if (%get_domain) {
            my $domid = (split)[0];
            next unless $get_domain{$domid};
        }
        my $domain;
        if ( eval { $domain = SCOP::Lite::Domain->new($_) } ) {
            $self->{domain_count}++;
            if ( $key eq 'both' || $key eq 'sid' ) {
                $self->{sid}->{ $domain->sid } = $domain;
            }
            if ( $key eq 'both' || $key eq 'sunid' ) {
                $self->{sunid}->{ $domain->sunid } = $domain;
            }
        }
        else {
            print $@;
            croak "Error reading domains from $file\n";
        }
    }

    return $self;
}

sub domain {
    my ( $self, $key ) = @_;
    if ( $key =~ /^(\d+)$/ ) {
        return $self->sunid->{$key};
    }
    elsif ( $key =~ /^[gsd]\d\w{3}[\w\.]\w/ ) {
        return $self->sid->{$key};
    }
    else {
        croak "Not a SCOP sid or sunid: $key";
    }
}

sub domains {
    my $self = shift;
    return [ values %{ $self->sid } ];
}

1;

=pod 

=head1 NAME

SCOP::Lite::Database

=head1 SYNOPSIS

Lightweight object for representing domains from a SCOP database release as
a collection of SCOP::Lite::Domain objects.

=head1 CONSTRUCTOR

$domain = SCOP::Lite->new(file => $file, domain_list => \@domain_list)

head2 $file 

SCOP classication file (.cla file).

head2 @domain_list 

(optional)

List of domains to be read from the SCOP file.  If only a subset of the SCOP
domains are required, using a domain list makes loading much faster.

If this argument is not specified, then all of the domains listed in the
classification file will be loaded.

=head1 METHODS

Methods for retrieving properties of the domain. Read the SCOP documentation
for details of what each property is.

sunid is the SCOP unique identifer (the integer identifier assigned to
every object in SCOP).

=head2 release 

SCOP release number

=head2 domain($sid)

=head2 domain($sunid)

Retrieve the SCOP domain with the specified sid or sunid.

=head2 domains

Retrieve a list of the domains that have been loaded.

=head1 DEPENDENCIES

Object::Tiny

=head1 AUTHOR

Tom Walsh (walshtp@gmail.com)

=head1 COPYRIGHT AND LICENSE

Copyright 2010 Tom Walsh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 1, or (at your option) any
later version.

=cut
