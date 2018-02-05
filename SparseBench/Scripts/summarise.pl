#!/usr/bin/perl

use strict;

my $mach = <STDIN>; chomp $mach; my $opt = <STDIN>; chomp $opt;
print "==== Performance summary for machine $mach, variant $opt ====\n";
my %speed; my ($problem,$component);

#
# read performance records and process them
# by maximising components over all problems
# in which they occur.
# each record is < problem.component top-rate asymptotic-rate >
#
while (<STDIN>) {
    chomp; s/ +/ /g; s/^ //;
    my ($id,$top,$asym) = split;
    ($problem,$component) = split(/\./,$id);

    # overall speed
    SetMax($top,"allmax");
    SetMax($asym,"allasm");

    # matrix vector product
    SetMax($asym,"mvp-reg") if $component =~ /mvp-reg/;
    SetMax($asym,"mvp-crs-symm")
	if $component =~ /mvp-crs/ && $problem =~ /cg/;
    SetMax($asym,"mvp-crs-full")
	if $component =~ /mvp-crs/ && $problem =~ /gmres/;

    # ilu preconditioning
    SetMax($asym,"ilu-reg") if $component =~ /ilu-reg/;
    SetMax($asym,"ilu-crs-symm")
	if $component =~ /ilu-crs/ && $problem =~ /cg/;
    SetMax($asym,"ilu-crs-full")
	if $component =~ /ilu-crs/ && $problem =~ /gmres/;

    # block jacobi preconditioning
    SetMax($asym,"bjac") if $component =~ /bjac/;

    # basic, storage-independent, operations
    SetMax($asym,"gvec") if $component =~ /vec/;
    SetMax($asym,"cg-vec")
	if $component =~ /vec/ && $problem =~ /cg/;
    SetMax($asym,"gmres-vec")
	if $component =~ /vec/ && $problem =~ /gmres/;
}

PrintMax("Maximum speed attained","allmax");
PrintMax("Regular mvp asymptotic","mvp-reg");
PrintMax("Symm-crs mvp asymptotic","mvp-crs-symm");
PrintMax("Full-crs mvp asymptotic","mvp-crs-full");
PrintMax("Block jacobi asymptotic","bjac");
PrintMax("Regular ilu asymptotic","ilu-reg");
PrintMax("Symm-crs ilu asymptotic","ilu-crs-symm");
PrintMax("Full-crs ilu asymptotic","ilu-crs-full");
PrintMax("Vectors ops general","gvec");
PrintMax("Vectors ops cg","cg-vec");
PrintMax("Vectors ops gmres","gmres-vec");

sub SetMax {
    my ($top,$component) = @_;
    return if $top==0;
    if ($top>$speed{$component}{maxval}) {
	$speed{$component}{maxval} = $top;
	$speed{$component}{maxat} = $problem;}
    if (!defined $speed{$component}{minval}
	|| $top<$speed{$component}{minval}) {
	$speed{$component}{minval} = $top;
	$speed{$component}{minat} = $problem;}}
sub PrintMax {
    my ($msg,$comp) = @_;
    my ($s,$ss) = (int($speed{$comp}{minval}),int($speed{$comp}{maxval}));
    print "$msg: ($speed{$comp}{minat}) $s -- ($speed{$comp}{maxat}) $ss\n";}
