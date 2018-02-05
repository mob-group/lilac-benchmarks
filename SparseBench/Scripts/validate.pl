#
# Perl script for comparing test results against reference results
#
my %ref;
open(REFERENCE,"reference_results") || (die "could not open reference file");
while (<REFERENCE>) {
    chomp;
    my ($y,$s,$f,$p,$i,$r) = split;
    $ref{"${y}_${s}_${f}_${p}"} = "$i $r";
}
close(REFERENCE);

my $test = 0; my $error = 0;
while (<STDIN>) {
    chomp;
    my ($y,$s,$f,$p,$i,$r) = split;
    my $ir = $ref{"${y}_${s}_${f}_${p}"};
    if (!$ir) {
	print "There is no reference data for $y $s $f $p\n";
	next;}
    my ($ii,$rr) = split ' ',$ir;
    if ($i!=$ii) {
	print STDERR "Test $y $s $f $p: wrong number of iterations ($i s/b $ii)\n";
	$error++;
    } elsif (abs(($r-$rr)/$rr)>.05) {
	print STDERR "Test $y $s $f $p: residual off by > 5% ($r s/b $rr)\n";
	$error++;
    }
}
if ($error==0) {
    print STDERR "All tests passed\n";
} else {
    print STDERR "Number of errors: $error\n";
}
