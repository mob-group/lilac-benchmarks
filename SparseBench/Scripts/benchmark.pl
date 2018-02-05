#
# Perl script for analysing and reporting benchmark results
#

my ($cg_prob,$gm_prob,$mf_prob,$mf_mem);
my $mf_max=0;
my $mv_reg=0; my $mv_crs=0; my $lu_reg=0; my $lu_crs=0;
my $cg_max=0; my $gm_max=0; 
my $bdir = <STDIN> ; chomp $bdir;
my $vector = <STDIN> ; chomp $vector;
{
    my ($x,$m,$s,$mt,$mm,$mp,$mv,$loc_max,$loc_mem,
	$asympt,$asympt_mvp,$asympt_prec,$g_e,$m_e,$p_e);
    my $prob = 0;
    while (<STDIN>) {
	my $x;
	chomp;
	$_ =~ s/  / /g; $_ =~ s/  / /g;
	if (/File:/) {
	    if ($prob ne 0) {
		($asympt,$asympt_mvp,$asympt_prec,$g_e,$m_e,$p_e)
		    = Asympt($prob);
		ProbPars($prob,$asympt,$asympt_mvp,$asympt_prec);
		OneProb($asympt,$asympt_mvp,$asympt_prec,$g_e,$m_e,$p_e,
			$loc_max,$loc_mem);
	    }
	    $overall{$prob} = {} ; $mvp{$prob} = {};
	    ($x,$prob) = split(/ /,$_);
	    $prob =~ s/.out//;
	    print STDOUT "Problem: $prob\n"; 
	    $loc_max=0;
	}
	if (/Size/) {
	    ($x,$s,$x,$m,$x,$x,$x,$mt,$x,$mm,$x,$mp,$x,$mv) = split(/ /,$_);
	    if ($mt > $mf_max) {$mf_max = $mt; $mf_mem = $m; $mf_prob = $prob;}
	    if ($mt > $loc_max) {$loc_max = $mt; $loc_mem = $m;}
	    $overall{$prob}{$m} = $mt;
	    $mvp{$prob}{$m} = $mm;
	    $prec{$prob}{$m} = $mp;
        }
    }
    ($asympt,$asympt_mvp,$asympt_prec,$g_e,$m_e,$p_e) = Asympt($prob);
    ProbPars($prob,$asympt,$asympt_mvp,$asympt_prec);
    OneProb($asympt,$asympt_mvp,$asympt_prec,$g_e,$m_e,$p_e,$loc_max,$loc_mem);
}

print "\n==== Performance summary ====\n";
print "Highest performance: $mf_max on problem $mf_prob; data size $mf_mem Mb.\n";
print "Highest asymptotic performance\n";
print ".. whole cg   : $cg_max on problem $cg_prob.\n";
print ".. whole gmres: $gm_max on problem $gm_prob.\n";
print ".. regular mvp: $mv_reg.\n";
print ".. crs mvp    : $mv_crs.\n";
print ".. regular ilu: $lu_reg.\n";
print ".. crs ilu    : $lu_crs.\n";
print "============================\n\n";
sub Asympt {
    my $prob = @_[0];
    ($asympt_overall_a,$asympt_overall_b,$overall_err)
	= AsymptAspect("overall",$prob);
    ($asympt_mvp_a,$asympt_mvp_b,$mvp_err) = AsymptAspect("mvp",$prob);
    ($asympt_prec_a,$asympt_prec_b,$prec_err) = AsymptAspect("prec",$prob);
    return ($asympt_overall_a,$asympt_mvp_a,$asympt_prec_a,
	    $overall_err,$mvp_err,$prec_err,
	    $asympt_overall_b,$asympt_mvp_b,$asympt_prec_b);
}
sub AsymptAspect {
    my ($asp,$prob) = @_;
    my @k = keys %{$ {$asp}{$prob}}; my $k=$#k+1;
    open(ASYMPT,"| $bdir/Scripts/lsq") || die "could not open pipe to lsq\n";
    print ASYMPT "$prob\n$asp\n$k\n$vector\n";
    for my $i (@k) {
	print ASYMPT "$i ${$asp}{$prob}{$i}\n";
    }
    close ASYMPT;
    open(ASYMPT,"< asympt.log") || die "no asympt log file?!\n";
    my $asympt = <ASYMPT>; chomp $asympt;
    $asympt =~ s/  / /g; $asympt =~ s/  / /g; $asympt =~ s/^ *//;
    my ($a,$b,$e) = split(/ /,$asympt);
    close ASYMPT;
    return ($a,$b,$e);
}
sub OneProb {
    my ($asympt,$asympt_mvp,$asympt_prec,$g_e,$m_e,$p_e,
	$loc_max,$loc_mem) = @_;
    print ".. maximum performance   : $loc_max for data set $loc_mem Mb\n";
    print ".. asymptotic overall performance: $asympt";
    if ($g_e) {print "*\n";} else {print "\n";}
    print ".. asymptotic mvp     performance: $asympt_mvp";
    if ($m_e) {print "*\n";} else {print "\n";}
    print ".. asymptotic prec    performance: $asympt_prec";
    if ($p_e) {print "*\n";} else {print "\n";}
    if ($g_e+$m_e+$p_e) {
	print "   * : This measurement was based on insufficient data points.\n";
	print "       Run a few more large test sizes if you want accurate data.\n";
    }
}
sub ProbPars {
    my ($prob,$mt,$mm,$mp) = @_;
    if ($prob=~/cg/) {
	if ($mt>$cg_max) {$cg_max=$mt; $cg_prob=$prob;}
    } else {
	if ($mt>$gm_max) {$gm_max=$mt; $gm_prob=$prob;}
    }
    if ($prob=~/reg/) {
	if ($mm>$mv_reg) {$mv_reg=$mm;}
	if ($prob=~/ilu/ && $mp>$lu_reg) {$lu_reg=$mp;}
    } else {
	if ($mm>$mv_crs) {$mv_crs=$mm;}
	if ($prob=~/ilu/ && $mp>$lu_crs) {$lu_crs=$mp;}
    }
}
