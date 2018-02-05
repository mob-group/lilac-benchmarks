#!/usr/bin/perl

use strict;

if ($#ARGV<0) {
    print "ERROR analyze.pl called without directory argument\n"; exit;
} else {
    my $dir = $ARGV[0]; my ($mvp_max,$pre_max,$vec_max,$all_max);
    my %mvp; my %pre; my %vec; my %all;
    opendir(DIR,$dir) || die "Could not open directory $dir\n";
    my @allfiles = grep /out/, readdir DIR;
    #
    # Go through all files and store speed of components,
    # indexed by method/storage/preconditioner/size
    #
    for my $file (@allfiles) {
	my ($meth,$stor,$prec,$size) = split(/\-|\./,$file);
	open(FILE,"$dir/$file") || die "Could not open file $file\n";
	#print "Now processing file $file\n";
	my ($m,$p,$v,$o,$msize) = (0,0,0,0,0);
	while (<FILE>) {
	    chomp; s/ +/ /g; s/^ //;
	    my ($w1,$w2,$w3,$w4) = split(/ /,$_);
	    # calculate size in .1Mbyte increments; that should
	    # be enough not to get duplicates.
	    if (/Integer/) {$msize = $w4;}
	    if (/Real/) {$msize = $msize+$w4; $size=int(10*$msize)/10;}
	    if ($w3=~/\*/) {$m=0;$p=0;$v=0;$o=0;next;}
	    if ($m==1) {$mvp{$meth}{$stor}{$prec}{$size} = $w3; $m=0;}
	    $m = 1 if /Matrix/;
	    if ($p==1) {$pre{$meth}{$stor}{$prec}{$size} = $w3; $p=0;}
	    $p = 1 if /Preconditioner/;
	    if ($v==1) {$vec{$meth}{$stor}{$prec}{$size} = $w3; $v=0;}
	    $v = 1 if /Vector/;
	    if ($o==1) {$all{$meth}{$stor}{$prec}{$size} = $w3; $o=0;}
	    $o = 1 if /Overall/;
	} close FILE or die "Could not close file\n";
    } closedir DIR or die "Could not close directory\n";
    #
    # Now generate plots for the components
    #
    GeneratePlot("all","","","",\%all);
    GeneratePlot("cg",   "cg",   "","",\%all);
    GeneratePlot("gmres","gmres","","",\%all);
    GeneratePlot("mvp","","","",\%mvp);
    GeneratePlot("mvp-reg","","reg","",\%mvp);
    GeneratePlot("mvp-crs","","crs","",\%mvp);
    GeneratePlot("ilu-reg","","reg","ilu",\%pre);
    GeneratePlot("ilu-crs","","crs","ilu",\%pre);
}

sub GeneratePlot {
    my ($graph,$meth_match,$stor_match,$prec_match,$dataref) = @_;
    open PLOT,">.plot/$graph.plot"
	or die "Could not open plot file $graph.plot for write\n";
    print ">> Generating plot for [$graph]\n";
    print PLOT
	"set terminal postscript eps\nset multiplot\nset origin 0,0\n";
    my %speed; my ($minx,$maxx,$miny,$maxy) = (0,0,0,0);
    #
    # condense from meth/stor/prec keys to one compound problem key
    #
    for my $meth ( keys %$dataref ) {
	next if $meth_match ne "" && $meth !~ $meth_match;
	for my $stor (keys %{$$dataref{$meth}}) {
	    next if $stor_match ne "" && $stor !~ $stor_match;
	    for my $prec (keys %{$$dataref{$meth}{$stor}}) {
		next if $prec_match ne "" && $prec !~ $prec_match;
		for my $size (keys %{$$dataref{$meth}{$stor}{$prec}}) {
		    $speed{"$meth-$prec-$stor"}{$size} =
			$$dataref{$meth}{$stor}{$prec}{$size};}}}}
    #
    # for each problem, write the graph sorted to a data file
    #
    for my $prob (keys %speed) {
	my @x = sort {$a<=>$b} keys %{$speed{$prob}};
	if ($minx==0 || $x[0]<$minx) {$minx = $x[0];}
	if ($x[$#x]>$maxx) {$maxx = $x[$#x];}
	open DAT,">.plot/$prob-$graph.dat"
	    or die "Could not open data file plot/$prob-$graph.dat\n";
	for my $size (@x) {
	    my $speed = $speed{$prob}{$size};
	    print DAT "$size $speed\n";
	    if ($miny==0 || $speed<$miny) {$miny = $speed;}
	    if ($speed>$maxy) {$maxy = $speed;}}
	$speed{$prob}{max} = $speed{$prob}{$x[$#x]};
	close DAT or die "Could not close data file plot/$prob.dat\n";}
    #
    # generate the plot file that will graph the data files
    #
    $maxx = 1.1*$maxx; $maxy = 1.1*$maxy;
    print PLOT "set xrange [0:$maxx]\nset yrange [0:$maxy]\n";
    my $n=1;
    for my $prob (sort {$speed{$b}{max}<=>$speed{$a}{max}} keys %speed) {
	my $cury = $maxy*(15-$n)/15; my $curx = .8*$maxx;
	print PLOT "set key $curx,$cury\n";
	print PLOT
	    "plot \"$prob-$graph.dat\" smooth bezier notitle with lines $n\n";
	print PLOT
	    "plot \"$prob-$graph.dat\"  title \"$prob\" with points $n\n";
	$n++;}
    close PLOT or die "Could not close plot file plot/$graph.plot\n";}
