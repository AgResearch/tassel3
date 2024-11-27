#!/usr/bin/perl -w

use strict;
my $top = '.';
my $libdir = "$top/lib";
opendir (DIR, "$libdir") || die "Could not open $libdir\n";
my @list = readdir(DIR);

my @fl = ();
foreach my $fn(@list){
   if ("$fn" =~ m/\.jar$/){
      push(@fl, "$libdir\/$fn");
   }
}
push(@fl, "$top/dist/sTASSELGBS.jar");
my $CP = join(":", @fl);
my @args = @ARGV;

exec("java -classpath '$CP' -Xms512m -Xmx512m net.maizegenetics.pipeline.TasselPipeline @args") or die("failed to exec java");
