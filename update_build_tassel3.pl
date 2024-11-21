#! /usr/bin/perl -w

use strict;

system "cp dist/sTASSELGBS.jar ../tassel3.0_standalone/sTASSEL.jar";


my $versionComment3 = `./run_pipeline_non_verbose.pl -versionComment`;
$versionComment3 =~ s/\s+$//;
my $versionTag3 = `./run_pipeline_non_verbose.pl -versionTag`;
$versionTag3 =~ s/\s+$//;

if (!defined $versionComment3 || length $versionComment3 == 0) {
   print "versionComment3 is not defined\n";
   exit;
}

if (!defined $versionTag3 || length $versionTag3 == 0) {
   print "versionTag3 is not defined\n";
   exit;
}

my @args = @ARGV;
if ($args[0] eq "-test") {
   print "Test Mode...\n\n";
   print "Tassel 3 Comment: " . $versionComment3 . "\n";
   print "Tassel 3 Tag: " . $versionTag3 . "\n";
   system "cd ../tassel3.0_standalone; git commit -a -m'$versionComment3' --dry-run";
}
elsif ($args[0] eq "-commit") {
   print "Commit Mode...\n\n";
   print "Tassel 3 Comment: " . $versionComment3 . "\n";
   print "Tassel 3 Tag: " . $versionTag3 . "\n";
   system "cd ../tassel3.0_standalone; git commit -a -m'$versionComment3'";
   system "cd ../tassel3.0_standalone; git tag -m'$versionComment3' $versionTag3";
   system "cd ../tassel3.0_standalone; git push";
   system "cd ../tassel3.0_standalone; git push --tags";
}
