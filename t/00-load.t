#!perl -T

use Test::More tests => 7;

BEGIN {
    use_ok( 'BioStudio::Basic' ) || print "Bail out!\n";
    use_ok( 'BioStudio::GFF3' ) || print "Bail out!\n";
    use_ok( 'BioStudio::MySQL' ) || print "Bail out!\n";
    use_ok( 'BioStudio::GBrowse' ) || print "Bail out!\n";
    use_ok( 'BioStudio::Diff' ) || print "Bail out!\n";
    use_ok( 'BioStudio::Foswiki' ) || print "Bail out!\n";
    use_ok( 'BioStudio::BLAST' ) || print "Bail out!\n";
}

diag( "Testing BioStudio::Basic $BioStudio::Basic::VERSION, Perl $], $^X" );
