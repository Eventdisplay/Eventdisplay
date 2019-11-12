#!/usr/bin/perl

#called by make_combined_lped_file.sh: combines low gain ped files for all telescopes/runs into single file. 
#storage is by FADC board and channel, so board swaps are handled automatically

use strict;
use warnings;

unless($#ARGV == 1)
{
    #print $#ARGV,"\n";
    print "usage: mutate_lped.pl <run> <calibdir>\n";
    exit(-1);
}


my $lpedrun = $ARGV[0];
my $calibdir = $ARGV[1];

print "Mutating $lpedrun\n";

# PERL MODULE WE WILL BE USING
use DBI;

# MySQL CONFIG VARIABLES
# CONFIG VARIABLES
my $platform = 'mysql';
my $database = 'VERITAS';
my $host = 'remus.ucsc.edu';
my $user = 'readonly';
my $pw = '';
#DATA SOURCE NAME
my $dsn = "DBI:$platform:$database:$host";
# PERL MYSQL CONNECT
my $connect = DBI->connect($dsn, $user, $pw) || die "Could not connect to database: $DBI::errstr";


# first: get the date
my $date;
my $run = substr($lpedrun, 0, 5);
my $statement = qq{select db_start_time from tblRun_Info where run_id = $run };
 
my $query_handle = $connect->prepare($statement); 
$query_handle->execute() or die "\n ($DBI::err): $DBI::errstr\n"; 
$query_handle->bind_columns(undef,\$date);
die "Could not get date for run $run.\n"unless ($query_handle->fetch()) ;

print "Start date for run $run: $date\n";


#now: get FADC crate/slot/channel per pixel for each telescope
$statement = qq{select telescope_id, fadc_crate, fadc_slot, fadc_channel, channel_id from tblFADC_Channel_Relation where db_start_time < '$date' and ( db_end_time is NULL or db_end_time > '$date');};

my $t_tel;
my $t_crate;
my $t_slot;
my $t_channel;
my $t_channel_id;

my @fadc_crate;
my @fadc_slot;
my @fadc_channel;
my @slot_info;

$query_handle = $connect->prepare($statement); 
$query_handle->execute() or die "\n ($DBI::err): $DBI::errstr\n"; 
$query_handle->bind_columns(undef,\$t_tel, \$t_crate, \$t_slot, \$t_channel, \$t_channel_id);
while( $query_handle->fetch() ) {
	$fadc_crate[ $t_tel ][ $t_channel_id ] = $t_crate;
	$fadc_slot[ $t_tel ][ $t_channel_id ] = $t_slot;
	$fadc_channel[ $t_tel ][ $t_channel_id ] = $t_channel;
	$slot_info[ $t_tel ][ $t_channel_id ] = 1;
}

#now: get FADC module per tel/crate/slot.
$statement = qq{select telescope_id, fadc_crate, fadc_slot, fadc_id from tblFADC_Slot_Relation where db_start_time < '$date' and ( db_end_time is NULL or db_end_time > '$date');};

print $statement, "\n";

my $t_fadc;

my @fadc_module;
my @module_info;

my $nline;

$query_handle = $connect->prepare($statement); 
$query_handle->execute() or die "\n ($DBI::err): $DBI::errstr\n"; 
$query_handle->bind_columns(undef,\$t_tel, \$t_crate, \$t_slot, \$t_fadc );
while( $query_handle->fetch() ) {
	$fadc_module[ $t_tel ][ $t_crate ][ $t_slot ] = $t_fadc;
	$module_info[ $t_tel ][ $t_crate ][ $t_slot ] = 1;
	#print "$t_tel $t_crate $t_slot $t_fadc\n";
	$nline++;
}

print $nline, "\n";

#output file
my $outfilename = "$calibdir/lpedfiles/$lpedrun.lped";
open OUT, "+>$outfilename";

for( my $i_tel=1; $i_tel<5; $i_tel++) {
	my $infilename = "$calibdir/Tel_".$i_tel."/$lpedrun.lped";
	print "$infilename\n";
	open IN, $infilename or next;
	while (my $line = <IN>) {
		my @numbers = split ' ', $line;
		my $tel = $numbers[0]; 
		my $chan = $numbers[1];
		# print "$tel $chan\n";
		next unless ( $slot_info[ $tel ][ $chan ]==1 && $module_info[ $tel ][ $fadc_crate[ $tel ][ $chan ] ][ $fadc_slot[$tel][$chan] ] == 1 );
		my $module = $fadc_module[ $tel ][ $fadc_crate[ $tel ][ $chan ] ][ $fadc_slot[$tel][$chan] ] ;
		my $channel= $fadc_channel[$tel ][ $chan ];
		print OUT "$module $channel $lpedrun ", $line;
		# print "$module $channel $lpedrun ", $line;
	}
	close( IN ) ;

}



close( OUT );


exit;
