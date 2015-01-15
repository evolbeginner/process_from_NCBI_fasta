#! /usr/bin/perl

# last updated on 2014-05-08

###################################################################################
=pod

=head1 basic info

function: select sequences from a list of organisms and rename the sequences according to your request given the sequences from BLAST result

author:   Sishuo Wang from the University of British Columbia

feedback: If you have any questions and/or suggestions, you are very welcome to write e-mail to sishuowang@hotmail.ca or tomassonwss@gmail.com. Your feedback is highly appreciated.

=head1 usage

perl process_from_blast_seq usage.pl <--blast_seq_file=>

Options:

--rank=            name of taxonomic levels that are required to be present in the sequence title
                   legal names include superkingdom, kingdom (not recommended), phylum, class, order, family,
genus, species, full
                   e.g. for Saccharomyces cerevisiae S288C, the corresponding names are
                        Eukaryota, Fungi, Ascomycota,
                        Saccharomycetes, Saccharomycetales, Saccharomycetaceae,
                        Saccharomyces, Saccharomyces cerevisiae, Saccharomyces cerevisiae S288C

--species_info_file=the file containing taxonomic information for multiple species.
                   B<If this argument is specified, only sequences from the species
                   included in this reference file will be outputted.>
                   A reference file is given as example/reference_list_check_superkindom3.misspelling

--accn             if Accession number needs to be present in the sequence title

--gi               if GI needs to be present in the sequence title

--taxa             if taxa name needs to be present in the sequence title

--no_space         if space is replaced by _ in the taxa name
                   e.g. Saccharomyces cerevisiae S288C will be converted to
                        Saccharomyces_cerevisiae_S288C

--connector        connector of ACCN, GI and taxa name
                   default: -
                   e.g. for gene with GI 6319790, the sequence title will be
                        Saccharomyces cerevisiae-6319790-NP_009871.2

--taxa_name_included= select sequences based on specific taxa names

--taxa_include=    taxa which needs to be included (those not belonging to specified taxa will be excluded)
                   format: taxonomic level, taxonomic name
                   e.g. if --taxa_include superkingdom,Bacteria is specified,
                        only bacterial sequences will be included

--taxa_exclude=    taxa which needs to be excluded
                   format: taxonomic level, taxonomic name

--taxa_limit_file= the file indicating the taxa that will be included and/or excluded
                   an example is given as 'taxa_limit_file_example'

=cut

###################################################################################
use 5.010;
use strict;
use Getopt::Long;
use Cwd;
use autodie;

my (@rank, %rank_hash, $blast_seq_file, $species_info_file, $all, 
	$is_no_space, $taxa_limit_file, $is_no_subspecies_name, $is_help);
my $connector = '-';
my (%seq_title_order, %taxa_names_included, $taxa_names_included, @taxa_include, @taxa_exclude, %taxa_limit);
# ------------------------ #
my (%swi_hash, @swi_array);
my ($search_taxa_path, $search_taxa_prog, $search_taxa_db);
my (%key_order_hash, %taxa_rank_info, @seq_array);

&show_usage() if not @ARGV;
%seq_title_order = %{&get_seq_title_order(\@ARGV)};
#foreach(keys %seq_title_order){print $_."\t".$seq_title_order{$_}."\n";}
GetOptions(
	'rank=s'                =>	\@rank,
	'blast_seq_file=s'      =>	\$blast_seq_file,
	'species_info_file=s'   =>	\$species_info_file,
	'taxa|TAXA!'            =>	\$swi_hash{taxa},
	'accn|ACCN!'            =>	\$swi_hash{ACCN},
	'gi|GI!'                =>	\$swi_hash{GI},
	'all!'                  =>	\$all,
	'connector=s'           =>	\$connector,
	'no_space!'             =>	\$is_no_space,
    'taxa_name_included=s'  =>  \$taxa_names_included,
	'taxa_include=s'        =>	\@taxa_include,
	'taxa_exclude=s'        =>	\@taxa_exclude,
	'taxa_limit_file=s'	    =>	\$taxa_limit_file,
	'no_subspecies_name!'	=>	\$is_no_subspecies_name,
	'h|help!'               =>	\$is_help,
) || die "param error!\n";

map {delete $swi_hash{$_} if not $swi_hash{$_}} keys %swi_hash;
if ($all){ @swi_hash{qw(taxa ACCN GI)} = (1)x3}
@rank_hash{map {split /[,]/} @rank} = (1) x scalar(@rank);
&show_usage if ($is_help);

my $pwd = getcwd();
$search_taxa_path = "$pwd/tools" if not $search_taxa_path;
$search_taxa_prog  = "$search_taxa_path/search_taxa_info.pl";
$search_taxa_db   = "$search_taxa_path/taxa_info.db";

my @key_order = qw(full species genus family order class phylum kingdom superkingdom);
foreach (0..$#key_order){
	$key_order_hash{$key_order[$_]} = $#key_order-$_;
}

foreach (split(/[,]/,$taxa_names_included)){
    $taxa_names_included{$_}=1;
}


########################################################################
&get_taxa_limit(\@taxa_include, 'include');
&get_taxa_limit(\@taxa_exclude, 'exclude');
&read_taxa_limit_file($taxa_limit_file) if $taxa_limit_file;

%taxa_rank_info = %{&read_species_info_file($species_info_file, \%taxa_limit)} if $species_info_file;

@seq_array = @{&read_blast_seq_file($blast_seq_file)};
if ($swi_hash{'taxa'} and (not $species_info_file)){
	map {&get_taxa_rank_info_using_script($_->{'taxa'})} @seq_array;
}

=cut
@seq_array = grep {
	my $species_name='';
	($species_name) = $_->{taxa} =~ /(\S+\s\S+)/ if $is_no_subspecies_name;
	exists $taxa_rank_info{$_->{taxa}} or exists $taxa_rank_info{$species_name}
	} @seq_array if $species_info_file;
=cut


if ($species_info_file){
	for my $index (0..$#seq_array){
		my $href = $seq_array[$index];
		my ($species_name) = $href->{taxa} =~ /(\S+\s\S+)/ if $is_no_subspecies_name;
		if (exists $taxa_rank_info{$href->{taxa}}){
			next;
		}
		elsif (exists $taxa_rank_info{$species_name}){
			$seq_array[$index]{taxa} = $species_name;
			next;
		}
		delete $seq_array[$index];
	}
}


@seq_array = sort {$taxa_rank_info{$a->{taxa}}{superkingdom} cmp $taxa_rank_info{$b->{taxa}}{superkingdom} or
                   $taxa_rank_info{$a->{taxa}}{kingdom}      cmp $taxa_rank_info{$b->{taxa}}{kingdom}      or
                   $taxa_rank_info{$a->{taxa}}{phylum}       cmp $taxa_rank_info{$b->{taxa}}{phylum}       or
                   $taxa_rank_info{$a->{taxa}}{class}        cmp $taxa_rank_info{$b->{taxa}}{class}        or
                   $taxa_rank_info{$a->{taxa}}{order}        cmp $taxa_rank_info{$b->{taxa}}{order}        or
                   $taxa_rank_info{$a->{taxa}}{family}       cmp $taxa_rank_info{$b->{taxa}}{family}       or
                   $taxa_rank_info{$a->{taxa}}{genus}        cmp $taxa_rank_info{$b->{taxa}}{genus}        or
                   $taxa_rank_info{$a->{taxa}}{species}      cmp $taxa_rank_info{$b->{taxa}}{species}      or
                   $taxa_rank_info{$a->{taxa}}{full}         cmp $taxa_rank_info{$b->{taxa}}{full}      
                  } @seq_array;


foreach my $seq_info_href (@seq_array){
	my ($taxa_info_name, $full_taxa_name, $taxa_name, @tmp_taxa_array, $seq_title);
	$full_taxa_name = $seq_info_href->{taxa};
    #print $full_taxa_name."\n";next;
	my ($species_name) = $full_taxa_name =~ /(\S+\s\S+)/;
	$taxa_name = exists $taxa_rank_info{$full_taxa_name} ? $full_taxa_name : $species_name;
    do {next if ! exists $taxa_names_included{$taxa_name}} if %taxa_names_included;

	for (sort {$key_order_hash{$a} cmp $key_order_hash{$b}}
	grep {exists $rank_hash{$_}}
	keys %{$taxa_rank_info{$taxa_name}}){
        push @tmp_taxa_array, $taxa_rank_info{$taxa_name}{$_};
	}

	$taxa_info_name = join ('-', @tmp_taxa_array);
	#print "taxa-$taxa_name"."\n" if exists $taxa_rank_info{'Dictyostelium discoideum AX4'};
	$seq_info_href->{taxa} = $taxa_info_name;
	$seq_title = join ($connector, @$seq_info_href{grep {$_ ne 'seq' and $seq_info_href->{$_} ne ''} sort {$seq_title_order{$a} <=> $seq_title_order{$b}} keys %swi_hash});
    #print $seq_title."\n";
	$seq_title =~ s/ /_/g if $is_no_space;
    if ($seq_title !~ /NotFound/){
        do {print ">$seq_title\n"; print $seq_info_href->{seq}."\n";}
    }
}


########################################################################
sub read_taxa_limit_file{
	my ($mode);
	my ($taxa_limit_file) = @_;
	open(my $IN, '<', $taxa_limit_file);
	while(my $line = <$IN>){
		chomp($line);
		given ($line){
			when (/^TAXA_INCLUDE$/) {$mode='include'}
			when (/^TAXA_EXCLUDE$/) {$mode='exclude'}
		}
		&get_taxa_limit([$line], $mode);
	}
}

sub get_taxa_limit{
	my ($taxa_limit_aref, $mode) = @_;
	for (@$taxa_limit_aref){
		my ($rank_name, $taxonomy_name) = split /[,]/;
		$rank_name = lc($rank_name);
		$taxonomy_name = ucfirst($taxonomy_name);
		$taxa_limit{$mode}{$rank_name}{$taxonomy_name}=1;
	}
}

sub get_seq_title_order{
	my (%seq_title_order, $seq_title_order_k);
	my ($argv_aref) = @_;
	foreach(@$argv_aref){
		given($_){
			when(/^\-+(?:GI|gi)$/)	{$seq_title_order{'GI'}  =++$seq_title_order_k}
			when(/^\-+(?:accn|ACCN)$/)	{$seq_title_order{'ACCN'}=++$seq_title_order_k}
			when(/^\-+(taxa|TAXA)$/)	{$seq_title_order{'taxa'}=++$seq_title_order_k}
		}
	}
	return(\%seq_title_order)
}

sub get_taxa_rank_info_using_script{
	my ($taxa_name_new) = @_;
    #my ($species_name_new) = $taxa_name_new =~ /(\S+\s\S+)/;
	my $rank_line = `perl $search_taxa_prog --taxa "$taxa_name_new" --db_file "$search_taxa_db" --rank all --concise --Bioperl`;
	chomp ($rank_line);
	%taxa_rank_info = %{&generate_rank_info($rank_line)};
}

sub read_species_info_file{
	my (%taxa_rank_info);
	my ($species_info_file, $taxa_limit_href) = @_;
	open(IN, '<', $species_info_file);
	while(<IN>){
		chomp;
		%taxa_rank_info = %{&generate_rank_info($_, \%taxa_rank_info)};
	}
	return (\%taxa_rank_info);
}

sub generate_rank_info{
	my ($line, $taxa_rank_info_href) = @_;
	my %taxa_rank_info = defined $taxa_rank_info_href ? %$taxa_rank_info_href : %taxa_rank_info;
	my @line = split /\t/, $line;
	my $full_name = $line[0];
	my $species_name = $line[8];
	push @line, $full_name;
	for my $taxa_name ($full_name, $species_name){
		for my $mode (keys %taxa_limit){
			for my $rank_name (map {split /[,]/} @key_order){
				my $taxonomy_name = $line[$key_order_hash{$rank_name}+1];
				next if ! exists $taxa_limit{$mode}{$rank_name};
				if ($mode eq 'include'){
					goto endOfFunc if ! exists $taxa_limit{$mode}{$rank_name}{$taxonomy_name};
				}
				if ($mode eq 'exclude'){
					goto endOfFunc if   exists $taxa_limit{$mode}{$rank_name}{$taxonomy_name};
				}
			}
		}
		foreach (@key_order){
			$taxa_rank_info{$taxa_name}{$_} = $line[$key_order_hash{$_}+1];
		}
	}
	endOfFunc:
	return (\%taxa_rank_info);
}

sub read_blast_seq_file{
	my (@seq_array, $num);
	my ($GI, $ACCN, $taxa_name, $readseq_href);
	my ($blast_seq_file) = @_;
	open (my $IN, '<', $blast_seq_file) || die "blast_seq_file $blast_seq_file cannot be opened";
	while(my $line=<$IN>){
		next if $line =~ /^$/;
		chomp($line);
		# if (/^> gi \| (\d+) \| ref \| ([A-Z]P_[^\|]+) \| [^\[\]]+ \[([^\[\[]+)\]$/x){
		if ($line =~ /^> gi \| (\d+) \| [^\|]+ \| ([^\|]+) \| [^\[\]]+ \[([^\[\[]+)\]/x){
			($GI, $ACCN, $taxa_name) = ($1,$2,$3);
			my ($seq);
			push @seq_array, {};
			$num=$#seq_array;
			$seq_array[$num]{GI}=$GI;
			$seq_array[$num]{ACCN}=$ACCN;
			$seq_array[$num]{taxa}=$taxa_name;
		}
		else{
			$seq_array[$num]{seq} .= $line;
		}
	}
	return (\@seq_array);
}

sub show_usage{
	system "perldoc $0";
	exit 0;
}

