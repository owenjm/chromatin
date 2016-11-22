#!/usr/bin/perl -w

# chromatin_transitions.pl
# Copyright © 2015-16 Owen Marshall

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 
# USA

use strict;
use File::Basename;
use 5.010;

my @cli;
my @in_files;

my %vars = (
	'nb' => "~/genomes/chromatin_dm6_nbs.dissected/current/",
	'imm' => "~/genomes/chromatin_dm6_immature_neurons/revised/",
	'mat' => "~/genomes/chromatin_dm6_mature_neurons/revised/",
	'arrow_width' => 15,
	'rscript' => "~/Dropbox/R/chromatin.transitions.final.vii.r",
	'process_dir' => "",
	'no_default_trans' => "",
);

my %vars_details = (
	'nb' => 'nb dir',
	'imm' => 'immature neuron dir',
	'mat' => 'mature neuron dir',
	'arrow_width' => 'arrow width',
	'rscript' => 'Rscript for transitions to call',
	'process_dir' => 'Process all txt lists of genes in this folder',
	'no_default_trans' => 'Only process the process_dir files',
);

# Home directory
my $HOME = (getpwuid($<))[7];

process_cli();

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $date = ($year+1900)."-".($mon+1)."-$mday.$hour:$min:$sec.arrow.width=$vars{'arrow_width'}";
switchDir($date);

my $nb_exp = expressed_genes($vars{'nb'});
my $imm_exp = expressed_genes($vars{'imm'});
my $mat_exp = expressed_genes($vars{'mat'});

unless ($vars{'no_default_trans'}) {
	# All (pc-separate)
	switchDir("all");
		plot_trans(NAME=>'All');
		chdir('..');
		
	switchDir("all.expressed");
		my $expr =  union($nb_exp,$imm_exp,$mat_exp);
		writegenes(@$expr);
		plot_trans(GENES=>'tmp.genes.txt',NAME=>"all.expressed");
		chdir('..');
	
	switchDir("all.turned-on");
		my $expr2 =  intersect(
						excl(gene_list($vars{'nb'}."all.analysed.genes.txt"), $nb_exp ),
						union($imm_exp,$mat_exp)
						);
		writegenes(@$expr2);
		plot_trans(GENES=>'tmp.genes.txt',NAME=>"all.turned-on");
		chdir('..');
	
	# PcG only
	switchDir("pc-only");
		my $pc_genes = union(genes('PcG',$vars{'nb'}), genes('PcG',$vars{'imm'}), genes('PcG',$vars{'mat'}));
		writegenes(@$pc_genes);
		plot_trans(GENES=>'tmp.genes.txt',NAME=>'PcG_only');
		chdir('..');
	
	# TrxG to Repressed PcG separate
	switchDir("trxG-to-repressed.pcg-separate");
		my $trxGenesSep = intersect(
								 genes('TrxG*permissive',$vars{'nb'}),
								 union(
									  genes('HP1',$vars{'mat'}),
									  genes('PcG*rep',$vars{'mat'}),
									  genes('TrxG*rep',$vars{'mat'}),
									  genes('Null*rep',$vars{'mat'})
									)
								 );
		
		writegenes(@$trxGenesSep);
		plot_trans(GENES=>'tmp.genes.txt',NAME=>'TrxG_to_Repressive');
		chdir('..');
}

# process a list of files in folders
if ($vars{'process_dir'}) {
	
	my @file_list = glob($vars{'process_dir'}."/*.txt");
	foreach my $file (@file_list) {
		my ($fhead) = fileparse($file, qr/\.[^.]*/);
		
		switchDir($fhead);
		plot_trans(GENES=>$file,NAME=>$fhead);
		chdir('..');
		
		switchDir("$fhead.expressed");
		my $expr = intersect(gene_list($file),
							 union($nb_exp,$imm_exp,$mat_exp)
							 );
		writegenes(@$expr);
		plot_trans(GENES=>'tmp.genes.txt',NAME=>"$fhead.expressed");
		chdir('..');
		
		switchDir("$fhead.turned-on");
		my $expr2 = intersect(gene_list($file),
							 intersect(
								excl(gene_list($vars{'nb'}."all.analysed.genes.txt"), $nb_exp ),
								union($imm_exp,$mat_exp)
								)
							 );
		writegenes(@$expr2);
		plot_trans(GENES=>'tmp.genes.txt',NAME=>"$fhead.turned-on");
		chdir('..');
	}
}

chdir('..');

print STDOUT "All done.\n\n";
exit 0;


sub switchDir {
	my $dir = shift;
	mkdir($dir);
	chdir($dir);
}

sub genes {
	my ($p, $f) = @_;
	my @files = glob("$f/simp.vit*$p*genes.txt");
	
	my @list;
	foreach my $file (@files) {
		open(my $fh, '<', $file) or die "$!\n";
		foreach (<$fh>) {
			next if m/^#/;
			chomp;
			s/\s//g;
			push @list, $_;
		}
		close $fh;
	}
	
	return \@list;
}

sub expressed_genes {
	my ($f) = @_;
	
	print "$f\n";
	my ($file) = glob("$f/*expressed.genes.txt");
	
	my @list;
	open(my $fh, '<', $file) or die "$!\n";
	foreach (<$fh>) {
		next if m/^#/;
		chomp;
		s/\s//g;
		push @list, $_;
	}
	close $fh;
	
	return \@list;
}

sub gene_list {
	my ($file) = @_;
	($file) = glob($file);
	
	my @list;
	open(my $fh, '<', $file) or die "$!\n";
	foreach (<$fh>) {
		next if m/^#/;
		chomp;
		s/\s//g;
		push @list, $_;
	}
	close $fh;
	
	return \@list;
}

sub writegenes {
	my @names = @_;
	open(my $fh, '>', 'tmp.genes.txt') || die;
	print $fh "$_\n" foreach @names;
	close $fh;
}

sub union {
	my @arefs = @_;
	my @list;
	foreach my $ar (@arefs) {
		push @list, @$ar;
	}
	
	my @ret = uniq(@list);
	return \@ret;
}

sub uniq {
	return unless @_;
	my %seen;
	return grep { !$seen{$_}++ } @_;
}

sub intersect {
	my @arefs = @_;
	my %list;
	foreach my $ar (@arefs) {
		map { $list{$_}++ } @$ar;
	}
		
	my @intersect = grep {$list{$_} == @arefs} keys %list;
	return(\@intersect)
}

sub excl {
	# in first but not in second
	my @arefs = @_;
	my %list;
	
	map { $list{$_}++ } @{$arefs[0]};
	map { delete($list{$_}) } @{$arefs[1]};
	
	my @ret = keys %list;
	return \@ret;
}

sub plot_trans {
	my (%opts) = @_;
	print STDERR "Now plotting $opts{NAME} ...\n";
	
	my $genes_opt = $opts{GENES} || '';

	my $nb = $vars{'nb'};
	my $imm = $vars{'imm'};
	my $mat = $vars{'mat'};
		
	$genes_opt &&= "--genes=$genes_opt";
	
	my $rout = `Rscript $vars{'rscript'} --arrow.width=$vars{'arrow_width'} $genes_opt --n1=$nb --n2=$imm --n3=$mat`;
	
	`pdfjam *nb*pdf *.mat*pdf --nup 2x1 --landscape --outfile $opts{NAME}.pdf`
}


sub process_cli {
	# CLI processing
	foreach (@ARGV) {
		if (/--(.*)=(.*)/) {
			unless (defined($vars{$1})) {
				print STDERR "Did not understand $_ ...\n";
				help();
			}
			my ($v, $opt) = ($1,$2);
			$opt =~ s/~/$HOME/;
			$vars{$v} = $opt;
			push @cli, "$v=$opt";
			next;
		} elsif (/--h[elp]*/) {
			help();
		} elsif (/--(.*)/) {
			# if no parameter is specified we assume it's a switch ...
			# (could be a bit nicer and check this is ok with a hash representing data type ...)
			if (defined($vars{$1})) {
				$vars{$1} = 1;
			} else {
				print STDERR "Did not understand $_ ...\n";
				help();
			}
			push @cli, "$1";
			next;
		}
		push @in_files, $_;
	}
}


sub help {
	print STDERR <<EOT;

chromatin_transitions.pl
	
Options:
EOT
	
	my $opt_len = 0;
	foreach (keys %vars) {
		my $l = length($_);
		#print "--> $_: $l\n";
		$opt_len = $l if $l > $opt_len;
	}
	
	$opt_len+=2;
	
	my $cols= `tput cols` || 80;
	
	my ($v, $val, $def, $def_format);
	my $help_format = "format STDOUT =\n"
		.' '.'^'.'<'x$opt_len . ' '. '^' . '<'x($cols-$opt_len-4) . "\n"
		.'$v, $def_format'."\n"
		.' '.'^'.'<'x$opt_len . '   '. '^' . '<'x($cols-$opt_len-6) . "~~\n"
		.'$v, $def_format'."\n"
		.".\n";
		
	eval $help_format;
	die $@ if $@;
	
	foreach my $k (sort (keys %vars)) {
		($v, $val, $def) = ($k, $vars{$k}, $vars_details{$k});
		$def||="";
		$def_format = $val ? "$def\n\r[Current value: $val]" : $def;
		$v = "--$v";
#		format =
# ^<<<<<<<<<<<<<<<<<<<< ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#$v, $def_format
# ^<<<<<<<<<<<<<<<<<<<<   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ~~
#$v, $def_format
#.

		write();
		
	}
	print STDOUT "\n";
	exit 1;
}