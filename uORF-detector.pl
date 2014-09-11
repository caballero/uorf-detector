#!/usr/bin/perl

use strict;
use warnings;

my $debug  = 0;
my %nuc2aa = (
'TCA' => 'S', # Serine
'TCC' => 'S', # Serine
'TCG' => 'S', # Serine
'TCT' => 'S', # Serine
'TTC' => 'F', # Phenylalanine
'TTT' => 'F', # Phenylalanine
'TTA' => 'L', # Leucine
'TTG' => 'L', # Leucine
'TAC' => 'Y', # Tyrosine
'TAT' => 'Y', # Tyrosine
'TAA' => 'X', # Stop
'TAG' => 'X', # Stop
'TGC' => 'C', # Cysteine
'TGT' => 'C', # Cysteine
'TGA' => 'X', # Stop
'TGG' => 'W', # Tryptophan
'CTA' => 'L', # Leucine
'CTC' => 'L', # Leucine
'CTG' => 'L', # Leucine
'CTT' => 'L', # Leucine
'CCA' => 'P', # Proline
'CAT' => 'H', # Histidine
'CAA' => 'Q', # Glutamine
'CAG' => 'Q', # Glutamine
'CGA' => 'R', # Arginine
'CGC' => 'R', # Arginine
'CGG' => 'R', # Arginine
'CGT' => 'R', # Arginine
'ATA' => 'I', # Isoleucine
'ATC' => 'I', # Isoleucine
'ATT' => 'I', # Isoleucine
'ATG' => 'M', # Methionine
'ACA' => 'T', # Threonine
'ACC' => 'T', # Threonine
'ACG' => 'T', # Threonine
'ACT' => 'T', # Threonine
'AAC' => 'N', # Asparagine
'AAT' => 'N', # Asparagine
'AAA' => 'K', # Lysine
'AAG' => 'K', # Lysine
'AGC' => 'S', # Serine
'AGT' => 'S', # Serine
'AGA' => 'R', # Arginine
'AGG' => 'R', # Arginine
'CCC' => 'P', # Proline
'CCG' => 'P', # Proline
'CCT' => 'P', # Proline
'CAC' => 'H', # Histidine
'GTA' => 'V', # Valine
'GTC' => 'V', # Valine
'GTG' => 'V', # Valine
'GTT' => 'V', # Valine
'GCA' => 'A', # Alanine
'GCC' => 'A', # Alanine
'GCG' => 'A', # Alanine
'GCT' => 'A', # Alanine
'GAC' => 'D', # Aspartic Acid
'GAT' => 'D', # Aspartic Acid
'GAA' => 'E', # Glutamic Acid
'GAG' => 'E', # Glutamic Acid
'GGA' => 'G', # Glycine
'GGC' => 'G', # Glycine
'GGG' => 'G', # Glycine
'GGT' => 'G'  # Glycine
);

$/ = "\n>"; # Fasta-slurp mode
while (<>) {
    s/>//g;
    my ($id, @seq)      = split (/\n/, $_);
    my $seq             = uc(join "", @seq);
    warn "ID=$id\nSEQ=$seq\n" if ($debug == 1);
    my @frames          = getFrames($seq);
    warn "F+1=$frames[0]\n" if ($debug == 1);
    warn "F+2=$frames[1]\n" if ($debug == 1);
    warn "F+3=$frames[2]\n" if ($debug == 1);
    warn "F-1=$frames[3]\n" if ($debug == 1);
    warn "F-2=$frames[4]\n" if ($debug == 1);
    warn "F-3=$frames[5]\n" if ($debug == 1);
    my $bestFrame       = getBestFrame(@frames);
    warn "BEST=$bestFrame\n" if ($debug == 1);
    my ($nf, $pos, $aa) = split (/:/, $bestFrame);
    next unless ($aa =~ /^M/);
    my $uorf            = checkUORF($nf, $pos, @frames);
    $uorf =~ tr/X/*/;
    warn "uORF=$uorf\n" if ($debug == 1);
    my $dir = 'F';
    if ($nf >= 3) {
        $nf -= 3;
        $dir = 'R';
    }
    $pos = (3 * $pos) + $nf;
    unless ($uorf eq "NO") {
        print "$id\t$uorf\t$pos\t$dir\t$aa\t$seq\n";
    }
}

sub getFrames {
    my $seq = shift @_;
    my $rc  = revcomp ($seq);
    my @frm = ();
    push  (@frm, translate($seq));             # +1
    push  (@frm, translate(substr ($seq, 1))); # +2
    push  (@frm, translate(substr ($seq, 2))); # +3
    push  (@frm, translate($rc));              # -1
    push  (@frm, translate(substr ($rc,  1))); # -2
    push  (@frm, translate(substr ($rc,  2))); # -3
    return @frm;
}

sub revcomp {
    my $seq = shift @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return reverse $seq;
}

sub translate {
    my $seq = shift @_;
    my $aa  = '';
    while ($seq) {
        my $codon = substr ($seq, 0, 3);
        last unless ((length $codon) == 3);
        if (defined $nuc2aa{ $codon }) { 
            $aa .= $nuc2aa{ $codon };
        } 
        else { 
            $aa .= '?'; 
        }
        substr ($seq, 0, 3) = '';
    }
    return $aa;
}

sub getBestFrame {
    my @frames  = @_;
    my $longest = -1;
    my $best    = '';
    my $nframe  =  0;
    foreach my $frame (@frames) {
        my @frag = split (/X/, $frame);
        my $pos  =  0;
        foreach my $frag (@frag) {
            my $len = length $frag;
            if ($len > $longest) {
                $longest = $len;
                $best    = "$nframe:$pos:$frag";
            }
            $pos += $len + 1;
        }
        $nframe++;
    }
    return $best;
}

sub checkUORF {
    my $nframe     = shift @_;
    my $maxpos     = shift @_;
    my @frames     = @_;
    my $uorf       = "NO";
    my @checkFrame = ();
    if ($nframe <= 2) { # +1, +2, +3
        @checkFrame = (0, 1, 2);
    }
    else {              # -1, -2, -3
        @checkFrame = (3, 4, 5);
    }
    
    foreach my $nf (@checkFrame) {
        my $frame   = $frames[$nf];
        my $upframe = substr ($frame, 0, $maxpos);
        while ($upframe =~ /M[^X]+?X/g) {
            my $orf = $&;
            my $len = length $orf;
            $uorf .= ":$orf" if ($len >= 5 and $len <= 50);
        }
    }
    if ($uorf =~ /:/) {
        $uorf =~ s/^NO://;
    }
    return $uorf;
}
