#!/usr/bin/perl -w

use strict;

my %H=();

# INPUT: space separated
# 16059402 rs36180591 A G 22_16059402_G_A
# 16059425 rs1462228124 G A,T 22_16059425_A_G
# 16059484 rs1367478724 C T 22_16059484_T_C


sub common_prefix_length{
    my $s1=shift;
    my $s2=shift;
    my $m=length($s1)>length($s2) ? length($s2) : length($s1);
    my $L=0;
    for (my $i=0;$i<$m;$i++){
	if (substr($s1,$i,1) eq substr($s2,$i,1)){
	    $L+=1;
	}
	else{
	    last;
	}
    }
    return $L;
}

sub common_prefix_length{
    my $s1=shift;
    my $s2=shift;
    my $m=length($s1)>length($s2) ? length($s2) : length($s1);
    my $L=0;
    for (my $i=0;$i<$m;$i++){
	if (substr($s1,$i,1) eq substr($s2,$i,1)){
	    $L+=1;
	}
	else{
	    last;
	}
    }
    return $L;
}

sub reduce_strings{
    my $s1=shift;
    my $s2=shift;

    my $l1=common_prefix_length($s1,$s2);
    $s1=substr($s1,$l1);
    $s2=substr($s2,$l1);

    if (length($s1)==0){
	return ("-",$s2);
    }elsif(length($s2)==0){
	return ($s1,"-");
    }
    
    my $l2=common_prefix_length(scalar reverse($s1),scalar reverse($s2));
    $s1=substr($s1,0,length($s1)-$l2);
    $s2=substr($s2,0,length($s2)-$l2);
    
    if (length($s1)==0){
	return ("-",$s2);
    }elsif(length($s2)==0){
	return ($s1,"-");
    }else{
	return ($s1,$s2);
    }
}

while(<STDIN>){
    chomp;
    my @a=split(/ /);
    my @b=split(/_/,$a[4]);
    my $a1=$b[2];
    my $a2=$b[3];
    my $ref=$a[2];
    my $rs=$a[1];
    my $id=$a[4];
    my $id2=$id;
    @alts=split(/,/,$a[3]);
    
    # remove common prefix/suffix
    my ($s1,$s2)=reduce_strings($a1,$a2);

    for my $a (@alts){
	my ($t1,$t2)=reduce_strings($ref,$a);
	if ((($s1 eq $t1) && ($s2 eq $t2)) || (($s1 eq $t2) && ($s2 eq $t1))){
	    $id2=$rs;
	    last;
	}

    }

    $H{$id}=$id2;
}

$,="\t";
foreach my $id (keys %H){
    print $id,$H{$id};
}
