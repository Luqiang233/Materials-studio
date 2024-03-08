#!perl

# Purpose: Calculate statistics (Length Angle Min/Max/Mean/Probability Distribution Function/Number) on the Hydrogen Bonds
#          in an MD trajectory. Apllicable to bothe AIMD and ff-MD. Put the results in a Study Tables.

use strict;
use warnings;
use MaterialsScript qw(:all);

my $hex0="0a204e4f54453a205479706520696e20746865207873642066696c65206e616d65206669727374210a0a"; 
my $hex1="63726561746564206279204855415355414e2d4d535f59414e472e0a";
my $hex2="436f6e7461637420576563686174204944202731333030353432373136302720666f72206d6f726520696e666f2e0a";
my $str1 = pack('H*', $hex0);
my $str2 = pack('H*', $hex1);
my $str3 = pack('H*', $hex2);
print $str2;
print $str3;
print $str1;

my $doc= $Documents{"file name.xtd"};

#Create a new Study Table for the results
my $statsDoc = Documents->New("HBondLengthStats.std");
my $statsDocA = Documents->New("HBondAngleStats.std");
my $evolution = Documents->New("HBondNumEvolution.std");
$evolution->ColumnHeading("A") = "Frame";
$evolution->ColumnHeading("B") = "NumHBonds";
my $simplified = Documents->New("HBonds.std");
my $PDF = Documents->New("PDF_HBondLength.std");
$PDF->ColumnHeading("A") = "Length";
$PDF->ColumnHeading("B") = "Probability";
my $PDFA = Documents->New("PDF_HBondAngle.std");
$PDFA->ColumnHeading("A") = "Angle";
$PDFA->ColumnHeading("B") = "Probability";



#Initialize variables for the stats calculations
my $totalLength = 0;
my $minLength   = 99999.9; #Arbitrary Big Num
my $maxLength   = 0;
my $row = 0;
my @Lengths;

my $totalAngle = 0;
my $minAngle   = 99999.9; #Arbitrary Big Num
my $maxAngle   = 0;
my $rowa = 0;
my @Angles;

my $totnum=0;
my $minnum=10000000; #Arbitrary Big Num
my $maxnum=0;
my @speciesList;
my $NumCl=0;
my $NumO=0;
my $NumF=0;
my $NumN=0;
my $NumS=0;
my $NumH=0;
@speciesList = ("Cl");#("N", "O", "S", "F", "Cl");
my $DonorMol = "C2 H6 O";

# To be changed, set bins
my $l2min=1;
my $l2max=2.4; # H bond Max Length
my $l2grid=0.02;
my $l2n=int(($l2max-$l2min)/$l2grid+0.001);
my @numlength=(0) x ($l2max/$l2grid+1);

my $a2min=0;
my $a2max=180;
my $a2grid=1;
my $a2n=int(($a2max-$a2min)/$a2grid+0.001);
my @numangle=(0) x ($a2max/$a2grid+1);

Tools->BondCalculation->ResetSettings;
Tools->BondCalculation->ChangeSettings(Settings(MaxHydrogenAcceptorDistance => "$l2max"));

COUNTDONOR();
sub COUNTDONOR{
    my $atoms = $doc->UnitCell->Atoms;
    foreach my $atom (@$atoms) {
        if ($atom->ElementSymbol eq "Cl"){
            $NumCl++;
        }
        if ($atom->ElementSymbol eq "N"){
            $NumN++;
        }
        if ($atom->ElementSymbol eq "F"){
            $NumF++;
        }
        if ($atom->ElementSymbol eq "S"){
            $NumS++;
        }
        if ($atom->ElementSymbol eq "O"){
            $NumO++;
        }
        if ($atom->ElementSymbol eq "H"){
            $NumH++;
        }
    }
    print "Number of atoms \n";
    print "O $NumO Cl $NumCl F $NumF N $NumN S $NumS H $NumH\n";

    #my $NumO = $doc->UnitCell->Sets("O")->Atoms->Count;
    #my $NumH = $doc->UnitCell->Sets("H")->Atoms->Count;
}

#DONORACCEPTER(); # HBonds formed by specified element, comment it out if not needed !!!!!!!!!!
sub DONORACCEPTER{
    #@speciesList = ("Cl");#("N", "O", "S", "F", "Cl");
    #Tools->BondCalculation->HBonds->ClearDonors;
    Tools->BondCalculation->HBonds->ClearAcceptors; 
    
    foreach my $element (@speciesList){
        Tools->BondCalculation->HBonds->AddDonor($element);
        Tools->BondCalculation->HBonds->AddAcceptor($element);
    }
    #Tools->BondCalculation->HBonds->Calculate($tmpxsd);

}



# sub main()
for (my $i=2; $i<=$doc->trajectory->NumFrames; ++$i){
    $doc->Trajectory->CurrentFrame = $i;
    my $tmpxsd= Documents->New("tmp_$i.xsd");
    $tmpxsd->CopyFrom($doc);
    $tmpxsd->CalculateHBonds;
    my $hbonds=$tmpxsd->UnitCell->HydrogenBonds;
    
    @Lengths=HBondStat($hbonds);
    @Angles=HBondAngleStat($hbonds);
    
    my $numbonds=$hbonds->Count;
    $totnum += $numbonds;
    $evolution->Cell($i-1, 0) = "$i";
    $evolution->Cell($i-1, 1) = "$numbonds";
    if($numbonds < $minnum) {
        $minnum = $numbonds;
    }
    elsif($numbonds > $maxnum) {
        $maxnum = $numbonds;
    }
    $tmpxsd->Delete;
}
#print $totnum;

my $numframes=$doc->trajectory->NumFrames;
$evolution->Cell($numframes,0)="AverageNumBonds";
$evolution->Cell($numframes,1)=$totnum/($numframes-1);
$evolution->Cell($numframes,2)="NumBondsMin";
$evolution->Cell($numframes,3)=$minnum;
$evolution->Cell($numframes,4)="NumBondsMax";
$evolution->Cell($numframes,5)=$maxnum;

$simplified->Cell(0,0)="AverageNumBonds";
$simplified->Cell(0,1)=$totnum/($numframes-1);
$simplified->Cell(0,2)="NumBondsMin";
$simplified->Cell(0,3)=$minnum;
$simplified->Cell(0,4)="NumBondsMax";
$simplified->Cell(0,5)=$maxnum;
$simplified->Cell(0,6)="CoordNumofH";
$simplified->Cell(0,7)=$totnum/($numframes-1)/$NumH;
$simplified->Cell(0,8)="CoordNumofO";
$simplified->Cell(0,9)=$totnum/($numframes-1)/$NumO;

$simplified->Cell(1,0)="AverageLength";
$simplified->Cell(1,1)=$Lengths[0];
$simplified->Cell(1,2)="MinLength";
$simplified->Cell(1,3)=$Lengths[1];
$simplified->Cell(1,4)="MaxLength";
$simplified->Cell(1,5)=$Lengths[2];

$simplified->Cell(2,0)="AverageAngle";
$simplified->Cell(2,1)=$Angles[0];
$simplified->Cell(2,2)="MinAngle";
$simplified->Cell(2,3)=$Angles[1];
$simplified->Cell(2,4)="MaxAngle";
$simplified->Cell(2,5)=$Angles[2];

#PRINTBINS;
#sub PDF{
    for (my $j= 0; $j<$l2n; ++$j){ 
        my $l2=$l2min + $l2grid*$j;
        my $temp2=$numlength[$j]/$row/0.02;
        $PDF->Cell($j,0)="$l2";
        $PDF->Cell($j,1)="$temp2";              
    }
    
    for (my $j= 0; $j<$a2n; ++$j){ 
        my $a2=$a2min + $a2grid*$j;
        my $temp3=$numangle[$j]/$rowa/0.02;
        $PDFA->Cell($j,0)="$a2";
        $PDFA->Cell($j,1)="$temp3";              
    }

#}

sub HBondStat {
    my ($hbonds) = @_;
    foreach my $hbond (@$hbonds) {
        #Output the bond length for each HBond
        $statsDoc->Cell($row, 0) = "HBond $row";
        my $l=$statsDoc->Cell($row, 1) = $hbond->Length;
        TOBINS($l);
        
        #Update the statistics information
        $totalLength += $hbond->Length;

        if($hbond->Length < $minLength) {
            $minLength = $hbond->Length;
        }
        elsif($hbond->Length > $maxLength) {
            $maxLength = $hbond->Length;
        }

        ++$row;
    }

    #printout the overall statistics
    $statsDoc->Cell($row, 0) = "Average";
    $statsDoc->Cell($row, 1) = $totalLength/$row;
    $statsDoc->Cell($row, 2) = "Min";
    $statsDoc->Cell($row, 3) = $minLength;
    $statsDoc->Cell($row, 4) = "Max";
    $statsDoc->Cell($row, 5) = $maxLength;
    return($totalLength/$row, $minLength, $maxLength)
}

sub HBondAngleStat {
    my ($hbonds) = @_;
    foreach my $hbond (@$hbonds) {
        #Output the bond length for each HBond
        $statsDocA->Cell($rowa, 0) = "HBond $rowa";
        my $angle=$statsDocA->Cell($rowa, 1) = $hbond->HBondAngle;
        TOBINSAngle($angle);
        #Update the statistics information
        $totalAngle += $hbond->HBondAngle;

        if($hbond->HBondAngle < $minAngle) {
            $minAngle = $hbond->HBondAngle;
        }
        elsif($hbond->HBondAngle > $maxAngle) {
            $maxAngle = $hbond->HBondAngle;
        }

        ++$rowa;
    }

    #printout the overall statistics
    $statsDoc->Cell($rowa, 0) = "Average";
    $statsDoc->Cell($rowa, 1) = $totalAngle/$rowa;
    $statsDoc->Cell($rowa, 2) = "Min";
    $statsDoc->Cell($rowa, 3) = $minAngle;
    $statsDoc->Cell($rowa, 4) = "Max";
    $statsDoc->Cell($rowa, 5) = $maxAngle;
    return($totalAngle/$rowa, $minAngle, $maxAngle)
}



sub TOBINS {
    my ($l) = @_;
    my $temp=int(($l-$l2min)/$l2grid+0.5);
    $numlength[$temp] += 1;
    #print "$temp $numlength[$temp]\n";
}

sub TOBINSAngle {
    my ($angle) = @_;
    my $temp4=int(($angle-$a2min)/$a2grid+0.5);
    $numangle[$temp4] += 1;
    #print "$temp4 $numlength[$temp]\n";
}

        