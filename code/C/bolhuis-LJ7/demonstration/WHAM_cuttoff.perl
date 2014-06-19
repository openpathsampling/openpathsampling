#!/usr/bin/perl

print "hello\n";

while (@ARGV) {
    $file = shift @ARGV;
    push(@textfiles,$file) if (-T $file && -e $file);
}

$firstround =1;
$nsets =20;
$isim=0;
foreach $textfile (@textfiles)
{
    $i=0;
    open(FILE,$textfile);
    for($j=0;$j<$nsets;$j++) {
	$max=0;
	for($i=0;$i<1000;$i++) {
	    $_=<FILE>;
	    @q=split(" ",$_);
	    $lambda[$isim][$i]=$q[0];
	    $hist[$isim][$i] =$q[1];
	    if ($q[1]> $max ) {
		$max =$q[1];
	    }
	}
	print "max = $max\n";
	for($i=0;$i<1000;$i++) {
	    if ($hist[$isim][$i] < 0.1* $max) {
		$hist[$isim][$i]=0;
	    }
	    $nt[$isim]+= $hist[$isim][$i];
	}
	
	$_=<FILE>;

	print "check if rule is zero".$_."check\n";
	print "$isim $nt[$isim]\n";	
	if($nt[$isim]> 0) {
	    $isim++;
	}
	if ($firstround==1) {
	  if ($isim==1) { 
	    $isim=0;
	    $firstround =0;
	    $nt[$isim]=0;
	  }
	}

	$ndata = 1000;;
    }
    close FILE;
}
    
$nsim=$isim;
print "nonzero replica's $nsim\n";
#$ndata=20000;

$ii=0;
open(FILE,"zvalues.dat");
while (<FILE>) {
    @q=split(" ",$_);
    $lnZ[$ii++] = $q[0];
}
close FILE;

#for($isim=0;$isim<$nsim;$isim++) {
#    $lnZ[$isim] =-1.5*$isim;
#}

for($index=0;$index<1000;$index++) {
    $sumhist[$index]=0;
    for ($k=0;$k<$nsim;$k++) {
	$sumhist[$index] += $hist[$k][$index];
    }
}

$dif=100;
$iteration=0;
while (($dif >0.00001) && ($iteration <1000)) {
    
    for($isim=0;$isim<$nsim;$isim++) {
	$newZi=0;	
	for($index=0;$index<1000;$index++) {
	    if ($hist[$isim][$index]>0) {
		$sum=0;
		for ($k=0;$k<$nsim;$k++) {
		    if ($hist[$k][$index]>0){
			$weight = 1;
		    } else {
			$weight =0;
		    }
		    $sum += $weight* $nt[$k]/exp($lnZ[$k]);
		}
		$newZi += $sumhist[$index]/$sum;
	    }
	}
	
#	print "new Zi =$newZi ";
	$lnZnew[$isim] = log($newZi);
#	print "new lnZi =$lnZnew[$isim]\n";
    }
    
    $dif=0;
    for ($k=0;$k<$nsim;$k++) {
#	print "new lnZi =$lnZnew[$k]\n";
	$dif += abs($lnZ[$k]- $lnZnew[$k]);
	$lnZ[$k]= $lnZnew[$k]-$lnZnew[0];
#	print "new lnZi =$lnZnew[$k], new it =$lnZ[$k]\n";
    }
    print "dif =$dif\n";
    $iteration++;
}	


#now create all  histograms

for ($k=0;$k<$nsim;$k++) {
   for($index=0;$index<1000;$index++) {
 	
       $finalhist[$k][$index] = $hist[$k][$index] *exp($lnZ[$k])/$nt[$k];
#	print " $lambda[$k][$index]  $finalhist[$k][$index] \n";
	}
}

open(OUTPUT,">allhists.dat");

for ($k=0;$k<$nsim;$k++) {
	for($index=0;$index<1000;$index++) {	

	print OUTPUT " $lambda[$k][$index]  $finalhist[$k][$index] \n";
	
    }
    print OUTPUT "\n";
	
}
close OUTPUT;


open(OUTPUT,">allloghists.dat");

for ($k=0;$k<$nsim;$k++) {
    for($index=0;$index<1000;$index++) {	

	
	if ($finalhist[$k][$index] >0)  {
	    $finalhist[$k][$index] = log($finalhist[$k][$index]);
	} else {
	    $finalhist[$k][$index] = -10;
	}
	print OUTPUT " $lambda[$k][$index]  $finalhist[$k][$index] \n";
	
    }
    print OUTPUT "\n";
	
}
close OUTPUT;


#now create final  histogram
for($index=0;$index<1000;$index++) {
    
    
    $sum=0;
    for ($k=0;$k<$nsim;$k++) {
	if ($hist[$k][$index]>0){
	    $weight = 1;
	} else {
	    $weight =0;
	}
	$sum += $weight* $nt[$k]/exp($lnZ[$k]); 	
	
    }
    if ($sum >0){
	$combinedhist[$index] = $sumhist[$index]/$sum;
    } else {
	print $integrand;
	$combinedhist[$index] =0;
    }
}



open(OUTPUT,">hist.dat");

$norm=0;
for($index=0;$index<1000;$index++) {	
    if($combinedhist[$index]>$norm) {
	$norm = $combinedhist[$index];
    }
}

for($index=0;$index<1000;$index++) {	
    $combinedhist[$index] /= $norm;
    print OUTPUT " $lambda[0][$index]  $combinedhist[$index] \n";
    
}
close OUTPUT;


open(OUTPUT,">lnhist.dat");

for($index=0;$index<1000;$index++) {	
    if ($combinedhist[$index] >0)  {
	$logh = log($combinedhist[$index]);
	print OUTPUT " $lambda[0][$index]  $logh \n";
    } 
}
close OUTPUT;

open(OUTPUT,">eqlnhist.dat");
$norm = $combinedhist[100];

for($index=0;$index<1000;$index++) {	
    $combinedhist[$index] /= $norm;
  if ($combinedhist[$index] >0)  {
	$logh = log($combinedhist[$index]);
	print OUTPUT " $lambda[0][$index]  $logh \n";
    } 

}






close OUTPUT;

$ii=0;

open(FILE,">zvalues.dat");
for($isim=0;$isim<$nsim;$isim++) {
    print FILE "$lnZ[$isim] \n";
}
close FILE;


open(OUTPUT,">maxheights.dat");

for ($k=0;$k<$nsim;$k++) {
    $maxh[$k]=0;
    for($index=0;$index<1000;$index++) {	
	if ($hist[$k][$index]> $maxh[$k] ) {
	    $maxh[$k] = $hist[$k][$index];
	    $maxl[$k] =$lambda[0][$index];
	}
    }
    print OUTPUT "$maxl[$k]  $maxh[$k] \n";
    
}
	
	
close OUTPUT
