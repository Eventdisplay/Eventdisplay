#!/bin/awk -f
BEGIN {}  # begin section 
($3=="swhigh" ) { 
	lmult[$2,$4,$6]=$8;
 }        
END{
	sw[6]=sw[10]=sw[16]=sw[18]=1;
	v["V5"]=v["V6"]=1;

	for (isw in sw) {
		for (iv in v) {
			for (tel=1; tel<5; tel++) {
				D2=tel iv;
				ivs= iv=="V5" ? "0\t63372" : "63373\t999999"
				printf "* LOWGAINMULTIPLIER_SUM\t%d\t%s\t%d\t%d\t4\t18\t", tel, ivs, meth, isw;
				for (swl=4; swl<19; swl++) printf "%.2f\t", lmult[D2, isw, swl];
				printf "\n";
			}
			print "";
		}
		print "";
	}
}    


# Tel 1V5 swhigh 6 swlow 4 lmult 9.87923

# * LOWGAINMULTIPLIER_SUM 1 0	63372	18 4 18  13.69   10.06   8.35    7.38    6.76    6.35    6.07    5.85    5.68    5.55    5.44    5.35    5.29    5.26    5.25
