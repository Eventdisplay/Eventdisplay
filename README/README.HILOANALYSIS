hiloAnalysis 

	-a	Inner HiLo Run Number
	-b	Outer HiLo Run Number
	-t	Telescope Number 1-4
	-s	Sum Window
	-d	Directory where the DST files are

Example:
./hiloAnalysis -a 64035 -b 64036 -t 1 -s 18 -d $VERITAS_USER_DATA_DIR/analysis/EVD400_DST/

---------------------------------------------------------------------------------------

1. To make the DST files use:

VTS.EVNDISP.sub_make_DST.sh <runlist> [pedestal calculation (default=1=on)] [Sumation Window]

2. Then use the above hiloAnalysis to create a root file with the HiLo Multipliers.

3. To Look inside the root files produced open the sharedlibrary and use:

	VHiLoTools::getMeanRMS( string ifile, int tel, int sumwindow, bool bPlot )

		ifile		output from hiloAnalysis
		tel		telescope number 1 - 4
		sumwindow	Summation Window
		bPlot		set to TRUE to make plot


4. To merge the Files for all the summation windows use:

	HiLoTools::mergeHiLoFiles( string infile1, string infile2, int SW, int TEL )

		Puts the trees from infile2 into infile1 for Summation window SW and Telescope TEL

5. Then put these files into: $VERITAS_EVNDISP_AUX_DIR/Calibration/Tel_X/

For example: 
	if you are looking at runs 12345 and 12346 give them the filename 1234546.lmult.root
	and add the file to $VERITAS_EVNDISP_AUX_DIR/Calibration/calibrationlist.LowGain.dat


