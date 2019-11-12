/*! class VLikelihoodFitter
    fit and plot spectral data using likelihood methods


*/
#include "VLikelihoodFitter.h"






	bool VLikelihoodFitter::initialize()
	{


		// cout << "New Starting\n";
		fRandom = new TRandom3();
		// cout << "TRandom\n";
		fRandom->SetSeed(0);
		// cout << "TRandom\n";



		// Getting the Effective Areas
		// cout << "Getting Effective Areas\n";
		fMeanEffectiveAreaRec = getEffectiveAreasRec();
		fMeanEffectiveAreaMC = getEffectiveAreasMC();


		int i_nRec = int(3/fBinWidth);
		double *i_fRecBins = new double[i_nRec + 1];
		for (int i =0; i <= i_nRec; i ++)
		{
			i_fRecBins[i] = -1.1 + i*fBinWidth;
		}

		int i_nMC = int(3/fBinWidth);
		double *i_fMCBins = new double[i_nMC+1];
		for (int i = 0; i <= i_nMC; i ++)
		{
			i_fMCBins[i] =  -1.1 + i*fBinWidth;
		}



		setBinningRecMC(i_nRec, i_fRecBins,  i_nMC, i_fMCBins);
		// setBinningRecMC(i_nRec, i_fRecBins,  i_nMC, i_fMCBins);

		return true;

	}


	void VLikelihoodFitter::setBinWidth( double i_BinWidth )
	{
		fBinWidth = i_BinWidth;



	}


	void VLikelihoodFitter::setBinningRecMC(int i_nRec, double *i_fRecBins, int i_nMC, double *i_fMCBins)
	{



		cout << "Rebinning Function\n";
		nRec = i_nRec;
		fRecBins = i_fRecBins;
		nMC = i_nMC;
		fMCBins = i_fMCBins;
		fRecBinCentres.clear();
		fMCBinCentres.clear();
		fEnergyBias.clear();

		for (int i = 0; i < nRec; i++)
		{
			fRecBinCentres.push_back( fRecBins[i] + 0.5*(fRecBins[i+1] - fRecBins[i] ) );
		}

		for (int i = 0; i < nMC; i++)
		{
			fMCBinCentres.push_back( fMCBins[i] + 0.5*(fMCBins[i+1] - fMCBins[i] ) );
		}
		// Reset data vectors
		fOnRebinnedHistograms.clear();
		fOffRebinnedHistograms.clear();
		fResponseMatrixRebinned.clear();
		fOnCounts.clear();
		fOffCounts.clear();
		fTotalOn.clear();
		fTotalOff.clear();
		fLastOn.clear();
		fLastOff.clear();

		// Get Raw Histograms
		cout << "Getting Raw Data\n";
		cout << "Getting On\n";
		vector <TH1D*> i_hOnRaw = getCountingHistogramOnRaw();
		cout << "Getting Off\n";
		vector <TH1D*> i_hOffRaw = getCountingHistogramOffRaw();
		cout << "Getting Response Matrix\n";
		vector <TH2D*> i_hResponseMatrixRaw = getResponseMatrixRaw();
		vector <double> i_RunBias;

		// cout << "Rebinning Raw Data\n";
		cout << "Rebinning Raw\n";
		for( unsigned int i = 0; i < i_hOnRaw.size(); i++ )
		{


			// Simulating Observations
			// TH1D *simOn = new TH1D("simOn", "simOn", nRec , fRecBins );
			// TH1D *simOff = new TH1D("simOff", "simOff", nRec , fRecBins );
			sprintf( hname, "Rebinned On Counts %d", fRunList[i].runnumber);
			// TH1D *tmpOn = (TH1D*) i_hOnRaw[i]->Rebin(nRec  , hname ,fRecBins);
			// hname = "Rebinned Off Counts " + std::to_string(fRunList[i].runnumber);
			sprintf( hname, "Rebinned Off Counts %d", fRunList[i].runnumber);

			// TH1D *tmpOff = (TH1D*) i_hOffRaw[i]->Rebin(nRec  , hname ,fRecBins);
			i_RunBias.clear();
			// TF1 *fRand = new TF1("fRand","1 - TMath::Exp(1*x)/10000000000000", 1,10);
			// Rebinning Counts
			// cout << i << " Getting On" <<endl;
			sprintf( hname, "Rebinned On Counts %d", fRunList[i].runnumber);
			// sprintf( hname, "Rebinned Off Counts %d", fRunList[i].runnumber);

			fOnRebinnedHistograms.push_back((TH1D*) i_hOnRaw[i]->Rebin(nRec  , hname,fRecBins));
			// for (int k = 1; k <= tmpOn->GetSize(); k++)
			// {
			// 	simOn->SetBinContent(k, fRandom->Poisson(tmpOn->GetBinContent(k) * fRand->Eval(TMath::Power(10.0, simOn->GetBinCenter(k))) ));
			// 	simOff->SetBinContent(k, fRandom->Poisson(tmpOff->GetBinContent(k)) );
			// 	// cout << simOn->GetBinContent(k) << " " <<simOff->GetBinContent(k) << endl;
			// }
			// fOnRebinnedHistograms.push_back((TH1D*) simOn->Clone());
			// cout << i << " Getting Off" <<endl;

			sprintf( hname, "Rebinned Off Counts %d", fRunList[i].runnumber);
			// hname = "Rebinned Off Counts " + std::to_string(fRunList[i].runnumber);
			fOffRebinnedHistograms.push_back((TH1D*) i_hOffRaw[i]->Rebin(nRec , hname ,fRecBins));
			// fOffRebinnedHistograms.push_back((TH1D*) simOff->Clone());

			// Calcualating totals
			fTotalOn.push_back(fOnRebinnedHistograms[i]->Integral(1, nRec));
			fTotalOff.push_back(fOffRebinnedHistograms[i]->Integral(1, nRec));
			fLastOn.push_back(fOnRebinnedHistograms[i]->GetBinCenter(fOnRebinnedHistograms[i]->FindLastBinAbove(0.9)));
			fLastOff.push_back(fOffRebinnedHistograms[i]->GetBinCenter(fOffRebinnedHistograms[i]->FindLastBinAbove(0.9)));
			// Rebinning Response Matrix
			// cout << i << " Getting Response Matrix" <<endl;

			// hname = "Rebinned Response Matrix " + std::to_string(fRunList[i].runnumber);
			sprintf( hname, "Rebinned Response Matrix %d", fRunList[i].runnumber);

			TH2D* i_htmp2D = new TH2D(hname, hname , nRec  , fRecBins, nMC  , fMCBins );

			TAxis *i_xAxis = i_hResponseMatrixRaw[i]->GetXaxis();
			TAxis *i_yAxis = i_hResponseMatrixRaw[i]->GetYaxis();

			// Looping over and filling histogram
			// cout << "Rebinning Response Matrix\n";

			for (int j = 1; j <= i_xAxis->GetNbins(); j++)
			{
				for (int k = 1; k <= i_yAxis->GetNbins() ; k++)
				{
					i_htmp2D->Fill(i_xAxis->GetBinCenter(j), i_yAxis->GetBinCenter(k),  i_hResponseMatrixRaw[i]->GetBinContent(j,k) );
					// cout << k << " "  << i_htmp2D->GetBinContent(j,k) << " " << i_ktmp << endl;
					// if (i_ktmp != 0)
					// {
					// 	cout << k << " " << i_ktmp << endl;

					// }
				}
			}

			// // Normailising
			// double i_jTot = 0;
			//
			// for (int j = 1; j <= i_htmp2D->GetYaxis()->GetNbins(); j++)
			// {
			// 	i_jTot = 0;
			// 	for (int k = 1; k <= i_htmp2D->GetXaxis()->GetNbins(); k++)
			// 	{
			// 		i_jTot += i_htmp2D->GetBinContent(k,j);
			// 	}
			//
			// 	for (int k = 1; k <= i_htmp2D->GetXaxis()->GetNbins(); k++)
			// 	{
			// 		if (i_jTot == 0){continue;}
			// 		i_htmp2D->SetBinContent(k, j, i_htmp2D->GetBinContent(k,j)/i_jTot);
			// 	}
			//
			// }

			VHistogramUtilities::normalizeTH2D_y( i_htmp2D );

			// Determining Bias
			for (int j = 1; j <= i_htmp2D->GetYaxis()->GetNbins(); j++)
			{
				TF1 *i_fFit = new TF1("i_fFit", "gaus", -2,2.5);
				TH1D *i_slice = (TH1D*)i_htmp2D->ProjectionX("i_slice", j,j);
				if(i_slice->GetEntries() == 0)
				{
					i_RunBias.push_back(1.0);
					continue;
				}

				i_slice->Fit("i_fFit","q0");
				// cout << fRunList[i].runnumber << " "
				// 	 << i_htmp2D->GetYaxis()->GetBinCenter(j) << " "
				// 	 << i_fFit->GetParameter(1) << " "
				// 	 <<  (i_htmp2D->GetYaxis()->GetBinCenter(j) - i_fFit->GetParameter(1))/i_htmp2D->GetYaxis()->GetBinCenter(j) << endl;
				// i_RunBias.push_back(
				double tmp = (TMath::Power(10.,i_fFit->GetParameter(1)) - TMath::Power(10.,i_htmp2D->GetYaxis()->GetBinCenter(j))) * (TMath::Power(10.,i_fFit->GetParameter(1)) - TMath::Power(10.,i_htmp2D->GetYaxis()->GetBinCenter(j))) ;
				tmp = TMath::Sqrt(tmp);
				tmp /= TMath::Power(10.,i_htmp2D->GetYaxis()->GetBinCenter(j));
				i_RunBias.push_back(abs(tmp));
				// i_RunBias.push_back(0.2);
			}
			// cout << "Filling Vectors\n";

			// Filling vector

			i_htmp2D->GetXaxis()->SetTitle("Energy_{rec} [TeV]");
			i_htmp2D->GetYaxis()->SetTitle("Energy_{mc} [TeV]");
			fResponseMatrixRebinned.push_back(i_htmp2D);
			fEnergyBias.push_back(i_RunBias);

			// delete simOn;
			// delete simOff;
		}

		// Getting counts from histograms
		fOffCounts = getCounts(fOffRebinnedHistograms);
		fOnCounts = getCounts(fOnRebinnedHistograms);


	}



	vector <TGraphAsymmErrors*> VLikelihoodFitter::getEffectiveAreasMC()
	{
		vector <TGraphAsymmErrors*> iVtemp ;
		cout << "Getting mean effective area for MC energy space" << endl;
		for( unsigned int i = 0; i < fRunList.size(); i++ )
		{
			cout << "Getting gMeanEffectiveAreaMC for run number: " << fRunList[i].runnumber << endl;
			sprintf(hname,"gMeanEffectiveAreaMC");
			cout << "Getting gMeanEffectiveAreaMC for run number: " << fRunList[i].runnumber << endl;

			TGraphAsymmErrors* i_gMeanEffectiveAreaMC = ( TGraphAsymmErrors* )getHistogram( hname, fRunList[i].runnumber, "EffectiveAreas", -9999 )->Clone();
			cout << "Successfully cloned gMeanEffectiveAreaMC for run number: " << fRunList[i].runnumber << endl;
			sprintf( hname, "gMeanEffectiveAreaMC_%d", fRunList[i].runnumber);

			// hname = hname + "_" + std::to_string(fRunList[i].runnumber);
			i_gMeanEffectiveAreaMC->SetTitle(hname);
			cout << "Successfully renamed gMeanEffectiveAreaMC for run number: " << fRunList[i].runnumber << " , " << hname <<  endl;

			iVtemp.push_back(i_gMeanEffectiveAreaMC);
			cout << "Successfully added gMeanEffectiveAreaMC to vector. " << endl;

		}
		cout << "Reading gMeanEffectiveAreaMCs completed " << endl;

		return iVtemp;

	}

	vector <TGraphAsymmErrors*> VLikelihoodFitter::getEffectiveAreasRec()
	{
		vector <TGraphAsymmErrors*> iVtemp ;


		for( unsigned int i = 0; i < fRunList.size(); i++ )
		{
			sprintf(hname,"gMeanEffectiveArea");
			TGraphAsymmErrors* i_gMeanEffectiveArea = ( TGraphAsymmErrors* )getHistogram( hname, fRunList[i].runnumber, "EffectiveAreas", -9999 )->Clone();
			sprintf( hname, "gMeanEffectiveArea_%d", fRunList[i].runnumber);

			// hname = hname + "_" + std::to_string(fRunList[i].runnumber);
			i_gMeanEffectiveArea->SetTitle(hname);
			iVtemp.push_back(i_gMeanEffectiveArea);
		}

		return iVtemp;

	}

	vector <TH1D*> VLikelihoodFitter::getCountingHistogramOnRaw()
	{
		vector <TH1D*> iVTemp;


		for( unsigned int i = 0; i < fRunList.size(); i++ )
		{

			sprintf(hname,"herecCounts_on");
			TH1D* i_hErecCountsOn = ( TH1D* )getHistogram( hname, fRunList[i].runnumber, "energyHistograms", -9999 );
			sprintf( hname, "herecCounts_on_%d", fRunList[i].runnumber);

			// hname = hname + "_" + std::to_string(fRunList[i].runnumber);
			i_hErecCountsOn->SetTitle(hname);
			iVTemp.push_back(i_hErecCountsOn);
		}
		return iVTemp;
	}

	vector <TH1D*> VLikelihoodFitter::getCountingHistogramOffRaw()
	{
		vector <TH1D*> iVTemp;


		for( unsigned int i = 0; i < fRunList.size(); i++ )
		{

			sprintf(hname,"herecCounts_off");
			TH1D* i_hErecCountsOff = ( TH1D* )getHistogram( hname, fRunList[i].runnumber, "energyHistograms", -9999 );
			// hname = hname + "_" + std::to_string(fRunList[i].runnumber);
			sprintf( hname, "herecCounts_off_%d", fRunList[i].runnumber);

			i_hErecCountsOff->SetTitle(hname);
			iVTemp.push_back(i_hErecCountsOff);
		}
		return iVTemp;
	}



	vector <TH2D*> VLikelihoodFitter::getResponseMatrixRaw()
	{
		vector <TH2D*>iVtemp ;


		for( unsigned int i = 0; i < fRunList.size(); i++ )
		{
			cout << "Getting Effective Areas for Run " << fRunList[i].runnumber << endl;
			sprintf(hname, "hMeanResponseMatrix");
			TH2D* i_hResponseMatrix = ( TH2D* )getHistogram( hname, fRunList[i].runnumber, "EffectiveAreas", -9999 )->Clone();
			sprintf( hname, "hMeanResponseMatrix_%d", fRunList[i].runnumber);

			// hname = hname + "_" + std::to_string(fRunList[i].runnumber);
			i_hResponseMatrix->SetTitle(hname);
			iVtemp.push_back(i_hResponseMatrix);


		}

		return iVtemp;


	}


	void VLikelihoodFitter::printRunInfo()
	{
		double i_total_tOn = 0;
		double i_total_On = 0;
		double i_total_Off = 0;

		cout << "Entry #\tRun\tLivetime\tTotal ON\tTotal Off\tAlpha\n";
		for( unsigned int i = 0; i < fRunList.size(); i++ )
		{
			cout << i << "\t" << fRunList[i].runnumber << "\t" << fRunList[i].tOn * fRunList[i].deadTimeFraction <<
			 "\t\t" << fTotalOn[i] << "\t\t" << fTotalOff[i] << "\t\t" << fRunList[i].alpha << "\t\t" << fRunList[i].MJD <<
			 "\t\t" << TMath::Power(10,fLastOn[i]) << "\t\t" << TMath::Power(10,fLastOff[i]) <<  endl;
			i_total_tOn += fRunList[i].tOn * fRunList[i].deadTimeFraction;
			i_total_On += fTotalOn[i];
			i_total_Off += fTotalOff[i];
		}
		cout << "Total\t\t" << i_total_tOn << "\t\t" << i_total_On <<"\t\t" << i_total_Off << endl;
	}



	vector < vector <double> > VLikelihoodFitter::getCounts(vector <TH1D*> i_hTemp)
	{
		vector < vector <double> > i_vTemp;
		vLastOnCount.clear();
		int tmpLastBin = 0;
		for ( unsigned int i = 0; i < i_hTemp.size(); i++)
		{
			vector <double> i_vRunCounts;
			for (int j = 0; j < i_hTemp[i]->GetXaxis()->GetNbins(); j++ )
			{
				if (i_hTemp[i]->GetBinContent(j+1) > 1e9)
				{
					cout << "Error \t " << i << " " << j << endl;
					i_vRunCounts.push_back( 0 + 1.e-10);
					continue;
				}

				if (i_hTemp[i]->GetBinContent(j+1) >= 1)
				{
					tmpLastBin = j;
				}
				i_vRunCounts.push_back( i_hTemp[i]->GetBinContent(j+1) + 1.e-9 );
			}
			vLastOnCount.push_back(tmpLastBin);
			i_vTemp.push_back( i_vRunCounts );
		}

		return i_vTemp;
	}





	void VLikelihoodFitter::setModel(int i_ID = 0)
	{

		// Note for spectral weighting of bins, the spectral index must be the 2nd parameter
		// i.e. [1] = Spectral index
		if (i_ID == 0)
		{
			sprintf( hname, "[0]*TMath::Power(TMath::Power(10.0,x) / %f, [1])", E_Norm);
			fModel = new TF1("fModel", hname,-1.0,5.0);
			fModelID = i_ID;
			nParms = 2;
		}

		// Power Law with Exponential Cut off
		else  if (i_ID == 1)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1] ) * TMath::Exp( -1. * TMath::Power( 10, x ) / [2] )", E_Norm);
			fModel = new TF1("fModel", hname,-1.0,5.0);
			fModelID = i_ID;
			nParms = 3;
		}


		// Curved Spectrum
		else if (i_ID == 2)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x )/ %f, [1]+[2]*TMath::Power( 10, x ) )", E_Norm);
			fModel = new TF1("fModel", hname,-1.0,5.0);
			fModelID = i_ID;
			nParms = 3;

		}




		// Log-parabola
		else if (i_ID == 3)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1]+[2]*TMath::Log( TMath::Power( 10, x ) ) )", E_Norm);
			fModel = new TF1("fModel", hname,-1.0,5.0);
			fModelID = i_ID;
			nParms = 3;

		}

		// EBL Corrected
		else if (i_ID == -999)
		{
			if (!fModel_int)
			{
				cout << "VLikelihoodFitter::setModel Error \n"
					 << "You must call VLikelihoodFitter::setIntrinsicModel(i_ID) first!\n"
					 << "Only use this if you know what you are doing!\n";
				fModel = 0;
				return;
			}
			fModel = new TF1("fModel",this,&VLikelihoodFitter::calculateIntrinsicSpectrum,-1,5, nParms);
			// nParms already set by
			// nParms = 3;

		}

		else
		{
			cout << "Model " << i_ID << " not found!\n";
		}
	}


	vector < vector <double> >  VLikelihoodFitter::getModelPredictedExcess(const double *parms, double i_eThresh = -999)
	{
		// cout << "Getting Predicted Excess\n";
		fModel->SetParameters(parms);

		vector < vector <double> > i_vModel;
		vector <double> i_TrueCenters;

		// Getting Spectrally weighted bin centres
		for (int i = 0; i < nMC; i++)
		{
			i_TrueCenters.push_back( VMathsandFunctions::getSpectralWeightedMeanEnergy(fMCBins[i], fMCBins[i+1], parms[1]) );
		}


		double i_ModelElement = 0;
		double i_EffetiveAreaElement = 0;
		double i_MatrixElement = 0;
		double i_ConversionElement = 1.e4;
		double i_EnergyElement = 0;
		double i_TOn = 0;

		// Current Bin
		int i_J_Rec_Cur = 0;
		int i_L_MC_Cur = 0;
		// Looping over each run
		for (unsigned int i = 0; i < fRunList.size(); i++)
		{

			TAxis *xpoints = fResponseMatrixRebinned[i]->GetXaxis();
			TAxis *ypoints = fResponseMatrixRebinned[i]->GetYaxis();
			vector <double> i_vTmp;

			// Looping over each Rec Energy bin
			for (int j = 0; j < nRec; j++)
			{
				// if (fRecBinCentres[j] < i_eThresh)
				// {
				// 	continue;
				// }
				i_vTmp.push_back(0);
				i_J_Rec_Cur = xpoints->FindBin(fRecBinCentres[j]);
				// if (fRecBinCentres[j] < fLowEnergyRange){continue;}
				// if (fRecBinCentres[j] > fHigherEnergyRange){continue;}
				// Looping over each MC Energy bin
				for (int l = 0; l < nMC; l++)
				{
					i_L_MC_Cur = ypoints->FindBin(fMCBinCentres[l]);

					// Lower range on sanity
					if (fEnergyBias[i][l] > fThresholdBias){continue;}
					// if (fRecBinCentres[j] < -0.9 ){continue;}

					i_ModelElement = fModel->Eval(fMCBinCentres[l]);
					// i_ModelElement = fModel->Integral(fMCBins[l], fMCBins[l+1]);

					i_EffetiveAreaElement = fMeanEffectiveAreaMC[i]->Eval(fMCBinCentres[l]);
					// i_EffetiveAreaElement = fMeanEffectiveAreaMC[i]->Integral(fMCBins[l], fMCBins[l+1]);

					i_MatrixElement = fResponseMatrixRebinned[i]->GetBinContent(i_J_Rec_Cur, i_L_MC_Cur);
					i_TOn = fRunList[i].tOn * fRunList[i].deadTimeFraction  ;

					// i_EnergyElement = 1.0;
					i_EnergyElement = pow(10.0,fMCBins[l+1]) - pow(10.0,fMCBins[l]);
					// i_EnergyElement = fMCBins[l+1] - fMCBins[l];

					// cout << i_EnergyElement << " " << i_ModelElement << " " << i_EffetiveAreaElement << " " << i_MatrixElement << " " << i_TOn << " " << endl;

					// if (i_MatrixElement > 1.0 )
					// {
					// 	cout << "Horrible Problem " << i_MatrixElement << endl;
					// }
					i_vTmp[j] +=   i_ModelElement * i_EffetiveAreaElement * i_MatrixElement * i_ConversionElement *  i_TOn  * i_EnergyElement ;


				}

			}
			i_vModel.push_back(i_vTmp);

		}
		return i_vModel;
	}





	vector < vector <double> >  VLikelihoodFitter::getModelPredictedOff(const double *parms)
	{

		// cout << "Getting Predicted Off\n";

		vector < vector <double> > i_OffMLE;
		vector < vector <double> > i_myModel = getModelPredictedExcess(parms);
		// cout <<" Getting Size "<< vRunNum.size() << endl;
		for (unsigned int i = 0; i < fRunList.size(); i++)
		{
			vector <double> i_vTmp;

			for (unsigned int j = 0; j < i_myModel[i].size(); j++)
			{

				double i_a = 0;
				double i_b = 0;
				i_a = fRunList[i].alpha *(fOnCounts[i][j] + fOffCounts[i][j]) - (fRunList[i].alpha  + 1)*i_myModel[i][j];
				i_b = 1.0/(2.0*fRunList[i].alpha *(fRunList[i].alpha  + 1.0)) * (i_a + sqrt(pow(i_a,2.0) +4.0 * fRunList[i].alpha * (fRunList[i].alpha + 1.0) * fOffCounts[i][j]*i_myModel[i][j] ));
				// To advoid future log(0) problems
				i_b += 1.e-12;
				i_vTmp.push_back(i_b);
			}

			i_OffMLE.push_back(i_vTmp);
		}

		return i_OffMLE;
	}

// 	double VLikelihoodFitter::getLogLi(const double *parms)
// 	{
// 		// cout << "Getting getLogLi\n";

// 		double LogLi = 0;
// 		// double Alpha = 1.0/8.0;
// 		vector < vector <double> > i_myModel = getModelPredictedExcess( parms );
// 		vector < vector <double> > i_myvOffMLE = getModelPredictedOff( parms );

// 		double a, b, c, d;
// 		double tmpNBinsFit = 0;
// 		bool i_bInc = true;
// 		nBinsFit = 0;

// 		for (int i = 0; i < fRunList.size(); i++){

// 			// nBinsFit = 0;
// 			// cout << vRunNum[i] << endl;

// 			if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
// 			{
// 				continue;
// 			}

// 			i_bInc = true;

// 			for (int e = 0; e < fExcludeMJD.size(); e++)
// 			{
// 				if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
// 				{
// 					i_bInc = false;
// 					break;
// 				}
// 			}

// 			for (int e = 0; e < fExcludeRun.size(); e++)
// 			{
// 				if (fRunList[i].runnumber == fExcludeRun[e])
// 				{
// 					i_bInc = false;
// 					break;
// 				}
// 			}

// 			for (int j = 0; j < nMC; j ++)
// 			{

// 				if (i_bInc == false){continue;}
// 				if (fRecBinCentres[j] < fFitMin){continue;}
// 				if (fRecBinCentres[j] > fFitMax) {continue;}
// 				// if (j >= vLastOnCount[i] ){continue;}

// 				nBinsFit += 1;
// 				if (fOnCounts[i][j] > 0)
// 				{
// 					a = fOnCounts[i][j]*log(i_myModel[i][j] + fRunList[i].alpha*i_myvOffMLE[i][j]);
// 				}

// 				else
// 				{
// 					a = 0;
// 				}


// 				if ( fOffCounts[i][j] > 0 )
// 				{
// 					b = fOffCounts[i][j]*log(i_myvOffMLE[i][j]);
// 				}

// 				else
// 				{
// 					b = 0;

// 				}
// 				// Getting c term
// 				c = -1.0*(fRunList[i].alpha + 1.0)*i_myvOffMLE[i][j];

// 				// Getting d term
// 				d = -i_myModel[i][j];

// 				LogLi = LogLi + a + b + c + d;
// 			}
// 			// if (nBinsFit > tmpNBinsFit)
// 			// {
// 			//   tmpNBinsFit = nBinsFit;
// 			// }
// 		}
// 		// cout << "Exiting GetLogLi\n";
// 		// nBinsFit ;
// 		return -1 * LogLi;
// 	}

// double VLikelihoodFitter::getLogLiTotal(const double *parms)
// 	{
// 		// cout << "Getting getLogLiTotal\n";

// 		double LogLi = 0;
// 		// double Alpha = 1.0/8.0;
// 		// cout << "Getting Model Counts" << endl;
// 		vector < vector <double> > i_myModel = getModelPredictedExcess( parms );
// 		vector < vector <double> > i_myvOffMLE = getModelPredictedOff( parms );

// 		double a, b, c, d;
// 		double tmpNBinsFit = 0;
// 		bool i_bInc = true;
// 		nBinsFit = 0;

// 		vector <double> i_total_On( nMC );
// 		// i_total_On.assign(nMC+1, 0);
// 		vector <double> i_total_Off( nMC );
// 		vector <double> i_total_Model( nMC );
// 		vector <double> i_total_ModelOff( nMC );
// 		// double i_total_On[nMC];

// 		double i_mean_alpha = 0;
// 		int count = 0;
// 		// cout << "Looping " << endl;
// 		for (int i = 0; i < fRunList.size(); i++){
// 			// cout << i << endl;
// 			nBinsFit = 0;

// 			if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
// 			{
// 				continue;
// 			}

// 			i_bInc = true;

// 			for (int e = 0; e < fExcludeMJD.size(); e++)
// 			{
// 				if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
// 				{
// 					i_bInc = false;
// 					break;
// 				}
// 			}

// 			for (int e = 0; e < fExcludeRun.size(); e++)
// 			{
// 				if (fRunList[i].runnumber == fExcludeRun[e])
// 				{
// 					i_bInc = false;
// 					break;
// 				}
// 			}

// 			for (int j = 0; j < nMC; j ++)
// 			{

// 				if (i_bInc == false){continue;}
// 				if (fRecBinCentres[j] < fFitMin){continue;}
// 				if (fRecBinCentres[j] > fFitMax) {continue;}
// 				// if (j >= vLastOnCount[i] ){continue;}
// 				// i_total_On[j] += 1;
// 				// cout << i << " , " << j << endl;
// 				i_total_On[j] += fOnCounts[i][j];
// 				i_total_Off[j] += fOffCounts[i][j];
// 				i_total_Model[j] += i_myModel[i][j];
// 				i_total_ModelOff[j] += i_myvOffMLE[i][j];

// 			}
// 			count += 1;
// 			i_mean_alpha += fRunList[i].alpha;

// 		}
// 		i_mean_alpha /= count;

// 		// cout << "Looping over counts" << endl;

// 		for (int i = 0; i < i_total_On.size(); i++ )
// 		{
// 			if (fRecBinCentres[i] < fFitMin){continue;}
// 			if (fRecBinCentres[i] > fFitMax) {continue;}
// 			nBinsFit += 1;
// 			// cout << i << " : a" << endl;
// 			if (i_total_On[i] > 0)
// 			{
// 				a = i_total_On[i]*log(i_total_Model[i] + i_mean_alpha*i_total_ModelOff[i]);
// 			}

// 			else
// 			{
// 				a = 0;
// 			}
// 			// cout << "a: " << a << endl;
// 			// cout << "b" << endl;

// 			if ( i_total_Off[i] > 0 )
// 			{
// 				b = i_total_Off[i]*log(i_total_ModelOff[i]);
// 			}

// 			else
// 			{
// 				b = 0;

// 			}
// 			// Getting c term
// 			// cout << "c" << endl;

// 			c = -1.0*(i_mean_alpha + 1.0)*i_total_ModelOff[i];

// 			// Getting d term
// 			// cout << "d" << endl;

// 			d = -1.0*i_total_Model[i];
// 			// cout << "adding" << endl;
// 			LogLi +=  a + b + c + d;
// 			// cout << "Done" << endl;
// 		}

// 		// cout << "Finished: " << -1 * LogLi << endl;
// 		// cout << LogLi << endl;
// 		// LogLi = -1*LogLi;
// 		// cout << "Exiting GetLogLi\n";
// 		// nBinsFit ;
// 		// double r = LogLi;
// 		// cout << "Returning" <<endl;
// 		return -1*LogLi;
// 	}

	// double VLikelihoodFitter::getLogL0Total(const double *parms)
	// {
	// 	// cout << "Getting getLogL0Total\n";

	// 	double LogLi = 0;
	// 	// double Alpha = 1.0/8.0;
	// 	// vector < vector <double> > i_myModel = getModelPredictedExcess( parms );
	// 	// vector < vector <double> > i_myvOffMLE = getModelPredictedOff( parms );

	// 	double a, b, c, d;
	// 	double tmpNBinsFit = 0;
	// 	bool i_bInc = true;
	// 	nBinsFit = 0;

	// 	vector <double> i_total_On(fOnCounts[0].size());
	// 	vector <double> i_total_Off(fOnCounts[0].size());
	// 	// vector <double> i_total_Model(i_myModel[0].size());
	// 	// vector <double> i_total_ModelOff(i_myModel[0].size());
	// 	double i_mean_alpha = 0;
	// 	int count = 0;

	// 	for (int i = 0; i < fRunList.size(); i++){

	// 		// nBinsFit = 0;
	// 		// cout << vRunNum[i] << endl;

	// 		if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
	// 		{
	// 			continue;
	// 		}

	// 		i_bInc = true;

	// 		for (int e = 0; e < fExcludeMJD.size(); e++)
	// 		{
	// 			if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
	// 			{
	// 				i_bInc = false;
	// 				break;
	// 			}
	// 		}
	// 		for (int e = 0; e < fExcludeRun.size(); e++)
	// 		{
	// 			if (fRunList[i].runnumber == fExcludeRun[e])
	// 			{
	// 				i_bInc = false;
	// 				break;
	// 			}
	// 		}

	// 		for (int j = 0; j < nMC; j ++)
	// 		{

	// 			if (i_bInc == false){continue;}
	// 			if (fRecBinCentres[j] < fFitMin){continue;}
	// 			if (fRecBinCentres[j] > fFitMax) {continue;}
	// 			// if (j >= vLastOnCount[i] ){continue;}
	// 			i_total_On[j] += fOnCounts[i][j];
	// 			i_total_Off[j] += fOffCounts[i][j];
	// 			// i_total_Model[i] += i_myModel[i][j];
	// 			// i_total_ModelOff[i] += i_myvOffMLE[i][j];

	// 		}
	// 		count += 1;
	// 		i_mean_alpha += fRunList[i].alpha;

	// 	}
	// 	i_mean_alpha /= count;

	// 	for (int i = 0; i < i_total_On.size(); i++ ){
	// 		nBinsFit += 1;

	// 		if (fRecBinCentres[i] < fFitMin){continue;}
	// 		if (fRecBinCentres[i] > fFitMax) {continue;}

	// 		if (i_total_On[i] > 0)
	// 		{
	// 			a = i_total_On[i]*log(i_total_On[i]);
	// 		}

	// 		else
	// 		{
	// 			a = 0;
	// 		}


	// 		if ( i_total_Off[i] > 0 )
	// 		{
	// 			b = i_total_Off[i]*log(i_total_Off[i]);
	// 		}

	// 		else
	// 		{
	// 			b = 0;

	// 		}
	// 		// Getting c term
	// 		c = -1.0*i_total_On[i];

	// 		// Getting d term
	// 		d = -1.0*i_total_Off[i];

	// 		LogLi = LogLi + a + b + c + d;
	// 	}

	// 	// cout << "Exiting GetLogLi\n";
	// 	// nBinsFit ;
	// 	return -1 * LogLi;
	// }



	// Function to optimise the likelihood function
	// Returns a TF1* with the best fit parameters
	TF1* VLikelihoodFitter::getLikelihoodFit()
	{


		//nParms = 2;
		ROOT::Math::Minimizer* i_min =
		   ROOT::Math::Factory::CreateMinimizer("Minuit", "Minos");

		// // set tolerance , etc...
		i_min->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
		i_min->SetMaxIterations(100000);  // for GSL
		i_min->SetTolerance(0.001);
		i_min->SetPrintLevel(1);
		i_min->SetErrorDef(0.5);


		ROOT::Math::Functor i_fitfunction(this,&VLikelihoodFitter::getLogLi,nParms);
		double *step = new double[nParms];
		double *variable = new double[nParms];

		i_min->SetFunction(i_fitfunction);


		if (fModelID == 0)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			variable[0] = 1.e-12;
			variable[1] = -2.5;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);

		}


		if (fModelID == 1)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.01;
			variable[0] = 1.e-13;
			variable[1] = -1.5;
			variable[2] = 4.0;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);
			i_min->SetVariable(2,"E_CutOff", variable[2], step[2]);

		}


		if (fModelID == 2)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.01;
			variable[0] = 1.e-13;
			variable[1] = -1.5;
			variable[2] = -0.01;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);
			i_min->SetVariable(2,"Beta",variable[2], step[2]);
		}


		if (fModelID == 3)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.001;
			variable[0] = 1.e-12;
			variable[1] = -1.5;
			variable[2] = 0.1;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Alpha",variable[1], step[1]);
			i_min->SetVariable(2,"Beta",variable[2], step[2]);
		}

		if (fModelID == 4)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.001;
			variable[0] = 1.e-12;
			variable[1] = -2.5;
			variable[2] = 1.0;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);
			i_min->SetLimitedVariable(2,"Alpha",variable[2], step[2], 0., 2.0);
		}

		// do the minimization
		i_min->Minimize();
		const double *xs = i_min->X();
		i_min->Hesse();
		double *i_FitValues = new double[nParms];
		const double *i_Errors = i_min->Errors();

		for (int i = 0; i < nParms; i++)
		{
			i_FitValues[i] = xs[i];
		}


		fModel->SetParameters(xs);
		fModel->SetParErrors(i_Errors);

		TF1 *i_BestFit = 0;

		// Removing exp(tau) component
		// if (fModel_intID != -999)
		// {
		// 	cout << "True doing EBL analysis" << endl;

		// 	for (int i = 0; i < nParms; i++)
		// 	{
		// 		fModel_int->SetParameter(i, fModel->GetParameter(i));
		// 		fModel_int->SetParError(i, fModel->GetParError(i));
		// 	}
		// 	fModel_int->SetParameter(0, fModel->GetParameter(0)*TMath::Exp(fEBLOpacityGraph->Eval( E_Norm )) );
		// 	fModel_int->SetParError(0, fModel->GetParError(0)*TMath::Exp(fEBLOpacityGraph->Eval( E_Norm )) );

		// 	cout << scientific << fModel_int->GetParameter(0) << " " << fModel_int->GetParError(0) << " " << E_Norm << " " << TMath::Exp(fEBLOpacityGraph->Eval( E_Norm )) << endl;
		// 	cout << scientific << fModel->GetParameter(0) << " " << fModel->GetParError(0) << endl;
		// 	i_BestFit = (TF1*)fModel_int->Clone();

		// }

		// else
		// {
		// 	// cout << "False Doing normal Analysis" << endl;
		// 	i_BestFit = (TF1*)fModel->Clone();
		// }


		// Removing exp(tau) component
		if (fModel_intID != -999)
		{
			cout << "True doing EBL analysis" << endl;

			for (int i = 0; i < nParms; i++)
			{
				fModel_int->SetParameter(i, fModel->GetParameter(i));
				fModel_int->SetParError(i, fModel->GetParError(i));
			}
			i_BestFit = (TF1*)fModel_int->Clone();

		}
		else
		{
			i_BestFit = (TF1*)fModel->Clone();
		}

		double ErrUp = 0;
		double ErrLow = 0;

		for (int i = 0; i < nParms; i++)
		{
			bool bErrors = i_min->GetMinosError(i, ErrLow, ErrUp);
			cout << scientific << "Variabile " << i << " Error Status: " << bErrors << ", E_Low = " << ErrLow << ", E_Up = " << ErrUp << endl;
		}


		unsigned int nSteps = 100;
		Index_Scan = new double[nSteps];
		Index_Scan_Likelihood = new double[nSteps];
		Norm_Scan = new double[nSteps];
		Norm_Scan_Likelihood = new double[nSteps];


		i_min->Scan(0, nSteps, Norm_Scan, Norm_Scan_Likelihood, 1.e-13,1.e-11);
		i_min->Scan(1, nSteps, Index_Scan, Index_Scan_Likelihood, -3.5,-1.5);

		// Getting Confidince interval for power law





		// if (fModel_intID == 0)
		// {
		// 	cout << "True doing EBL analysis" << endl;
		// 	for (int i = 0; i < nParms; i++)
		// 	{
		// 		fModel_int->SetParameter(i, fModel->GetParameter(i));
		// 		fModel_int->SetParError(i, fModel->GetParError(i));
		// 	}
		// }

		// else
		// {
		// 	cout << "False Doing normal Analysis" << endl;
		// }

		// double i_energy = fFitMin;
		// vector<double> i_verr;
		// vector<double> i_venergy;
		// vector<double> i_vFluxMin;
		// vector<double> i_vFluxMax;
		// for (int j = 0 ; j < i_count; j++)
		// {

		// 	double x = TMath::Power(10,i_energy) ;
		// 	// double i_err = i_min->CovMatrix(0,0) * TMath::Power(x,i_FitValues[1]) * TMath::Power(x,i_FitValues[1]);
		// 	// // cout << i_err << "\t" ;
		// 	// i_err += i_min->CovMatrix(1,1) *  (i_FitValues[0] * TMath::Power(x ,i_FitValues[1])*log(x )) *(i_FitValues[0] * TMath::Power(x ,i_FitValues[1])*log(x ))  ;
		// 	// // cout << i_err << "\t" ;

		// 	// i_err += 2. * i_min->CovMatrix(0,1)  *  TMath::Power(x,i_FitValues[1]) * (i_FitValues[0] * TMath::Power(x ,i_FitValues[1])*log(x)) ;
		// 	// // cout << i_err << "\t" ;

		// 	// i_err = sqrt(i_err);
		// 	// cout << i_err << "\t" ;



		// 	double i_err = i_Errors[0] * i_Errors[0] / i_FitValues[0] / i_FitValues[0] ;
		// 	i_err += log(x/E_Norm) * log(x/E_Norm) * i_Errors[1] * i_Errors[1];
		// 	i_err += -2*log(x/E_Norm) / i_FitValues[0] * i_min->CovMatrix(0,1) ;
		// 	i_err = i_BestFit->Eval(i_energy) * sqrt(i_err);
		// 	// cout << i_err << "\t";
		// 	i_verr.push_back(i_err);
		// 	i_venergy.push_back(i_energy);
		// 	// double i_err = 0;
		// 	// double* i_energy = new double[2];
		// 	// i_energy &(fFitMin + j * delE);

		// 	// Getting error calculation using TF1::GradientPar

		// 	// double* p = new double;
		// 	// *p = fFitMin + j * delE;
		// 	// for (int k  = 0; k < nParms; k++)
		// 	// {
		// 	// 	for (int l = 0; l < nParms; l++)
		// 	// 	{
		// 	// 		i_err += fModel->GradientPar(k,p) * fModel->GradientPar(l,p) * i_min->CovMatrix(k,l) ;
		// 	// 	}

		// 	// }

		// 	// cout << "Setting Point: " << j << " " << i_energy << " " << i_BestFit->Eval(i_energy) << " " << i_err << endl;
		// 	fConfidinceInterval->SetPoint(j, i_energy, i_BestFit->Eval(i_energy) );
		// 	fConfidinceInterval->SetPointError(j , 0 , i_err );
		// 	// cout << j << " " << i_vFluxMax[j] << " " << fModel->Eval(i_venergy[j]) << " " << i_vFluxMin[j] << endl;

		// 	// if (fModel_intID == 0)
		// 	// {
		// 	// 	i_vFluxMin.push_back(fModel_int->Eval(i_energy) - i_err);
		// 	// 	i_vFluxMax.push_back(fModel_int->Eval(i_energy) + i_err);
		// 	// }

		// 	// else
		// 	// {
		// 	// 	i_vFluxMin.push_back(fModel->Eval(i_energy) - i_err);
		// 	// 	i_vFluxMax.push_back(fModel->Eval(i_energy) + i_err);
		// 	// }

		// 	// cout << fModel->Eval(i_energy) << "\t" << fModel->Eval(i_energy) + i_err << "\t" << fModel->Eval(i_energy) - i_err << endl;

		// 		i_energy += delE;
		// }


		cout << "Binned Likelihood Fit: \n" ;
		cout << "Parameter \t Best Fit \t Error\n";
		for ( int i = 0; i < nParms; i++)
		{
			cout << scientific << "[" << i << "]\t\t" << " " << i_min->VariableName(i) << " \t\t" << i_BestFit->GetParameter(i) << "\t\t" << i_BestFit->GetParError(i) << endl;
		}
		cout << endl;

		cout << "\n\nPrinting Covariance Matrix\n";
		double *i_covmat = new double [nParms*nParms];
		for (int i =0; i < nParms; i++)
		{
			for (int j = 0; j < nParms; j++)
			{
				i_covmat[ i*nParms + j ] = i_min->CovMatrix(i,j);
				cout << scientific << i_min->CovMatrix(i,j) << "\t" ;
			}
			cout << endl;
		}

		cout << "\n";
		// GetChi2(xs);
		getChi2 (xs);
		cout << "\n";
		cout << "Calculating Total Chi^2\n";
		getChi2Total (xs);
		cout << "\n";

		// fModel->Draw();
		// fModel->GetXaxis()->
		// gConfidinceInterval = GetConfidinceInterval();

		// int i_count = 10000;
		// double delE = (fFitMax - fFitMin)/i_count ;
		if (fConfidinceInterval)
		{
			delete fConfidinceInterval;
		}
		fConfidinceInterval = calculateConfidinceInterval(i_covmat, i_BestFit, fModelID, nParms);

		// fConfidinceInterval = new TGraphErrors(i_count);
		float *i_flux = getIntergralFlux(fFitMin, fFitMax , i_BestFit, true);

		cout << "Intergral Flux:\n";
		cout << "F (" << TMath::Power(10,fFitMin) << " TeV < E < " << TMath::Power(10,fFitMax) << ") = " << i_flux[0] << "+/-" << i_flux[1] << " [Photons/cm^2/s] \n";
		cout << "F (" << TMath::Power(10,fFitMin) << " TeV < E < " << TMath::Power(10,fFitMax) << ") = " << i_flux[2] << "+/-" << i_flux[3] << " [Crab] \n";


		// Getting Decorrelation Energy
		double E_d = E_Norm * TMath::Exp( -i_min->CovMatrix(0,1) / xs[0] / i_Errors[1] / i_Errors[1] );

		cout << "Printing Decorrelation Energy (Assuming a Power Law Model, consider reapplying the fit.):\nE_d : " << E_d << endl;

		// if (fModel_intID == 0 )
		// {
		// 	cout << "True" << endl;
		// 	fModel_int->SetChisquare(fModel->GetChisquare() );
		// 	fModel_int->SetNDF(fModel->GetNDF() );
		// 	return fModel_int;
		// }
		i_BestFit->SetChisquare(fModel->GetChisquare() );
		i_BestFit->SetNDF(fModel->GetNDF() );
		return i_BestFit;
	}




TF1* VLikelihoodFitter::getLikelihoodFitTotal()
	{

		//nParms = 2;
		ROOT::Math::Minimizer* i_min =
		   ROOT::Math::Factory::CreateMinimizer("Minuit", "Minos");

		// // set tolerance , etc...
		i_min->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
		i_min->SetMaxIterations(100000);  // for GSL
		i_min->SetTolerance(0.001);
		i_min->SetPrintLevel(1);
		i_min->SetErrorDef(0.5);


		ROOT::Math::Functor i_fitfunction(this,&VLikelihoodFitter::getLogLiTotal,nParms);
		double *step = new double[nParms];
		double *variable = new double[nParms];

		i_min->SetFunction(i_fitfunction);


		if (fModelID == 0)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			variable[0] = 1.e-12;
			variable[1] = -2.5;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);

		}


		if (fModelID == 1)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.01;
			variable[0] = 1.e-13;
			variable[1] = -1.5;
			variable[2] = 4.0;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);
			i_min->SetVariable(2,"E_CutOff", variable[2], step[2]);

		}


		if (fModelID == 2)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.01;
			variable[0] = 1.e-13;
			variable[1] = -1.5;
			variable[2] = -0.01;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);
			i_min->SetVariable(2,"Beta",variable[2], step[2]);
		}


		if (fModelID == 3)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.001;
			variable[0] = 1.e-12;
			variable[1] = -1.5;
			variable[2] = 0.1;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Alpha",variable[1], step[1]);
			i_min->SetVariable(2,"Beta",variable[2], step[2]);
		}

		if (fModelID == 4)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.001;
			variable[0] = 1.e-12;
			variable[1] = -2.5;
			variable[2] = 1.0;

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);
			i_min->SetLimitedVariable(2,"Alpha",variable[2], step[2], 0., 2.0);
		}


		// do the minimization
		i_min->Minimize();
		const double *xs = i_min->X();
		i_min->Hesse();
		double *i_FitValues = new double[nParms];
		const double *i_Errors = i_min->Errors();

		for (int i = 0; i < nParms; i++)
		{
			i_FitValues[i] = xs[i];
		}


		fModel->SetParameters(xs);
		fModel->SetParErrors(i_Errors);
		TF1 *i_BestFit = 0;

		// Removing exp(tau) component
		if (fModel_intID != -999)
		{
			cout << "True doing EBL analysis" << endl;

			for (int i = 0; i < nParms; i++)
			{
				fModel_int->SetParameter(i, fModel->GetParameter(i));
				fModel_int->SetParError(i, fModel->GetParError(i));
			}
			i_BestFit = (TF1*)fModel_int->Clone();

		}

		else
		{
			// cout << "False Doing normal Analysis" << endl;
			i_BestFit = (TF1*)fModel->Clone();
		}


		double ErrUp = 0;
		double ErrLow = 0;

		for (int i = 0; i < nParms; i++)
		{
			bool bErrors = i_min->GetMinosError(i, ErrLow, ErrUp);
			cout << "Variabile " << i << " Error Status: " << bErrors << ", E_Low = " << ErrLow << ", E_Up = " << ErrUp << endl;
		}


		unsigned int nSteps = 100;
		Index_Scan = new double[nSteps];
		Index_Scan_Likelihood = new double[nSteps];
		Norm_Scan = new double[nSteps];
		Norm_Scan_Likelihood = new double[nSteps];


		i_min->Scan(0, nSteps, Norm_Scan, Norm_Scan_Likelihood, 1.e-13,1.e-11);
		i_min->Scan(1, nSteps, Index_Scan, Index_Scan_Likelihood, -3.5,-1.5);

		// Getting Confidince interval for power law






		// double i_energy = fFitMin;
		// vector<double> i_verr;
		// vector<double> i_venergy;
		// vector<double> i_vFluxMin;
		// vector<double> i_vFluxMax;
		// for (int j = 0 ; j < i_count; j++)
		// {

		// 	double x = TMath::Power(10,i_energy) ;
		// 	// double i_err = i_min->CovMatrix(0,0) * TMath::Power(x,i_FitValues[1]) * TMath::Power(x,i_FitValues[1]);
		// 	// // cout << i_err << "\t" ;
		// 	// i_err += i_min->CovMatrix(1,1) *  (i_FitValues[0] * TMath::Power(x ,i_FitValues[1])*log(x )) *(i_FitValues[0] * TMath::Power(x ,i_FitValues[1])*log(x ))  ;
		// 	// // cout << i_err << "\t" ;

		// 	// i_err += 2. * i_min->CovMatrix(0,1)  *  TMath::Power(x,i_FitValues[1]) * (i_FitValues[0] * TMath::Power(x ,i_FitValues[1])*log(x)) ;
		// 	// // cout << i_err << "\t" ;

		// 	// i_err = sqrt(i_err);
		// 	// cout << i_err << "\t" ;



		// 	double i_err = i_Errors[0] * i_Errors[0] / i_FitValues[0] / i_FitValues[0] ;
		// 	i_err += log(x/E_Norm) * log(x/E_Norm) * i_Errors[1] * i_Errors[1];
		// 	i_err += -2*log(x/E_Norm) / i_FitValues[0] * i_min->CovMatrix(0,1) ;
		// 	i_err = i_BestFit->Eval(i_energy) * sqrt(i_err);
		// 	// cout << i_err << "\t";
		// 	i_verr.push_back(i_err);
		// 	i_venergy.push_back(i_energy);
		// 	// double i_err = 0;
		// 	// double* i_energy = new double[2];
		// 	// i_energy &(fFitMin + j * delE);

		// 	// Getting error calculation using TF1::GradientPar

		// 	// double* p = new double;
		// 	// *p = fFitMin + j * delE;
		// 	// for (int k  = 0; k < nParms; k++)
		// 	// {
		// 	// 	for (int l = 0; l < nParms; l++)
		// 	// 	{
		// 	// 		i_err += fModel->GradientPar(k,p) * fModel->GradientPar(l,p) * i_min->CovMatrix(k,l) ;
		// 	// 	}

		// 	// }

		// 	// cout << "Setting Point: " << j << " " << i_energy << " " << i_BestFit->Eval(i_energy) << " " << i_err << endl;
		// 	fConfidinceInterval->SetPoint(j, i_energy, i_BestFit->Eval(i_energy) );
		// 	fConfidinceInterval->SetPointError(j , 0 , i_err );
		// 	// cout << j << " " << i_vFluxMax[j] << " " << fModel->Eval(i_venergy[j]) << " " << i_vFluxMin[j] << endl;

		// 	// if (fModel_intID == 0)
		// 	// {
		// 	// 	i_vFluxMin.push_back(fModel_int->Eval(i_energy) - i_err);
		// 	// 	i_vFluxMax.push_back(fModel_int->Eval(i_energy) + i_err);
		// 	// }

		// 	// else
		// 	// {
		// 	// 	i_vFluxMin.push_back(fModel->Eval(i_energy) - i_err);
		// 	// 	i_vFluxMax.push_back(fModel->Eval(i_energy) + i_err);
		// 	// }

		// 	// cout << fModel->Eval(i_energy) << "\t" << fModel->Eval(i_energy) + i_err << "\t" << fModel->Eval(i_energy) - i_err << endl;

		// 		i_energy += delE;
		// }


		cout << "Binned Likelihood Fit: \n" ;
		cout << "Parameter \t Best Fit \t Error\n";
		for ( int i = 0; i < nParms; i++)
		{
			cout << "[" << i << "]\t\t" << " " << i_min->VariableName(i) << " \t\t" << xs[i] << "\t\t" << i_Errors[i] << endl;
		}
		cout << endl;

		cout << "\n\nPrinting Covariance Matrix\n";
		double *i_covmat = new double [nParms*nParms];
		for (int i =0; i < nParms; i++)
		{
			for (int j = 0; j < nParms; j++)
			{
				i_covmat[ i*nParms + j ] = i_min->CovMatrix(i,j);
				cout << i_min->CovMatrix(i,j) << "\t" ;
			}
			cout << endl;
		}

		cout << "\n";
		// GetChi2(xs);
		getChi2 (xs);
		cout << "\n";
		cout << "Calculating Total Chi^2\n";
		getChi2Total (xs);
		cout << "\n";

		// fModel->Draw();
		// fModel->GetXaxis()->
		// gConfidinceInterval = GetConfidinceInterval();

		// int i_count = 10000;
		// double delE = (fFitMax - fFitMin)/i_count ;
		if (fConfidinceInterval)
		{
			delete fConfidinceInterval;
		}


		fConfidinceInterval = calculateConfidinceInterval(i_covmat, i_BestFit, fModelID, nParms);

		// fConfidinceInterval = new TGraphErrors(i_count);
		float *i_flux = getIntergralFlux(fFitMin, fFitMax , i_BestFit, true);

		cout << "Intergral Flux:\n";
		cout << "F (" << TMath::Power(10,fFitMin) << " TeV < E < " << TMath::Power(10,fFitMax) << ") = " << i_flux[0] << "+/-" << i_flux[1] << " [Photons/cm^2/s] \n";
		cout << "F (" << TMath::Power(10,fFitMin) << " TeV < E < " << TMath::Power(10,fFitMax) << ") = " << i_flux[2] << "+/-" << i_flux[3] << " [Crab] \n";


		// if (fModel_intID == 0 )
		// {
		// 	cout << "True" << endl;
		// 	fModel_int->SetChisquare(fModel->GetChisquare() );
		// 	fModel_int->SetNDF(fModel->GetNDF() );
		// 	return fModel_int;
		// }
		i_BestFit->SetChisquare(fModel->GetChisquare() );
		i_BestFit->SetNDF(fModel->GetNDF() );
		return i_BestFit;
	}

	// double VLikelihoodFitter::getChi2Total (const double *parms)
	// {

	// 	// double *parms2 = new double[nParms];
	// 	// parms2[0] = parms[0];
	// 	// parms2[1] = parms[1];
	// 	cout << "Getting getLogLi " << -1*getLogLi(parms) << endl;
	// 	double L = -1*getLogLiTotal(parms);
	// 	cout << "Got getLogLiTotal" << endl;
	// 	cout << "Getting getLogL0Total" << endl;
	// 	double L_0 = -1*getLogL0Total(parms);
	// 	cout << "Got getLogL0Total" << endl;

	// 	double i_TotalTime = 0;
	// 	for (int i = 0; i < fRunList.size(); i++)
	// 	{
	// 		i_TotalTime += fRunList[i].tOn * fRunList[i].deadTimeFraction;
	// 	}

	// 	cout << "Total Live time (s): " << i_TotalTime << endl
	// 		 << "Number of Bins: " << nBinsFit << endl
	// 		 << "Number of Fitting Parameters: " << nParms << endl
	// 		 << "Number of Degrees of Freedom: "<< nBinsFit - nParms << endl
	// 		 << "Chi^2 :" << -2*(L - L_0) << endl
	// 		 << "Chi^2/NDF : " << -2*(L - L_0) *(1.0/(nBinsFit - nParms)) << endl
	// 		 << "L : " << L << endl
	// 	     << "L_0 : " << L_0 << endl;

	// 	//      fModel->SetChisquare(-2*(L - L_0));
	// 	//      fModel->SetNDF((nBinsFit - nParms));


	// 	// return -2;
	// 	return -2*(L -1 *(L_0));
	// }





	// Function to get the runwise results
	int VLikelihoodFitter::getRunWiseFitInfo(int i_runnum, TF1* i_fit)
	{

		// Things I want to see in here:
		// L
		// L0
		// Chi^2
		// On, off, model on, model off

		// Checking if within MJD range

		if (fRunList[i_runnum].MJD < fMJD_Min ||  fRunList[i_runnum].MJD > fMJD_Max )
		{
			return -1;
		}

		// Checking if date is excluded
		for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
		{
			if ( fRunList[i_runnum].MJD > fExcludeMJD[e][0] &&  fRunList[i_runnum].MJD < fExcludeMJD[e][1] )
			{
				return -1;
			}
		}

		const double *parms = i_fit->GetParameters();
		double i_lRun = -1*getRunwiseLogL(i_runnum, parms);
		double i_l0Run = -1*getRunwiseLogL0(i_runnum, parms);
		double i_RunNBins = nBinsFit_runwise - nParms;
		// double i_RunOn = getRunwiseOnCounts(i_runnum);
		// double i_RunOff = getRunwiseOffCounts(i_runnum);
		// double i_RunModExcess = getRunwiseModelExcessCounts(i_runnum, parms);
		// double i_RunModOff = getRunwiseModelExcessCounts(i_runnum, parms);
		double i_Chi2Run = -2*(i_lRun - i_l0Run);


		cout << i_runnum << " " << fRunList[i_runnum].runnumber <<  " " << i_lRun << " " << i_l0Run << " " << i_Chi2Run << " " << i_RunNBins - nParms << " "  << i_Chi2Run/(i_RunNBins - nParms) << endl;
		     // << i_RunOn << " " << i_RunOff << " " << i_RunModExcess << " " << i_RunModOff << endl;

		return 0;

	}


	void VLikelihoodFitter::getRunwiseFitStatus(TF1 *i_fit)
	{
		// cout << i_runnum << " " << fRunList[i_runnum].runnumber <<  " " << i_lRun << " " << i_l0Run << " " << i_Chi2Run << " " << i_RunNBins << " "  << i_Chi2Run/i_RunNBins << endl;
		cout << "Run\tlog(L)\tlog(L_0)\tChi^2\tNDF\treduced Chi^2" << endl;
		for (unsigned int i = 0; i < fRunList.size(); i++)
		{
			getRunWiseFitInfo(i, i_fit);
		}
	}


	double VLikelihoodFitter::getRunwiseLogL(int i_run, const double *parms)
	{
		vector < vector <double> > i_myModel = getModelPredictedExcess( parms );
		vector < vector <double> > i_myvOffMLE = getModelPredictedOff( parms );
		double a, b, c, d;
		// double tmpNBinsFit = 0;
		bool i_bInc = true;
		nBinsFit_runwise = 0;
		double LogLi = 0;



		// Checking if run is excluded
		for (unsigned int e = 0; e < fExcludeRun.size(); e++)
		{
			if (fRunList[i_run].runnumber  == fExcludeRun[e])
			{
				i_bInc = false;
			}
		}

		// Checking if within MJD range
		if (fRunList[i_run].MJD < fMJD_Min ||  fRunList[i_run].MJD > fMJD_Max )
		{
			return -1;
		}

		// Checking if date is excluded
		for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
		{
			if ( fRunList[i_run].MJD > fExcludeMJD[e][0] &&  fRunList[i_run].MJD < fExcludeMJD[e][1] )
			{
				i_bInc = false;
				break;
			}
		}

		for (int j = 0; j < nMC; j ++)
		{

			if (i_bInc == false){continue;}
			if (fRecBinCentres[j] < fFitMin){continue;}
			if (fRecBinCentres[j] > fFitMax) {continue;}
			// if (j >= vLastOnCount[i] ){continue;}

			nBinsFit_runwise += 1;
			if (fOnCounts[i_run][j] > 0)
			{
				a = fOnCounts[i_run][j]*log(i_myModel[i_run][j] + fRunList[i_run].alpha*i_myvOffMLE[i_run][j]);
			}

			else
			{
				a = 0;
			}


			if ( fOffCounts[i_run][j] > 0 )
			{
				b = fOffCounts[i_run][j]*log(i_myvOffMLE[i_run][j]);
			}

			else
			{
				b = 0;

			}
			// Getting c term
			c = -1.0*(fRunList[i_run].alpha + 1.0)*i_myvOffMLE[i_run][j];

			// Getting d term
			d = -i_myModel[i_run][j];

			LogLi = LogLi + a + b + c + d;

		}

		return -1 * LogLi;

	}

	double VLikelihoodFitter::getRunwiseLogL0(int i_run, const double *parms)
	{

		double LogLi = 0;
		vector < vector <double> > i_myModel = getModelPredictedExcess( parms );
		vector < vector <double> > i_myvOffMLE = getModelPredictedOff( parms );

		double a, b, c, d;
		bool i_bInc = false;

		i_bInc = true;

		// Checking if run is excluded
		for (unsigned int e = 0; e < fExcludeRun.size(); e++)
		{
			if (fRunList[i_run].runnumber  == fExcludeRun[e])
			{
				i_bInc = false;
			}
		}

		// Checking if within MJD range
		if (fRunList[i_run].MJD < fMJD_Min ||  fRunList[i_run].MJD > fMJD_Max )
		{
			return -1;
		}

		// Checking if date is excluded
		for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
		{
			if ( fRunList[i_run].MJD > fExcludeMJD[e][0] &&  fRunList[i_run].MJD < fExcludeMJD[e][1] )
			{
				i_bInc = false;
				break;
			}
		}


		for (int j = 0; j < nMC; j ++)
		{

			if (i_bInc == false){continue;}
			if (fRecBinCentres[j] < fFitMin){continue;}
			if (fRecBinCentres[j] > fFitMax) {continue;}
			// if (j >= vLastOnCount[i] ){continue;}
			// if (fEnergyBias[i][j] > fThresholdBias){continue;}

			// Getting a term
			// if (fOnCounts[i][j] == 0 && fOffCounts[i][j] == 0){continue;}
			// nBinsFit += 1;




			if (fOnCounts[i_run][j] < 1)
			{
				a = 0;
			}
			else
			{
				a = fOnCounts[i_run][j] * log( fOnCounts[i_run][j] );

			}

			if (fOffCounts[i_run][j] < 1)
			{
				b = 0;
			}
			else
			{
				b = fOffCounts[i_run][j] * log( fOffCounts[i_run][j] );

			}

			c = -1*fOnCounts[i_run][j];
			d = -1*fOffCounts[i_run][j];

			LogLi += a + b + c + d;
		}


		return -1 * LogLi;
	}

	// VLikelihoodFitter::getRunwiseNBins()
	// {

	// }

	// VLikelihoodFitter::getRunwiseOnCounts()
	// {

	// }

	// VLikelihoodFitter::getRunwiseOffCounts()
	// {

	// }

	// VLikelihoodFitter::getRunwiseModelExcessCounts()
	// {

	// }

	// vector <double> VLikelihoodFitter::getRunwiseModelExcessCounts(int i_run, const double *parms)
	// {
	// 			// cout << "Getting Predicted Excess\n";
	// 	fModel->SetParameters(parms);

	// 	vector <double> i_TrueCenters;

	// 	// Getting Spectrally weighted bin centres
	// 	for (int i = 0; i < nMC; i++)
	// 	{
	// 		i_TrueCenters.push_back( VMathsandFunctions::getSpectralWeightedMeanEnergy(fMCBins[i], fMCBins[i+1], parms[1]) );
	// 	}


	// 	double i_ModelElement = 0;
	// 	double i_EffetiveAreaElement = 0;
	// 	double i_MatrixElement = 0;
	// 	double i_ConversionElement = 1.e4;
	// 	double i_EnergyElement = 0;
	// 	double i_TOn = 0;

	// 	// Current Bin
	// 	int i_J_Rec_Cur = 0;
	// 	int i_L_MC_Cur = 0;

	// 	TAxis *xpoints = fResponseMatrixRebinned[i_run]->GetXaxis();
	// 	TAxis *ypoints = fResponseMatrixRebinned[i_run]->GetYaxis();
	// 	vector <double> i_vTmp;

	// 	// Looping over each Rec Energy bin
	// 	for (int j = 0; j < nRec; j++)
	// 	{
	// 		if (fRecBinCentres[j] < i_eThresh)
	// 		{
	// 			continue;
	// 		}
	// 		i_vTmp.push_back(0);
	// 		i_J_Rec_Cur = xpoints->FindBin(fRecBinCentres[j]);

	// 		for (int l = 0; l < nMC; l++)
	// 		{
	// 			i_L_MC_Cur = ypoints->FindBin(fMCBinCentres[l]);

	// 			// Lower range on sanity
	// 			if (fEnergyBias[i_run][l] > fThresholdBias){continue;}

	// 			i_ModelElement = fModel->Eval(fMCBinCentres[l]);

	// 			i_EffetiveAreaElement = fMeanEffectiveAreaMC[i_run]->Eval(fMCBinCentres[l]);

	// 			i_MatrixElement = fResponseMatrixRebinned[i_run]->GetBinContent(i_J_Rec_Cur, i_L_MC_Cur);
	// 			i_TOn = fRunList[i_run].tOn * fRunList[i_run].deadTimeFraction  ;

	// 			i_EnergyElement = pow(10.0,fMCBins[l+1]) - pow(10.0,fMCBins[l]);
	// 			i_vTmp[j] +=   i_ModelElement * i_EffetiveAreaElement * i_MatrixElement * i_ConversionElement *  i_TOn  * i_EnergyElement ;
	// 		}

	// 	}
	// 	return i_vTmp;

	// }


	// // To get L0 we need to evaluate the Loglikelihood function at our expected values. I.e. the observed coutns are the predicted ones
	// double VLikelihoodFitter::getLogL0 (const double *parms)
	// {


	// 	// cout << "Getting getLogLi\n";

	// 	double LogLi = 0;
	// 	// double Alpha = 1.0/8.0;
	// 	// vector < vector <double> > i_myModel = getModelPredictedExcess( parms );
	// 	// vector < vector <double> > i_myvOffMLE = getModelPredictedOff( parms );

	// 	double a, b, c, d;
	// 	double tmpNBinsFit = 0;
	// 	bool i_bInc = true;
	// 	nBinsFit = 0;

	// 	for (int i = 0; i < fRunList.size(); i++){


	// 		if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
	// 		{
	// 			continue;
	// 		}

	// 		i_bInc = true;

	// 		for (int e = 0; e < fExcludeMJD.size(); e++)
	// 		{
	// 			if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
	// 			{
	// 				i_bInc = false;
	// 				break;
	// 			}
	// 		}

	// 		for (int e = 0; e < fExcludeRun.size(); e++)
	// 		{
	// 			if (fRunList[i].runnumber == fExcludeRun[e])
	// 			{
	// 				i_bInc = false;
	// 				break;
	// 			}
	// 		}

	// 		for (int j = 0; j < nMC; j ++)
	// 		{

	// 			if (i_bInc == false){continue;}
	// 			if (fRecBinCentres[j] < fFitMin){continue;}
	// 			if (fRecBinCentres[j] > fFitMax) {continue;}
	// 			// if (j >= vLastOnCount[i] ){continue;}
	// 			// if (fEnergyBias[i][j] > fThresholdBias){continue;}

	// 			// Getting a term
	// 			// if (fOnCounts[i][j] == 0 && fOffCounts[i][j] == 0){continue;}
	// 			nBinsFit += 1;




	// 			// if (fOnCounts[i][j] < 1)
	// 			// {
	// 			// 	a = 0;
	// 			// }
	// 			// else
	// 			// {
	// 			// 	a = fOnCounts[i][j] * log( fOnCounts[i][j] );

	// 			// }

	// 			if (fOnCounts[i][j] > 0)
	// 			{
	// 				a = fOnCounts[i][j] * log( fOnCounts[i][j] );
	// 			}

	// 			else
	// 			{
	// 				a = 0;
	// 			}

	// 			// if (fOffCounts[i][j] < 1)
	// 			// {
	// 			// 	b = 0;
	// 			// }
	// 			// else
	// 			// {
	// 			// 	b = fOffCounts[i][j] * log( fOffCounts[i][j] );

	// 			// }


	// 			if ( fOffCounts[i][j] > 0 )
	// 			{
	// 				b = fOffCounts[i][j] * log( fOffCounts[i][j] );
	// 			}

	// 			else
	// 			{
	// 				b = 0;

	// 			}

	// 			c = -1*fOnCounts[i][j];
	// 			d = -1*fOffCounts[i][j];

	// 			LogLi = LogLi + a + b + c + d;
	// 		}
	// 	}

	// 	return -1 * LogLi;

	// }



	double VLikelihoodFitter::getLogLi(const double *parms)
	{

		double LogLi = 0;
		vector < vector <double> > i_myModel = getModelPredictedExcess( parms );
		vector < vector <double> > i_myvOffMLE = getModelPredictedOff( parms );

		double a, b, c, d;
		bool i_bInc = true;
		nBinsFit = 0;

		for (unsigned int i = 0; i < fRunList.size(); i++){


			// Excluding based on MJD Min/Max
			if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
			{
				continue;
			}

			i_bInc = true;
			// Looking for Excluded MJD datas
			for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
			{
				if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
				{
					i_bInc = false;
					break;
				}
			}
			// Looking for excluded Runs
			for (unsigned int e = 0; e < fExcludeRun.size(); e++)
			{
				if (fRunList[i].runnumber == fExcludeRun[e])
				{
					i_bInc = false;
					break;
				}
			}

			// Looping over data bins
			for (int j = 0; j < nMC; j ++)
			{
				// If run is excluded
				if (i_bInc == false){continue;}
				// Fit min/max
				if (fRecBinCentres[j] < fFitMin){continue;}
				if (fRecBinCentres[j] > fFitMax) {continue;}
				// if (j >= vLastOnCount[i] ){continue;}

				nBinsFit += 1;

				if (fOnCounts[i][j] > 0)
				{
					a = fOnCounts[i][j]*log(i_myModel[i][j] + fRunList[i].alpha*i_myvOffMLE[i][j]);
				}

				else
				{
					a = 0;
				}


				if ( fOffCounts[i][j] > 0 )
				{
					b = fOffCounts[i][j]*log(i_myvOffMLE[i][j]);
				}

				else
				{
					b = 0;

				}
				// Getting c term
				c = -1.0*(fRunList[i].alpha + 1.0)*i_myvOffMLE[i][j];

				// Getting d term
				d = -i_myModel[i][j];

				LogLi = LogLi + a + b + c + d;
			}

		}
		return -1 * LogLi;
	}



	double VLikelihoodFitter::getLogL0 (const double *parms)
	{

		double LogLi = 0;


		double a, b, c, d;
		bool i_bInc = true;
		nBinsFit = 0;

		for (unsigned int i = 0; i < fRunList.size(); i++){

			// Checking if within accepted MJD
			if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
			{
				continue;
			}

			i_bInc = true;

			// Checking if within excluded MJD
			for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
			{
				if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
				{
					i_bInc = false;
					break;
				}
			}
			// checking if excluded run
			for (unsigned int e = 0; e < fExcludeRun.size(); e++)
			{
				if (fRunList[i].runnumber == fExcludeRun[e])
				{
					i_bInc = false;
					break;
				}
			}

			// Looping over data bins
			for (int j = 0; j < nMC; j ++)
			{

				// If excluded
				if (i_bInc == false){continue;}
				// Checking if within min and max
				if (fRecBinCentres[j] < fFitMin){continue;}
				if (fRecBinCentres[j] > fFitMax) {continue;}

				// if (j >= vLastOnCount[i] ){continue;}
				// if (fEnergyBias[i][j] > fThresholdBias){continue;}

				nBinsFit += 1;

				if (fOnCounts[i][j] > 0)
				{
					a = fOnCounts[i][j] * log( fOnCounts[i][j] );
				}

				else
				{
					a = 0;
				}



				if ( fOffCounts[i][j] > 0 )
				{
					b = fOffCounts[i][j] * log( fOffCounts[i][j] );
				}

				else
				{
					b = 0;

				}

				c = -1*fOnCounts[i][j];
				d = -1*fOffCounts[i][j];

				LogLi = LogLi + a + b + c + d;
			}
		}

		return -1 * LogLi;

	}

	// Function to return the likelihood ratio. -2ln(L/L_0) which is equivilant to a chi2 statistic
	double VLikelihoodFitter::getChi2 (const double *parms)
	{


		double L = -1*getLogLi(parms);
		double L_0 = -1*getLogL0(parms);

		double i_TotalTime = 0;
		for (unsigned int i = 0; i < fRunList.size(); i++)
		{
			i_TotalTime += fRunList[i].tOn * fRunList[i].deadTimeFraction;
		}

		cout << "Total Live time (s): " << i_TotalTime << endl
			 << "Number of Bins: " << nBinsFit << endl
		     << "Number of Fitting Parameters: " << nParms << endl
		     << "Number of Degrees of Freedom: "<< nBinsFit - nParms << endl
		     << "Chi^2 :" << -2*(L - L_0) << endl
		     << "Chi^2/NDF : " << -2*(L - L_0) *(1.0/(nBinsFit - nParms)) << endl
		     << "L : " << L << endl
		     << "L_0 : " << L_0 << endl;

		// Setting Model Parameters
		fModel->SetChisquare(-2*(L - L_0));
		fModel->SetNDF((nBinsFit - nParms));



		return -2*(L -1 *(L_0));
	}



	double VLikelihoodFitter::getLogLiTotal(const double *parms)
	{

		vector < vector <double> > i_myModel = getModelPredictedExcess( parms );
		vector < vector <double> > i_myModelOff = getModelPredictedOff( parms );

		double LogLi = 0;
		vector <double> i_total_On(fOnCounts[0].size());
		vector <double> i_total_Off(fOnCounts[0].size());
		vector <double> i_total_Model(fOnCounts[0].size());
		vector <double> i_total_ModelOff(fOnCounts[0].size());
		double i_mean_alpha = 0;
		int i_counter = 0;

		double a, b, c, d;
		bool i_bInc = true;
		nBinsFit_Total = 0;

		for (unsigned int i = 0; i < fRunList.size(); i++){


			// Excluding based on MJD Min/Max
			if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
			{
				continue;
			}

			i_bInc = true;
			// Looking for Excluded MJD datas
			for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
			{
				if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
				{
					i_bInc = false;
					break;
				}
			}
			// Looking for excluded Runs
			for (unsigned int e = 0; e < fExcludeRun.size(); e++)
			{
				if (fRunList[i].runnumber == fExcludeRun[e])
				{
					i_bInc = false;
					break;
				}
			}

			if (i_bInc == false){continue;}


			for (int j = 0; j < nMC; j++)
			{
				i_total_On[j] += fOnCounts[i][j];
				i_total_Off[j] += fOffCounts[i][j];
				i_total_Model[j] += i_myModel[i][j];
				i_total_ModelOff[j] += i_myModelOff[i][j];

				// Run average Alpha (should be time averaged..... TODO LIST!)
				// cout << i_total_On[j] << " " << i_total_Off[j] << " " << i_total_Model[j] << " " << i_total_ModelOff[j] << endl;
			}
			i_mean_alpha += fRunList[i].alpha;

			i_counter +=1 ;
		}

		i_mean_alpha /= i_counter ;
		// cout << "Mean Alpha: " << i_mean_alpha << endl;

		// for (int j = 0; j < nMC; j ++)
		// {
		// 	// cout << i_total_On[j] << " " << i_total_Off[j] << " " << i_total_Model[j] << " " << i_total_ModelOff[j] <<  endl;
		// }
		// Looping over data bins
		for (int j = 0; j < nMC; j ++)
		{
			// If run is excluded
			// Fit min/max
			if (fRecBinCentres[j] < fFitMin){continue;}
			if (fRecBinCentres[j] > fFitMax) {continue;}
			// if (j >= vLastOnCount[i] ){continue;}

			nBinsFit_Total++;

			if (i_total_On[j] >= 1)
			{
				a = i_total_On[j]*log(i_total_Model[j] + i_mean_alpha*i_total_ModelOff[j]);
			}

			else
			{
				a = 0;
			}



			if ( i_total_Off[j] >= 1 )
			{
				b = i_total_Off[j]*log(i_total_ModelOff[j]);
			}

			else
			{
				b = 0;

			}
			// Getting c term
			c = -1.0*(i_mean_alpha + 1.0)*i_total_ModelOff[j];

			// Getting d term
			d = -i_total_Model[j];

			// cout << "a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
			LogLi = LogLi + a + b + c + d;
		}

		// cout << "LogLi : " << LogLi  << endl;
		return -1 * LogLi;
	}




	double VLikelihoodFitter::getLogL0Total (const double *parms)
	{

		double LogLi = 0;

		vector <double> i_total_On(fOnCounts[0].size());
		vector <double> i_total_Off(fOnCounts[0].size());

		double a, b, c, d;
		bool i_bInc = true;
		double i_mean_alpha = 0;
		int i_counter = 0;


		for (unsigned int i = 0; i < fRunList.size(); i++){

			// Checking if within accepted MJD
			if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
			{
				continue;
			}

			i_bInc = true;

			// Checking if within excluded MJD
			for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
			{
				if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
				{
					i_bInc = false;
					break;
				}
			}
			// checking if excluded run
			for (unsigned int e = 0; e < fExcludeRun.size(); e++)
			{
				if (fRunList[i].runnumber == fExcludeRun[e])
				{
					i_bInc = false;
					break;
				}
			}

			if (i_bInc == false){continue;}


			for (int j = 0; j < nMC; j++)
			{
				i_total_On[j] += fOnCounts[i][j];
				i_total_Off[j] += fOffCounts[i][j];

				// Run average Alpha (should be time averaged..... TODO LIST!)
				i_mean_alpha += fRunList[i].alpha;
				// cout << i_total_On[j] << " " << i_total_Off[j] <<  endl;


			}
			i_counter += 1;
		}

		i_mean_alpha /= i_counter ;


		// // Printing Counts in Bins:
		// for (int j = 0; j < nMC; j ++)
		// {
		// 	cout << i_total_On[j] << " " << i_total_Off[j] <<  endl;
		// }
		// Looping over data bins
		for (int j = 0; j < nMC; j ++)
		{

			// If excluded
			// Checking if within min and max
			if (fRecBinCentres[j] < fFitMin){continue;}
			if (fRecBinCentres[j] > fFitMax) {continue;}

			// if (j >= vLastOnCount[i] ){continue;}
			// if (fEnergyBias[i][j] > fThresholdBias){continue;}


			if (i_total_On[j] >= 1)
			{
				a = i_total_On[j] * log( i_total_On[j] );
			}

			else
			{
				a = 0;
			}



			if ( i_total_Off[j] >= 1 )
			{
				b = i_total_Off[j] * log( i_total_Off[j] );
			}

			else
			{
				b = 0;

			}

			c = -1*int(i_total_On[j]);
			d = -1*int(i_total_Off[j]);
			// cout << "a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
			LogLi = LogLi + a + b + c + d;
		}


		return -1 * LogLi;

	}


	double VLikelihoodFitter::getChi2Total (const double *parms)
	{

		// double *parms2 = new double[nParms];
		// parms2[0] = parms[0];
		// parms2[1] = parms[1];
		// cout << "Getting getLogLi " << -1*getLogLi(parms) << endl;
		double L = -1*getLogLiTotal(parms);
		cout << "Got getLogLiTotal: " << L << endl;
		cout << "Getting getLogL0Total" << endl;
		double L_0 = -1*getLogL0Total(parms);
		cout << "Got getLogL0Total" << endl;

		double i_TotalTime = 0;
		for (unsigned int i = 0; i < fRunList.size(); i++)
		{
			i_TotalTime += fRunList[i].tOn * fRunList[i].deadTimeFraction;
		}

		cout << "Total Live time (s): " << i_TotalTime << endl
			 << "Number of Bins: " << nBinsFit_Total << endl
			 << "Number of Fitting Parameters: " << nParms << endl
			 << "Number of Degrees of Freedom: "<< nBinsFit_Total - nParms << endl
			 << "Chi^2 :" << -2*(L - L_0) << endl
			 << "Chi^2/NDF : " << -2*(L - L_0) *(1.0/(nBinsFit_Total - nParms)) << endl
			 << "L : " << L << endl
		     << "L_0 : " << L_0 << endl;

		//      fModel->SetChisquare(-2*(L - L_0));
		//      fModel->SetNDF((nBinsFit - nParms));


		// return -2;
		return -2*(L -1 *(L_0));
	}




	double VLikelihoodFitter::getNDF()
	{
		return (nBinsFit - nParms);
	}


	// Plotting Energy Bias for a given run
	TCanvas *VLikelihoodFitter::plotEnergyBias(int i_Entry)
	{


		TCanvas *cEnergyBias = new TCanvas();

		vector <double> i_bias = fEnergyBias[i_Entry];
		vector <double> i_Energy = fMCBinCentres;




		TGraphAsymmErrors *gBias = new TGraphAsymmErrors(i_Energy.size(), &(i_Energy[0]), &(i_bias[0]));
		sprintf( hname, "Energy Bias for %d", fRunList[i_Entry].runnumber);

		gBias->SetTitle(hname);
		gBias->GetXaxis()->SetTitle("Log(E_{MC}) (TeV)");
		gBias->GetYaxis()->SetTitle("Energy Bias (E_{MC} - #bar{E_{Rec}})/E_{MC}");
		gBias->SetMarkerSize(1.0);
		gBias->SetMarkerStyle(8);
		gBias->SetMarkerColor(kBlack);
		gBias->Draw("APL");


		TLine *lThresh = new TLine(i_Energy[0], fThresholdBias, i_Energy[i_Energy.size() -1], fThresholdBias);
		lThresh->SetLineColor(kRed);
		lThresh->Draw("SAME");

		return cEnergyBias;
	}


	// EBL Stuff

	// Set the source redshift
	void VLikelihoodFitter::setRedShift(double i_Z)
	{
		Z = i_Z;
	}

	// Set the EBL Opacity
	int VLikelihoodFitter::setEBLOpacity(TGraph2D *i_OpacityTable)
	{
		TH2D *hOPTable = (TH2D*) i_OpacityTable->GetHistogram();

		int i_ZBinNum ;

		// Currently not interpolated between Redshift
		i_ZBinNum = hOPTable->GetYaxis()->FindBin(Z);

		// Get Slice at the redshift

		sprintf(hname, "Projection Z = %f" , Z);
		TH1D *hTmp = hOPTable->ProjectionX(hname, i_ZBinNum, i_ZBinNum +1);


		vector <double> xValues, yValues;
		for (int i = 1 ; i <= hTmp->GetNbinsX(); i++)
		{
			xValues.push_back( hTmp->GetBinCenter(i) );
			yValues.push_back( hTmp->GetBinContent(i) );
		}
		// Creating a TGraph to allow quick interpolating (using TGraph::Eval(X))
		fEBLOpacityGraph = new TGraph(xValues.size(), &(xValues[0]),&(yValues[0]));


		if (!fEBLOpacityGraph)
		{
			return 1;
		}

		return 0;
	}

	void VLikelihoodFitter::setEBLOpacity (TGraph* i_EBLOpacity)
	{
		fEBLOpacityGraph = (TGraph*)i_EBLOpacity->Clone();
	}


	// Fit fuction that takes into account EBL attenuation
	double VLikelihoodFitter::calculateIntrinsicSpectrum(Double_t *x, Double_t *parm)
	{
		for (int i = 0; i < nParms; i++)
		{
			// Skipping case of scaled EBL power law model
			if ( i == nParms -1 && fModelID == 4 )
			{
				continue;
			}
			fModel_int->SetParameter(i,parm[i]);
		}
		// Evaluating Note Spectral Model is in log and EBL lookup table is in linear space
		if (fModelID == 4 )
		{
			cout  << scientific << parm[nParms -1] << " "  << x[0] << " " << -parm[nParms -1]*fEBLOpacityGraph->Eval(TMath::Power(10,x[0])) << " " << fModel_int->Eval(x[0])*TMath::Exp(-parm[nParms -1]*fEBLOpacityGraph->Eval(TMath::Power(10,x[0]))) << endl;
			return fModel_int->Eval(x[0])*TMath::Exp(-parm[nParms -1]*fEBLOpacityGraph->Eval(TMath::Power(10,x[0])));

		}

		else
		{
			return fModel_int->Eval(x[0])*TMath::Exp(-1*fEBLOpacityGraph->Eval(TMath::Power(10,x[0])));

		}
	}


	// Setting the intrinsic Source Model
	void VLikelihoodFitter::setIntrinsicModel(int i_ID = 0)
	{

		// Note for spectral weighting of bins, the spectral index must be the 2nd parameter
		// i.e. [1] = Spectral index
		fModel_intID = i_ID;
		if (i_ID == 0)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1] )", E_Norm);
			fModel_int = new TF1("fModel_int",hname,-1.0,5.0);
			nParms = 2;
			fModelID = 0;
		}

		// Power Law with Exponential Cut off
		else  if (i_ID == 1)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f , [1] ) * TMath::Exp( -1. * TMath::Power( 10, x ) / [2] )", E_Norm);
			fModel_int = new TF1("fModel_int",hname,-1,5);
			nParms = 3;
			fModelID = 1;

		}


		// Curved Spectrum
		else if (i_ID == 2)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1]+[2]*TMath::Power( 10, x )  )", E_Norm);
			fModel_int = new TF1("fModel_int",hname,-1,5);
			nParms = 3;
			fModelID = 2;

		}

		else if (i_ID == 3)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1]+[2]*TMath::Log( TMath::Power( 10, x ) ) )", E_Norm);
			fModel_int = new TF1("fModel", hname,-1.0,5.0);
			fModelID = i_ID;
			nParms = 3;

		}


		// Scaled Alpha Power Law Spectrum
		else if (i_ID == 4)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1] )", E_Norm);
			fModel_int = new TF1("fModel_int",hname,-1.0,5.0);
			nParms = 3;
			fModelID = 4;

		}

		else
		{
			cout << "Model " << i_ID << " not found!\n";
		}

		// Setting fit model to EBL corrected model
		setModel(-999);
	}


	TGraphErrors* VLikelihoodFitter::calculateConfidinceInterval( double* i_covmat, TF1 *i_fitfunction, int i_model, int i_nparms)
	{

		cout << "Model Details: " << endl << i_model << " " << fModelID << " " << fModel_intID << endl;
		int i_nsteps = 1000;
		double i_start = fFitMin;
		double i_stop = fFitMax;

		double i_step = (i_stop - i_start) / i_nsteps;

		TGraphErrors *i_ConfidinceInterval = new TGraphErrors(i_nsteps);

		double *i_deriv = new double [i_nparms];

		for ( int i = 0; i < i_nsteps; i++ )
		{
			double i_energy = i_start + i*i_step ;
			double i_flux_err = 0;

			for ( int j = 0; j < i_nparms ; j++ )
			{
				for ( int k = 0; k < i_nparms; k++ )
				{

					// Power Law Model
					if (i_model == 0 )
					{
						// dF/dN
						// (E/E_0)^gamma
						i_deriv[0] = TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) );
						// dF/dGamma
						// -1 * N (E/E_0)^gamma log(E/E_0)
						i_deriv[1] = i_fitfunction->GetParameter(0) * TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) ) * TMath::Log(TMath::Power(10., i_energy) / E_Norm + 1.e-9);
					}

					// Power Law with Exp Cutoff
					if (i_model == 1  )
					{
						// dF/dN
						// (E/E_0)^gamma * exp(-E/E_C)
						i_deriv[0] = TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) ) * TMath::Exp(-1*TMath::Power(10., i_energy) / i_fitfunction->GetParameter(2) );

						// dF/dGamma
						// -1 * N (E/E_0)^gamma log(E/E0) * exp(-E/E_C)
						i_deriv[1] = i_fitfunction->GetParameter(0) * TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) );
						i_deriv[1] *= TMath::Log(TMath::Power(10., i_energy) / E_Norm + 1.e-9) * TMath::Exp(-1*TMath::Power(10., i_energy) / i_fitfunction->GetParameter(2) );

						// dF/dE_cutoff
						// N (E/E_0)^gamma *(-1)*exp(-E/E_C) * ( E / E_C^2)
						i_deriv[2] = i_fitfunction->GetParameter(0) * TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) );
						i_deriv[2] *=   TMath::Exp(TMath::Power(10., i_energy) / i_fitfunction->GetParameter(2) ) * TMath::Power(10., i_energy) / i_fitfunction->GetParameter(2) / i_fitfunction->GetParameter(2);

					}

					// Curved Power Law
					if (i_model == 2 )
					{
						// dF/dN
						// (E/E_0)^(alpha + beta *E)
						i_deriv[0] = TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) +  i_fitfunction->GetParameter(2) * TMath::Power(10., i_energy) );

						// dF/dalpha
						// N (E/E_0)^(alpha + beta*E) * log(E/E_0) *(-1)
						i_deriv[1] = i_fitfunction->GetParameter(0) * TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) +  i_fitfunction->GetParameter(2) * TMath::Power(10., i_energy) );
						i_deriv[1] *= TMath::Log( TMath::Power(10., i_energy) / E_Norm + 1.e-9);

						// dF/dBeta
						// N (E/E_0) ^(alpha + beta*E) * E *log(E/E_0)
						i_deriv[2] = i_fitfunction->GetParameter(0) * TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) +  i_fitfunction->GetParameter(2) * TMath::Power(10., i_energy) );
						i_deriv[2] *= TMath::Power(10., i_energy) * TMath::Log( TMath::Power(10., i_energy) / E_Norm + 1.e-9);
					}


					// log parabola
					if (i_model == 3 )
					{
						// dF/dN
						// (E/E_0)^(alpha + beta * log(E) )
						i_deriv[0] = TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) +  i_fitfunction->GetParameter(2) * TMath::Log( TMath::Power(10., i_energy) + 1.e-9) );

						// dF/dalpha
						// N (E/E_0) ^(alpha + beta * log(E) ) * -1 * log(E/E_0)
						i_deriv[1] = i_fitfunction->GetParameter(0) * TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) +  i_fitfunction->GetParameter(2) * TMath::Log( TMath::Power(10., i_energy) + 1.e-9) );
						i_deriv[1] *= TMath::Log( TMath::Power(10., i_energy) / E_Norm + 1.e-9);

						// dF/dBeta
						// N * (E/E_0)^(alpha + beta * log(E) ) * log(E) * log(E/E_0)
						i_deriv[2] = i_fitfunction->GetParameter(0) * TMath::Power( TMath::Power(10., i_energy) / E_Norm , i_fitfunction->GetParameter(1) +  i_fitfunction->GetParameter(2) * TMath::Log( TMath::Power(10., i_energy) + 1.e-9) );
						i_deriv[2] *= TMath::Log(TMath::Power(10., i_energy) + 1.e-9) * TMath::Log( TMath::Power(10., i_energy) / E_Norm + 1.e-9);

					}

					i_flux_err += i_covmat[ j*i_nparms+k ] * i_deriv[j] * i_deriv[k];

				}

			}




			double i_flux = i_fitfunction->Eval(i_energy);

			if (fModel_intID == -999)
			{
				i_ConfidinceInterval->SetPoint(i, i_energy, i_flux);
				i_ConfidinceInterval->SetPointError(i, 0, TMath::Sqrt(i_flux_err) );
			}

			else
			{
				i_ConfidinceInterval->SetPoint(i, i_energy, i_flux  );
				i_ConfidinceInterval->SetPointError(i, 0, TMath::Sqrt(i_flux_err) ) ;// * TMath::Exp(fEBLOpacityGraph->Eval( TMath::Power(10,i_energy) ) ) );

				// cout << i <<  " " <<  i_energy << " " << TMath::Sqrt(i_flux_err) << " " << fEBLOpacityGraph->Eval( TMath::Power(10,i_energy) ) << " " << TMath::Exp(fEBLOpacityGraph->Eval( TMath::Power(10,i_energy) ) ) << endl;
				// i_ConfidinceInterval->SetPoint(i, i_energy, i_flux);
				// i_ConfidinceInterval->SetPointError(i, 0, TMath::Sqrt(i_flux_err) );
			}


		}

		delete []i_deriv;

		return i_ConfidinceInterval;
	}




	void VLikelihoodFitter::printCountsInBins()
	{
		vector <double> i_TotalCountsOn(fOnCounts[0].size()	);
		vector <double> i_TotalCountsOff(fOnCounts[0].size());


		for (unsigned int i = 0; i < fRunList.size(); i++)
		{
			for (unsigned int j = 0; j < fOnCounts[i].size(); j++)
			{
				i_TotalCountsOn[j] += fOnCounts[i][j];
				i_TotalCountsOff[j] += fOffCounts[i][j];
			}

		}


		for (unsigned int i = 0 ; i < i_TotalCountsOn.size(); i++)
		{
			cout << TMath::Power(10.0,fRecBinCentres[i])	 << " " << i_TotalCountsOn[i] << " " << i_TotalCountsOff[i] << endl;
		}

	}




	// Function to apply fit to only the energy bin of interest
	float* VLikelihoodFitter::getSpectralPoint( double *parms , double BinMin, double BinMax, double iE_Norm, TF1* iBestFit)
	{


		setFitMinMax(BinMin, BinMax);


		if (fModelID == 0)
		{
			sprintf( hname, "[0]*TMath::Power(TMath::Power(10.0,x) / %f, [1])", iE_Norm);
			fModel = new TF1("fModel", hname,-1.0,5.0);
			nParms = 2;
		}

		// Power Law with Exponential Cut off
		else  if (fModelID == 1)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1] ) * TMath::Exp( -1. * TMath::Power( 10, x ) / [2] )", iE_Norm);
			fModel = new TF1("fModel", hname,-1.0,5.0);
			nParms = 3;
		}


		// Curved Spectrum
		else if (fModelID == 2)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x )/ %f, [1]+[2]*TMath::Power( 10, x ) )", iE_Norm);
			fModel = new TF1("fModel", hname,-1.0,5.0);
			nParms = 3;

		}


		// Log-parabola
		else if (fModelID == 3)
		{
			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1]+[2]*TMath::Log( TMath::Power( 10, x ) ) )", iE_Norm);
			fModel = new TF1("fModel", hname,-1.0,5.0);
			nParms = 3;

		}




		//nParms = 2;
		ROOT::Math::Minimizer* i_min =
		   ROOT::Math::Factory::CreateMinimizer("Minuit", "Minos");

		// // set tolerance , etc...
		i_min->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
		i_min->SetMaxIterations(100000);  // for GSL
		i_min->SetTolerance(0.001);
		i_min->SetPrintLevel(1);
		i_min->SetErrorDef(0.5);
		// // // create funciton wrapper for minmizer
		// // // a IMultiGenFunction type

		ROOT::Math::Functor i_fitfunction(this,&VLikelihoodFitter::getLogLi,nParms);
		// ROOT::Math::Functor i_fitfunction(this,nParms, "getLikelihoodFit");
		double *step = new double[nParms];
		double *variable = new double[nParms];

		i_min->SetFunction(i_fitfunction);


		if (fModelID == 0)
		{
			step[0] = 1.e-17;
			step[1] = 0.01;
			variable[0] = 1.e-16;
			variable[1] = parms[1];
			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);

			i_min->FixVariable(1);

		}


		if (fModelID == 1)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.01;
			variable[0] = 1.e-12;
			variable[1] = parms[1];
			variable[2] = parms[2];
			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);
			i_min->SetVariable(2,"E_CutOff", variable[2], step[2]);

			i_min->FixVariable(1);
			i_min->FixVariable(2);

		}


		if (fModelID == 2)
		{
			step[0] = 1.e-13;
			step[1] = 0.01;
			step[2] = 0.01;
			variable[0] = 1.e-12;
			variable[1] = parms[1];
			variable[2] = parms[2];

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Index",variable[1], step[1]);
			i_min->SetVariable(2,"Beta",variable[2], step[2]);

			i_min->FixVariable(1);
			i_min->FixVariable(2);
		}


		if (fModelID == 3)
		{
			step[0] = 1.e-14;
			step[1] = 0.01;
			step[2] = 0.001;
			variable[0] = 1.e-13;
			variable[1] = parms[1];
			variable[2] = parms[2];

			// Set the free variables to be minimized!
			i_min->SetVariable(0,"Norm",variable[0], step[0]);
			i_min->SetVariable(1,"Alpha",variable[1], step[1]);
			i_min->SetVariable(2,"Beta",variable[2], step[2]);

			i_min->FixVariable(1);
			i_min->FixVariable(2);

		}



		// Checking TS before applying fit
		cout << "Getting Total Counts for TS" << endl;
		float *iFluxPoint = new float[3];
		int i_count = 0;
		double i_onTotal = 0;
		double i_offTotal = 0;
		double i_alphaAvg = 0;
		bool i_bInc = true;
		for (unsigned int i = 0; i < fOnCounts.size(); i++)
		{
			i_bInc = true;

			for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
			{
				if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
				{
					i_bInc = false;
					break;
				}
			}

			if (i_bInc == false){continue;}

			for (unsigned int j = 0; j < fOnCounts[0].size(); j++ )
			{
				if (fRecBinCentres[j] < fFitMin){continue;}
				if (fRecBinCentres[j] > fFitMax) {continue;}
				i_onTotal += fOnCounts[i][j];
				i_offTotal += fOffCounts[i][j];
				i_alphaAvg += fRunList[i].alpha;
				i_count += 1;

			}
		}
		i_alphaAvg /= i_count;
		iFluxPoint[2] = getBinTS(i_onTotal, i_offTotal, i_alphaAvg);
		cout << "Test Statistic: " << i_onTotal << " " << i_offTotal << " " << i_alphaAvg << " " << iFluxPoint[2] << endl;


		if ( iFluxPoint[2] > 25. )
		{
					// // do the minimization
			i_min->Minimize();
			const double *i_FitValues = i_min->X();
			const double *i_Errors = i_min->Errors();

			// float i_predflux = iBestFit->Eval( TMath::Log10(iE_Norm) );

			iFluxPoint[0] = i_FitValues[0];
			iFluxPoint[1] = i_Errors[0];
			// i_min->Minimize();
			unsigned int i_steps = 10000;
			double i_norm[10000];
			double i_likelihood[10000];


			i_min->Scan(0, i_steps, i_norm, i_likelihood, iFluxPoint[0]- 2*iFluxPoint[1] ,iFluxPoint[0] + 2*iFluxPoint[1]);

			double max_likelihood = 999;

			for (unsigned int h = 0; h < i_steps; h++)
			{
				if (i_likelihood[h] < max_likelihood)
				{
					max_likelihood = i_likelihood[h];
				}
			}

			for ( unsigned int h  = 0; h < i_steps; h++)
			{
				i_likelihood[h] = max_likelihood - i_likelihood[h];
			}

			TGraph *i_limitSetter = new TGraph (i_steps, i_norm, i_likelihood);
			// TCanvas *cTest = new TCanvas();
			i_limitSetter->Draw("APL");
			// cTest->SetLogy();
		}


		else
		{
			// i_min->SetErrorDef(2.71);

			float i_predflux = iBestFit->Eval( TMath::Log10(iE_Norm) );
			cout << "Predicted Flux : " << i_predflux << endl;
			iFluxPoint[0] = i_predflux;
			iFluxPoint[1] = 0;
			i_min->Minimize();
			unsigned int i_steps = 10000;
			double i_norm[10000];
			double i_likelihood[10000];


			i_min->Scan(0, i_steps, i_norm, i_likelihood, i_predflux,1.e-7);

			double max_likelihood = 999;

			for (unsigned int h = 0; h < i_steps; h++)
			{
				if (i_likelihood[h] < max_likelihood)
				{
					max_likelihood = i_likelihood[h];
				}
			}

			// for ( int h  = 0; h < i_steps; h++)
			// {
			// 	i_likelihood[h] = max_likelihood - i_likelihood[h];
			// }

			TGraph *i_limitSetter = new TGraph (i_steps, i_norm, i_likelihood);
			TCanvas *cTest = new TCanvas();
			i_limitSetter->Draw("APL");
			cTest->SetLogy();

		}



		return iFluxPoint;
	}




// // Function to apply fit to only the energy bin of interest
// 	double* VLikelihoodFitter::getSpectralPointTotal( double *parms , double BinMin, double BinMax, double iE_Norm)
// 	{


// 		setFitMinMax(BinMin, BinMax);


// 		if (fModelID == 0)
// 		{
// 			sprintf( hname, "[0]*TMath::Power(TMath::Power(10.0,x) / %f, [1])", iE_Norm);
// 			fModel = new TF1("fModel", hname,-1.0,5.0);
// 			nParms = 2;
// 		}

// 		// Power Law with Exponential Cut off
// 		else  if (fModelID == 1)
// 		{
// 			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1] ) * TMath::Exp( -1. * TMath::Power( 10, x ) / [2] )", iE_Norm);
// 			fModel = new TF1("fModel", hname,-1.0,5.0);
// 			nParms = 3;
// 		}


// 		// Curved Spectrum
// 		else if (fModelID == 2)
// 		{
// 			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x )/ %f, [1]+[2]*TMath::Power( 10, x ) )", iE_Norm);
// 			fModel = new TF1("fModel", hname,-1.0,5.0);
// 			nParms = 3;

// 		}


// 		// Log-parabola
// 		else if (fModelID == 3)
// 		{
// 			sprintf( hname, "[0] * TMath::Power( TMath::Power( 10, x ) / %f, [1]+[2]*TMath::Log( TMath::Power( 10, x ) ) )", iE_Norm);
// 			fModel = new TF1("fModel", hname,-1.0,5.0);
// 			nParms = 3;

// 		}




// 		//nParms = 2;
// 		ROOT::Math::Minimizer* i_min =
// 		   ROOT::Math::Factory::CreateMinimizer("Minuit", "Minos");

// 		// // set tolerance , etc...
// 		i_min->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
// 		i_min->SetMaxIterations(100000);  // for GSL
// 		i_min->SetTolerance(0.001);
// 		i_min->SetPrintLevel(1);
// 		i_min->SetErrorDef(0.5);
// 		// // // create funciton wrapper for minmizer
// 		// // // a IMultiGenFunction type

// 		ROOT::Math::Functor i_fitfunction(this,&VLikelihoodFitter::getLogLiTotal,nParms);
// 		// ROOT::Math::Functor i_fitfunction(this,nParms, "getLikelihoodFit");
// 		double *step = new double[nParms];
// 		double *variable = new double[nParms];

// 		i_min->SetFunction(i_fitfunction);


// 		if (fModelID == 0)
// 		{
// 			step[0] = 1.e-17;
// 			step[1] = 0.01;
// 			variable[0] = 1.e-16;
// 			variable[1] = parms[1];
// 			// Set the free variables to be minimized!
// 			i_min->SetVariable(0,"Norm",variable[0], step[0]);
// 			i_min->SetVariable(1,"Index",variable[1], step[1]);

// 			i_min->FixVariable(1);

// 		}


// 		if (fModelID == 1)
// 		{
// 			step[0] = 1.e-13;
// 			step[1] = 0.01;
// 			step[2] = 0.01;
// 			variable[0] = 1.e-12;
// 			variable[1] = parms[1];
// 			variable[2] = parms[2];
// 			// Set the free variables to be minimized!
// 			i_min->SetVariable(0,"Norm",variable[0], step[0]);
// 			i_min->SetVariable(1,"Index",variable[1], step[1]);
// 			i_min->SetVariable(2,"E_CutOff", variable[2], step[2]);

// 			i_min->FixVariable(1);
// 			i_min->FixVariable(2);

// 		}


// 		if (fModelID == 2)
// 		{
// 			step[0] = 1.e-13;
// 			step[1] = 0.01;
// 			step[2] = 0.01;
// 			variable[0] = 1.e-12;
// 			variable[1] = parms[1];
// 			variable[2] = parms[2];

// 			// Set the free variables to be minimized!
// 			i_min->SetVariable(0,"Norm",variable[0], step[0]);
// 			i_min->SetVariable(1,"Index",variable[1], step[1]);
// 			i_min->SetVariable(2,"Beta",variable[2], step[2]);

// 			i_min->FixVariable(1);
// 			i_min->FixVariable(2);
// 		}


// 		if (fModelID == 3)
// 		{
// 			step[0] = 1.e-13;
// 			step[1] = 0.01;
// 			step[2] = 0.001;
// 			variable[0] = 1.e-12;
// 			variable[1] = parms[1];
// 			variable[2] = parms[2];

// 			// Set the free variables to be minimized!
// 			i_min->SetVariable(0,"Norm",variable[0], step[0]);
// 			i_min->SetVariable(1,"Alpha",variable[1], step[1]);
// 			i_min->SetVariable(2,"Beta",variable[2], step[2]);

// 			i_min->FixVariable(1);
// 			i_min->FixVariable(2);

// 		}



// 		// starting point



// 		// // do the minimization
// 		i_min->Minimize();
// 		const double *i_FitValues = i_min->X();
// 		const double *i_Errors = i_min->Errors();


// 		double *iFluxPoint = new double[2];
// 		iFluxPoint[0] = i_FitValues[0];
// 		iFluxPoint[1] = i_Errors[0];

// 		return iFluxPoint;
// 	}


	// Getting Energy Spectrum Points
	TGraphAsymmErrors* VLikelihoodFitter::getEnergySpectrum(TF1 *iBestFit)
	{


		vector <float> iFlux;
		vector <float> iFluxErr;
		vector <float> iEnergy;
		vector <float> iEnergyErr;
		vector <float> iBinTS;

		double iGlobalMin = fFitMin;
		double iGlobalMax = fFitMax;
		float *iFluxPoint;

		double iE_Norm;
		double *parms = iBestFit->GetParameters();
		cout << "\n\nGetting Energy Spectrum\n";
		cout << "Norm : " << parms[0] << endl;
		cout << "Gamma : " << parms[1] << endl;
		for (int i = 0 ; i < nMC; i++)
		{
			// cout << i << " " << fMCBins[i] << " - " << fMCBins[i+1] << " ; " <<  iGlobalMin << " - " << iGlobalMax << endl ;
			if (fMCBins[i] < iGlobalMin )
			{
				continue;
			}

			if (fMCBins[i] > iGlobalMax || fMCBins[i+1] > iGlobalMax )
			{
				continue;
			}

			cout << "Getting energy bin " << i << " : " << fMCBins[i]  << " - " <<  fMCBins[i+1] << endl;

			setFitMinMax(fMCBins[i] , fMCBins[i+1]);
			iE_Norm = VMathsandFunctions::getSpectralWeightedMeanEnergy(fMCBins[i], fMCBins[i+1], parms[1]);
			iE_Norm = TMath::Power(10.,iE_Norm);
			iFluxPoint = getSpectralPoint( parms , fMCBins[i], fMCBins[i+1], iE_Norm, iBestFit);
			cout << iE_Norm << " : " << iFluxPoint[0] << " +/- " << iFluxPoint[1] <<  " : " << iFluxPoint[2] <<  endl;

			iFlux.push_back(iFluxPoint[0]);
			iFluxErr.push_back(iFluxPoint[1]);
			iEnergy.push_back(TMath::Log10(iE_Norm));
			iEnergyErr.push_back(0);
			iBinTS.push_back(iFluxPoint[2]);

		}

		setFitMinMax(iGlobalMin, iGlobalMax);


		TGraphAsymmErrors * iEnergySpectrum = new TGraphAsymmErrors(iEnergy.size(), &(iEnergy[0]), &(iFlux[0]),
			&(iEnergyErr[0]), &(iEnergyErr[0]), &(iFluxErr[0]), &(iFluxErr[0]));

		cout << "E\tEMin\tEMax\tdN/dE\terr(dN/dE)\tTS\n";
		unsigned int i = 0;
		int j = 0;
		while (i < iEnergy.size())
		{
			if (fMCBins[j] < iGlobalMin )
			{
				j++;
				continue;
			}

			if (fMCBins[j] > iGlobalMax || fMCBins[j+1] > iGlobalMax )
			{
				j++;
				continue;
			}
			cout << TMath::Power(10.,iEnergy[i]) << "\t" << TMath::Power(10.,fMCBins[j]) << "\t"
			     << TMath::Power(10.,fMCBins[j+1]) << "\t" << iFlux[i] << "\t" << iFluxErr[i] << "\t" << iBinTS[i] << endl;
			i++;
			j++;
		}
		return iEnergySpectrum;
	}



	// // Fitting to the total energy spectrum
	// TGraphAsymmErrors* VLikelihoodFitter::getEnergySpectrumTotal(TF1 *iBestFit)
	// {


	// 	vector <double> iFlux;
	// 	vector <double> iFluxErr;
	// 	vector <double> iEnergy;
	// 	vector <double> iEnergyErr;

	// 	double iGlobalMin = fFitMin;
	// 	double iGlobalMax = fFitMax;
	// 	double *iFluxPoint;

	// 	double iE_Norm;
	// 	double *parms = iBestFit->GetParameters();
	// 	cout << "\n\nGetting Energy Spectrum\n";
	// 	cout << "Norm : " << parms[0] << endl;
	// 	cout << "Gamma : " << parms[1] << endl;
	// 	for (int i = 0 ; i < nMC; i++)
	// 	{
	// 		// cout << i << " " << fMCBins[i] << " - " << fMCBins[i+1] << " ; " <<  iGlobalMin << " - " << iGlobalMax << endl ;
	// 		if (fMCBins[i] < iGlobalMin )
	// 		{
	// 			continue;
	// 		}

	// 		if (fMCBins[i] > iGlobalMax || fMCBins[i+1] > iGlobalMax )
	// 		{
	// 			continue;
	// 		}

	// 		cout << "Getting energy bin " << i << " : " << fMCBins[i]  << " - " <<  fMCBins[i+1] << endl;

	// 		setFitMinMax(fMCBins[i] , fMCBins[i+1]);
	// 		iE_Norm = VMathsandFunctions::getSpectralWeightedMeanEnergy(fMCBins[i], fMCBins[i+1], parms[1]);
	// 		iE_Norm = TMath::Power(10.,iE_Norm);
	// 		iFluxPoint = getSpectralPointTotal( parms , fMCBins[i], fMCBins[i+1], iE_Norm);
	// 		cout << iE_Norm << " : " << iFluxPoint[0] << " +/- " << iFluxPoint[1] << endl;

	// 		iFlux.push_back(iFluxPoint[0]);
	// 		iFluxErr.push_back(iFluxPoint[1]);
	// 		iEnergy.push_back(iE_Norm);
	// 		iEnergyErr.push_back(0);

	// 	}

	// 	setFitMinMax(iGlobalMin, iGlobalMax);


	// 	TGraphAsymmErrors * iEnergySpectrum = new TGraphAsymmErrors(iEnergy.size(), &(iEnergy[0]), &(iFlux[0]), &(iEnergyErr[0]), &(iEnergyErr[0]), &(iFluxErr[0]), &(iFluxErr[0]));

	// 	cout << "E\tEMin\tEMax\tdN/dE\terr(dN/dE)\n";
	// 	int i = 0;
	// 	int j = 0;
	// 	while (i < iEnergy.size())
	// 	{
	// 		if (fMCBins[j] < iGlobalMin )
	// 		{
	// 			j++;
	// 			continue;
	// 		}

	// 		if (fMCBins[j] > iGlobalMax || fMCBins[j+1] > iGlobalMax )
	// 		{
	// 			j++;
	// 			continue;
	// 		}
	// 		cout << iEnergy[i] << "\t" << TMath::Power(10.,fMCBins[j]) << "\t" << TMath::Power(10.,fMCBins[j+1]) << "\t" << iFlux[i] << "\t" << iFluxErr[i] << endl;
	// 		i++;
	// 		j++;
	// 	}
	// 	return iEnergySpectrum;
	// }



	// Getting plots with total counts and residuals
	TCanvas* VLikelihoodFitter::getTotalCountsPlots()
	{

		double *i_Parms = new double[nParms];

		for (int i = 0; i < nParms; i ++)
		{
			i_Parms[i] = fModel->GetParameter(i);

		}

		vector < vector <double> > i_vOffMLE = getModelPredictedOff(i_Parms);
		vector < vector <double> > i_vModel = getModelPredictedExcess(i_Parms);

		vector < double > i_OnTotal( i_vOffMLE[0].size());
		vector < double > i_OffTotal( i_vOffMLE[0].size());
		vector < double > i_ExcessTotal( i_vOffMLE[0].size());
		vector < double > i_vOffMLETotal( i_vOffMLE[0].size());
		vector < double > i_vModelTotal( i_vOffMLE[0].size());
		vector < double > i_OffError( i_vOffMLE[0].size());
		vector < double > i_ModelOn( i_vOffMLE[0].size());

		vector < double > i_OnRes( i_vOffMLE[0].size());
		vector < double > i_OffRes( i_vOffMLE[0].size());
		vector < double > i_ExcessRes( i_vOffMLE[0].size());

		vector < double > i_OnResErr( i_vOffMLE[0].size());
		vector < double > i_OffResErr( i_vOffMLE[0].size());

		double i_max_On = 0;
		double i_max_Off = 0;
		bool i_bInc = true;


		for (unsigned int i = 0; i < fRunList.size(); i++ )
		{

			if (fRunList[i].MJD < fMJD_Min ||  fRunList[i].MJD > fMJD_Max )
			{
				continue;
			}

			i_bInc = true;

			for (unsigned int e = 0; e < fExcludeMJD.size(); e++)
			{
				if ( fRunList[i].MJD > fExcludeMJD[e][0] &&  fRunList[i].MJD < fExcludeMJD[e][1] )
				{
					i_bInc = false;
					break;
				}
			}

			if (i_bInc == false){continue;}

			for (unsigned int j = 0; j < i_vOffMLE[0].size(); j++)
			{

				if ( fOnCounts[i][j] > 1.e9 || fOffCounts[i][j] > 1.e9 || i_vOffMLE[i][j] > 1.e9 || i_vModel[i][j] > 1.e9 )
				{
					// cout << i << " " << j << " " << fOnCounts[i][j] << " " << fOffCounts[i][j] << " " << i_vOffMLE[i][j] << " " << i_vModel[i][j] << endl;
				}

				if (fEnergyBias[i][j] > fThresholdBias){continue;}

				i_OnTotal[j] += fOnCounts[i][j];
				i_OffTotal[j] += fOffCounts[i][j];
				i_ExcessTotal[j] += fOnCounts[i][j] -  fRunList[i].alpha * fOnCounts[i][j];
				i_vOffMLETotal[j] += i_vOffMLE[i][j];
				i_vModelTotal[j] += i_vModel[i][j];
				i_ModelOn[j] +=  i_vModel[i][j] + fRunList[i].alpha * i_vOffMLE[i][j];

			}
		}

		for ( unsigned int i = 0 ; i < i_vOffMLE[0].size(); i++)
		{
			i_OnRes [i] = (i_OnTotal[i] - i_ModelOn[i]) / i_ModelOn[i];

			i_OffRes [i] = 	( i_OffTotal[i] - i_vOffMLETotal[i] ) / i_vOffMLETotal[i];
			i_ExcessRes [i] = 	( i_ExcessTotal[i] - i_vModelTotal[i] ) / i_vModelTotal[i];

			if (i_OnTotal[i] < 1.e-5 )
			{
				i_OnResErr [i] = 0;
				i_OnRes [i] = 0;
			}
			else
			{
				i_OnResErr [i] = TMath::Sqrt(i_OnTotal[i]) / i_ModelOn[i] ;
				// i_OnResErr [i] = 2;


			}


			if (i_OffTotal[i] < 1.e-5 )
			{
				i_OffResErr [i] = 0;
				i_OffRes [i] = 0;
			}
			else
			{
				i_OffResErr [i] = TMath::Sqrt(i_OffTotal[i]) / i_vOffMLETotal[i] ;

			}

			// cout << i_OnTotal[i] << " " << i_ModelOn[i] << " " << i_OnRes [i] << " " << i_OnResErr[i] << endl;


		}

//		i_max_On = *max_element(std::begin(i_OnTotal), std::end(i_OnTotal));
//		i_max_Off = *max_element(std::begin(i_OffTotal), std::end(i_OffTotal));
		i_max_On = 1.e5;
		i_max_Off = 1.e5;

		TCanvas *i_cTemp = new TCanvas();
		i_cTemp->Divide(2,2);




		TF1 *fZero = new TF1("fZero", "0", -10,10);
		fZero->SetLineStyle(7);
		fZero->SetLineColor(kBlack);



		// On
		TPad *p = (TPad*)i_cTemp->cd(1);
		TGraphAsymmErrors *i_gOnCounts = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_OnTotal[0]) );
		TGraphAsymmErrors *i_gOnModel = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_ModelOn[0]) );


		TLine *line1_On = new TLine(fFitMin,0.1,fFitMin,10*i_max_On);
		TLine *line2_On = new TLine(fFitMax,0.1,fFitMax,10*i_max_On);
		TLine *line1_Off = new TLine(fFitMin,0.1,fFitMin,10*i_max_Off);
		TLine *line2_Off = new TLine(fFitMax,0.1,fFitMax,10*i_max_Off);
		line1_On->SetLineColor(kRed);
		line2_On->SetLineColor(kRed);
		line1_Off->SetLineColor(kRed);
		line2_Off->SetLineColor(kRed);

		i_gOnCounts->SetTitle("On Counts");
		i_gOnCounts->GetXaxis()->SetTitle("log(Energy) [TeV]");
		i_gOnCounts->GetYaxis()->SetTitle("Counts");
		i_gOnCounts->GetYaxis()->SetRangeUser(0.1,10*i_max_On);

		i_gOnCounts->SetMarkerStyle(8);
		i_gOnCounts->SetMarkerSize(1.0);
		i_gOnCounts->SetMarkerColor(kBlue);
		i_gOnCounts->Draw("AP");


		i_gOnModel->SetMarkerStyle(8);
		i_gOnModel->SetMarkerSize(1.0);
		i_gOnModel->SetMarkerColor(kGreen);
		i_gOnModel->Draw("SAME,P");


		line1_On->Draw("SAME");
		line2_On->Draw("SAME");
		p->SetLogy();



		// Off
		p = (TPad*)i_cTemp->cd(2);
		TGraphAsymmErrors *i_gOffCounts = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_OffTotal[0]) );
		TGraphAsymmErrors *i_gPredOffCounts = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_vOffMLETotal[0]) );

		i_gOffCounts->SetTitle("Off Counts");
		i_gOffCounts->GetXaxis()->SetTitle("log(Energy) [TeV]");
		i_gOffCounts->GetYaxis()->SetTitle("Counts");
		i_gOffCounts->GetYaxis()->SetRangeUser(0.1,10*i_max_Off);

		i_gOffCounts->SetMarkerStyle(8);
		i_gOffCounts->SetMarkerSize(1.0);
		i_gOffCounts->SetMarkerColor(kBlue);
		i_gOffCounts->Draw("AP");

		i_gPredOffCounts->SetMarkerStyle(8);
		i_gPredOffCounts->SetMarkerSize(1.0);
		i_gPredOffCounts->SetMarkerColor(kGreen);
		i_gPredOffCounts->Draw("SAME,P");

		line1_Off->Draw("SAME");
		line2_Off->Draw("SAME");
		p->SetLogy();



		TLine *line3= new TLine(fFitMin,-2,fFitMin,2);
		TLine *line4 = new TLine(fFitMax,-2,fFitMax,2);
		line3->SetLineColor(kRed);
		line4->SetLineColor(kRed);
		// On Residuals
		p = (TPad*)i_cTemp->cd(3);
		TGraphErrors *i_gOnRes = new TGraphErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_OnRes[0]), 0, &(i_OnResErr[0]));

		i_gOnRes->SetTitle("On Count Residuals");
		i_gOnRes->GetXaxis()->SetTitle("log(Energy) [TeV]");
		i_gOnRes->GetYaxis()->SetTitle("(Counts - Model) / Model");
		i_gOnRes->GetYaxis()->SetRangeUser(-2,2);

		i_gOnRes->SetMarkerStyle(8);
		i_gOnRes->SetMarkerSize(1.0);
		i_gOnRes->SetMarkerColor(kBlack);
		i_gOnRes->Draw("AP");
		line3->Draw("SAME");
		line4->Draw("SAME");
		fZero->Draw("SAME");



		// Off Residuals
		p = (TPad*)i_cTemp->cd(4);
		TGraphErrors *i_gOffRes = new TGraphErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_OffRes[0]), 0, &(i_OffResErr[0]));

		i_gOffRes->SetTitle("Off Count Residuals");
		i_gOffRes->GetXaxis()->SetTitle("log(Energy) [TeV]");
		i_gOffRes->GetYaxis()->SetTitle("(Counts - Model) / Model");
		i_gOffRes->GetYaxis()->SetRangeUser(-2,2);

		i_gOffRes->SetMarkerStyle(8);
		i_gOffRes->SetMarkerSize(1.0);
		i_gOffRes->SetMarkerColor(kBlack);
		i_gOffRes->Draw("AP");

		line3->Draw("SAME");
		line4->Draw("SAME");
		fZero->Draw("SAME");


		return i_cTemp;

	}



// Getting useful plots for runs
TCanvas* VLikelihoodFitter::getRunPlots(int i_Entry)
	{

		double *i_Parms = new double[nParms];
		for (int i = 0; i < nParms; i ++)
		{
			i_Parms[i] = fModel->GetParameter(i);

		}

		vector < vector <double> > i_vOffMLE = getModelPredictedOff(i_Parms);
		vector < vector <double> > i_vModel = getModelPredictedExcess(i_Parms);

		TCanvas *i_cTemp = new TCanvas();
		i_cTemp->Divide(3,1);

		TLine *line1 = new TLine(fFitMin,0.1,fFitMin,5000);
		TLine *line2 = new TLine(fFitMax,0.1,fFitMax,5000);
		line1->SetLineColor(kRed);
		line2->SetLineColor(kRed);

		// On
		TPad *p = (TPad*)i_cTemp->cd(1);
		TGraphAsymmErrors *i_gOnCounts = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(fOnCounts[i_Entry][0]) );
		// TGraphAsymmErrors *i_gOn

		i_gOnCounts->SetTitle("On Counts");
		i_gOnCounts->GetXaxis()->SetTitle("log(Energy) [TeV]");
		i_gOnCounts->GetYaxis()->SetTitle("Counts");
		i_gOnCounts->GetYaxis()->SetRangeUser(0.1,1.e5);

		i_gOnCounts->SetMarkerStyle(8);
		i_gOnCounts->SetMarkerSize(1.0);
		i_gOnCounts->SetMarkerColor(kRed);
		i_gOnCounts->Draw("AP");

		line1->Draw("SAME");
		line2->Draw("SAME");

		p->SetLogy();


		// Off
		p = (TPad*)i_cTemp->cd(2);
		TGraphAsymmErrors *i_gOffCounts = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(fOffCounts[i_Entry][0]) );
		TGraphAsymmErrors *i_gPredOffCounts = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_vOffMLE[i_Entry][0]) );

		i_gOffCounts->SetTitle("Off Counts");
		i_gOffCounts->GetXaxis()->SetTitle("log(Energy) [TeV]");
		i_gOffCounts->GetYaxis()->SetTitle("Counts");
		i_gOffCounts->SetMarkerStyle(8);
		i_gOffCounts->SetMarkerSize(1.0);
		i_gOffCounts->SetMarkerColor(kBlue);
		i_gOffCounts->GetYaxis()->SetRangeUser(0.1,1.e5);

		i_gOffCounts->Draw("AP");

		i_gPredOffCounts->SetMarkerStyle(8);
		i_gPredOffCounts->SetMarkerSize(1.0);
		i_gPredOffCounts->SetMarkerColor(kGreen);
		i_gPredOffCounts->Draw("SAME,P");

		line1->Draw("SAME");
		line2->Draw("SAME");
		p->SetLogy();

		// Excess
		p = (TPad*)i_cTemp->cd(3);
		vector < double > i_vExcess;
		for (unsigned int i = 0; i < fOffCounts[i_Entry].size(); i++)
		{
			i_vExcess.push_back(fOnCounts[i_Entry][i] - fRunList[i_Entry].alpha * fOffCounts[i_Entry][i]);
		}
		TGraphAsymmErrors *i_gExcess = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_vExcess[0]) );
		TGraphAsymmErrors *i_gPredExcess = new TGraphAsymmErrors(fRecBinCentres.size(), &(fRecBinCentres[0]), &(i_vModel[i_Entry][0]) );

		i_gExcess->SetTitle("Excess Counts");
		i_gExcess->GetXaxis()->SetTitle("log(Energy) [TeV]");
		i_gExcess->GetYaxis()->SetTitle("Counts");
		i_gExcess->GetYaxis()->SetRangeUser(0.1,1.e5);

		i_gExcess->SetMarkerStyle(8);
		i_gExcess->SetMarkerSize(1.0);
		i_gExcess->SetMarkerColor(kBlue);
		i_gExcess->Draw("AP");

		i_gPredExcess->SetMarkerStyle(8);
		i_gPredExcess->SetMarkerSize(1.0);
		i_gPredExcess->SetMarkerColor(kGreen);
		i_gPredExcess->Draw("SAME,P");

		line1->Draw("SAME");
		line2->Draw("SAME");

		p->SetLogy();


		return i_cTemp;

	}



// Dates to exclude from Spectral Fit
void VLikelihoodFitter::addExclusionDate(double i_MJDStart, double i_MJDStop)
{
	vector <double> tmp;
	tmp.push_back(i_MJDStart);
	tmp.push_back(i_MJDStop);

	fExcludeMJD.push_back(tmp);
}

float* VLikelihoodFitter::getIntergralFlux(double i_EMin, double i_EMax, TF1* i_Model, bool i_log)
{


	// Getting the ingtegral flux from the best fit model

	if (i_log)
	{
		i_EMin = TMath::Power(10,i_EMin);
		i_EMax = TMath::Power(10,i_EMax);
	}

	float *i_flux = new float[4];

	// Calculating flux from best fit model;
	i_flux[0] = i_Model->GetParameter(0) * (TMath::Power(i_EMax, i_Model->GetParameter(1) +1 ) - TMath::Power(i_EMin, i_Model->GetParameter(1) +1 ) ) /( i_Model->GetParameter(1) +1) / (TMath::Power(E_Norm,i_Model->GetParameter(1))) ;

	// Calculating the Error
	float a,b,c,d;
	// dN_int/dNo
	a = (TMath::Power(i_EMax, i_Model->GetParameter(1) +1 ) - TMath::Power(i_EMin, i_Model->GetParameter(1) +1 ) )/( i_Model->GetParameter(1) +1) / (TMath::Power(E_Norm,i_Model->GetParameter(1))) ;

	// dN_int/dGamma
	b = i_Model->GetParameter(0) * (TMath::Power(i_EMax, i_Model->GetParameter(1) +1 ) - TMath::Power(i_EMin, i_Model->GetParameter(1) +1 ) ) /( i_Model->GetParameter(1) +1) / ( i_Model->GetParameter(1) +1)  / (TMath::Power(E_Norm,i_Model->GetParameter(1))) ;

	c = i_Model->GetParameter(0) * (TMath::Power(i_EMax, i_Model->GetParameter(1) +1 ) * TMath::Log(i_EMax) - TMath::Power(i_EMin, i_Model->GetParameter(1) +1 ) * TMath::Log(i_EMin) ) /( i_Model->GetParameter(1) +1)  / (TMath::Power(E_Norm,i_Model->GetParameter(1))) ;

	d =  - i_Model->GetParameter(0) * (TMath::Power(E_Norm,-1*i_Model->GetParameter(1))) * TMath::Log(E_Norm) * (TMath::Power(i_EMax, i_Model->GetParameter(1) +1 ) - TMath::Power(i_EMin, i_Model->GetParameter(1) +1 ) ) /( i_Model->GetParameter(1) +1) ;



	// Calculating Error
	i_flux[1] = sqrt( a* a * i_Model->GetParError(0) * i_Model->GetParError(0) + (b + c + d) * (b + c + d) *  i_Model->GetParError(1) * i_Model->GetParError(1));



	i_flux[2] = getCrabFlux(i_flux[0],i_EMin, i_EMax);
	i_flux[3] = getCrabFlux(i_flux[1],i_EMin, i_EMax);

	return i_flux;
}


double VLikelihoodFitter::getCrabFlux( double iF, double i_EMin, double i_EMax, double i_Gamma)
{
	double i_N0 = 3.20e-11;
	double i_Crab = i_N0 * (TMath::Power(i_EMax, i_Gamma +1 ) - TMath::Power(i_EMin, i_Gamma +1 ) ) /( i_Gamma +1);

	return (iF/i_Crab);

}


double VLikelihoodFitter::getBinTS(double i_on, double i_off, double i_alpha)
{
	double i_a = 0;
	double i_b = 0;

	i_a = i_on * TMath::Log( (1. + i_alpha) *i_on / ( i_alpha * (i_on + i_off) ) );
	i_b = i_off * TMath::Log( (1. + i_alpha) *i_off / ( (i_on + i_off) ) );

	return 2.*( i_a + i_b );
}
