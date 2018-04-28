#ifndef PAIR_PURITY_FITTER_H
#define PAIR_PURITY_FITTER_H

#include "HistoAnalyzer.h"
#include "XmlFunction.h"
#include "XmlHistogram.h"
#include "CutCollection.h"
#include "FitConfidence.h"

#include "RooPlotLib.h"

#include "TRandom3.h"

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "vendor/loguru.h"

// globals for fit
TH1 * hPDFMuMu     = nullptr;
TH1 * hPDFPiPi     = nullptr;
TH1 * hPDFPiMu     = nullptr;

Double_t fitfun3( double *x, double *par ){

	double x0=x[0];
	
	double scale_sigsig  = par[0];
	double scale_bgbg    = par[1];
	double scale_bgsig   = par[2];
	
	int sigsig_bin   = hPDFMuMu->GetXaxis()->FindBin( x0 );
	int bgbg_bin     = hPDFPiPi->GetXaxis()->FindBin(x0);
	int bgsig_bin    = hPDFPiMu->GetXaxis()->FindBin(x0);

	float sigsig_val = hPDFMuMu->GetBinContent( sigsig_bin );
	float bgbg_val   = hPDFPiPi->GetBinContent( bgbg_bin );
	float bgsig_val  = hPDFPiMu->GetBinContent( bgsig_bin );

	return sigsig_val * scale_sigsig + bgbg_val * scale_bgbg + bgsig_val * scale_bgsig;
}

class PairPurityFitter : public HistoAnalyzer {
protected:

	string rpName;

	map<string, HistoBins> bins;


	bool export_img;
	TCanvas *can;
	TCanvas *can2;

	TH2 * nn_sig_pos, *nn_sig_neg;
	TH2 * nn_pi_pos, *nn_pi_neg;

	map<string, vector<TH1 *>> nn_pt;
	map<string, vector<TH2 *>> deltas;

	size_t nSamples;
	bool weight_by_pt  = false;
	bool weight_by_pid = false;

	bool use_kaon = false;

	string bg_source;
	string sig_source;

	TRandom3 r3;
	float purityCut = 0.0;

public:

	virtual void initialize(){
		HistoAnalyzer::initialize();

		r3.SetSeed( config.get<int>( "seed", 0 ) );

		export_img = config.getBool( "can:export", false );
		weight_by_pt = config.getBool( "weight:pt", false );
		weight_by_pid = config.getBool( "weight:pid", false );
		use_kaon = config.getBool( "fit:kaon", false );
		nSamples = config.getInt( "nSamples", 1000 );

		purityCut = config.get<float>( "purityCut", 1.2 );

		vector<string> paths = config.childrenOf( "bins" );
		for ( auto p : paths ){
			string tn = config.tagName(p);
			LOG_F(  INFO, "Loading Bins: %s", tn.c_str() );
			bins[ tn ].load( config, p );
			LOG_F( INFO, "%s", bins[tn].toString().c_str() );
		}

		// build pt projections

		book->cd();
		sig_source = "sig";
		build_pt_projections( this->sig_source );
		bg_source = "bg";
		build_pt_projections( this->bg_source );
		if ( use_kaon )
			build_pt_projections( "kaon" );

		if ( weight_by_pt ){
			build_mass_projections( "uls_delta_pt", "deltaPt", "deltaPt" );
			build_mass_projections( "ls_delta_pt", "deltaPt", "deltaPt" );
		}
		if ( weight_by_pid ){
			build_mass_projections( "uls_delta_pid", "minpid", "deltaPid" );
			build_mass_projections( "ls_delta_pid", "minpid", "deltaPid" );
		}


		// book->makeAll( config, "" )
	}

	void build_pt_projections( string hname ){
		LOG_SCOPE_FUNCTION( INFO );
		TH2 * h2rawpos = get<TH2>( hname + "_pos_mlp", "mc" );
		TH2 * h2rbpos = HistoBins::rebin2D( hname + "posrb", h2rawpos, bins["ptTemplate"], bins["pid"] );

		TH2 * h2rawneg = get<TH2>( hname + "_neg_mlp", "mc" );
		TH2 * h2rbneg = HistoBins::rebin2D( hname + "negrb", h2rawpos, bins["ptTemplate"], bins["pid"] );

		vector<TH1 *> projections;
		LOG_F( 9, "nBins = %d", h2rbpos->GetXaxis()->GetNbins() );
		for ( size_t i = 0; i <= h2rbpos->GetXaxis()->GetNbins(); i++ ){
			TH1 * h = h2rbpos->ProjectionY( TString::Format( "%s_pt_%lu", hname.c_str(), i ), i+1, i+1 );
			TH1 * hn = h2rbneg->ProjectionY( TString::Format( "%s_neg_pt_%lu", hname.c_str(), i ), i+1, i+1 );
			h->Add( hn );

			LOG_F( 9, "Integral pt[%lu] = %f", i, h->Integral() );
			projections.push_back( h );
		}
		nn_pt[ hname ] = projections;
	}

	void build_mass_projections( string hname, string binsx, string binsy ){
		LOG_SCOPE_FUNCTION( INFO );
		TH3 * h3raw = get<TH3>( hname, "pairPid" );
		TH3 * h3rb = HistoBins::rebin3D( hname + "rb", h3raw, bins[binsx], bins[binsy], bins["mass"] );

		vector<TH2 *> projections;
		for ( size_t i = 0; i <= h3rb->GetZaxis()->GetNbins(); i++ ){
			h3rb->GetZaxis()->SetRange( i, i );

			TH2 * h2 = (TH2*)h3rb->Project3D( "yx" );
			h2 = (TH2*)h2->Clone( TString::Format( "%s_mass_%lu", hname.c_str(), i ) );
			projections.push_back( h2 );
		}

		deltas[ hname ] = projections;
	}

	virtual void loop_on_mass(TH2 * h2pid, TH2 * h2pt, string prefix){
		assert( h2pid != nullptr );
		assert( h2pt != nullptr );

		for ( size_t i = 1; i <= h2pid->GetXaxis()->GetNbins(); i++ ){
			// if ( i < config.get<int>( "fit:mMin", 1 ) || i > config.getInt( "fit:mMax", 1 ) ) continue;

			string projName = prefix + "_pid_mass" + ts( (int) i);
			string projNamePt = prefix + "_pt_mass" + ts( (int) i);
			TH1 * hpid = h2pid->ProjectionY( projName.c_str(), i, i );
			TH1 * hpt  = h2pt->ProjectionY( projNamePt.c_str(), i, i );

			float m1 = h2pid->GetXaxis()->GetBinLowEdge( i );
			float m2 = h2pid->GetXaxis()->GetBinUpEdge( i );

			hpid->SetTitle( TString::Format( "%0.2f < M < %0.2f", m1, m2 ) );

			if ( hpid->Integral() <= 0 ) continue;

			TH1 * hpidraw = (TH1*)hpid->Clone( TString::Format( "%s_pid_m%lu", prefix.c_str(), i ) );

			hpid->Sumw2();
			hpid->Scale( 1.0 / hpid->Integral() );
			hpid->SetMinimum( 5e-4 );
			hpid->Draw( "" );
			gPad->SetLogy(1);
			// can->Print( rpName.c_str() );

			hpt->SetTitle( TString::Format( "%0.2f < M < %0.2f", m1, m2 ) );
			hpt->GetXaxis()->SetRangeUser( 0, 5.0 );
			hpt->Scale( 1.0 / hpt->Integral() );
			hpt->SetMinimum( 1e-4 );
			hpt->Draw( "" );
			gPad->SetLogy(1);
			// can->Print( rpName.c_str() );

			fit_pair_pid( hpid, hpt, prefix, i, m1, m2 );

			// break;

		} //loop on mass bins
	} // loop_on_mass


	/* Generates template pairPid shapes from MC
	 * Samples the input pT distribution and builds a weighted template for the given pair kinematics
	 *
	 */
	void generate_templates( TH1 * hpt, TH1 * hpipi, TH1 * hpimu, TH1 * hmumu, TH2 * hdeltaPid, TH2 * hdeltaPt ){
		LOG_SCOPE_FUNCTION( INFO );
		assert( hpt != nullptr );
		assert( nullptr != hpipi );
		assert( nullptr != hpimu );
		assert( nullptr != hmumu );
		assert( nullptr != hdeltaPid );
		assert( nullptr != hdeltaPt );

		if ( hpt->Integral() <= 0 ) {
			LOG_F( INFO, "pt is empty" );
			return;
		}
		assert( nn_pt.count( "bg" ) > 0 );
		if ( use_kaon ){
			assert( nn_pt.count( "kaon" ) > 0 );
		}
		assert( nn_pt.count( this->sig_source ) > 0 );

		vector<TH1 * > pi_pt  = nn_pt[ this->bg_source ];
		vector<TH1 * > k_pt  = nn_pt[ "kaon" ];
		vector<TH1 * > sig_pt = nn_pt[ this->sig_source ];

		for ( size_t i = 0; i < nSamples; i++ ){

			
			vector<TH1 *> *bg_pt = &pi_pt;


			float pt1 = hpt->GetRandom();
			float pt2 = hpt->GetRandom();

			if ( pt1 > 3 ) pt1 = 3;
			if ( pt2 > 3 ) pt2 = 3;
			
			int ipt1 = bins[ "ptTemplate" ].findBin( pt1 );
			int ipt2 = bins[ "ptTemplate" ].findBin( pt2 );
			if ( ipt1 < 0 || ipt1 >= sig_pt.size() ) continue;
			if ( ipt2 < 0 || ipt2 >= sig_pt.size() ) continue;

			float IBg1  = (*bg_pt)[ipt1]->Integral();
			float IBg2  = (*bg_pt)[ipt2]->Integral();
			float ISig1 = sig_pt[ipt1]->Integral();
			float ISig2 = sig_pt[ipt2]->Integral();

			float rBG1  = -999;
			float rBG2  = -999;
			float rSig1 = -999;
			float rSig2 = -999;

			if ( IBg1 > 0 ){
				
				if ( r3.Uniform( 1.0 ) < 0.1 && k_pt[ipt1]->Integral() > 0 )
					bg_pt = &k_pt;
				else 
					bg_pt = &pi_pt;
				rBG1  = (*bg_pt)[ipt1]->GetRandom();
			}
			if ( IBg2 > 0 ){
				if ( r3.Uniform( 1.0 ) < 0.1 && k_pt[ipt2]->Integral() > 0 )
					bg_pt = &k_pt;
				else 
					bg_pt = &pi_pt;
				rBG2  = (*bg_pt)[ipt2]->GetRandom();
			}
			
			if ( ISig1 > 0 )
				rSig1 = sig_pt[ipt1]->GetRandom() - fabs(r3.Gaus( 0, 0.00 ) );
			if ( ISig2 > 0 )
				rSig2 = sig_pt[ipt2]->GetRandom() - fabs(r3.Gaus( 0, 0.00 ) );

			float pairPid_bgbg   = sqrt( pow( rBG1, 2 ) + pow( rBG2, 2) );
			float pairPid_bgsig  = sqrt( pow( rBG1, 2 ) + pow( rSig2, 2) );
			float pairPid_sigbg  = sqrt( pow( rSig1, 2 ) + pow( rBG2, 2) );
			float pairPid_sigsig = sqrt( pow( rSig1, 2 ) + pow( rSig2, 2) );

			if ( IBg1 > 0 && IBg2 > 0 )
				hpipi->Fill( pairPid_bgbg );
			if ( ISig1 > 0 && ISig2 > 0 )
				hmumu->Fill( pairPid_sigsig );
			if ( IBg1 > 0 && ISig2 > 0 )
				hpimu->Fill( pairPid_bgsig );
			if ( IBg2 > 0 && ISig1 > 0 )
				hpimu->Fill( pairPid_sigbg );
		}
	} // generate_templates

	void fit_pair_pid( TH1 * hpid, TH1 * hpt, string prefix, size_t im, float m1, float m2 ){
		LOG_SCOPE_FUNCTION( INFO );
		LOG_F( INFO, "Fitting %s[%lu], %0.3f < M < %0.3f", prefix.c_str(), im, m1, m2 );
		assert( hpid != nullptr );

		TH1 * hpipi = nullptr, *hpimu = nullptr, *hmumu = nullptr;
		TH2 * hdeltaPid = nullptr, *hdeltaPt = nullptr;
		hpipi   = new TH1F( TString::Format( "template_%s_pipi_m%lu", prefix.c_str(), im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		// hkk   = new TH1F( TString::Format( "template_%s_kk_m%lu", prefix.c_str(), im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		hpimu  = new TH1F( TString::Format( "template_%s_pimu_m%lu", prefix.c_str(), im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		hmumu = new TH1F( TString::Format( "template_%s_mumu_m%lu", prefix.c_str(), im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		hdeltaPid  = new TH2F( TString::Format( "template_%s_deltaPid_m%lu", prefix.c_str(), im ), "", bins["pid"].nBins(), bins["pid"].getBins().data(), bins["deltaPid"].nBins(), bins["deltaPid"].getBins().data() );
		hdeltaPt   = new TH2F( TString::Format( "template_%s_deltaPt_m%lu", prefix.c_str(), im ), "", bins["pt"].nBins(), bins["pt"].getBins().data(), bins["deltaPt"].nBins(), bins["deltaPt"].getBins().data() );
		RooPlotLib rpl;
		rpl.style( hpipi ).set( config, "style.bgbg" );
		rpl.style( hpimu ).set( config, "style.bgsig" );
		rpl.style( hmumu ).set( config, "style.sigsig" );

		// first pass builds the correlations
		generate_templates( hpt, hpipi, hpimu, hmumu, hdeltaPid, hdeltaPt );

		/* First Pass Fitting
		 * use the vanilla mixed templates without any correlation information
		 */
		if ( hmumu->Integral() <= 0 ) return;
		if ( hpimu->Integral() <= 0 ) return;
		if ( hpipi->Integral() <= 0 ) return;

		hmumu->Scale( 1.0 / hmumu->Integral() );
		hpipi->Scale( 1.0 / hpipi->Integral() );
		hpimu->Scale( 1.0 / hpimu->Integral() );

		hPDFMuMu = hmumu;
		hPDFPiMu  = hpimu;
		hPDFPiPi   = hpipi;

		TF1 * ff = nullptr;
		ff = new TF1( "ff", fitfun3, -1, 2, 3 );
		ff->SetParNames( "sigsig", "bgbg", "bgsig" );
		ff->SetParLimits( 0, 1e-4, 1 );
		ff->SetParLimits( 1, 1e-4, 1 );
		ff->SetParLimits( 2, 1e-4, 1 );
		ff->SetParameters( 0.1, 0.1, 0.01 );

		ff->SetNpx( 1000 );
		ff->SetLineColor(kBlack);
		ff->SetLineWidth(2);

		string fitOpt = config[ "fit:opt" ] + "S";
		LOG_F( INFO, "Fitting (opt=%s) in range (%f, %f)", fitOpt.c_str(), config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
		TFitResultPtr fitPointer = hpid->Fit( "ff", fitOpt.c_str(), "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );

		report_fit( hpid, ff, prefix, im, m1, m2 );


		//  report raw purity
		float total = ff->GetParameter( "sigsig" ) + ff->GetParameter( "bgsig" ) + ff->GetParameter( "bgbg" );
		
		book->get( prefix + "_rawPurityMuMu" )->SetBinContent( im, ff->GetParameter( "sigsig" ) / total );
		book->get( prefix + "_rawPurityMuMu" )->SetBinError( im, ff->GetParError( ff->GetParNumber( "sigsig" ) ) );
		book->get( prefix + "_rawPurityPiMu" )->SetBinContent( im, ff->GetParameter( "bgsig" ) / total );
		book->get( prefix + "_rawPurityPiMu" )->SetBinError( im, ff->GetParError( ff->GetParNumber( "bgsig" ) ) );
		book->get( prefix + "_rawPurityPiPi" )->SetBinContent( im, ff->GetParameter( "bgbg" ) / total );
		book->get( prefix + "_rawPurityPiPi" )->SetBinError( im, ff->GetParError( ff->GetParNumber( "bgbg" ) ) );

		// report cut purity
		int iBin = hPDFMuMu->GetXaxis()->FindBin( purityCut );
		double eMuMu = 0, ePiMu = 0, ePiPi = 0;
		total = hPDFMuMu->IntegralAndError( iBin, -1, eMuMu ) + hPDFPiMu->IntegralAndError( iBin, -1, ePiMu ) + hPDFPiPi->IntegralAndError( iBin, -1, ePiPi );

		book->get( prefix + "_purityMuMu" )->SetBinContent( im, hPDFMuMu->Integral( iBin, -1 ) / total );
		// book->get( prefix + "_purityMuMu" )->SetBinError( im,  sqrt(hPDFMuMu->Integral( iBin, -1 )) );
		book->get( prefix + "_purityPiMu" )->SetBinContent( im, hPDFPiMu->Integral( iBin, -1 ) / total );
		// book->get( prefix + "_purityPiMu" )->SetBinError( im, sqrt(hPDFPiMu->Integral( iBin, -1 )) );
		book->get( prefix + "_purityPiPi" )->SetBinContent( im, hPDFPiPi->Integral( iBin, -1 ) / total );
		// book->get( prefix + "_purityPiPi" )->SetBinError( im, sqrt(hPDFPiPi->Integral( iBin, -1 ) ) );

		book->get( prefix + "_chi2ndf" )->SetBinContent( im, ff->GetChisquare() / (float)ff->GetNDF() );

		delete ff;
	} // fit pair_pid

	void report_fit( TH1 * hpid, TF1 * ff, string prefix, size_t im, float m1, float m2 ) {
		LOG_SCOPE_FUNCTION( INFO );
		assert( nullptr != hpid );
		assert( nullptr != ff );
		RooPlotLib rpl;
		rpl.style( hpid ).set( config, "style.data" ).draw();
		gPad->SetLogy(1);

		hPDFMuMu->Scale( ff->GetParameter( "sigsig" ) );
		hPDFMuMu->Draw("same");

		hPDFPiMu->Scale( ff->GetParameter( "bgsig" ) );
		hPDFPiMu->Draw("same");

		hPDFPiPi->Scale( ff->GetParameter( "bgbg" ) );
		hPDFPiPi->Draw("same");

		TH1 * hSum = (TH1*) hPDFPiPi->Clone( TString::Format( "template_%s_sum_m%lu", prefix.c_str(), im ) );
		hSum->Add( hPDFPiMu );
		hSum->Add( hPDFMuMu );

		TH1 * hSumBg = (TH1*) hPDFPiPi->Clone( TString::Format( "template_%s_sumbg_m%lu", prefix.c_str(), im ) );
		hSumBg->Add( hPDFPiMu );

		rpl.style( hSum ).set( config, "style.sum" );

		hSum->Draw("same");

		TLatex tl;
		tl.SetTextSize( 12.0 / 360.0 );
		tl.DrawLatexNDC( 0.20, 0.85, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", ff->GetChisquare(), ff->GetNDF(), ff->GetChisquare() / (float)ff->GetNDF() ) );
		tl.DrawLatexNDC( 0.20, 0.8, TString::Format("sigsig = %0.4f #pm %0.4f", ff->GetParameter( "sigsig" ), ff->GetParError( ff->GetParNumber( "sigsig" ) )) );
		tl.DrawLatexNDC( 0.20, 0.75, TString::Format("bgbg = %0.4f #pm %0.4f", ff->GetParameter( "bgbg" ), ff->GetParError( ff->GetParNumber( "bgbg" ) )) );
		tl.DrawLatexNDC( 0.20, 0.70, TString::Format("bgsig = %0.4f #pm %0.4f", ff->GetParameter( "bgsig" ), ff->GetParError( ff->GetParNumber( "bgsig" ) )) );

		tl.SetTextSize( 18.0 / 360.0 );
		tl.DrawLatexNDC( 0.75, 0.90, TString::Format( "%0.2f < M < %0.2f", m1, m2 ) );
		if ( "uls" == prefix )
			tl.DrawLatexNDC( 0.75, 0.86, "unlike-sign" );
		if ( "ls" == prefix )
			tl.DrawLatexNDC( 0.75, 0.86, "like-sign" );
		

		TLegend * leg = new TLegend( 0.2, 0.88, 0.6, 0.95 );
		leg->SetBorderSize( 0 );
		leg->SetNColumns( 5 );
		leg->AddEntry( hpid, "data" );
		leg->AddEntry( hPDFMuMu, "sig+sig" );
		leg->AddEntry( hPDFPiPi, "bg+bg" );
		leg->AddEntry( hPDFPiMu, "bg+sig" );
		leg->AddEntry( hSum, "Fit" );
		leg->Draw();

		can->Print( rpName.c_str() );
	}

	float fweight( TH1 * h1, float x, float y ){
		TH2 * hh = dynamic_cast<TH2*>(h1);
		if ( nullptr == hh )
			return 0;

		int iBin1 = hh->GetXaxis()->FindBin( x );
		int iBin2 = hh->GetYaxis()->FindBin( y );
		float v = hh->GetBinContent( iBin1, iBin2 );
		if ( v <= 0 )
			return 0;
		return v;
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION( INFO );

		book->cd();
		book->makeAll( config, nodePath + ".histograms" );


		can = new TCanvas( "can", "can", config.get<int>( "can:w", 500 ), config.get<int>( "can:h", 500 ) );
		can2 = new TCanvas( "can2", "can2", config.get<int>( "can:w", 500 ), config.get<int>( "can:h2", 500 ) );
		can->SetTopMargin( 0.05 );
		can->SetRightMargin( 0.11 );
		can->SetBottomMargin( 0.11 );
		can->SetLeftMargin( 0.15 );

		can2->SetTopMargin( 0.05 );
		can2->SetRightMargin( 0.01 );
		can2->SetBottomMargin( 0.11 );
		can2->SetLeftMargin( 0.15 );


		rpName = config[nodePath + ".output.Report:url" ];

		can->Print( (rpName+"[").c_str() );
		can->cd();
		
		gStyle->SetOptFit(111);
		gStyle->SetOptStat(0);

		RooPlotLib rpl;


		TH2 * hulsPairPid = get<TH2>( "uls_pair_pid_vs_mass", "pairPid" );
		TH2 * hulsPt      = get<TH2>( "uls_pt_vs_mass", "pairPid" );

		TH2 * hlsPairPid = get<TH2>( "ls_pair_pid_vs_mass", "pairPid" );
		TH2 * hlsPt      = get<TH2>( "ls_pt_vs_mass", "pairPid" );

		TH2 * hulsPairPidrb = HistoBins::rebin2D( "hulsPairPidrb", hulsPairPid, bins["mass"], bins["pairPid"] );
		TH2 * hulsPtrb 		= HistoBins::rebin2D( "hulsPtrb", hulsPt, bins["mass"], bins["pt"] );

		TH2 * hlsPairPidrb = HistoBins::rebin2D( "hlsPairPidrb", hlsPairPid, bins["mass"], bins["pairPid"] );
		TH2 * hlsPtrb 		= HistoBins::rebin2D( "hlsPtrb", hlsPt, bins["mass"], bins["pt"] );

		hulsPairPidrb->Draw("colz");
		can->Print( rpName.c_str() );
		can->SetRightMargin( 0.01 );

		loop_on_mass( hulsPairPidrb, hulsPtrb, "uls" );
		loop_on_mass( hlsPairPidrb, hlsPtrb, "ls" );

		gPad->SetLogy(0);
		book->get( "uls_chi2ndf" )->Draw(  );
		can->Print( rpName.c_str() );
		book->get( "ls_chi2ndf" )->Draw(  );
		can->Print( rpName.c_str() );
		book->get( "uls_purityMuMu" )->Draw(  );
		can->Print( rpName.c_str() );
		book->get( "ls_purityMuMu" )->Draw(  );
		can->Print( rpName.c_str() );



		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}

		can->Print( (rpName+"]").c_str()  );
	} // make

};



#endif