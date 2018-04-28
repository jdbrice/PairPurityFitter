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
TH1 * hPDFSigSig     = nullptr;
TH1 * hPDFBgBg       = nullptr;
TH1 * hPDFBgSig      = nullptr;

Double_t fitfun3( double *x, double *par ){

	double x0=x[0];
	
	double scale_sigsig  = par[0];
	double scale_bgbg    = par[1];
	double scale_bgsig   = par[2];
	
	int sigsig_bin   = hPDFSigSig->GetXaxis()->FindBin( x0 );
	int bgbg_bin     = hPDFBgBg->GetXaxis()->FindBin(x0);
	int bgsig_bin    = hPDFBgSig->GetXaxis()->FindBin(x0);

	float sigsig_val = hPDFSigSig->GetBinContent( sigsig_bin );
	float bgbg_val   = hPDFBgBg->GetBinContent( bgbg_bin );
	float bgsig_val  = hPDFBgSig->GetBinContent( bgsig_bin );

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

	string bg_source;
	string sig_source;

public:

	virtual void initialize(){
		HistoAnalyzer::initialize();

		export_img = config.getBool( "can:export", false );
		weight_by_pt = config.getBool( "weight:pt", false );
		weight_by_pid = config.getBool( "weight:pid", false );
		nSamples = config.getInt( "nSamples", 1000 );

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

		if ( weight_by_pt ){
			build_mass_projections( "uls_delta_pt", "deltaPt", "deltaPt" );
			build_mass_projections( "ls_delta_pt", "deltaPt", "deltaPt" );
		}
		if ( weight_by_pid ){
			build_mass_projections( "uls_delta_pid", "minpid", "deltaPid" );
			build_mass_projections( "ls_delta_pid", "minpid", "deltaPid" );
		}
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
			string projName = prefix + "_pid_mass" + ts( (int) i);
			string projNamePt = prefix + "_pt_mass" + ts( (int) i);
			TH1 * hpid = h2pid->ProjectionY( projName.c_str(), i, i );
			TH1 * hpt  = h2pt->ProjectionY( projNamePt.c_str(), i, i );

			float m1 = h2pid->GetXaxis()->GetBinLowEdge( i );
			float m2 = h2pid->GetXaxis()->GetBinUpEdge( i );

			hpid->SetTitle( TString::Format( "%0.2f < M < %0.2f", m1, m2 ) );

			if ( hpid->Integral() <= 0 ) continue;

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
	void generate_templates( TH1 * hpt, TH1 * hbgbg, TH1 * hbgsig, TH1 * hsigsig, TH2 * hdeltaPid, TH2 * hdeltaPt, TH2 * hDataDeltaPid = nullptr, TH2 * hDataDeltaPt = nullptr, float Ibgbg = 1.0, float Ibgsig = 1.0, float Isigsig = 1.0 ){
		LOG_SCOPE_FUNCTION( INFO );
		assert( hpt != nullptr );
		assert( nullptr != hbgbg );
		assert( nullptr != hbgsig );
		assert( nullptr != hsigsig );
		assert( nullptr != hdeltaPid );
		assert( nullptr != hdeltaPt );

		if ( Ibgbg   <= 0 ) Ibgbg   = 1.0;
		if ( Ibgsig  <= 0 ) Ibgsig  = 1.0;
		if ( Isigsig <= 0 ) Isigsig = 1.0;

		if ( hpt->Integral() <= 0 ) {
			LOG_F( INFO, "pt is empty" );
			return;
		}
		assert( nn_pt.count( this->bg_source ) > 0 );
		assert( nn_pt.count( this->sig_source ) > 0 );

		bool weightPt  = false;
		bool weightPid = false;
		if ( nullptr != hDataDeltaPt )
			weightPt = true;
		if ( nullptr != hDataDeltaPid )
			weightPid = true;

		vector<TH1 * > bg_pt  = nn_pt[ this->bg_source ];
		vector<TH1 * > sig_pt = nn_pt[ this->sig_source ];

		for ( size_t i = 0; i < nSamples; i++ ){

			float pt1 = hpt->GetRandom();
			float pt2 = hpt->GetRandom();
			
			int ipt1 = bins[ "ptTemplate" ].findBin( pt1 );
			int ipt2 = bins[ "ptTemplate" ].findBin( pt2 );
			if ( ipt1 < 0 || ipt1 >= sig_pt.size() ) continue;
			if ( ipt2 < 0 || ipt2 >= sig_pt.size() ) continue;

			if ( false == weightPt )
				hdeltaPt->Fill( TMath::Min( pt1, pt2 ), fabs( pt1 - pt2 ) );

			float IBg1  = bg_pt[ipt1]->Integral();
			float IBg2  = bg_pt[ipt2]->Integral();
			float ISig1 = sig_pt[ipt1]->Integral();
			float ISig2 = sig_pt[ipt2]->Integral();

			float rBG1  = -999;
			float rBG2  = -999;
			float rSig1 = -999;
			float rSig2 = -999;

			if ( IBg1 > 0 )
				rBG1  = bg_pt[ipt1]->GetRandom();
			if ( IBg2 > 0 )
				rBG2  = bg_pt[ipt2]->GetRandom();
			
			if ( ISig1 > 0 )
				rSig1 = sig_pt[ipt1]->GetRandom();
			if ( ISig2 > 0 )
				rSig2 = sig_pt[ipt2]->GetRandom();

			float pairPid_bgbg   = sqrt( pow( rBG1, 2 ) + pow( rBG2, 2) );
			float pairPid_bgsig  = sqrt( pow( rBG1, 2 ) + pow( rSig2, 2) );
			float pairPid_sigbg  = sqrt( pow( rSig1, 2 ) + pow( rBG2, 2) );
			float pairPid_sigsig = sqrt( pow( rSig1, 2 ) + pow( rSig2, 2) );


			if ( false == weightPid ){
				hdeltaPid->Fill( TMath::Min( rBG1, rBG2 ), fabs( rBG1 - rBG2 ), Ibgbg );
				hdeltaPid->Fill( TMath::Min( rBG1, rSig2 ), fabs( rBG1 - rSig2 ), 0.5 * Ibgsig );
				hdeltaPid->Fill( TMath::Min( rSig1, rBG2 ), fabs( rSig1 - rBG2 ), 0.5 * Ibgsig );
				hdeltaPid->Fill( TMath::Min( rSig1, rSig2 ), fabs( rSig1 - rSig2 ), Isigsig );
			}

			float wpt = 1.0;
			if ( weightPt ){
				float minPt = TMath::Min( pt1, pt2 );
				float diffPt = fabs( pt1 - pt2 );
				float mixedW = fweight( hdeltaPt, minPt, diffPt );
				float sameW = fweight( hDataDeltaPt, minPt, diffPt );
				if ( mixedW > 0 && sameW > 0){
					wpt = sameW / mixedW;
				}
				
			}
			float pairPid_bgbg_wpid = 1.0;
			float pairPid_bgsig_wpid = 1.0;
			float pairPid_sigbg_wpid = 1.0;
			float pairPid_sigsig_wpid = 1.0;
			if ( weightPid ){
				// bg bg
				float minPid = TMath::Min( rBG1, rBG2 );
				float diffPid = fabs( rBG1 - rBG2 );
				float mixedW = fweight( hdeltaPid, minPid, diffPid );
				float sameW = fweight( hDataDeltaPid, minPid, diffPid );
				if ( mixedW > 0 ){
					pairPid_bgbg_wpid = sameW / mixedW;
				}

				// bg sig
				minPid = TMath::Min( rBG1, rSig2 );
				diffPid = fabs( rBG1 - rSig2 );
				mixedW = fweight( hdeltaPid, minPid, diffPid );
				sameW = fweight( hDataDeltaPid, minPid, diffPid );
				if ( mixedW > 0 ){
					pairPid_bgsig_wpid = sameW / mixedW;
				}

				// sig bg
				minPid = TMath::Min( rSig1, rBG2 );
				diffPid = fabs( rSig1 - rBG2 );
				mixedW = fweight( hdeltaPid, minPid, diffPid );
				sameW = fweight( hDataDeltaPid, minPid, diffPid );
				if ( mixedW > 0 ){
					pairPid_sigbg_wpid = sameW / mixedW;
				}

				// sig sig
				minPid = TMath::Min( rSig1, rSig2 );
				diffPid = fabs( rSig1 - rSig2 );
				mixedW = fweight( hdeltaPid, minPid, diffPid );
				sameW = fweight( hDataDeltaPid, minPid, diffPid );
				if ( mixedW > 0 ){
					pairPid_sigsig_wpid = sameW / mixedW;
				}
			}

			if ( IBg1 > 0 && IBg2 > 0 )
				hbgbg->Fill( pairPid_bgbg, wpt * pairPid_bgbg_wpid * Ibgbg );
			if ( ISig1 > 0 && ISig2 > 0 )
				hsigsig->Fill( pairPid_sigsig, wpt * pairPid_sigsig_wpid * Isigsig );
			if ( IBg1 > 0 && ISig2 > 0 )
				hbgsig->Fill( pairPid_bgsig );
			if ( IBg2 > 0 && ISig1 > 0 )
				hbgsig->Fill( pairPid_sigbg );
		}
	} // generate_templates

	void fit_pair_pid( TH1 * hpid, TH1 * hpt, string prefix, size_t im, float m1, float m2 ){
		LOG_SCOPE_FUNCTION( INFO );
		LOG_F( INFO, "Fitting %s[%lu], %0.3f < M < %0.3f", prefix.c_str(), im, m1, m2 );
		assert( hpid != nullptr );

		TH1 * hbgbg = nullptr, *hbgsig = nullptr, *hsigsig = nullptr;
		TH2 * hdeltaPid = nullptr, *hdeltaPt = nullptr;
		hbgbg   = new TH1F( TString::Format( "template_%s_bgbg_m%lu", prefix.c_str(), im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		hbgsig  = new TH1F( TString::Format( "template_%s_bgsig_m%lu", prefix.c_str(), im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		hsigsig = new TH1F( TString::Format( "template_%s_sigsig_m%lu", prefix.c_str(), im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		hdeltaPid  = new TH2F( TString::Format( "template_%s_deltaPid_m%lu", prefix.c_str(), im ), "", bins["pid"].nBins(), bins["pid"].getBins().data(), bins["deltaPid"].nBins(), bins["deltaPid"].getBins().data() );
		hdeltaPt   = new TH2F( TString::Format( "template_%s_deltaPt_m%lu", prefix.c_str(), im ), "", bins["pt"].nBins(), bins["pt"].getBins().data(), bins["deltaPt"].nBins(), bins["deltaPt"].getBins().data() );
		RooPlotLib rpl;
		rpl.style( hbgbg ).set( config, "style.bgbg" );
		rpl.style( hbgsig ).set( config, "style.bgsig" );
		rpl.style( hsigsig ).set( config, "style.sigsig" );

		// first pass builds the correlations
		generate_templates( hpt, hbgbg, hbgsig, hsigsig, hdeltaPid, hdeltaPt );

		/* First Pass Fitting
		 * use the vanilla mixed templates without any correlation information
		 */
		if ( hsigsig->Integral() <= 0 ) return;
		if ( hbgsig->Integral() <= 0 ) return;
		if ( hbgbg->Integral() <= 0 ) return;

		hsigsig->Scale( 1.0 / hsigsig->Integral() );
		hbgbg->Scale( 1.0 / hbgbg->Integral() );
		hbgsig->Scale( 1.0 / hbgsig->Integral() );

		hPDFSigSig = hsigsig;
		hPDFBgSig  = hbgsig;
		hPDFBgBg   = hbgbg;

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
		delete ff;

		// if ( weight_by_pt || weight_by_pid ){

		// 	generate_templates( hpt, hbgbg, hbgsig, hsigsig, hdeltaPid, hdeltaPt, nullptr, nullptr, hbgbg->Integral(), hbgsig->Integral(), hsigsig->Integral() );

		// 	if ( hsigsig->Integral() <= 0 ) return;
		// 	if ( hbgsig->Integral() <= 0 ) return;
		// 	if ( hbgbg->Integral() <= 0 ) return;

		// 	hsigsig->Scale( 1.0 / hsigsig->Integral() );
		// 	hbgbg->Scale( 1.0 / hbgbg->Integral() );
		// 	hbgsig->Scale( 1.0 / hbgsig->Integral() );

		// 	hPDFSigSig = hsigsig;
		// 	hPDFBgSig  = hbgsig;
		// 	hPDFBgBg   = hbgbg;

		// 	TF1 * ff = nullptr;
		// 	ff = new TF1( "ff", fitfun3, -1, 2, 3 );
		// 	ff->SetParNames( "sigsig", "bgbg", "bgsig" );
		// 	ff->SetParLimits( 0, 1e-4, 1 );
		// 	ff->SetParLimits( 1, 1e-4, 1 );
		// 	ff->SetParLimits( 2, 1e-4, 1 );
		// 	ff->SetParameters( 0.1, 0.1, 0.01 );

		// 	ff->SetNpx( 1000 );
		// 	ff->SetLineColor(kBlack);
		// 	ff->SetLineWidth(2);

		// 	string fitOpt = config[ "fit:opt" ] + "S";
		// 	LOG_F( INFO, "Fitting (opt=%s) in range (%f, %f)", fitOpt.c_str(), config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
		// 	TFitResultPtr fitPointer = hpid->Fit( "ff", fitOpt.c_str(), "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );

		// 	report_fit( hpid, ff, prefix, im, m1, m2 );
		// 	delete ff;
		// }

		// NOTES:
		// what we should do
		// instead of filling a single delta, separate by bgbg, bgsig, sigsig
		// 1. fit to the vanilla mixed templates
		// 2. use the first pass yields to make a properly weighted delta = Y_bgbg * delta_bgbg + Y_bgsig * delta_bgsig + Y_sigsig * delta_sigsig
		// 3. generate new templates with weighting
		// 4. Do final fit

		if ( weight_by_pt || weight_by_pid ){
			// second pass applies correlation weighting

			TH2 * hdpt = nullptr, *hdpid = nullptr;
			if ( deltas.count( prefix + "_delta_pt" ) > 0 && im < deltas[ prefix + "_delta_pt" ].size() ){
				hdpt = deltas[ prefix + "_delta_pt" ][im];
				hdpt->Scale( 1.0 / hdpt->Integral() );
			}
			if ( deltas.count( prefix + "_delta_pid" ) > 0 && im < deltas[ prefix + "_delta_pid" ].size() ){
				hdpid = deltas[ prefix + "_delta_pid" ][im];
				hdpid->Scale( 1.0 / hdpid->Integral() );
			}

			hdeltaPid->Scale(1.0 / hdeltaPid->Integral() );
			hdeltaPt->Scale(1.0 / hdeltaPt->Integral() );

			hbgbg->Reset();
			hbgsig->Reset();
			hsigsig->Reset();
			generate_templates( hpt, hbgbg, hbgsig, hsigsig, hdeltaPid, hdeltaPt, hdpid, hdpt );

			/* Second Pass Fitting 
			 * Using the weighted templates
			 */
			if ( hsigsig->Integral() <= 0 ) return;
			if ( hbgsig->Integral() <= 0 ) return;
			if ( hbgbg->Integral() <= 0 ) return;

			hsigsig->Scale( 1.0 / hsigsig->Integral() );
			hbgbg->Scale( 1.0 / hbgbg->Integral() );
			hbgsig->Scale( 1.0 / hbgsig->Integral() );

			hPDFSigSig = hsigsig;
			hPDFBgSig  = hbgsig;
			hPDFBgBg   = hbgbg;

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

			LOG_F( INFO, "Fitting 2nd Pass (opt=%s) in range (%f, %f)", fitOpt.c_str(), config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
			fitPointer = hpid->Fit( "ff", fitOpt.c_str(), "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );

			report_fit( hpid, ff, prefix, im, m1, m2 );
			delete ff;
		} // weighted second pass
	} // fit pair_pid

	void report_fit( TH1 * hpid, TF1 * ff, string prefix, size_t im, float m1, float m2 ) {
		LOG_SCOPE_FUNCTION( INFO );
		assert( nullptr != hpid );
		assert( nullptr != ff );
		RooPlotLib rpl;
		rpl.style( hpid ).set( config, "style.data" ).draw();
		gPad->SetLogy(1);

		hPDFSigSig->Scale( ff->GetParameter( "sigsig" ) );
		hPDFSigSig->Draw("same");

		hPDFBgSig->Scale( ff->GetParameter( "bgsig" ) );
		hPDFBgSig->Draw("same");

		hPDFBgBg->Scale( ff->GetParameter( "bgbg" ) );
		hPDFBgBg->Draw("same");

		TH1 * hSum = (TH1*) hPDFBgBg->Clone( TString::Format( "template_sum_m%lu", im ) );
		hSum->Add( hPDFBgSig );
		hSum->Add( hPDFSigSig );

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
		leg->AddEntry( hPDFSigSig, "sig+sig" );
		leg->AddEntry( hPDFBgBg, "bg+bg" );
		leg->AddEntry( hPDFBgSig, "bg+sig" );
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
		// loop_on_mass( hlsPairPidrb, hlsPtrb, "ls" );

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}

		can->Print( (rpName+"]").c_str()  );
	} // make

};



#endif