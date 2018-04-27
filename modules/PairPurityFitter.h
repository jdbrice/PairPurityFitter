#ifndef PAIR_PURITY_FITTER_H
#define PAIR_PURITY_FITTER_H

#include "HistoAnalyzer.h"
#include "XmlFunction.h"
#include "XmlHistogram.h"
#include "CutCollection.h"

#include "RooPlotLib.h"

#include "TRandom3.h"

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TLatex.h"

#include "vendor/loguru.h"

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

	size_t nSamples;

public:

	virtual void initialize(){
		HistoAnalyzer::initialize();

		export_img = config.getBool( "can:export", false );
		nSamples = config.getInt( "nSamples", 1000 );

		bins["pairPid"].load( config, "bins.pairPid" );
		bins["mass"].load( config, "bins.mass" );
		bins["pt"].load( config, "bins.pt" );
		bins["ptTemplate"].load( config, "bins.ptTemplate" );
		bins["pid"].load( config, "bins.pid" );

		// build pt projections

		book->cd();
		build_pt_projections( "sig_pos" );
		build_pt_projections( "bg_pos" );


	}

	void build_pt_projections( string hname ){

		TH2 * h2raw = get<TH2>( hname + "_mlp", "mc" );
		TH2 * h2rb = HistoBins::rebin2D( hname + "rb", h2raw, bins["ptTemplate"], bins["pid"] );

		vector<TH1 *> projections;
		for ( size_t i = 0; i <= h2rb->GetXaxis()->GetNbins(); i++ ){
			TH1 * h = h2rb->ProjectionY( TString::Format( "%s_pt_%lu", hname.c_str(), i ), i, i );
			projections.push_back( h );
		}
		nn_pt[ hname ] = projections;
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

		} //loop on mass bins
	} // loop_on_mass


	// NEED
	// pt correlation 
	// pid correlation
	void generate_templates( TH1 * hpt, TH1 * hbgbg, TH1 * hbgsig, TH1 * hsigsig, TH2 * hdelta ){
		assert( hpt != nullptr );
		if ( hpt->Integral() <= 0 ) {
			LOG_F( INFO, "pt is empty" );
			return;
		}
		assert( nn_pt.count( "bg_pos" ) > 0 );
		assert( nn_pt.count( "sig_pos" ) > 0 );

		vector<TH1 * > bg_pt  = nn_pt["bg_pos"];
		vector<TH1 * > sig_pt = nn_pt["sig_pos"];

		for ( size_t i = 0; i < nSamples; i++ ){

			float pt1 = hpt->GetRandom();
			float pt2 = hpt->GetRandom();
			
			int ipt1 = bins[ "ptTemplate" ].findBin( pt1 );
			int ipt2 = bins[ "ptTemplate" ].findBin( pt2 );
			if ( ipt1 < 0 || ipt1 >= sig_pt.size() ) continue;
			if ( ipt2 < 0 || ipt2 >= sig_pt.size() ) continue;


			float IBg1  = nn_pt[ "bg_pos" ][ipt1]->Integral();
			float IBg2  = nn_pt[ "bg_pos" ][ipt2]->Integral();
			float ISig1 = nn_pt[ "sig_pos" ][ipt1]->Integral();
			float ISig2 = nn_pt[ "sig_pos" ][ipt2]->Integral();

			float rBG1  = -999;
			float rBG2  = -999;
			float rSig1 = -999;
			float rSig2 = -999;

			if ( IBg1 > 0 )
				rBG1  = nn_pt[ "bg_pos" ][ipt1]->GetRandom();
			if ( IBg2 > 0 )
				rBG2  = nn_pt[ "bg_pos" ][ipt2]->GetRandom();
			
			if ( ISig1 > 0 )
				rSig1 = nn_pt[ "sig_pos" ][ipt1]->GetRandom();
			if ( ISig2 > 0 )
				rSig2 = nn_pt[ "sig_pos" ][ipt2]->GetRandom();
			

			float pairPid_bgbg   = sqrt( pow( rBG1, 2 ) + pow( rBG2, 2) );
			float pairPid_bgsig  = sqrt( pow( rBG1, 2 ) + pow( rSig2, 2) );
			float pairPid_sigbg  = sqrt( pow( rSig1, 2 ) + pow( rBG2, 2) );
			float pairPid_sigsig = sqrt( pow( rSig1, 2 ) + pow( rSig2, 2) );

			if ( IBg1 > 0 && IBg2 > 0 )
				hbgbg->Fill( pairPid_bgbg );
			if ( ISig1 > 0 && ISig2 > 0 )
				hsigsig->Fill( pairPid_sigsig );
			if ( IBg1 > 0 && ISig2 > 0 )
				hbgsig->Fill( pairPid_bgsig );
			if ( IBg2 > 0 && ISig1 > 0 )
				hbgsig->Fill( pairPid_sigbg );
		}
	} // generate_templates

	void fit_pair_pid( TH1 * hpid, TH1 * hpt, string prefix, size_t im, float m1, float m2 ){
		assert( hpid != nullptr );

		TH1 * hbgbg = nullptr, *hbgsig = nullptr, *hsigsig = nullptr;
		hbgbg   = new TH1F( TString::Format( "template_bgbg_m%lu", im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		hbgsig  = new TH1F( TString::Format( "template_bgsig_m%lu", im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		hsigsig = new TH1F( TString::Format( "template_sigsig_m%lu", im ), "", bins["pairPid"].nBins(), bins["pairPid"].getBins().data() );
		RooPlotLib rpl;
		rpl.style( hbgbg ).set( config, "style.bgbg" );
		rpl.style( hbgsig ).set( config, "style.bgsig" );
		rpl.style( hsigsig ).set( config, "style.sigsig" );

		generate_templates( hpt, hbgbg, hbgsig, hsigsig );

		if ( hsigsig->Integral() <= 0 ) return;
		if ( hbgsig->Integral() <= 0 ) return;
		if ( hbgbg->Integral() <= 0 ) return;

		rpl.style( hpid ).set( config, "style.data" ).draw();
		// hpid->Draw( "pe" );
		gPad->SetLogy(1);

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

		string fitOpt = config[ "fit:opt" ];
		LOG_F( INFO, "Fitting (opt=%s) in range (%f, %f)", fitOpt.c_str(), config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
		hpid->Fit( "ff", fitOpt.c_str(), "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );

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



	} // fit pair_pid


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

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}

		can->Print( (rpName+"]").c_str()  );

	} // make

};



#endif