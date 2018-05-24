#ifndef CUT_OPTIMIZATION_MAKER_H
#define CUT_OPTIMIZATION_MAKER_H

#include "HistoAnalyzer.h"
#include "RooPlotLib.h"

class CutOptimizationMaker : public HistoAnalyzer {
protected: 
	map<string, HistoBins> bins;
	TCanvas * can;
	string rpName;

	bool cumulative;
public:
	CutOptimizationMaker() {}
	~CutOptimizationMaker() {}

	virtual void initialize(){
		HistoAnalyzer::initialize();

		// Load the histogram bins
		vector<string> paths = config.childrenOf( "bins" );
		for ( auto p : paths ){
			string tn = config.tagName(p);
			LOG_F(  INFO, "Loading Bins: %s", tn.c_str() );
			bins[ tn ].load( config, p );
			LOG_F( INFO, "%s", bins[tn].toString().c_str() );
		}

		rpName = config[nodePath + ".output.Report:url" ];

		cumulative = config.get<bool>(nodePath + ":cumulative", true );

	}


	virtual void make(){
		LOG_SCOPE_FUNCTION(INFO);
		book->cd();
		book->makeAll( config, nodePath + ".histograms" );

		RooPlotLib rpl;
		rpl.link( book ) ;

		can = new TCanvas( "can", "can", config.get<int>( "can:w", 500 ), config.get<int>( "can:h", 500 ) );
		
		can->SetTopMargin( 0.05 );
		can->SetRightMargin( 0.01 );
		can->SetBottomMargin( 0.11 );
		can->SetLeftMargin( 0.11 );

		can->Print( (rpName+"[").c_str() );

		gStyle->SetOptStat(0);

		assert( bins.count( "mass" ) > 0 ) ;
		string sign = "uls";
		for ( size_t i = 0 ;i < bins["mass"].nBins(); i++ ){

			string n = string(TString::Format( "template_%s_%s_m%lu", sign.c_str(), "pipi", i+1 ));
			LOG_F( INFO, "n=%s", n.c_str() );
			TH1 * hpipi  = get<TH1>( n, "pairPid" );
			TH1 * hpimu  = get<TH1>( string(TString::Format( "template_%s_%s_m%lu", sign.c_str(), "pimu", i+1 )), "pairPid" );
			TH1 * hmumu  = get<TH1>( string(TString::Format( "template_%s_%s_m%lu", sign.c_str(), "mumu", i+1 )), "pairPid" );
			TH1 * hsum   = get<TH1>( string(TString::Format( "template_%s_%s_m%lu", sign.c_str(), "sum", i+1 )), "pairPid" );
			TH1 * hsumbg = get<TH1>( string(TString::Format( "template_%s_%s_m%lu", sign.c_str(), "sumbg", i+1 )), "pairPid" );

			TH1 * hdata     = get<TH1>( string(TString::Format( "%s_pid_m%lu", sign.c_str(), i+1 )), "pairPid" );
			TH1 * hdatanorm = get<TH1>( string(TString::Format( "%s_pid_m%lu", sign.c_str(), i+1 )), "pairPid" );
			

			if ( nullptr == hpipi ) continue;
			if ( nullptr == hpimu ) continue;
			if ( nullptr == hmumu ) continue;
			if ( nullptr == hsum ) continue;
			if ( nullptr == hsumbg ) continue;
			if ( nullptr == hdata ) continue;
			if ( nullptr == hdatanorm ) continue;

			float Idata = hdata->Integral();
			hpipi->Scale( Idata );
			hpimu->Scale( Idata );
			hmumu->Scale( Idata );
			hsum->Scale( Idata );
			hsumbg->Scale( Idata );


			float m1 = bins["mass"].getBins()[i];
			float m2 = bins["mass"].getBins()[i+1];

			LOG_F( INFO, "" );

			LOG_F( INFO, "hpipi=%p", hpipi );
			LOG_F( INFO, "hpimu=%p", hpimu );
			LOG_F( INFO, "hmumu=%p", hmumu );
			LOG_F( INFO, "hsum=%p", hsum );
			LOG_F( INFO, "hsumbg=%p", hsumbg );

			TH1 * heff_sig      = (TH1*)hmumu->Clone( TString::Format( "eff_%s_%s_m%lu", sign.c_str(), "sig", i ) );
			TH1 * hpurity_sig      = (TH1*)hmumu->Clone( TString::Format( "purity_%s_%s_m%lu", sign.c_str(), "sig", i ) );
			TH1 * hpurity_bg       = (TH1*)hmumu->Clone( TString::Format( "purity_%s_%s_m%lu", sign.c_str(), "bg", i ) );
			TH1 * hpurity_sigbg    = (TH1*)hmumu->Clone( TString::Format( "purity_%s_%s_m%lu", sign.c_str(), "sigbg", i ) );
			TH1 * hsignificance    = (TH1*)hmumu->Clone( TString::Format( "significance_%s_m%lu", sign.c_str(), i ) );
			TH1 * hmaxsignificance = (TH1*)hmumu->Clone( TString::Format( "maxsignificance_%s_m%lu", sign.c_str(), i ) );
			heff_sig->Reset();
			hpurity_sig->Reset();
			hpurity_bg->Reset();
			hpurity_sigbg->Reset();
			hsignificance->Reset();
			hmaxsignificance->Reset();

			float max_significance = 0;
			int imax_significance = 1;

			LOG_F( INFO, "hpurity_sig->GetXaxis()->GetNbins()=%d", hpurity_sig->GetXaxis()->GetNbins() );
			for ( int j = 1; j < hpurity_sig->GetXaxis()->GetNbins(); j ++ ){

				float S    = hmumu->Integral( j, -1 );
				float Stot = hmumu->Integral();
				float B    = hsumbg->Integral( j, -1 );

				if ( false == cumulative  ) {
					S    = hmumu->Integral( j, j );
					B    = hsumbg->Integral( j, j );
				}

				float psig = S / ( S+B );
				if ( psig == psig )
					hpurity_sig->SetBinContent( j, psig );

				float pbg = B / ( S+B );
				if ( pbg == pbg )
					hpurity_bg->SetBinContent( j, pbg );

				if ( (S/B) == (S/B) )
					hpurity_sigbg->SetBinContent( j, (S / B) * 0.1 );
				else {
					LOG_F( INFO, "NAN %d", j );
				}

				if ( S/Stot == S/Stot )
					heff_sig->SetBinContent( j, (S / Stot) );


				float significance = S / sqrt( S + B ) * 0.01;
				if ( significance == significance )
					hsignificance->SetBinContent( j, significance  );
				if ( significance > max_significance ){
					max_significance = significance ;
					imax_significance = j;
				}
			} // loop on j

			LOG_F( INFO, "[%d]max=%f", imax_significance, max_significance );
			hmaxsignificance->SetBinContent( imax_significance, max_significance );

			book->get( "optimal_cut" )->SetBinContent( i+1,  hpurity_sig->GetBinLowEdge(imax_significance) );

			rpl.style( hpurity_sig ).set( config, "style.sig" );
			rpl.style( hpurity_bg ).set( config, "style.bg" );
			rpl.style( hpurity_sigbg ).set( config, "style.sigbg" );
			rpl.style( hsignificance ).set( config, "style.significance" );
			rpl.style( hmaxsignificance ).set( config, "style.maxsignificance" );
			rpl.style( heff_sig ).set( config, "style.eff" );
			
			hpurity_sig->Draw("C");
			heff_sig->Draw("same C");
			// hpurity_bg->Draw("same C");
			hpurity_sigbg->Draw("same C");
			hsignificance->Draw("same C");
			hmaxsignificance->Draw("same pe");

			TLatex tl;
			tl.SetTextSize( 16.0 / 360.0 );
			tl.DrawLatexNDC( 0.70, 0.90, TString::Format("%0.2f < M < %0.2f (GeV/c^{2})", m1, m2 ) );


			TLegend * leg = new TLegend( 0.2, 0.73, 0.6, 0.95 );
			leg->SetTextSize( 16.0 / 360.0 );
			leg->SetNColumns( 2 );
			leg->SetBorderSize( 0 );
			leg->AddEntry( hpurity_sig, "signal purity" );
			leg->AddEntry( heff_sig, "signal efficiency" );
			leg->AddEntry( hpurity_sigbg, "S / B #times 0.1" );
			leg->AddEntry( hsignificance, "S / sqrt( S + B ) #times 0.01" );
			leg->AddEntry( hmaxsignificance, "Max( S / sqrt( S + B ) ) #times 0.01", "p" );
			leg->Draw();


			can->Print( rpName.c_str() );


		}

		rpl.style( "optimal_cut" ).set(config, "style.optimal_cut");
		book->get( "optimal_cut" )->Draw();


		can->Print( rpName.c_str() );

		can->Print( (rpName+"]").c_str() );

	}
protected:
	
};


#endif