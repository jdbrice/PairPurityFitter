#ifndef RATIO_MAKER_H
#define RATIO_MAKER_H

#include "HistoAnalyzer.h"
#include "RooPlotLib.h"

class RatioMaker : public HistoAnalyzer
{
protected: 
	map<string, HistoBins> bins;
	TCanvas * can;
	string rpName;
public:
	RatioMaker() {}
	~RatioMaker() {}

	virtual void initialize(){
		HistoAnalyzer::initialize();

		vector<string> paths = config.childrenOf( "bins" );
		for ( auto p : paths ){
			string tn = config.tagName(p);
			LOG_F(  INFO, "Loading Bins: %s", tn.c_str() );
			bins[ tn ].load( config, p );
			LOG_F( INFO, "%s", bins[tn].toString().c_str() );
		}

		rpName = config[nodePath + ".output.Report:url" ];

	}


	virtual void make(){
		LOG_SCOPE_FUNCTION(INFO);
		book->cd();
		book->makeAll( config, nodePath + ".histograms" );


		can = new TCanvas( "can", "can", config.get<int>( "can:w", 500 ), config.get<int>( "can:h", 500 ) );
		
		can->SetTopMargin( 0.05 );
		can->SetRightMargin( 0.01 );
		can->SetBottomMargin( 0.11 );
		can->SetLeftMargin( 0.15 );

		can->Print( (rpName+"[").c_str() );

		gStyle->SetOptStat(0);

		assert( bins.count( "mass" ) > 0 ) ;
		for ( size_t i = 0 ;i < bins["mass"].nBins(); i++ ){



			TH1 * huls = get<TH1>( "template_uls_sumbg_m" + ts( (int) i+1 ), "pairPid" );
			TH1 * hls = get<TH1>( "template_ls_sumbg_m" + ts( (int) i+1 ), "pairPid" );

			LOG_F( INFO, "huls=%p", huls );
			LOG_F( INFO, "hls=%p", hls );

			if ( nullptr == huls || nullptr == hls ) continue;

			float m1 = bins["mass"].getBins()[i];
			float m2 = bins["mass"].getBins()[i+1];

			TH1 * hdatauls = get<TH1>( "uls_pid_m" + ts((int)i+1), "pairPid" );
			TH1 * hdatals = get<TH1>( "ls_pid_m" + ts((int)i+1), "pairPid" );

			huls->Scale( hdatauls->Integral() );
			hls->Scale( hdatals->Integral() );

			TH1 * hrdata = (TH1*)hdatauls->Clone( TString::Format( "ratio_data_m%lu", i+1 ) );
			TH1 * hrtemp = (TH1*)huls->Clone( TString::Format( "ratio_template_m%lu", i+1 ) );

			RooPlotLib rpl;
			hrdata->Divide( hdatals );
			hrtemp->Divide( hls );
			hrtemp->SetLineColor(kRed);
			hrtemp->GetXaxis()->SetRangeUser( 0, 1.42 );
			hrdata->GetYaxis()->SetRangeUser( 0.0, 5.0 );

			rpl.style( hrdata ).set( config, "style.data" );
			rpl.style( hrtemp ).set( config, "style.template" );
			

			hrdata->Draw();
			hrtemp->Draw("same");

			TLatex tl;
			tl.SetTextSize( 16.0 / 360.0 );
			tl.DrawLatexNDC( 0.70, 0.90, TString::Format("%0.2f < M < %0.2f", m1, m2 ) );

			TLegend * leg = new TLegend( 0.2, 0.88, 0.6, 0.95 );
			leg->SetNColumns( 2 );
			leg->SetBorderSize( 0 );
			leg->AddEntry( hrdata, "data" );
			leg->AddEntry( hrtemp, "templates" );
			leg->Draw();

			can->Print( rpName.c_str() );


		}

		can->Print( (rpName+"]").c_str() );

	}
protected:
	
};


#endif