#ifndef RATIO_MAKER_H
#define RATIO_MAKER_H

#include "HistoAnalyzer.h"
#include "RooPlotLib.h"
#include "XmlRange.h"

class RatioMaker : public HistoAnalyzer
{
protected: 
	map<string, HistoBins> bins;
	TCanvas * can;
	string rpName;

	XmlRange pidRange;


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

		pidRange.loadConfig( config, "pidRange" );
		LOG_F( INFO, "%s", pidRange.toString().c_str() );

	}

	float count_signal( TH1 * huls, TH1 * hls, float cf, float &error ){

		assert( huls && hls );

		TAxis * ax = huls->GetXaxis();
		int ix1 = ax->FindBin( pidRange.min );
		int ix2 = ax->FindBin( pidRange.max );

		double eS = 0;
		double eB = 0;
		float rS = huls->IntegralAndError( ix1, ix2, eS );
		float rB = hls->IntegralAndError( ix1, ix2, eB );

		error = eS + eB*cf;
		float S = rS - ( rB * cf );
		return S;
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION(INFO);
		book->cd();
		book->makeAll( config, nodePath + ".histograms" );
		RooPlotLib rpl;
		rpl.link( book );

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

			hdatauls->Sumw2();
			hdatals->Sumw2();

			huls->Scale( hdatauls->Integral() );
			hls->Scale( hdatals->Integral() );

			huls->Sumw2();
			hls->Sumw2();

			huls = HistoBins::rebin1D( "rb_template_uls_m" + ts(i+1), huls, bins["pairPid"] );
			hls  = HistoBins::rebin1D( "rb_template_uls_m" + ts(i+1), hls, bins["pairPid"] );

			hdatauls = HistoBins::rebin1D( "rb_data_uls_m" + ts(i+1), hdatauls, bins["pairPid"] );
			hdatals  = HistoBins::rebin1D( "rb_data_uls_m" + ts(i+1), hdatals, bins["pairPid"] );

			TH1 * hrdata = (TH1*)hdatauls->Clone( TString::Format( "ratio_data_m%lu", i+1 ) );
			TH1 * hrtemp = (TH1*)huls->Clone( TString::Format( "ratio_template_m%lu", i+1 ) );

			
			hrdata->Divide( hdatals );
			hrtemp->Divide( hls );
			hrtemp->SetLineColor(kRed);
			hrtemp->GetXaxis()->SetRangeUser( 0, 1.42 );
			hrdata->GetYaxis()->SetRangeUser( 0.0, 5.0 );

			// set the bin errors on the template ratio equal to those on the data
			// so that the fit is weighted somewhat properly

			for ( size_t j = 1; j < hrdata->GetXaxis()->GetNbins() + 1; j++ ){
				hrtemp->SetBinError( j, hdatauls->GetBinError( j ) );
				// hrtemp->SetBinError( j, 0.000001 );
			}

			rpl.style( hrdata ).set( config, "style.data" );
			rpl.style( hrtemp ).set( config, "style.template" );
			

			hrdata->Draw("hp");
			hrtemp->Draw("same hist");

			TLatex tl;
			tl.SetTextSize( 16.0 / 360.0 );
			tl.DrawLatexNDC( 0.70, 0.90, TString::Format("%0.2f < M < %0.2f", m1, m2 ) );

			TLegend * leg = new TLegend( 0.2, 0.88, 0.6, 0.95 );
			leg->SetNColumns( 2 );
			leg->SetBorderSize( 0 );
			leg->AddEntry( hrdata, "data" );
			leg->AddEntry( hrtemp, "templates" );
			leg->Draw();

			TF1 * fpol0 = new TF1( "fpol0", "pol0" );
			hrtemp->Fit( fpol0, "NRE", "", pidRange.min, pidRange.max );
			float Rfactor = fpol0->GetParameter(0);
			fpol0->SetRange( pidRange.min, pidRange.max );
			rpl.style( fpol0 ).set( config, "style.fit" ).draw("same");


			tl.DrawLatexNDC( 0.70, 0.85, TString::Format("R=%f", Rfactor ) );
			can->Print( rpName.c_str() );

			LOG_F( INFO, "E = %f", fpol0->GetParameter(0) );
			book->get("R_mass")->SetBinContent( i+1, Rfactor );
			book->get("R_mass")->SetBinError( i+1, 0.001 );

			float error = 0;
			book->get( "N_mass" )->SetBinContent( i+1, count_signal( hdatauls, hdatals, Rfactor, error ) );
			book->get( "N_mass" )->SetBinError( i+1, error );


		}


		rpl.style( "R_mass" ).set( config, "style.R_mass" ).draw();
		can->Print( rpName.c_str() );


		TH1 * h = book->get("N_mass");
		TH1 * hrb = HistoBins::rebin1D( "N_mass_rb", h, bins["sigmass"] );
		hrb->Scale( 1.0 * h->GetBinWidth( 3 ), "width" );
		rpl.style( "N_mass" ).set( config, "style.N_mass" ).draw();
		rpl.style( hrb ).set( config, "style.N_mass_rb" ).draw("same");
		gPad->SetLogy(1);
		can->Print( rpName.c_str() );

		can->Print( (rpName+"]").c_str() );

	}
protected:
	
};


#endif