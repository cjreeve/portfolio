//sfa7
// group of functions for loading and analysing the data produced by vmelist2root
#define NDETS 24
 
#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCutG.h"

#include "TFormula.h"
#include "TRandom.h"
#include "TFunction.h"
#include "TMethodCall.h"
#include "TObjString.h"
#include "TError.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TColor.h"
#include "TDiamond.h"

#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
	
#include "TVector.h"
#include "TMatrix.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstring>
#include <vector>

using namespace std;

struct event{
Int_t	 mod;					// module number - bit useless
Int_t	 trig_det;				// the detector number that had the highest energy of the last adc module... bit useless
Double_t hi_tim;				// high presision time. Multiply by 8E-6 to convert to seconds. Cycles every few hours?
Double_t time;					// time of event in unix time
Double_t energy[NDETS];			// energy of detectors
Int_t	 status[4];				// temperature1, temperature2, humidity, flowmeter
Int_t	 flag;					// event flag. Reject if flag==1;
Int_t	 veto;	    			// used to store TTL flag.
Int_t    runno;					// the squential run number that the event belongs to
Int_t    groupno; 				// number to which the calibration group the event belongs. starts at 0
};
#define EVTREE2buffer sprintf(buffer,\
	"mod/I:trig_det/I:hi_tim/D:time/D:energy[%d]/D:status[4]/I:flag/I:veto/I:runno/I:groupno/I",\
	NDETS);

struct run_stats{
Int_t	 firstevno; 			// event number of the first event in the run
Int_t    starttime;				// time of the start of the run in unix time
Int_t    groupno;				// the group that the run belongs to
Int_t    flag[NDETS];			// flags used to reject runs. reject if flag[det-1]==1
Int_t	 flaggedtime[NDETS];	// used for small cuts... not yet used.
Int_t    nevents;				// number of saved events per run
Int_t	 nvetoevents;			// the number of HV veto events
Int_t	 neventsall;			// actual number of events in daq file. Does not including nvetoedevents. Deadtime = (neventsall+nevetoedevents)* deadTimePerEvent
Int_t	 det_on[NDETS];			// whether the detector is on for this run. This number is set by the daq. Use data if det_on==1
Int_t	 hv[NDETS];				// the high voltage setting reported by the daq
Int_t	 amp[NDETS];			// the amplifacation setting (0-127) reported by the daq
Int_t	 thr[NDETS];			// the thresholds setting (0-127) reported by the daq
Int_t	 runlen;				// the run length in seconds (usually 3600s unless short runs are accepted when compiling data)
Double_t lo_thresh[NDETS];		// the threshold above which the data can be analysied (keV)
};
#define RSTREE2buffer sprintf(buffer,\
	"firstevno/I:starttime/I:groupno/I:flag[%d]/I:flaggedtime[%d]/I:nevents/I:nvetoedevents/I:neventsall/I:det_on[%d]/I:hv[%d]/I:amp[%d]/I:thr[%d]/I:runlen/I:lo_thresh[%d]/D",\
	NDETS,NDETS,NDETS,NDETS,NDETS,NDETS,NDETS);

struct group_stats{
Int_t	 firstrunno; 			// run number of first run in group
Int_t    starttime[NDETS];		// start time of group (for each detector)  // not much use
Int_t    nruns;					// number of runs in a group
Int_t    runs_toobig;			// number of runs rejected due to file size cut
Int_t    det_on[NDETS];			// specifies if a detector is on or off. This is set by the user in group_stats
Double_t hi_thresh[NDETS];		// upper energy threshold of data
Double_t lo_thresh[NDETS];		// lower threshold can be set high by user in groups_stats. This is then moved over to rs.lo_thresh in sfaplot.
Double_t cal_eqn[4][NDETS]; 	// ={slope,error,constant,error}
Double_t res_eqn[4][NDETS];		// ={slope,error,constant,error}
Int_t    group_name[50];		// name of the group. Replace each number with its corespondic ASCI character to reviel.
Int_t    flag;					// to flag whole group and groups sellections. Useful for comparing results from one group against another.
Int_t	 nevents;				// number of events in the group
};
#define GSTREE2buffer sprintf(buffer,\
	"firstrunno/I:starttime[%d]/I:nruns/I:runs_toobig/I:det_on[%d]/I:hi_thresh[%d]/D:lo_thresh[%d]/D:cal_eqn[4][%d]/D:res_eqn[4][%d]/D:group_name[50]/I:flag/I:nevents/I",\
	NDETS,NDETS,NDETS,NDETS,NDETS,NDETS);

struct eppstruct{
	vector <Double_t> epp;			// number of events per period
	vector <Double_t> time;			// time of start of period
	vector <Int_t>    runno;		// sequential run number of first run in period
	vector <Int_t>    used_run_cnt; // number of runs used in period
} epp;

struct CutStats{
	Double_t mean;
	Double_t error;
	Double_t probX2;
	Int_t availruns;
	Int_t extracut;
} *cs;


// #ifdef __CINT__
// #pragma link C++ class std::vector<Int_t>+;
// #pragma link C++ class std::vector<Double_t>+;
// #endif


#define SKIPFLAGS()																				\
		while((gs[ev.groupno].flag || !gs[ev.groupno].det_on[detn-1]) && i<nentries){		\
			i += gs[ev.groupno].nevents;														\
			evtree->GetEvent(i);																\
		}																						\
																								\
		while( (!rs[ev.runno].det_on || (use_flags[1]!='_' && rs[ev.runno].flag[detn-1]==(use_flags[1]=='0'))) && i<nentries){	\
			i += rs[ev.runno].nevents; 														\
			evtree->GetEvent(i);																\
		}																						\


// Global variables
event ev;
run_stats *rs;
group_stats *gs;
Double_t gbl_deadtimeperevent=1.2E-4;  // this is the approx deadtime per event Uncertainty? What impact will this have? Check!!!
Double_t samescalef=1;
Long64_t nentries;
Int_t ngroups, nruns;
Int_t ndetectors = NDETS;
Int_t *flag;
Int_t *BFlag;
TCanvas *c1;
TGraph  *g_waterfall;
TGraph *g_handle, *g_handle2;
TH1D   *h_handle, *h_handle2;
TH2D   *h2_handle, *h2_handle2;
Char_t *fname;
TFile *ssfile;
TTree *evtree;
Int_t fit_locut =0, fit_hicut = 0; // used for output of fit_poisson
Int_t gblperiod =1; // global user period, ie the number of runs per period
Int_t gbl_verbose=1;
Int_t gbl_draw=1;
TF1 *gblfit; // used for fit_poisson. needs to be global if I want to be able to get information of fit externally.
int gblpause=0;
Int_t gbldetx, gbldety;
vector <Double_t> gblresults;
// vector <Double_t> *gblcuts = new vector <Double_t> [NDETS];
// vector <Double_t> *newcuts = new vector <Double_t> [NDETS];	



// decleration of some functions
Double_t SSgetdeadtime(Int_t det);
Double_t SSgetusedtime(Int_t detn, Int_t incl_flaggedtime=1, Int_t incl_deadtime=1);
//Int_t getentries(Int_t group=-1, Int_t return_output=0); outof date
void llfn_getgroupname(int groupno);
void SShelp();	
void llfn_makewhite(TCanvas *c1);
void SSgroups_set(Char_t use_groups[20]=" ");
Int_t SSgroups_get();
void SSfit_poisson(Double_t mean =-1, TH1D *hist1 = h_handle, Int_t min=0, Int_t max=100);
void SSfit_gaus(Double_t mean =-1., Double_t sigma = -0.10, Double_t scale = -1., TH1D *hist1 = h_handle);
void SSfit_poly();
void llfn_poly1(Double_t *a, Double_t x1, Double_t y1, Double_t x2, Double_t y2);
//subfunctions
char *llfn_getname(FILE *fd);
Int_t llfn_time2runno(Double_t time);
Double_t llfn_time2unixtime(Double_t time);
void llfn_runno2runname(Int_t starttime, char *runlist);
Double_t SSgetuptime(Int_t detn);
Double_t SSgetflaggedtime(Int_t detn);
Double_t SSgetusedtime(Int_t detn, Int_t incl_flaggedtime, Int_t incl_deadtime);
Double_t SSgetusedtimeR(Int_t detn, Double_t Emin, Double_t Emax, Double_t deadtimeperevent=1.2E-4);
Double_t SSgetdeadtime(Int_t det);
Int_t SSgetdatetimedif(int i=0, int f=nentries-1);
void SSgetResEqn(Int_t group, Int_t det);
void SSaddleggend(char title1[],char title2[],TH1D *h_handle1, TH1D *h_handle2);
void SSprint(Int_t opt = 0, char name[64]="c1", TCanvas *c1=c1);
void llfn_setstyles(char OptStat[16]="e");
void SScheckorder();
// Double_t SSgetcuttime(Int_t detn);
// void llfn_addcuts();
int llfn_NOT(int val, int sign);
Int_t llfn_Decimal2Bit(Int_t decimal, Int_t bitN, Int_t Nbits);


/// ################################################ FUNDERMENTAL FUNCTIONS ##############################################


	/// ---------------------------------------------- low level functions -----------------------------------------------
	

/// ############################################ END OF FUNDERMENTAL FUNCTIONS ###########################################	


/// #################################################### MAIN FUNCTIONS ##################################################





/// #################################################### OLD FUNCTIONS ##################################################


// calibration functuion threshold number to channel number
// TF1 SSth2ch("th2ch","26.6*x+24+5",0,3000);
// TF1 SSch2th("ch2th","1./26.6*x-24/26.6",0,127);
TF1 SSth2ch("th2ch","26.6*x+24+26",0,3000);
TF1 SSch2th("ch2th","1./26.6*x-50/26.6",0,127);
Double_t SSe2ch(Double_t E, Int_t det, Int_t group=1) {
	return E/gs[group-1].cal_eqn[0][det-1] - gs[group-1].cal_eqn[2][det-1]/gs[group-1].cal_eqn[0][det-1];
}
Double_t SSch2e(Double_t Ch, Int_t det, Int_t group=1) {
	return gs[group-1].cal_eqn[0][det-1]*Ch + gs[group-1].cal_eqn[2][det-1];
}




Double_t SSget_coincident_usedtime(Int_t *dets) {
	// this function gets the used time for periods when all of the "dets" are not flagged.
	// The function only takes into account runs so far
	
	int i,j;
	TH1D *f1 = new TH1D("f1","f1",nruns,0,nruns);
	TH1D *f2 = new TH1D("f2","f2",nruns,0,nruns);
	
	// put the used runs of dets[0] into a histogram
	for(j=0;j<nruns;++j) if(rs[j].flag[dets[0]] && gs[rs[j].groupno].flag) f1->Fill(j,1);
	for(i=1;dets[i];++i) {
		// put the used runs of the next dets[] into a histogram
		for(j=0;j<nruns;++j) if(rs[j].flag[dets[i]] && !gs[rs[j].groupno].flag) f2->Fill(j,1);
		// f2 to f1.
		f1->Add(f2);
	}
	
	Int_t time=0;
	for(j=0;j<nruns;++j) if(f1->GetBinContent(j+1)==0) ++time;
	return Double_t(time);
}
		
void SScut_plot_coincedencesHelp() {
	cout << "This function serves two porposes. \n"
	     << "It can plot coinidences and or can flag coincidences of event across detectors.\n\n"
		 << "Specify the detectors to the analysed by creating an integer array of the required\n"
		 << "   detector numbers setting the last one in the list to 0 which marks the end of the\n"
		 << "   list similar to the \\0 character in a string\n"
		 << "eg: $ Int_t dets[]={2,3,5,7,8,0}; // selects the current LHS detectors only\n\n";
		 
	cout << "flagevents = 0         - Do not flag coincident events\n"
		 << "           = 1         - Do flag coincident events\n"
		 << "           = 2, 3, etc - Flag only events were there 2, 3 etc coincidences\n\n";
		 
	cout << "detn       = 0,        - Coincidences between all detectors are used\n"
		 << "           = det#      - Coincidences between det# (replace with number) and all other\n"
		 << "                         detectors is used\n\n";
		 
	cout << "plot_layer = 0         - Search for coincidences between all layers of detectors\n"
		 << "                         (2D graph is not drawn) \n"
		 << "           = layer#    - Search for coincidences only within specified layer \n"
		 << "                         A 2D colour plot of all the coincidence counts is also plotted\n"
		 << "                         where the detectors are numbered as they are in Gran Sasso\n"
		 << "                         and the detector between coincidences are search, is marked\n\n";
		 
	cout << "Emin       = 0         - Do not restrict lower energy range\n"
		 << "           = # (>0)    - Restrict search above specified energy (keV)\n"
		 << "                         NB The lower thresholds are always restricted by rs.lo_thresh\n\n";
		 
	cout << "Emax       = 0         - Do not restrict uppper energy rangea\n"
		 << "           = # (>0)    - restrict search below specfied energy (keV)\n\n";
		 
	cout << "Do not use the remaining options for now\n\n";
	
	cout << "So what do the plots mean? I'll update that section soon!\n";
}

void SScut_plot_coincedences(Int_t *dets, Int_t flagevents=0, Char_t use_flags[3] = "000", Int_t detn=0, Int_t plot_layer=0, Double_t Emin=0, Double_t Emax=0, Int_t useDetnDefaultThresh = 0, Double_t cutsecsbefore=0., Double_t cutsecsafter=0.){
	//This function flags events that are coincident with another event in the list of searched detectors "dets"
	// NB if flag events > 1 then only coincidences of greater than flagevents are flagged.
	// if detn==0 then coincidences between all 'dets' are used. Else coincidences between detn is only used.
	// plots are only made when 'flagevents' is not zero.
	// plot_layer: set to 0 to plot all layers (2D graph is not drawn) otherwise set the the layer that should be plotted.
	// if dets=0 then using defaults:{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,0};, to use 'dets' create an array of detectors to be used
	// set useDetnDefaultThresh: specify if lo_thresh (0) of the run or the trigger threshold (1) is used
	

	int i, j, k, cnt;
	char buffer[50];
	TDiamond *diamond=0;
	Double_t cur_time;
	vector <Int_t> coincidence_counts(NDETS,0);
	vector <Int_t> coincidence_counts_total(NDETS,0);
// 	for(j=0;j<NDETS;j++) newcuts[j].clear();
	Int_t default_dets[]={ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
						 /*  17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
	                      33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
						   49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,*/  0};
	if(!dets) dets = default_dets;

	int plot_energies=0;//!flagevents;
	int count=3;
	if(plot_energies) gblresults.clear();
	Double_t depossitionE[NDETS];
	for(i=0;i<NDETS;++i) depossitionE[i]=0;
// 	cout << "hi1\n";
	//Int_t *Nhits = new Int_t[NDETS];
	//for(i=0;i<NDETS;i++) Nhits[i] = 0;
	h_handle = new TH1D("Coincidences","Coincidences",NDETS-1,1.,NDETS);
	if(Emin>0) sprintf(buffer,"Number depositions between %1.0lf and %1.0lfkeV",Emin,Emax);
	else sprintf(buffer,"Number depositions between");
	h_handle->GetXaxis()->SetTitle(buffer);
	h2_handle = new TH2D("CoinArray","Coincidence Array",4,0.5,4.5,4,0.5,4.5);

Double_t lower_thresh = 0;

int cnt2=0, cnt3=0, cnt4=0;
	for(cnt=0,i=0,evtree->GetEvent(0);i<nentries;++i,evtree->GetEvent(i)){
		for(j=0;dets[j];++j) {
			if(useDetnDefaultThresh) lower_thresh = SSch2e(SSth2ch.Eval(rs[ev.runno].thr[dets[j]-1]),dets[j],ev.groupno+1);
			else lower_thresh = rs[ev.runno].lo_thresh[dets[j]-1];
			
			if( 
				((ev.energy[dets[j]-1] > Emin && Emax ==0) || (ev.energy[dets[j]-1] > Emin && ev.energy[dets[j]-1] < Emax)) &&
				
				// check if within allowed thresholds
				(ev.energy[dets[j]-1] > lower_thresh && (Emax == 0 || (Emax < gs[ev.groupno].hi_thresh[dets[j]-1]))) &&

				rs[ev.runno].det_on[dets[j]-1] &&
				(use_flags[0]=='_' || !rs[ev.runno].flag[dets[j]-1]) &&
				gs[ev.groupno].det_on[dets[j]-1] &&
				!gs[ev.groupno].flag &&
					
				// old 12/03/2008 &&(use_flags==0 || (use_flags==1 && !flag[i] && !rs[ev.runno].flag[dets[j]-1])) ) {
				(use_flags[0]=='_' || ((use_flags[0]=='0') && !flag[i]) || ((use_flags[0]=='1') && flag[i]) || (use_flags[0]==flag[i])) &&
				(use_flags[2]=='_' || !llfn_Decimal2Bit(BFlag[i],detn,NDETS)==(use_flags[2]=='0'))
			) {
// 				cout << cnt << endl;
				++cnt; //check for coincidences
				++coincidence_counts[dets[j]-1];
				if(plot_energies) depossitionE[dets[j]-1]=ev.energy[dets[j]-1];
			}
		}
		
		if(cnt>1) ++cnt2; // is cnt2 not used?
		
		if(detn!=0 && coincidence_counts[detn-1]==0){ // check that detn fired, reject if it did not.
			cnt=0;
			for(k=0;k<NDETS;k++) coincidence_counts[k]=0;
		}
		else for(k=0;k<NDETS;k++) coincidence_counts_total[k]+=coincidence_counts[k];
	
		h_handle->Fill(cnt);
		if(plot_energies) {
			if(cnt==count) {
				gblresults.push_back(i);
				for(j=0;j<NDETS;++j) {
					if(depossitionE[j]>0) {
						gblresults.push_back(j+1);
						gblresults.push_back(depossitionE[j]);
					}
				}
				gblresults.push_back(-1);
			}
			for(j=0;j<NDETS;++j) depossitionE[j] = 0;
		}
		
		if(flagevents > 0 && cnt>flagevents) {
			cur_time = ev.time;
			if(cur_time-cutsecsbefore <ev.time ) { // go back "cutsecsbefore" seconds
				for(;(cur_time-cutsecsbefore)<ev.time && i>=0;i--,evtree->GetEvent(i)); 
				i++;
				evtree->GetEvent(i);
			}
			++cnt3;
			for(;(cur_time+cutsecsafter)>=ev.time && i<nentries;i++,evtree->GetEvent(i)) {flag[i]=1, ++cnt4;};
			--i;
		}
		//delete counters
		cnt=0;
		for(k=0;k<NDETS;k++) coincidence_counts[k]=0;
	}
// 	cout << "hi5\n";
	if(flagevents==0) {
		c1 = new TCanvas("c1","c1",10,10,600,400);
		c1->SetLogy();
		c1->cd();
		h_handle->Draw();
		int k_start=(plot_layer-1)*16;
		if(plot_layer>0) {
			TCanvas *c2 = new TCanvas("c2","c2",10,10,600,400);
			c2->cd();
			for(i=4,k=k_start;i>0 && k<(k_start+16);--i) {
				for(j=1;j<5 && k<(k_start+16);++j) {
					if(k==detn-1) { 
						h2_handle->Fill(i,j,0);
						if(int((detn-1)/16)+1==plot_layer) diamond = new TDiamond(i-0.25,j-0.25,i+0.25,j+0.25);
					}
					else h2_handle->Fill(i,j,coincidence_counts_total[k]);
					
					k++;	
				}
			}
			
			h2_handle->Draw("colz");

			if(int((detn-1)/16)+1==plot_layer) {
				diamond->SetFillColor(1);
				diamond->SetTextAlign(22);
				diamond->SetTextSize(0);
				diamond->SetBorderSize(0);			
				diamond->Draw();
			}
			
			
		}
	}
// 	cout << "hi6\n";
	if(plot_energies) {
		TCanvas *c3 = new TCanvas("c3","c3",10,10,600,400);c3->cd();
		TH2D *h2_handle3= new TH2D("h3"," ",4,0.5,4.5,4,0.5,4.5);
		for(i=1;i<int(gblresults.size());++i) {
			while(i< int(gblresults.size()) && gblresults[i]!=-1) {
				h2_handle3->Fill( -1*(int(gblresults[i])-16)/4+1, (int(gblresults[i])%4!=0)*int(gblresults[i])%4+(int(gblresults[i])%4==0)*4, gblresults[i+1]);
				cout << "det: " << gblresults[i] << "  filling " << -1*(int(gblresults[i])-16)/4+1 << " " << (int(gblresults[i])%4!=0)*int(gblresults[i])%4+(int(gblresults[i])%4==0)*4 << " " << gblresults[i+1] << endl;
				i+=2;
			}
			cout << endl;
			++i;
			h2_handle3->Draw("colz");
// 			c3->Update();
//			gPad->WaitPrimitive();
			sprintf(buffer,"coin%d.gif",i);
			c3->Print(buffer);
			h2_handle3->Reset();
			
		}
	}
				
		
	
		
// 	llfn_addcuts(); 
	cout << cnt2 << " " << cnt3 << " " << cnt4 << endl;
	cout << "NB in the event of a call for a flag, all consequative events with the same system time are flagged.\n";
	cout << "To remove this feature the hi_time must be compared. This is essential for shorter time cuts anyway\n";
// 	cout << ".........\nNB this function does not work well across groups. (does not do range safety checks)\n";
	cout << "\n\nUsed time = " << SSget_coincident_usedtime(dets) << " hours\n\n";

	
}


// because I don't trust SScut_plot_coincidences I'm writing a new simple function...
// the conclusion was that it does not work either. Basically this remove all the events from about 2 groups
void SScut_coincidences() {
// this function flags all data that is coincident with other events above the trigger threshold
// the only problem in doing this is that it is possible that a noise event that has been cut out with the cleaning is included in this
// search and could be coincident with a real event, which is then rejected.
	int i, j, cnt;
	for(i=0,evtree->GetEvent(0);i<nentries;++i,evtree->GetEvent(i)){
		cnt = 0;
		for(j=0;j<NDETS;++j) {
			if( rs[ev.runno].det_on[j] &&
// 				gs[ev.groupno].det_on[j] &&
// 				!gs[ev.groupno].flag &&
				(ev.energy[j] > SSch2e(SSth2ch.Eval(rs[ev.runno].thr[j]),j+1,ev.groupno+1))
			) ++cnt;
		}
		if(cnt>1) flag[i]=1;
	}
}


// NB itoa can be used to convert an integer to binay string
//#include <string.h>

int llfn_Binary2Decimal(char *bin){

	int dec = 0, fac = 1;
	char *s = bin + strlen(bin) -1;
	
	while (s >= bin) {
		if (*s == '1') {
			dec += fac;
		}
		fac *= 2;
		s--;
	}
	
	return dec;
}


void SSflag_greatest(Int_t plot_triggered=0){
	int i, j;
	char *expandedBFlag = new char[NDETS+3];
	for(i=0;i<NDETS;i++) expandedBFlag[i] = '1'; expandedBFlag[NDETS] = '\0';
	Double_t Egreatest;
	Int_t DetNo = 0;
	
	h_handle = new TH1D("Triggered","Triggered",NDETS,0.,NDETS);
	
	for(i=0,evtree->GetEvent(0);i<nentries;i++,evtree->GetEvent(i)){
		Egreatest = 0.;
		for(j=0;j<NDETS;j++){ 
			if(ev.energy[j] > Egreatest){ 
				Egreatest = ev.energy[j];
				DetNo = j;
			}
		}
		if(plot_triggered) h_handle->Fill(DetNo);
		expandedBFlag[DetNo]='0';
		BFlag[i] = llfn_Binary2Decimal(expandedBFlag);
		expandedBFlag[DetNo]='1';
// 		cout << DetNo << " " << BFlag[i] << " " << llfn_Decimal2Bit(BFlag[i],detn,NDETS)<< endl;
	}
	
	if(plot_triggered) h_handle->Draw();
	delete expandedBFlag;
}
			




// TH1F *th1_thresholds;

void SSlocate_Thresholds(){
	int i, j, k, largest_thresh=0, first_run=0, numberOfRuns=0;
	for(i=0;i<ngroups;i++){
		numberOfRuns += gs[rs[first_run].groupno].nruns;
// 		cout << "group " << i << endl;
		for(j=0;j<NDETS;j++){
			//locate largest threshold used in group
// 			cout << "locating largest threshold of det " << j << endl;
			for(k=first_run,largest_thresh=0; k<numberOfRuns; ++k) if(rs[k].thr[j]>rs[largest_thresh].thr[j]) largest_thresh=k;
			//set thresholds		
// 			cout << "setting threshold\n";
// cout << "det = " << j+1 << " largest = " << largest_thresh << " ch = " << SSth2ch.Eval(rs[largest_thresh].thr[j]) << endl;
			gs[i].lo_thresh[j] = gs[i].cal_eqn[0][j]*(SSth2ch.Eval(rs[largest_thresh].thr[j])) + gs[i].cal_eqn[2][j];
		}
// 		cout << "setting first_run to " << first_run << endl;
		first_run += gs[i].nruns;
	}
	
		
}

void SSplot_settings_thr(char help='H'){
	help = help;
	cout << "SSplot_settings_thr plots the threshold settings as a function of time\n";
	cout << "To select plotting formate enter mode number:\n";
	cout << "0 ..... threshold numbers are plotted\n";
	cout << "1 ..... thresholds in channels are plotted\n";
	cout << "2 ..... thresholds in energy are plotted\n";
}
	
void SSplot_settings_thr(Int_t detn, int mode =0, int drawPlot=1) {
	Double_t *run_time = new Double_t[nruns];
	Double_t *run_thr = new Double_t[nruns];
	for(int run_counter=0;run_counter<nruns;++run_counter) {
		run_thr[run_counter] = rs[run_counter].thr[detn-1];
		if(mode==1) run_thr[run_counter] = SSth2ch.Eval(rs[run_counter].thr[detn-1]);
		if(mode==2) run_thr[run_counter] = gs[rs[run_counter].groupno].cal_eqn[0][detn-1]*(SSth2ch.Eval(rs[run_counter].thr[detn-1])) + gs[rs[run_counter].groupno].cal_eqn[2][detn-1];
		run_time[run_counter] = rs[run_counter].starttime;/*cout << run_thr[run_counter] << " " << run_time[run_counter] << endl;*/
	}
	g_handle = new TGraph(nruns,run_time,run_thr);
// 	g_handle->SetStats(0);
	if(drawPlot) {
		g_handle->GetXaxis()->SetTimeDisplay(1);
		g_handle->GetXaxis()->SetTimeFormat("");
		g_handle->GetXaxis()->SetNdivisions(-510);

		g_handle->GetXaxis()->SetTitle("Date");
		g_handle->GetYaxis()->SetTitle("Threshold Number");
		if(mode==1) g_handle->GetYaxis()->SetTitle("Threshold (Channels)");
		if(mode==2) g_handle->GetYaxis()->SetTitle("Threshold (keV)");
	
		g_handle->Draw("APl");
	}
}


void SSplot_settings_userthr(Int_t detn, Char_t plot_opt[10]="same") {
	Double_t time_scale=gblperiod;
	
	if(time_scale==0) h_handle = new TH1D("thresholds","run thresholds",(rs[nruns-1].starttime-rs[0].starttime)/3600,rs[0].starttime,rs[nruns-1].starttime);
	else h_handle = new TH1D("thresholds","run thresholds",(rs[nruns-1].starttime-rs[0].starttime)/3600,0,(rs[nruns-1].starttime-rs[0].starttime)/3600*3600./time_scale);
	
	for(int run_counter=0;run_counter<nruns;++run_counter) {
		if(time_scale==0) h_handle->Fill(rs[run_counter].starttime,rs[run_counter].lo_thresh[detn-1]);
		else h_handle->Fill(Double_t(rs[run_counter].starttime-rs[0].starttime)/time_scale,rs[run_counter].lo_thresh[detn-1]);
	}

	if(time_scale==0) {
		h_handle->GetXaxis()->SetTimeDisplay(1);
		h_handle->GetXaxis()->SetTimeFormat("");
		h_handle->GetXaxis()->SetNdivisions(-510);
		h_handle->GetXaxis()->SetTitle("Date");
	}
	else h_handle->GetXaxis()->SetTitle("Time");
		
	h_handle->GetYaxis()->SetTitle("Threshold (keV)");
	h_handle->SetLineColor(kBlue);
	h_handle->Draw(plot_opt);
		

}

void SSplot_settings_hv(Int_t detn, int drawPlot=1) {
	Double_t *run_time = new Double_t[nruns];
	Double_t *setting = new Double_t[nruns];
	for(int run_counter=0;run_counter<nruns;++run_counter) {
		setting[run_counter] = rs[run_counter].hv[detn-1];
		run_time[run_counter] = rs[run_counter].starttime;/*cout << run_thr[run_counter] << " " << run_time[run_counter] << endl;*/
	}
	g_handle = new TGraph(nruns,run_time,setting);
// 	g_handle->SetStats(0);
	c1 = new TCanvas("settings", "hv settings",14,30,700,500);
	if(drawPlot) {
		g_handle->GetXaxis()->SetTimeDisplay(1);
		g_handle->GetXaxis()->SetTimeFormat("");
		g_handle->GetXaxis()->SetNdivisions(-510);

		g_handle->GetXaxis()->SetTitle("Date");
		g_handle->GetYaxis()->SetTitle("Voltage (V)");
	
		g_handle->Draw("APl");
	}
}

void SSlocate_NoiseThresholds(int detn, double tolerance = 0.1){
	int i, j;
	Int_t cur_group=-1, cur_run=-1, lo_thresh=0;
	int nbins=0, last_nruns=0;
	int *fraction_ok_regions=NULL;
	
	// function of typical approximate spectrum Normalised to Counts/10keV/hr
	TF1 SSGausFunc("gausfunc","1.11*TMath::Gaus(x,110,98)",0,500);
	
	TH1D * spec=NULL;
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
// 		SKIPFLAGS();  //skip flagged groups and runs (defined above)
// cout << "h1\n";
		if(ev.groupno!=cur_group){
			// test and set thresholds
			if(fraction_ok_regions) {
				for(j=nbins-1;j>=0;j--) {
					if(double(fraction_ok_regions[j])/last_nruns < tolerance && gs[cur_group].lo_thresh[detn-1] < spec->GetBinCenter(j+1) + 5) {
						gs[cur_group].lo_thresh[detn-1] = spec->GetBinCenter(j+1) + 5;
						break;
					}
				}
			}
			for(j=0;j<nbins;j++) printf("%3.2lf ",double(fraction_ok_regions[j])/last_nruns); cout << endl;
			for(j=0;j<nbins;j++) printf("%3.0lf ",spec->GetBinCenter(j+1)); cout << endl;
// 		cout << "h2\n";
			cur_group=ev.groupno;
			if(int(gs[ev.groupno].lo_thresh[detn-1])%10) lo_thresh=Int_t(gs[ev.groupno].lo_thresh[detn-1])+10;
			else lo_thresh=Int_t(gs[ev.groupno].lo_thresh[detn-1]);
			if(lo_thresh<=400) nbins = (400-lo_thresh)/10	;
			else nbins = 1;
			last_nruns = gs[ev.groupno].nruns;
			
			if(fraction_ok_regions) delete fraction_ok_regions;
// 			cout << nbins << endl;
			fraction_ok_regions = new int[nbins];
// 			cout << "h6\n";
			for(j=0;j<nbins;j++) fraction_ok_regions[j] = 0;
		}
		if(ev.runno!=cur_run) {
// 		cout << "h3\n";
			cur_run = ev.runno;
			// count ok regions if counts are less than 4 sigmas away from standard spectrum 
			if(spec) for(j=0;j<nbins;j++) if(spec->GetBinContent(j+1) < SSGausFunc.Eval(spec->GetBinCenter(j+1)) + 4.*sqrt(SSGausFunc.Eval(spec->GetBinCenter(j+1)))) fraction_ok_regions[j]++;
// 			if(spec) cout << spec->GetBinContent(3)  << " " << GausFunc.Eval(spec->GetBinCenter(3)) + 4.*sqrt(GausFunc.Eval(spec->GetBinCenter(3))) << endl;
// cout << "h4\n";
			// create new test spectrum			
			if(spec) spec->Delete();
// 			cout << "h5\n";
			spec = new TH1D("spec","spec",nbins,lo_thresh,400);

			
		}
		spec->Fill(ev.energy[detn-1]);

	
		
		
	}
	
	if(fraction_ok_regions) 
				for(j=nbins-1;j>=0;j--) 
					if(double(fraction_ok_regions[j])/last_nruns < tolerance) 
						if(gs[cur_group].lo_thresh[detn-1] > spec->GetBinCenter(j+1) + 10) gs[cur_group].lo_thresh[detn-1] = spec->GetBinCenter(j+1) + 5;
				for(j=0;j<nbins;j++) printf("%1.2lf ",double(fraction_ok_regions[j])/last_nruns); cout << endl;
				for(j=0;j<nbins;j++) printf("%3.0lf ",spec->GetBinCenter(j+1)); cout << endl;	
}

void SSplot_spec2(char help='H'){
	help = help;
	cout << "Use plot_spec2 to plot the spectrum of detector \"detn\" using the lower threshold rs.lo_thresh\n";
	cout << "  where the life time is calculated for each keV bin. The spectrum is plotted with h_handle. \n";
	cout << "  eg you can use h_handle to redraw (h_handle->Draw()). The lifetime is saved to h_handle2.\n\n";
	cout << "rebin: define binning. NB histogramed is rescaled so that it is still in bins/keV\n";
	cout << "use_flags:  flag options: _ = ignore, 0 = do not use flagged data, 1 = only use flagged data\n";
	cout << "		     1st number is ev flags. That is the flag fag for each event (blind to detector number)\n";
	cout << "			 2nd number is run flags. An hour of data is rejected if flagged.\n";
	cout << "			 3rd number is Binary ev flag. With this flag the specific detector numbers can be flagged\n";
	cout << "sumw2:   set this to draw spectrum with errors\n";
	cout << "rawdata: set this if you want the scale in Counts/keV. NB rebinned data is also in /keV.\n";
	cout << "reset_used_groups: set this to reset groups sellection with  \"groups_set();\"\n";
}		
void SSplot_spec2(Int_t detn, Int_t rebin = 1, Char_t use_flags[3] = "000", Int_t use_sumw2=0, Int_t rawdata=0, Double_t undercut=10.){ 
// This function plots a spectrum like plot_spec() but uses thresholds from each run and life times per keV bin to 
// generate the spectrum. This way the full data can be used.
// undercut is the amount bellow the upper threshold that the data is used and avoids inclusion of saturation events.
	Char_t buffer[512];
	
	if(gbl_draw) {
		c1 = new TCanvas("cspec", "cspec",14,30,700,500);
		c1->SetLogy();
	}

	int i, j, k, j_start;
	int specUpperRange = 12000;
	
	// initialise plot ranges	
	sprintf(buffer,"D%d",detn);
	h_handle = new TH1D("hspec",buffer,specUpperRange,0.,specUpperRange);
	sprintf(buffer,"D%d time",detn);
	h_handle2 = new TH1D("spec_time",buffer,specUpperRange,0.,specUpperRange);
	
	i=0;
	evtree->GetEvent(i);
	
	for(i=0,j=0,j_start=0,evtree->GetEvent(j);i<nruns;++i) { //loop over runs
		
		// if run or group is flagged, skip to next run (set j to i+=event in run
		// skip run if...
		if( !gs[rs[i].groupno].det_on[detn-1] || 														// skip off detectors for whole group
			gs[rs[i].groupno].flag || 																// skip flagged groups
			!rs[i].det_on[detn-1] || 																			// skip off detectors
			(use_flags[1]!='_' && rs[i].flag[detn-1]==(use_flags[1]=='0'))								// skip flagged runs
			) {	
			
			j = j_start+rs[i].nevents;
			j_start = j;
			// maybe check that it is still in sync (ie that the previous event if of the previous run
			// idealy it would go like this: j = rs[i+1].evnumber which would take it to the first event of the i+1'th run
			continue;
		}
		
		for(evtree->GetEvent(j);j<j_start+rs[i].nevents;++j,evtree->GetEvent(j)) { 						// loop over events in run  (could have been j<rs[i+1].evnumber)
			// check if event is flagged
			if( int(rs[i].lo_thresh[detn-1]) < ev.energy[detn-1] &&										// Fill if above lower run threshold
				(gs[rs[i].groupno].hi_thresh[detn-1]-undercut > ev.energy[detn-1]) &&					// Fill if below upper group threshold (so far there is no upper run threshold as it is unlikely to be used)
				(use_flags[0]=='_' || !flag[j]==(use_flags[0]=='0')) &&									// Fill if using evflags and is not evflagged
				(use_flags[2]=='_' || !llfn_Decimal2Bit(BFlag[j],detn,NDETS)==(use_flags[2]=='0')) ) {  // Fill if using evBitFlagged and is not evBitFlagged (evBitFlag is a binary number which is decomposed into a flag for each detector)
			
				// Add event to spectrum
				h_handle->Fill(ev.energy[detn-1]+(gRandom->Uniform()-0.5));
			}
		}
		// add lifetime in seconds to runningtime histogram
		for(k=int(rs[i].lo_thresh[detn-1]);Double_t(k)<gs[rs[i].groupno].hi_thresh[detn-1]-undercut && k<specUpperRange;++k)
			h_handle2->Fill(k,3600.-(rs[i].neventsall+rs[i].nvetoevents)*gbl_deadtimeperevent); // 
		
		j = j_start+rs[i].nevents;
		j_start = j;
		evtree->GetEvent(j);
		//if(i+1!=ev.runno) cout << "Warning, out of sync!\n";
	}
	
	//export?
	
	// Add error bars are required
	if(use_sumw2) {
		h_handle->Sumw2();
		h_handle2->Sumw2();
	}
	// normalise / rebin
	// before 28/02/2008 a detector mass of 5.78g was assumed so any plots must be rescaled accordingly.
	// an average mass of 6.5g is a better approximation. although idealy the precise mass should be used.
	if(!rawdata) h_handle->Divide(h_handle,h_handle2,3600.*24./6.5e-3);
	//rebin 
	h_handle->Rebin(rebin); // I check, these functions work as predicted. Only point is that the error on 1 and and 0 is sqrt(1).
	h_handle->Scale(1./rebin); 
	// label and plot if required
	if(gs[0].cal_eqn[0][detn-1]!=1.) h_handle->GetXaxis()->SetTitle("Energy (keV)");
	else h_handle->GetXaxis()->SetTitle("Channels");
	if(rawdata && gs[0].cal_eqn[0][detn-1]!=1.) h_handle->GetYaxis()->SetTitle("Counts/keV");
	else if(rawdata && gs[0].cal_eqn[0][detn-1]==1.) h_handle->GetYaxis()->SetTitle("Counts/Channel");
	else if(!rawdata && gs[0].cal_eqn[0][detn-1]!=1.) h_handle->GetYaxis()->SetTitle("Counts/keV/kg/day");
	else if(!rawdata && gs[0].cal_eqn[0][detn-1]==1.) h_handle->GetYaxis()->SetTitle("Counts/Channel/kg/day");
	if(gbl_draw) h_handle->Draw();

}

void SSplot_spec(char help='H'){


	help = help;
	cout << "Use plot_spec to plot the spectra of detector number \"detn\"\n\n";
	cout << "rebin: define binning. NB histogramed is rescaled so that it is still in bins/keV\n";
	cout << "Emin:  Lower energy on x-axis for spectrum. Group used only if between gs[i].hi_thresh[detn-1] and gs[i].lo_thresh[detn-1]\n";
	cout << "Emax:  Higher energy boundery. Set to zero to use maximum but note that different groups mays have different end points.\n";
	cout << "use_flags:  flag options: _ = ignore, 0 = do not use flagged data, 1 = only use flagged data\n";
	cout << "		     1st number is ev flags. That is the flag fag for each event (blind to detector number)\n";
	cout << "			 2nd number is run flags. An hour of data is rejected if flagged.\n";
	cout << "			 3rd number is Binary ev flag. With this flag the specific detector numbers can be flagged\n";
	cout << "rawdata: set this if you want the scale in Counts/keV. NB rebin data is also in /keV.\n";
	cout << "sumw2:   set this to draw spectrum with errors\n";
	cout << "close_file: set this to save spectrum to spec.C, instead of plotting\n";
	cout << "reset_used_groups: set this to reset groups sellection with  \"groups_set();\"\n";
}

void SSplot_spec(Int_t detn, Int_t rebin = 1, Double_t Emin = 120, Double_t Emax = 10000, Char_t use_flags[3] = "000", Int_t rawdata=0, Int_t use_sumw2=0, Int_t closefile=0, Int_t reset_used_groups =1){
	Char_t buffer[512];
	//TROOT simple("simple","nTuple");
	gROOT->SetStyle("Plain");
	//sprintf(fname,"%s.root", argv[1]);
// 	TFile *hfile = new TFile("hist.root","RECREATE","hist",9);	
	TFile *hfile = new TFile("hist.root","UPDATE","hist",9);
	
	if(gbl_draw) {
		c1 = new TCanvas("spec", "spec",14,30,700,500);
		c1->SetLogy();
	}
// 	llfn_makewhite(spec);
	int i, j;
	
	int *used_group = new int[ngroups]; 
 	for(i=0;i<ngroups;i++) used_group[i] = 0;
	Double_t used_time =0;
	
	sprintf(buffer,"D%d",detn);
	TH1D *h_spec = new TH1D("spec",buffer,Int_t(Emax-Emin),Emin,Emax);
	
	if(reset_used_groups) SSgroups_set(); 
	
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		//for(j=1;j<=24;j++) cout << llfn_Decimal2Bit(BFlag[i],j,24); cout  << endl;
		SKIPFLAGS();  //skip flagged groups and runs (defined above)
		if((/*gs[ev.groupno].lo_thresh[detn-1] <= Emin && */gs[ev.groupno].hi_thresh[detn-1] >= Emax) // don't include groups that have smaller ranged (this is now redundent now that rs.lo_thresh has been introduced)
			&& (use_flags[0]=='_' || !flag[i]==(use_flags[0]=='0'))
			&& (use_flags[2]=='_' || !llfn_Decimal2Bit(BFlag[i],detn,NDETS)==(use_flags[2]=='0')) 
			&& (Emin < ev.energy[detn-1])
			/// NB I just added the following line and that the life time is calcalated with the new function...
			&& (Emin > rs[ev.runno].lo_thresh[detn-1])){ 
				h_spec->Fill(ev.energy[detn-1]+(gRandom->Uniform()-0.5));
				used_group[ev.groupno] = 1;  //could use a more efficient way then this as 
		}
	}
	//flag skipped groups
	for(i=0;i<ngroups;i++) if(!used_group[i]) gs[i].flag = 1;

	if(use_sumw2) h_spec->Sumw2();
	h_spec->Rebin(rebin);
	h_spec->Scale(1./rebin);

	//Get Running Time.....
	used_time = SSgetusedtimeR(detn,Emin,Emax);
	
	if(use_flags[1]!='0') cout << "WARNING: don't trust the running time\n";
	// Scale to g/keV/day
	if(!rawdata) {
// 		h_spec->Scale(3600.*24./5.78e-3/used_time);
// 		h_spec->SetBinContent(0,5.78e-3*used_time/3600./24.);  // changed 28/02/2008
		h_spec->Scale(3600.*24./6.5e-3/used_time);
		h_spec->SetBinContent(0,6.5e-3*used_time/3600./24.);
	}
	else h_spec->SetBinContent(0,used_time);
	//h_spec->SetEntries(used_time);
// 	sprintf(buffer,"D%d %1.2f hours. (Deadtime of %1.2f hours subtracted)",detn,used_time/3600.,getdeadtime(detn)/3600);
	sprintf(buffer,"D%d %1.2f hours.",detn,used_time/3600.);
	h_spec->SetTitle(buffer);

	if(gs[0].cal_eqn[0][detn-1]!=1.) h_spec->GetXaxis()->SetTitle("Energy (keV)");
	else h_spec->GetXaxis()->SetTitle("Channels");
	if(rawdata && gs[0].cal_eqn[0][detn-1]!=1.) h_spec->GetYaxis()->SetTitle("Counts/keV");
	else if(rawdata && gs[0].cal_eqn[0][detn-1]==1.) h_spec->GetYaxis()->SetTitle("Counts/Channel");
	else if(!rawdata && gs[0].cal_eqn[0][detn-1]!=1.) h_spec->GetYaxis()->SetTitle("Counts/keV/kg/day");
	else if(!rawdata && gs[0].cal_eqn[0][detn-1]==1.) h_spec->GetYaxis()->SetTitle("Counts/Channel/kg/day");
	if(gbl_draw) h_spec->Draw();
	h_handle = h_spec;
	if(gblpause) gPad->WaitPrimitive();

	hfile->Write();
 	if(closefile) hfile->Close();
	if(gbl_verbose){
		// print out which groups were used and total timed used. 
		cout << "The following groups were used for this plot:\n";
		j = SSgroups_get();
		cout << j << "/" << ngroups << " groups used.  Total " << used_time/3600 << " hours of data.\n\n";
		cout << "includes deadtime = " << SSgetdeadtime(detn)/3600 << " hours. (" << SSgetdeadtime(detn)/used_time*100 << "%)\n";	
	}
}

// use detn = 0 to analyse all detectors
void SSplot_epp(Int_t detn, Double_t period = 1E-3, Int_t flag_finds = 0, Double_t Emin=150, Double_t Emax=0, Double_t range_Xmax = 100){
// if(time_period>=3600) {
// 	cout << "This function is only for period searches

// plot events per given time period. If trig_filter = 0 then applies to all detectors. 
// Energy range only works for detn != 0 currantly. This could be changed to "if any detector detectors energy > 0" with a fors loop and break statements.
//!!!!  WARNING: if search periods of greater than 1 hour then the periods will not be complete before a group (or set) is complete! ie, every 


	c1 = new TCanvas("c1", "c1",14,30,700,500);
	c1->SetLogy();
	int i, j, k, l;
	int eventCounter;
	Double_t period_start, period_start_hitime, skip_steps;//, hi_time_start;
	Char_t buffer[100];
	
	// save the event number of multiple events so that they can be flagged
	vector <int> flagged_event_log;
	
	// ensure that 3600 is devisable by period
	if(period>3600.) period = 3600.;
	period = 3600./int(3600./period);
	
	sprintf(buffer,"Events per %f seconds. D%d",period,detn);
	h_handle = new TH1D("Events per hour",buffer,Int_t(range_Xmax),0,range_Xmax);
	h_handle->GetXaxis()->SetTitle("number of events");
	h_handle->GetYaxis()->SetTitle("frequecy");
	
	
	for(i=0;i<ngroups;++i) {  // loop over groups
		// skip group if flagged
		if( !gs[i].det_on[detn-1] || gs[i].flag /*|| Emax > gs[i].hi_thresh[detn-1] */) continue;
		for(j=gs[i].firstrunno;j<nruns && j<gs[i+1].firstrunno;++j) { // loop over runs
			// skip runs if flagged
			if(	!rs[j].det_on[detn-1] || rs[j].flag[detn-1] /*|| Emin < rs[j].lo_thresh[detn-1]*//*SSch2e(SSth2ch.Eval(rs[j].thr[detn-1]),detn)*/ ) continue;
			// initialise everything... counting is not done accross runs with this tool
			eventCounter = 0;
			flagged_event_log.clear();
			period_start = rs[j].starttime;
			evtree->GetEvent(rs[j].firstevno);
			
			period_start_hitime = ev.hi_tim + (rs[j].starttime - ev.time)/8e-6; // set to the value the hi_time would have had at start of run
			
			for(k=rs[j].firstevno;k<nentries && k<rs[j+1].firstevno;++k,evtree->GetEvent(k)) { // loop over event in run
				if( flag[k] || 
					llfn_Decimal2Bit(BFlag[k],detn,NDETS) ||
					ev.energy[detn-1] < Emin ||
					ev.energy[detn-1] > Emax ||
					ev.energy[detn-1] < rs[j].lo_thresh[detn-1] ||
					ev.energy[detn-1] > gs[i].hi_thresh[detn-1] ||
					ev.energy[detn-1] < SSth2ch.Eval(rs[j].thr[detn-1])
					) continue;
					
				if(period >= 1) { // use unix time
					if(ev.time - period_start > period) {  
						h_handle->Fill(eventCounter);
						if(flag_finds && eventCounter>flag_finds) for(l=0;l<int(flagged_event_log.size());++l) flag[flagged_event_log[l]]='p';
						eventCounter = 0;
						flagged_event_log.clear();
						period_start += period;
						--k; // count this event in next period
					}
					else {
						++eventCounter;
						flagged_event_log.push_back(k);
					}
				}
				else { // time stamp
					if((ev.hi_tim - period_start_hitime)*8e-6 > period) {
						if(eventCounter>0) {
							h_handle->Fill(eventCounter);
							if(flag_finds && int(flagged_event_log.size())>flag_finds) for(l=0;l<int(flagged_event_log.size());++l) flag[flagged_event_log[l]]='p';
							flagged_event_log.clear();
							eventCounter = 0;
							period_start_hitime += period/8e-6;
						}
						else {  // it is important not to recheck un filled regions to save time...
							skip_steps = (ev.hi_tim - period_start_hitime)*8e-6/period;
							for(l=0;l<skip_steps;++l) h_handle->Fill(0);
							period_start_hitime += skip_steps*period/8e-6;
						}
						--k; // count this event in next period
					}
					else if((ev.hi_tim - period_start_hitime)*8e-6 < 0) { // check if counter has restarted
						period_start_hitime = period_start_hitime - pow(2.,32) -1; // hopefully this sets the timer to the new time
					}
					else {
						++eventCounter;
						flagged_event_log.push_back(k);
					}
				}
			}
		}
	}
	
	h_handle->Draw();

// check what happens if not many events and time_start is always before the last event. Do the the events run away from time_start?
}


void cut_highenergyevents(Int_t detn, Int_t Emin, Double_t safetime){
	int i;
	Double_t evtime;

	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		if(ev.energy[detn-1] > Emin){
			evtime=ev.time;
			//go back safetime seconds
			for(;i>=0 && TMath::Abs(ev.time-evtime)<safetime;i--,evtree->GetEvent(i));
			// flag data over this period
			for(;i<nentries && (ev.time-evtime)<safetime;i++,evtree->GetEvent(i)) flag[i]=1;
		}
	}
}

void SSplot_eph(Int_t detn, Int_t period = gblperiod, Char_t use_flags[3] = "000", Double_t range_Xmax = 100, Double_t Emin=150, Double_t Emax=0, Int_t plottimeline=0, Int_t autoselectgroups=1, Double_t samescale = 0, Int_t ploterrors=1, Int_t pause = 0, Int_t logYaxis = 0){
// plot events per given time period.
// detn = 0 is no longer supported
//!!!!  WARNING: if search periods of greater than 1 hour then the periods will not be complete before a group (or set) is complete ans so will not be plotted.
//               NB counts across skipped runs. Maybe this is not a problem and maybe even ok to counts accross runs?
//!!! WARNING fit_poisson does not work with normalised plot and ploterrors=0
	if(pause) gPad->WaitPrimitive();
// 	if(!c1) c1 = new TCanvas("c1", "c1",14,30,700,500);
	c1 = new TCanvas("c1", "c1",14,30,700,500);

		
	if(plottimeline==1){
		c1->Divide(1,2); 
		c1->cd(1);
	}		
	if(logYaxis) c1->SetLogy();
	gROOT->SetStyle("Plain");
	gStyle->SetOptFit(1111);
	//gStyle->SetOptStat(0);
	llfn_makewhite(c1); 
	int i=0, j=0, k=0, p_cnt = 0, p_cnt_buf=0;
	Int_t cur_group = 0, cur_run =0, run_cnt =0, enteredloop=0; 
	Int_t runno_buf=0;
	Double_t time_start;//, hi_time_start;
	Char_t buffer[100];
	gblperiod = period;
	
	
	sprintf(buffer,"D%d",detn);
	//if(h_handle) h_handle->Delete();
	h_handle = new TH1D("Events per hour",buffer,Int_t(range_Xmax),0,range_Xmax);
	sprintf(buffer,"Counts / %d hours(s) / (%d-%d)keV",period,int(Emin),int(Emax));
	h_handle->GetXaxis()->SetTitle(buffer);
	h_handle->GetYaxis()->SetTitle("frequecy (normalised)");
	
//-------------- code for the timeline plot ------------
	// estimate space needed and create the vectors to store the counts
	Int_t runningtime = 0;
	evtree->GetEvent(0); Double_t start_time = ev.time;
	for(i=0;i<ngroups;i++) runningtime+=nruns*3600;
	epp.epp.clear();
	epp.time.clear();
	epp.runno.clear();
	

//-----------------------------------------------------------
	

			
	time_start = rs[ev.runno].starttime; // this time should really be loaded from the start time which is save in the stats branch
	//hi_time_start = ev[i].hi_tim;
	for(i=0,j=0,evtree->GetEvent(i); i<nentries ;){
	
		//SKIPFLAGS();  //skip flagged groups and runs (defined above)
		//cout << "hi1\n";
		while( (gs[ev.groupno].flag || // do not plot if group is flagged
		(autoselectgroups && (/*gs[ev.groupno].lo_thresh[detn-1] > Emin ||*/ gs[ev.groupno].hi_thresh[detn-1] < Emax)) || // only plots if the data energy range falls within specified range
		!gs[ev.groupno].det_on[detn-1]) &&  
		i<nentries){
			if(gbl_verbose) cout << "skipping group " << ev.groupno+1 << "...\n";
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
			time_start = rs[ev.runno].starttime;
		}
		while( (	!rs[ev.runno].det_on || 														//Skip if det is off
					(use_flags[1]!='_' && rs[ev.runno].flag[detn-1]==(use_flags[1]=='0')) ||		//Skip if run is flagged
					(Emin < rs[ev.runno].lo_thresh[detn-1])					//Skip if runthreshold is too high
				) && i<nentries){	
			i += rs[ev.runno].nevents;
			evtree->GetEvent(i);
			time_start = rs[ev.runno].starttime;
		}
		if(ev.groupno > cur_group){ 
			cur_group = ev.groupno;
			p_cnt = 0;
			p_cnt_buf = 0;
			run_cnt = 0;  // do not count across groups
		}
		if(run_cnt==0){ 
			epp.time.push_back((ev.time - start_time)/3600/24); // timeline: record time of period start
			runno_buf = ev.runno;
		}
		for(cur_run=rs[ev.runno].starttime;rs[ev.runno].starttime==cur_run && i<nentries;i++,evtree->GetEvent(i)){ // loop over run
			enteredloop = 1;
			if((use_flags[0]=='_' || ((use_flags[0]=='0') && !flag[i]) || ((use_flags[0]=='1') && flag[i]) || int((use_flags[0])==flag[i]))
			&& (use_flags[2]=='_' || !llfn_Decimal2Bit(BFlag[i],detn,NDETS)==(use_flags[2]=='0')))
				if((ev.energy[detn-1] >= Emin && Emax == 0) || (ev.energy[detn-1] >= Emin && ev.energy[detn-1] <= Emax)) 
					p_cnt_buf++; // count the number of events in period
		}
		run_cnt++;
		p_cnt+=p_cnt_buf;
		p_cnt_buf=0;
		if(run_cnt==period && enteredloop){
			h_handle->Fill(p_cnt);	
			epp.epp.push_back(p_cnt);
			epp.runno.push_back(runno_buf+1);
			run_cnt=0;
			p_cnt=0;
			enteredloop = 0; // this prevents false zeros being entered.
			j++;
		}
	}
	if(ploterrors) h_handle->Sumw2();
	Double_t scale = 1./h_handle->Integral();
	if(samescale == 0) {
		h_handle->Scale(scale);
		samescalef = scale;
	}
	else {
		h_handle->Scale(samescalef);
		h_handle->SetLineColor(kRed);
		h_handle2=h_handle;
	}
	if(plottimeline!=-1) h_handle->Draw();
	
	

	
	
	
//-------------- code for the timeline plot ------------
	Double_t *yerrors = new Double_t[j];
	for(k=0;k<int(epp.epp.size());k++) yerrors[k]=1 + sqrt(epp.epp[k] + 0.75); // This appoximation to the mean is for Poisson distribution and obtained from STSDAS Help Pages "explain_errors" http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?explain_errors 
	if(g_handle2) g_handle2->Delete();
	TGraphErrors *g_handle2 = new TGraphErrors(epp.epp.size(),&epp.time[0],&epp.epp[0],0,yerrors);
	g_handle2->SetTitle("time line");
	g_handle2->GetXaxis()->SetTitle("Time (days)");
	//sprintf(buffer,"Counts/%d hour(s)",period);
	g_handle2->GetYaxis()->SetTitle(buffer);
	//TGraph *graph = new TGraph(k,time_stamp,period_cnt);
	if(plottimeline!=0){
		if(plottimeline!=-1) c1->cd(2);
		if(samescale){ 
			g_handle2->SetLineColor(kRed);
		}
		g_handle2->Draw("A*");
		g_handle=g_handle2;
	}
//-----------------------------------------------------------


c1->cd(1);
}

void SSdebug_nevents(){
	int i, pre, cur, nxt;
		for(i=0;i<nentries;){
			evtree->GetEvent(i);
			i += rs[ev.runno].nevents;
				evtree->GetEvent(i-1);
				pre=rs[ev.runno].starttime;
				//cout << "i -1 " << i-1 << endl;		

			evtree->GetEvent(i); 
			cur = rs[ev.runno].starttime;
			evtree->GetEvent(i+1); 
			nxt = rs[ev.runno].starttime;
			evtree->GetEvent(i);
			if(pre==cur){ 
				cout << i << " " << pre << " " << cur << " " << nxt << endl;
				//cin >> cur;
			}

		}
}

void SSplot_eppt(Int_t detn, Double_t Emin, Double_t Emax, int period=24, Char_t use_flags[3] = "000", Double_t axisScale=3600*24){
// plots counts per period against time
// this function differs from plot_pph in that does not borrow runs from the next hour but rescales the counts and errors accordingly if data is cut out

	c1 = new TCanvas("rates","rates",10,10,600,400);
	llfn_makewhite(c1);
	int i=0, j=0, k=0, p_cnt = 0, p_cnt_buf=0;
	Int_t cur_group = 0, cur_run =0, run_cnt =0, used_run_cnt = 0,enteredloop=0; 
	Int_t runno_buf=0;
	Double_t time_start;//, hi_time_start;
	Char_t buffer[100];
	gblperiod = period;

	epp.epp.clear();
	epp.time.clear();
	epp.runno.clear();
	epp.used_run_cnt.clear();
	
	evtree->GetEvent(0);
	Double_t start_time = ev.time;
	time_start = rs[ev.runno].starttime; // this time should really be loaded from the start time which is save in the stats branch
	//hi_time_start = ev[i].hi_tim;
	for(i=0,j=0,evtree->GetEvent(i); i<nentries ;){
	
		//SKIPFLAGS();  //skip flagged groups and runs (defined above)
		//cout << "hi1\n";
		while( 	(
					gs[ev.groupno].flag || // do not plot if group is flagged
					gs[ev.groupno].hi_thresh[detn-1] < Emax || // only plots if the data energy range falls within specified range
					!gs[ev.groupno].det_on[detn-1] 
				) &&  
				i<nentries 
				) {
			if(gbl_verbose) cout << "skipping group " << ev.groupno+1 << "...\n";
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
			time_start = rs[ev.runno].starttime;
		}
		if(ev.groupno > cur_group){ 
			cur_group = ev.groupno;
			p_cnt = 0;
			p_cnt_buf = 0;
			run_cnt = 0;  // do not count across groups
			used_run_cnt = 0;
		}
		if(run_cnt==0){ 
			time_start = rs[ev.runno].starttime;
			runno_buf = ev.runno;
		}
		if( ( 
				!rs[ev.runno].det_on || 
				(use_flags[1]!='_' && rs[ev.runno].flag[detn-1]==(use_flags[1]=='0')) ||
				rs[ev.runno].lo_thresh[detn-1] > Emin
					
			) && 
			i<nentries
			){	
			
			i += rs[ev.runno].nevents;
			evtree->GetEvent(i);
		}		
		else {
			for(cur_run=rs[ev.runno].starttime;rs[ev.runno].starttime==cur_run && i<nentries;i++,evtree->GetEvent(i)){ // loop over run
				enteredloop = 1;
				if((use_flags[0]=='_' || ((use_flags[0]=='0') && !flag[i]) || ((use_flags[0]=='1') && flag[i]) || int((use_flags[0])==flag[i]))
				&& (use_flags[2]=='_' || !llfn_Decimal2Bit(BFlag[i],detn,NDETS)==(use_flags[2]=='0')))
					if((ev.energy[detn-1] >= Emin && Emax == 0) || (ev.energy[detn-1] >= Emin && ev.energy[detn-1] <= Emax)) 
						p_cnt_buf++; // count the number of events in period
			}
			used_run_cnt++;
		}
		run_cnt++;
		p_cnt+=p_cnt_buf;
		p_cnt_buf=0;
		if(run_cnt >= period && enteredloop){
			if(used_run_cnt){
				epp.epp.push_back(p_cnt);
				epp.used_run_cnt.push_back(used_run_cnt);
				if(axisScale>0) epp.time.push_back((time_start-start_time)/axisScale);
				else epp.time.push_back(time_start-25*3600*24*365);
			}
			run_cnt=0;
			used_run_cnt = 0;
			p_cnt=0;
			enteredloop = 0; // this prevents false zeros being entered.
			j++;
		}
	}
	

	
	//-------------- code for the timeline plot ------------
		Double_t *yerrors = new Double_t[epp.epp.size()];
		Double_t *y = new Double_t[epp.epp.size()];
		
		for(k=0;k<int(epp.epp.size());k++) {
			y[k] = double(epp.epp[k]*period)/epp.used_run_cnt[k];
			yerrors[k] = (1. + sqrt(epp.epp[k] + 0.75))*period/epp.used_run_cnt[k]; // This appoximation to the mean is for Poisson distribution and obtained from STSDAS Help Pages "explain_errors" http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?explain_errors 
		}
		TGraphErrors *g_handle2 = new TGraphErrors(epp.epp.size(),&epp.time[0],y,0,yerrors);
		sprintf(buffer,"D%d   time line",detn);
		g_handle2->SetTitle(buffer);
		g_handle2->GetXaxis()->SetTitle("Time (days)");
		sprintf(buffer,"Counts/%d-%d keV/%d hour(s)",int(Emin),int(Emax),period);
		g_handle2->GetYaxis()->SetTitle(buffer);
		//TGraph *graph = new TGraph(k,time_stamp,period_cnt);
		
		if(axisScale==0) {
			g_handle2->GetXaxis()->SetTimeDisplay(1);
			g_handle2->GetXaxis()->SetTimeFormat("");
			g_handle2->GetXaxis()->SetNdivisions(-510);
			g_handle2->GetXaxis()->SetTitle("Real time");
		}
		
		g_handle2->Draw("A*");
		g_handle=g_handle2;
		
}

int llfn_NOT(int val, int sign){
	if(sign>0) return val;
	else return !val;
}

void SSplot_waterfall(char help='H'){
	help = help;
	cout << "detn:................Detector Number\n"
		 << "markerstyle..........6 = *\n"
		 << "                     1 = .\n"
		 << "time_scale...........x-axis scale. set to 1 to plot in seconds, 60 = minutes, 3600 = hours etc.\n"
		 << "                                   set to 0 to plot in real time / date format \n"
		 << "use_flags............\"XYZ\"\n";
	cout << "                     X = events flags\n"
		 << "                     Y = run flags\n"
		 << "                     Z = binary event flag\n"
		 << "                     options for A,B,C: _, 0, 1, where\n"
		 << "                     _ = ignore flags\n"
		 << "                     0 = do not use flagged events\n"
		 << "                     1 = only use flagged events\n"
		 << "                     set use_flags = \"A\" to auto matically highlight flagged events\n"
		 << "     new===>         You can now set the event flag in the data to any number you like except \n"
		 << "                        95=int('_'), NULL=int('0'), SOH=int('1'), and 65 = int('A'). You can then\n"
		 << "                        only plot these events by setting parameter X to your chosen coresponding character.\n"
		 << "                        eg: if you flagged some events with flag[i]=int('f') then set use_flags=\"f00\n";
	cout << "overlay..............setting this draws events ontop of existing waterfall canvas\n"
		 << "colour...............set marker colour, eg 0=while, 1=black...\n";
}
void SSplot_waterfall(Int_t detn, Int_t markerstyle = 6, Int_t time_scale = 3600, Char_t use_flags[3] = "___", Int_t overlay = 0, Int_t colour=1, Int_t logy=1){
// plots energy vs time of events from detn
// ev_runflags: "__" = ignor ev and run flags
//				"0_" = subtract ev flags, don't use run flagged data
//				"10" = only plot ev flags, don't plot run flagged data

// 				"A1" = Auto mode 1
	gblperiod = time_scale;
	gbldety = detn;
	

	if(use_flags[0]=='A'){
		SSplot_waterfall(detn,markerstyle,time_scale,"___",0,2);
		SSplot_waterfall(detn,markerstyle,time_scale,"000",1,1);
		return;
	}
	
	Char_t buffer[100];
	int i, j;
	gROOT->SetStyle("Plain");
	Double_t *energy, *timex;
	energy = (Double_t *)malloc(sizeof(Double_t)*nentries);
	if(energy == NULL){
		cerr << "Out of Memory!" << endl;
		exit(-1);
	}
	timex = (Double_t *)malloc(sizeof(Double_t)*nentries);
	if(timex == NULL){
		cerr << "Out of Memory!" << endl;
		exit(-1);
	}

	
	//Double_t *cut_energy = new Double_t[nentries], *cut_time = new Double_t[nentries];
	Double_t start_time = 0.;
	if(time_scale > 0) {
		evtree->GetEvent(0); 
		start_time = rs[ev.runno].starttime;
	}
	
	for(i=0,evtree->GetEvent(i), j=0; i<nentries ;i++,evtree->GetEvent(i)){
		// skip flagged groups
		while((gs[ev.groupno].flag || !gs[ev.groupno].det_on[detn-1]) && i<nentries){
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
		}
		// skip flagged run: Skip flagged data if (flagegd && 0) || (!flagged && !0)
		while( (!rs[ev.runno].det_on[detn-1] || (use_flags[1]!='_' && rs[ev.runno].flag[detn-1]==(use_flags[1]=='0'))) && i<nentries){
			i += rs[ev.runno].nevents;
			evtree->GetEvent(i);
		}
		// check ev flags  
		if(	(use_flags[0]=='_' || ((use_flags[0]=='0') && !flag[i]) || ((use_flags[0]=='1') && flag[i]) || (use_flags[0]==flag[i])) && 
			(use_flags[2]=='_' || !llfn_Decimal2Bit(BFlag[i],detn,NDETS)==(use_flags[2]=='0'))
			){
			energy[j] = ev.energy[detn-1];
			if(time_scale==0) timex[j] = ev.time - 25*3600*24*365;
			else timex[j] = (ev.time - start_time)/time_scale;
			j++;
		}
	}
	
	
	c1 = new TCanvas("wf","wf",10,10,600*2,400);
	if(logy) c1->SetLogy();
	
	llfn_makewhite(c1);
	//sprintf(buffer,"Waterfall plot. D%d",trig_filter);
// 	if(g_waterfall != NULL) g_waterfall->Delete();
	TGraph *g_waterfall = new TGraph(j,timex,energy);
	sprintf(buffer,"Waterfall plot. ADC%d",detn);
	g_waterfall->SetTitle(buffer);
	if(time_scale==0) sprintf(buffer,"Real Time");
	else sprintf(buffer,"Time / %d seconds",time_scale);
	g_waterfall->GetXaxis()->SetTitle(buffer);
	if(gs[0].cal_eqn[0][detn-1]!=1.) g_waterfall->GetYaxis()->SetTitle("Energy (keV)");
	else g_waterfall->GetYaxis()->SetTitle("Channels");
	g_waterfall->GetYaxis()->SetTitleOffset(1.3);
	g_waterfall->SetMarkerStyle(markerstyle);
	g_waterfall->SetMarkerColor(4);
	
	// Options 0 -> set colour to black, 1 -> Set colour to red, 2 -> set colour to black and overlay, any other number set colour to blue
	if (overlay==1){ 
		g_waterfall->SetMarkerColor(colour);
		g_handle2 = g_waterfall;		
		g_handle->Draw("AP");
		g_handle2->Draw("P");
		if(time_scale==0) {
			g_handle->GetXaxis()->SetTimeDisplay(1);
			g_handle->GetXaxis()->SetTimeFormat("");
			g_handle->GetXaxis()->SetNdivisions(-510);
		}
	}
	if(overlay==0) {
		g_waterfall->SetMarkerColor(colour);
		g_handle = g_waterfall;
		g_handle->Draw("AP");
		if(time_scale==0) {
			g_handle->GetXaxis()->SetTimeDisplay(1);
			g_handle->GetXaxis()->SetTimeFormat("");
			g_handle->GetXaxis()->SetNdivisions(-510);
		}
		if(gblpause) gPad->WaitPrimitive();
	}
		
	SSplot_settings_userthr(detn);
		
	
	//TH1D *h_spec = new TH1D("wf_plot","waterfall plot",Int_t(time[j]-time[0]),time[0],time[j]);
	//for(i = 0; i<=j;i++) h_spec->Fill(time[i],energy[i]);
	//h_spec->Draw();
	
	
	delete[] energy;
	delete[] timex;
	//delete[] cut_energy;
	//delete[] cut_time;

}

void SSplot_waterfall_window_help() {
	cout << "This function plots the energy of the events on the y-axis and the time on the x-axis\n"
	     << "for the specified time window. Times are in unix time.\n";
}
	
void SSplot_waterfall_window(Int_t detn, Double_t t_start, Double_t t_end/*, Chart_t use_flags = "___"*/) {
	int i, j, k;
	vector <Double_t> time;
	vector <Double_t> energy;
// 	time.push_back(0);
// 	energy.push_back(0);
	
	Double_t offset=0;
	int inwindow=0;
	
	
	for(i=0;i<ngroups;++i) {  // loop over groups
		// skip group if flagged
		if( gs[i].starttime[detn-1] > t_end || 
			(ngroups>=i+1)*gs[i+1].starttime[detn-1]<t_start || 
			!gs[i].det_on[detn-1] || 
			gs[i].flag /*|| Emax > gs[i].hi_thresh[detn-1] */
		 	) continue;
		for(j=gs[i].firstrunno;j<nruns && j<gs[i+1].firstrunno;++j) { // loop over runs
			// skip runs if flagged
			if(	rs[j].starttime > t_end ||
				(nruns>=j+1)*rs[j+1].starttime<t_start ||
				!rs[j].det_on[detn-1] || 
				rs[j].flag[detn-1] 
				/*|| Emin < rs[j].lo_thresh[detn-1]*//*SSch2e(SSth2ch.Eval(rs[j].thr[detn-1]),detn)*/ 
				) continue;
		
			for(k=rs[j].firstevno;k<nentries && k<rs[j+1].firstevno;++k,evtree->GetEvent(k)) { // loop over event in run
				if( ev.time < t_start ||
					ev.time > t_end ||
					flag[k] || 
					llfn_Decimal2Bit(BFlag[k],detn,NDETS)
					) continue;
				
				if(ev.time-t_start > 10 ) {
					time.push_back(ev.time-t_start);
					
				}
				else {
					if(!inwindow) {
						inwindow=1;
						offset=ev.hi_tim*8e-6 - t_start;
					}
					time.push_back(ev.hi_tim*8e-6-t_start-offset);
				}
				energy.push_back(ev.energy[detn-1]);
			}
		}
	}	
	
// 	time.push_back(t_end-t_start);
// 	energy.push_back(0);
	
	g_handle = new TGraph(time.size(),&time[0],&energy[0]);	
	g_handle->SetMarkerStyle(6);
	g_handle->Draw("AP");
}
	

	

void SSplot_3Dwaterfall(Int_t detn, Int_t time_scale = 3600, Double_t Emax =0, Char_t use_flags[2] = "__"){
	int i, j;
	gROOT->SetStyle("Plain");
	Double_t *energy, *timex, *veto;
	energy = (Double_t *)malloc(sizeof(Double_t)*nentries);
	if(energy == NULL){
		cerr << "Out of Memory!" << endl;
		return;
	}
	timex = (Double_t *)malloc(sizeof(Double_t)*nentries);
	if(timex == NULL){
		cerr << "Out of Memory!" << endl;
		return;
	}
	veto = (Double_t *)malloc(sizeof(Double_t)*nentries);
	if(timex == NULL){
		cerr << "Out of Memory!" << endl;
		return;
	}	
	cout << "hi\n";
	//Double_t *cut_energy = new Double_t[nentries], *cut_time = new Double_t[nentries];
	Double_t start_time = 0.;
	if(time_scale == 0) time_scale = 1;  // by setting time_scale=0 the unix time of the event is shown
	else {evtree->GetEvent(0); start_time = rs[ev.runno].starttime;}
	
	for(i=0,evtree->GetEvent(i), j=0; i<nentries ;i++,evtree->GetEvent(i)){
		SKIPFLAGS();  //skip flagged groups and runs (defined above)
		if((use_flags == 0 || !flag[i]) && (Emax==0 || ev.energy[detn-1] < Emax) ){
			energy[j] = ev.energy[detn-1];
			timex[j] = (ev.time - start_time)/time_scale;
			veto[j] = ev.veto;
			j++;
		}
	}
	cout << "hi2\n";
	
	c1 = new TCanvas("wf","wf",10,10,600,400);
	llfn_makewhite(c1);
	TGraph2D *g2handle = new TGraph2D("g2d","g2d",j,timex,energy,veto);
// 	sprintf(buffer,"Time / %d seconds",time_scale);
// 	g2handle->GetXaxis()->SetTitle(buffer);
// 	g2handle->GetYaxis()->SetTitle("Energy (keV)");
// 	g2handle->GetZaxis()->SetTitle("veto");
		
	g2handle->Draw("AP");
	
}
 
/// xxx I think I should rewrite this function
void SSplot_timediffs(Double_t xmin = 0, Double_t xmax = 3500, Int_t nbins = 3500,  Int_t use_flags = 1, Double_t eph = 0, Int_t detn = 0, Double_t Emin = 0, Double_t Emax = 0){  // This could be overloaded with an option to fit theoretical curve. 
	// WARNING: function updated but not tested (22/04/06)
	// NB If just one point was flagged note that the function will skip onto next point and not count the time between the last unflagged point
	// eph is the the number events per hour in chosen range
	int i, i_log=0, last_was_flagged = 0; 
	int *used_group = new int[ngroups]; for(i=0;i<ngroups;i++) used_group[i] = 0;
	evtree->GetEvent(0); Double_t time_old = ev.time, hi_tim_old = ev.hi_tim; // this should really be checked to be an event that passes the if condition below.
	Double_t used_time =0;
	//Double_t total_flagged_time = 0;
	Char_t buffer[50];
	c1 = new TCanvas("c1", "c1",14,30,700,500); 
	c1->SetLogy();
	llfn_makewhite(c1);
	sprintf(buffer,"Time between events. D%d",detn);
	TH1D *h_td = new TH1D("hdt",buffer,nbins,xmin,xmax);
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		
		// SKIP FLAGS	
		while( (gs[ev.groupno].flag || 						// do not plot if group is flagged
			( Emax!=0 && (gs[ev.groupno].hi_thresh[detn-1] < Emax)) ||
			!gs[ev.groupno].det_on[detn-1]) &&  			// only plots if the data energy range falls within specified range
			i<nentries){
			
			if(gbl_verbose) cout << "skipping group " << ev.groupno+1 << "...\n";
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
			time_old = ev.time;
			i_log = i; 
			hi_tim_old = ev.hi_tim;			
		}
		while(i<nentries && 
				(	
					(use_flags && rs[ev.runno].flag[detn-1]) ||
				 	!rs[ev.runno].det_on[detn-1]
				)
			) {
				
			i += rs[ev.runno].nevents;
			evtree->GetEvent(i);
			time_old = ev.time;
			i_log = i; 
			hi_tim_old = ev.hi_tim;			
		}
// 		&& (use_flags[0]=='_' || !flag[i]==(use_flags[0]=='0'))
// 			&& (use_flags[2]=='_' || !llfn_Decimal2Bit(BFlag[i],detn,NDETS)==(use_flags[2]=='0')) 
				
		if( (use_flags == 0 || !flag[i]) && 
			( detn == 0 || (ev.energy[detn-1] > Emin && Emax == 0) || (ev.energy[detn-1] > Emin && ev.energy[detn-1] < Emax) ) &&
			(detn ==0 || (Emin > rs[ev.runno].lo_thresh[detn-1])) 
			){
			
			if(!last_was_flagged){
				if(ev.time-time_old > 10){ 
					h_td->Fill( (ev.time-time_old));
					used_group[ev.groupno] = 1;
				}
				else {
					h_td->Fill( (ev.hi_tim - hi_tim_old)*8E-6 );
					used_group[ev.groupno] = 1;
				}
			}
			else last_was_flagged = 0; 
			time_old = ev.time;
			i_log = i; 
			hi_tim_old = ev.hi_tim;
		}
		else if(flag[i]) last_was_flagged = 1;
	}
	
	// normalise   NEED TO UPDATE THIS
// 	for(i=0;i<ngroups;i++){
// 		if(use_flags && used_group[i]) 	used_time += gs[i].uptime[detn-1] - gs[i].flaggedtime[detn-1];
// 		else if(used_group[i]) used_time += gs[i].uptime[detn-1];
// 	}
	if(gbl_verbose) cout << "used time = " << used_time << endl;
	h_td->Sumw2();
	h_td->GetXaxis()->SetTitle("Event separation (Seconds)");
	h_td->GetYaxis()->SetTitle("Probability");
	//h_td->Draw();
	h_handle = h_td;
	
	 
	/// this is my old way of fitting. A better way is to normalise to one and then fit an exponential.
	if(eph != 0) { 
		Double_t td_rate = eph/3600;
		sprintf(buffer,"[0]*%f*exp(-1*%f*x)",td_rate, td_rate);
		TF1 *exp1 = new TF1("exp1",buffer,xmin,xmax);
		exp1->SetParName(0,"scale");
		exp1->SetParameter(0,1);
		exp1->SetLineColor(4);
		exp1->SetLineWidth(2);
		h_handle->Fit("exp1","0");
		cout << "   5  ProbChi      " << TMath::Prob(exp1->GetChisquare(),exp1->GetNDF()) << "\n";
		h_td->Scale(1./exp1->GetParameter(0));
		sprintf(buffer,"%f*exp(-1*%f*x)",td_rate, td_rate);
		TF1 *exp2 = new TF1("exp2",buffer,xmin,xmax);
		exp2->SetLineColor(4);
		exp2->SetLineWidth(2);
		h_td->Draw();	
		exp2->Draw("same");
		
	}
	else h_td->Draw();
}

void SSplot_energyratios(){
	gROOT->SetStyle("Plain");
	c1 = new TCanvas("cp","cp",10,10,600,400);
	llfn_makewhite(c1);	
	int i, j, k;
	char buffer[512];
	Double_t *energy, *det;
	energy = (Double_t *)malloc(sizeof(Double_t)*nentries);
	if(energy == NULL){
		cerr << "Out of Memory!" << endl;
		exit(-1);
	}
	det = (Double_t *)malloc(sizeof(Double_t)*nentries);
	if(det == NULL){
		cerr << "Out of Memory!" << endl;
		exit(-1);
	} 
	for(i=0, k=0, evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		if(flag[i]){ 
			for(j=1;j<=NDETS;j++){
				energy[k] = ev.energy[j-1];
				det[k] = j;
				k++;
			}
		}
	}
	
	g_handle = new TGraph(k,det,energy);
	sprintf(buffer,"Correlated energies of all detectors");
	g_handle->SetTitle(buffer);
	g_handle->GetXaxis()->SetTitle("Detector number");
	g_handle->GetYaxis()->SetTitle("Energy");
	g_handle->GetYaxis()->SetTitleOffset(1.3);
	g_handle->SetMarkerStyle(6);
	g_handle->SetMarkerColor(4);
	g_handle->Draw("AP");
}

void SSplot_correlations(Int_t detx, Int_t dety, Int_t offset = 0, Int_t use_flags = 0, Int_t overlay012 = 0, Int_t markerstyle=6, Double_t p1x = 0, Double_t p1y=0, Double_t p2x = 0, Double_t p2y=0, Double_t p3x=0, Double_t p3y=0, Double_t p4x=0, Double_t p4y=0){
	// plots the energy of one detector against the energy of another detector for the same adc reading.
	// offset: can be used to plot one event on a detector against a pre or proceeding event of the same or another detector.
	gROOT->SetStyle("Plain");
	c1 = new TCanvas("cp","cp",10,10,600,400);
	llfn_makewhite(c1);	
	int i;
	char buffer[512];
	Double_t *ex, *ey, a1[2],a2[2],a3[2],a4[2];
	
// 	ex = (Double_t *)malloc(sizeof(Double_t)*nentries);
	ex = new Double_t[nentries];
	if(ex == NULL){
		cerr << "Out of Memory!" << endl;
		exit(-1);
	}
// 	ey = (Double_t *)malloc(sizeof(Double_t)*nentries);
	ey = new Double_t[nentries];
	if(ey == NULL){
		cerr << "Out of Memory!" << endl;
		exit(-1);
	}
	
	// save used detector numbers to global buffers
	gbldetx = detx;
	gbldety = dety;
	
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		// skip flagged groups
		while((gs[ev.groupno].flag || !gs[ev.groupno].det_on[detx-1] || !gs[ev.groupno].det_on[dety-1]) && i<nentries){
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
		}
		
		if(use_flags == 0 || (use_flags==1 && !flag[i]) || (use_flags && flag[i]!=use_flags)){
			ex[i]=ev.energy[detx-1];
			if(!offset){
				if(i+offset < 0 || i+offset >=nentries){
					ey[i] = -1; 
					continue;
				}
			}
				evtree->GetEvent(i+offset);
			ey[i]=ev.energy[dety-1];
			if(!offset) evtree->GetEvent(i-offset);
			if(p1x>0){
				llfn_poly1(a1,p1x,p1y,p2x,p2y);
				llfn_poly1(a2,p2x,p2y,p3x,p3y);
				llfn_poly1(a3,p3x,p3y,p4x,p4y);
				llfn_poly1(a4,p4x,p4y,p1x,p1y);
				if (ey[i] <= a1[1]*ex[i]+a1[0] && ex[i] <= ey[i]/a2[1]-a2[0]/a2[1] && ey[i] >= a3[1]*ex[i]+a3[0] && ex[i] >= ey[i]/a4[1]-a4[0]/a4[1]) flag[i]=1;
			}
		}
	}
	cout << "Equations of lines: 1)top, 2)right, 3)bottom, 4)left\n";
	cout << "y1 = " << a1[1] << " * x1 + " << a1[0] << endl;
	cout << "x2 = " << a2[1] << " * y2 + " << a2[0] << endl;
	cout << "y3 = " << a3[1] << " * x3 + " << a3[0] << endl;
	cout << "x4 = " << a4[1] << " * y4 + " << a4[0] << endl;
	
	
	TGraph *g_h = new TGraph(nentries,ex,ey);
	sprintf(buffer,"Correlation of energies D%d vs D%d",dety,detx);
	g_h->SetTitle(buffer);
	if(gs[0].cal_eqn[0][detx-1]!=1.) sprintf(buffer,"Energy (keV)    D%d",detx);
	else sprintf(buffer,"Energy (Channels)    D%d",detx);
	g_h->GetXaxis()->SetTitle(buffer);
	if(gs[0].cal_eqn[0][dety-1]!=1.) sprintf(buffer,"Energy (keV)    D%d",dety);
	else sprintf(buffer,"Energy (Channels)    D%d",dety);	
	g_h->GetYaxis()->SetTitle(buffer);
	g_h->GetYaxis()->SetTitleOffset(1.3);
	g_h->SetMarkerStyle(markerstyle);
	g_h->SetMarkerColor(4);
	
	
	
	if (overlay012 == 2){ 
		g_h->SetMarkerColor(4);
// 		g_handle2->Delete();
		g_handle2 = g_h;
		g_handle->Draw("AP");
		g_handle2->Draw("P");
	}else{
		if(overlay012 == 1) g_h->SetMarkerColor(2);
		if(overlay012 == 0) g_h->SetMarkerColor(4);
// 		g_handle->Delete();
		g_handle = g_h;
		g_handle->Draw("AP");
	}	
	
	delete [] ex;
	delete [] ey;
}

void SSplot_aboveThresholdEvents(){
	int i, j;
	TH1D *thECounts = new TH1D("thECounts","Number of counts above detector thresholds",NDETS,1.,NDETS+1.);
	thECounts->GetXaxis()->SetTitle("Detector");
	thECounts->GetYaxis()->SetTitle("Counts per hour");
	for(i=0,evtree->GetEvent(0);i<nentries;++i,evtree->GetEvent(i)) {
		for(j=0;j<NDETS;++j) {
			if( !flag[j] && 
				!rs[ev.runno].flag[j] && 
				rs[ev.runno].det_on[j] &&
				gs[ev.groupno].det_on[j] &&
				!gs[ev.groupno].flag
				) {
				if(ev.energy[j] > SSch2e(SSth2ch.Eval(rs[ev.runno].thr[j]),j+1,ev.groupno)) thECounts->Fill(j+1);
			}
		}
	}
	
	//scale to counts / hour
	for(j=0;j<NDETS;++j) thECounts->SetBinContent(j+1,thECounts->GetBinContent(j+1)/(SSgetusedtimeR(j+1,0,0)+1E-5)*3600);
	
	thECounts->SetFillColor(2);
	thECounts->SetFillStyle(3001);
	thECounts->Draw();
	
			
	cout << "This function does not handle groups correctly\n";
	// NB this function can only act as a guide because sfalist2root only uses the latest current largest threshold and not the finaly largest
}


/// Searching functions		

void SSsearchc(Int_t detn, Double_t tperiod, Double_t Emin1=0, Double_t Emax1=0, Double_t Emin2=0, Double_t Emax2=0, Int_t use_flags = 1, Int_t verbose =1){ // search for conincidences
	int i, j, counted = 0, coincidences = 0, looking_for_first_event=1, i_log=0; 
	evtree->GetEvent(0); Double_t first_time = ev.time, first_hitime = ev.hi_tim, E1=0;
	Int_t detn_start = detn, detn_end = detn, use_deltaE=0;

FILE * s_fp=NULL;	
if(verbose==2){	
		if(!(s_fp = fopen("search.dat","w"))){
		cout << "Error creating file: search.dat" << endl;
		exit(1);
	}
}
	
	if(Emin1==-1){ 
		use_deltaE = 1; // to use E +- deltaE set Emin1 = -1 and Emax1 = E
		E1 = Emax1;
	}
	if(detn == 0){
		detn_start = 1;
		detn_end = ndetectors;
	}
	Double_t energy_log;
	printf("   E1      E2       dt      group\n");
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){ 

		while( 	
				(
					gs[ev.groupno].flag || // do not plot if run is flagged
					((gs[ev.groupno].hi_thresh[detn-1] < Emax1) && (gs[ev.groupno].hi_thresh[detn-1] < Emax2) ) ||
					!gs[ev.groupno].det_on[detn-1]
				) &&  // only plots if the data energy range falls within specified range
				i<nentries
			){
			cout << "skipping group " << ev.groupno+1 << "...\n";
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
			looking_for_first_event = 1;
		}
		while( 
				(
					(use_flags && rs[ev.runno].flag[detn-1]) ||
					(Emin1 < rs[ev.runno].lo_thresh[detn-1]) ||
					(Emin2 < rs[ev.runno].lo_thresh[detn-1]) ||
					!rs[ev.runno].det_on[detn-1]
				) && 
				i<nentries
			) {	
			i += rs[ev.runno].nevents;
			evtree->GetEvent(i);
			looking_for_first_event = 1;
		}	
	
		if((use_flags == 0 || !flag[i]) 
			){//&& counted != i ){  // check that the previous event has not already been used. 
			for(j=detn_start;j<=detn_end;j++){
				if(looking_for_first_event){
					if(use_deltaE){
						Emin1 = E1 - (gs[ev.groupno].res_eqn[0][j-1]*E1 + gs[ev.groupno].res_eqn[2][j-1]); 
						Emax1 = E1 + (gs[ev.groupno].res_eqn[0][j-1]*E1 + gs[ev.groupno].res_eqn[2][j-1]);
					}
					if ( (ev.energy[j-1] > Emin1 && Emax1 == 0) || (ev.energy[j-1] > Emin1 && ev.energy[j-1] < Emax1) ){
						looking_for_first_event = 0;
						first_time = ev.time;
						first_hitime = ev.hi_tim;
						i_log = i;
					}
					break;
				}
				else{
					if( ((ev.time - first_time) < tperiod && (ev.time - first_time) > 10) || ((ev.time - first_time) <= 10 && (ev.hi_tim - first_hitime)*8E-6 < tperiod && (ev.hi_tim - first_hitime)*8E-6 >=0) ){ // check if time between events is less than period. If normal system time difference is zero compare high accuracy time. the hi timer is used for times less then 10 seconds as this is well away from the limit of the clock
						if( (ev.energy[j-1] > Emin2 && Emax2 == 0) || (ev.energy[j-1] > Emin2 && ev.energy[j-1] < Emax2) ){
							coincidences++;
							looking_for_first_event = 1;
							counted = i; // might need to use flags so no further events are re-used. this is reset currantly so only rmembers the last one
							//ev[i].flag = 1; ev[i_log].flag = 1;
							evtree->GetEvent(i_log);energy_log=ev.energy[j-1]; evtree->GetEvent(i);
							if(verbose==1){ 
								printf("%7.2f,  %7.2f,  %1.2e, ", energy_log,ev.energy[j-1], (ev.hi_tim - first_hitime)*8E-6);
								llfn_getgroupname(ev.groupno); cout << endl;
							}
							if(verbose==2) fprintf(s_fp,"%7.2f  %7.2f  %1.2e\n",
								energy_log,ev.energy[j-1], (ev.hi_tim - first_hitime)*8E-6);
							if(verbose==3) printf("%7.2f, %1.2e\n",
								energy_log/ev.energy[j-1], (ev.hi_tim - first_hitime)*8E-6);	 
							//cout << "E1 = " << ev[i_log].energy[j-1] << " E2 = " << ev[i].energy[j-1] <<  " " << Emin2 
							//		<< " dt = " << (ev[i].hi_tim - first_hitime)*8E-6  << "s" 
							//		<< " group: "<< gs[ev[i_log].groupno].group_name << endl;
									//<< " " << first_hitime << " " << ev[i].hi_tim << endl;
							i = i_log; evtree->GetEvent(i);
							break; // so far this only checks if at least one detectors fires. To reject coincidences where more than one detector fires will need another condition. also might want to check if another detector fires using a difference energy range. possibly should overload input to do this.
						}
						else{ 
							if( detn != 0 || (detn == 0 && j == detn_end) ){
								looking_for_first_event = 1;
								i = i_log; evtree->GetEvent(i); // go back to start position
							}
						}
					}
					else{
						looking_for_first_event = 1;
						i = i_log; evtree->GetEvent(i);
						break;
					}
				}
			} 
		}
	}
	cout << "Number of coinidences found = " << coincidences << endl;
	if(verbose==2) fclose(s_fp);
} 

struct StructCoindence{
	Double_t dt;
	Double_t E1;
	Double_t E2;
};

void SSsearch_back(char help = 'H'){
	
	cout << "Options:\n";
	cout << "enter times in units of seconds\n";
	cout << "ploth = 0: Do not plot results \n";
	cout << "ploth = 1: Plot dt spectrum   (opt is not used)\n";
	cout << "ploth = 2: Plot energy spectrum\n";
	cout << "  opt = 1: Plot spectrum of E2\n";
	cout << "  opt = 2: Plot spectrum of E1\n";
	cout << "ploth = 3: Plot E2 vs E1      (opt is not used)\n";
	cout << "flagfinds: - set to zero if found events should not be flagged\n";
	cout << "              set to int('c') where c can be any character except _,0,1 or A.\n";
	cout << "                  to set flag to that charactors number\n\n";
	
	help = help;
}

Int_t SSsearch_back(Int_t detn, Double_t Emin2, Double_t Emax2, Double_t Emin1, Double_t Emax1, Double_t tmin=0, Double_t tmax = 1, Int_t nbins=10, Int_t use_flags=1, Int_t ploth=0, Int_t opt=3, Int_t flagfinds = 0, Int_t n=1, Double_t tscale = 1.e-3){
// This function searches for n events after first events and plots thier delta t's energy to one 2d histogram 
	Int_t i, cnt_evnt=0, event = 1, i_log=0, cnt=0;
	Double_t E1=0, E2=0, t0=0, dt;
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat("e");//gStyle->SetOptStat(10);//gStyle->SetOptStat(1000010); 	// set stat box options
	gStyle->SetOptFit(1111);			// set fitting options
	char buffer[100];
	
	int skipped_groups = 0;
	int skipped_runs = 0;
	
// 	StructCoindence coin;
// 	TFile *file = new TFile("search_back.root","UPDATE","search_back",9);
// 	sprintf(buffer,"det%d",detn);
// 	TTree *ttree = new TTree(buffer,buffer);
// 	ttree->Branch(buffer,&coin,"dt/D:E1/D:E2/D");
	
	//TH2D *th2d = new TH2D("th3d","th3d",1000,0,1,100,0,10000);
	if(ploth==1){ 
		sprintf(buffer,"D%d",detn);
		h_handle = new TH1D(buffer,buffer,nbins,1./tscale*tmin,tmax/tscale);
		sprintf(buffer,"Time between events (%e seconds)",tscale);
		h_handle->GetXaxis()->SetTitle(buffer);
		h_handle->GetYaxis()->SetTitle("Counts");
	}
	else if(ploth==2) { 	
		sprintf(buffer,"D%d",detn);
		Double_t upper_range; 
		if( Emin1*Emin2*Emax1*Emax2>0) upper_range = (Emax2>=Emax1)*Emax2+(Emax1>Emax2)*Emax1;
		else upper_range = 100;
		h_handle = new TH1D(buffer,buffer,Int_t(upper_range/nbins),0.,upper_range);
		h_handle->GetXaxis()->SetTitle("Energy");
		h_handle->GetYaxis()->SetTitle("Counts");
	}
	else if(ploth==3) {
		h2_handle = new TH2D("e1vse2","E2 vs E1",nbins,0.,Double_t(Emax1),nbins,0.,Double_t(Emax2));	
		h2_handle->GetYaxis()->SetTitle("E2");
		h2_handle->GetXaxis()->SetTitle("E1");
	}
	else if(ploth==4) {
		h2_handle = new TH2D("e2vsdt","E2 vs dt",nbins,0.,tmax/tscale,nbins,0.,Double_t(Emax2));	
		h2_handle->GetYaxis()->SetTitle("E2 (keV)");
		//h2_handle->GetXaxis()->SetTitle("dt (ms)");
		sprintf(buffer,"Time between events (%e seconds)",tscale);
		h2_handle->GetXaxis()->SetTitle(buffer);
	}	
	
	
	cout << " dt,  E2,  E1  i2, i1, t1\n";
	for(i=nentries-1,evtree->GetEvent(i); i>=0 ;i--,evtree->GetEvent(i)){
	
		while(  (gs[ev.groupno].flag ||  									// do not search if run is flagged
				(Emax2>0 && gs[ev.groupno].hi_thresh[detn-1] < Emax2) ||	// do not search if Emax2 is outside uppper range	 
				(Emax1>0 && gs[ev.groupno].hi_thresh[detn-1] < Emax1) ||	// do not search if Emax1 is outside upper range
				!gs[ev.groupno].det_on[detn-1]) 							// do not search if the detector is not to be used for this group
			) {  
			
			
			if(i<=0) break;
			
			++skipped_groups;
			cout << "skipping group " << ev.groupno+1 << "...\n";
			i -= gs[ev.groupno].nevents;
			evtree->GetEvent(i);
			event = 1;
			cnt_evnt = 0;
			
		}
		// not that during a long timing search if a flagged period intercepts the search the counter is reset to events with a flagged period are not included.
		while(	!rs[ev.runno].det_on[detn-1] ||						// Do not search if the detector is not to be used for this run
				(use_flags && rs[ev.runno].flag[detn-1]) || 			// Do not search if the run is flagged
				(Emin2>0 && rs[ev.runno].lo_thresh[detn-1] > Emin2) ||  				// Do not use if Emin1 and Emin1 are below the lower threshold for this run
				(Emin1>0 && rs[ev.runno].lo_thresh[detn-1] > Emin1 )
			) {
			
			if(i<=0) break;
			
			++skipped_runs;
			i -= rs[ev.runno].nevents;
			evtree->GetEvent(i);
			event = 1;
			cnt_evnt = 0;
		}	
		
		if(i<=0) break;
		
		if(event==1){
			if((ev.energy[detn-1] > Emin2 && ev.energy[detn-1] < Emax2) && 
				(use_flags==0 || (use_flags==1 && !flag[i]) || (use_flags>1 && use_flags==flag[i]))
				){
				E2 = ev.energy[detn-1];
				t0 = ev.hi_tim;
				event = 2;
				i_log = i;
				continue; // don't want to check for second event before moving back a step
			}
		}
		if(event==2){
			if(((ev.energy[detn-1] > Emin1 && Emax1 ==0) || (ev.energy[detn-1] > Emin1 && ev.energy[detn-1] < Emax1)) && 
				(use_flags==0 || (use_flags==1 && !flag[i]) || (use_flags>1 && use_flags==flag[i]))
				){
				E1 = ev.energy[detn-1];
				dt = (t0-ev.hi_tim)*8E-6; // convert to seconds
				//th2d->Fill(dt,E1);
				if(dt > tmin && dt < tmax){
					if(ploth==1) h_handle->Fill(dt/tscale);
					else if(ploth==2){
						if(opt==1) h_handle->Fill(E2);
						if(opt==2) h_handle->Fill(E1);
					}
					if(ploth==3) h2_handle->Fill(E1,E2);
					if(ploth==4 && opt==1) h2_handle->Fill(dt/tscale,E2);
					if(ploth==4 && opt==2) h2_handle->Fill(dt/tscale,E1);
					llfn_getgroupname(ev.groupno);
					printf("%8.2e %8.2lf %8.2lf %10d %10d %10.0lf\n",dt,E2,E1,i_log,i,ev.time);
// // // 					coin.dt = dt;
// // // 					coin.E2 = E2;
// // // 					coin.E1 = E1;
// // // 					ttree->Fill();
					//cout << dt << " " << E2 << " " << E1 << endl;
					cnt++;
					if(flagfinds){ //flag events
						if(flagfinds<0)flag[i] = abs(flagfinds);  
						if(flagfinds>0) flag[i_log] = abs(flagfinds);
					}	
					
					//if(E2>7000 && E2 < 7800) flag[i_log]=int('f');
					
				}
				//if(TMath::Abs(dt) < 1) cout << dt << " " << E2 << " " << E1 << endl;
				cnt_evnt++;
			}
		}
		if(cnt_evnt==n){
			event = 1;
			cnt_evnt = 0;
			i= i_log;
		}
	}
	if(ploth==1 || ploth==2){
		h_handle->Sumw2();
		h_handle->Draw();
	}
	else h2_handle->Draw("colz");
	
	printf( "Skipped groups:   %1.1lf%%\n",skipped_groups*100./ngroups);
	printf( "Skipped runss:    %1.1lf%%\n",skipped_runs*100./nruns);
	cout << "Number of events: " << cnt << endl;
	
// // // 	file->Write();
// // // 	//file->Close();
// // // 	ttree->Delete();	
	
	return cnt;
}


/// ########################################### sec:flag manipulating and reporting ########################################

void SScut_bursts( Int_t detn, Double_t sample_period = 500., Int_t nevents_in_sp = 20, Double_t Emin=150, Double_t Emax=4000, Double_t rest_period = 150., Int_t trig_filter = 0){
	int i=0, j, e_cnt = 0, i_log = 0;
	int det_min=detn-1;
	if(detn == 0) det_min = 0;		// so can take cut time away from uptime for different detectors
	
	Double_t time_start, last_burst_time = 0;
	time_start = rs[ev.runno].starttime; // this time should really be loaded from the start time which is save in the stats branch
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		//if(trig_filter == 0 || trig_filter == ev[i].trig_det){
		if(detn == 0 || (ev.energy[detn-1] > Emin && Emax == 0) || (ev.energy[detn-1] > Emin && ev.energy[detn-1] < Emax)){
			if( (ev.time - time_start) > sample_period ){
				if(e_cnt > nevents_in_sp){
					for(j=i_log,evtree->GetEvent(j);j<=i;j++,evtree->GetEvent(j)) if(trig_filter == 0 || trig_filter == ev.trig_det) flag[i] = 1;
					evtree->GetEvent(i); 
					last_burst_time = ev.time;
/// 					for(j=det_min;j<ndetectors;j++) gs[ev.groupno].flaggedtime[j] += Int_t(ev.time - time_start); // note how much time is cut out. NB I changed this from using trig_filter to detn to decide whether to count cut time on all of the detectors.
				}
				else if( (ev.time - last_burst_time) < rest_period) { // NB if the sample time is greater than the rest period this condition will never be used!
					for(j=i_log,evtree->GetEvent(j);j<=i;j++,evtree->GetEvent(j)) if(trig_filter == 0 || trig_filter == ev.trig_det)  flag[i] = 1; 
					evtree->GetEvent(i);
///					for(j=det_min;j<ndetectors;j++) gs[ev.groupno].flaggedtime[j] += Int_t(ev.time - time_start); 
				}
				e_cnt = 0;
				time_start = ev.time;
				i_log = i;
			}
			e_cnt++;
		}
	}
}

void SScut_datbursts(Int_t detn, Double_t sample_period=1., Int_t locut=0, Int_t hicut = 27, Double_t Emin=200, Double_t Emax=4000, Int_t verbose = 1){
// This function is as cut_burst but rejects events by hourly daq*dat basis and could report which daq*dat files have rejected. 
// if detn = 0 then all energies are used
	int i=0, j=0, e_cnt = 0, i_log=0, cnt=0;
	int det_min=detn-1;
	if(detn == 0) det_min = 0;		// so can take cut time away from uptime for different detectors
	sample_period *=3600.; 
	if(sample_period<3600) cout << "warning, this tool does not handle times less then 3600 seconds yet!!\n";
	
	evtree->GetEvent(i); Int_t time_start = rs[ev.runno].starttime, run_length = rs[ev.runno].runlen;
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		// skip flagged groups
		while(gs[ev.groupno].flag && i<nentries){
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
			time_start = rs[ev.runno].starttime;	
		}
		//if(time_start != rs[ev.runno].starttime) time_start != rs[ev.runno].starttime;
		if(detn == 0 || (ev.energy[detn-1] > Emin && Emax == 0) || (ev.energy[detn-1] > Emin && ev.energy[detn-1] < Emax)){ 
			if( (ev.time - time_start) >= sample_period){ 
				if(e_cnt < locut || e_cnt > hicut){
					for(j=0;j<int(sample_period/3600.);j++){
					cout << i << endl;
						i = i_log;
						evtree->GetEvent(i);
						rs[ev.runno].flag[detn-1] = 1;
							if(verbose==1){ cout << "Rejected run numbers: " << rs[ev.runno].starttime; verbose++;}  
							else if(verbose==2) cout << ", " << rs[ev.runno].starttime ;
						i+=rs[ev.runno].nevents;
					}
					cnt++;
				}
				else{
					for(i=i_log,evtree->GetEvent(i); rs[ev.runno].starttime==time_start ;i++,evtree->GetEvent(i)); // This will take us back to the run after the last starting point
					i_log=i; // the above for loop is the reason that this tool cannot handle sample times of less than 3600 seconds.  Mayb a bit inaficient to go back all the way just to find the beginning od the next run.
				}
				time_start = rs[ev.runno].starttime;  // (another reason that the tool does not handle times less than 3600 s.)
				run_length = rs[ev.runno].runlen;
				e_cnt = 0;
			}
			e_cnt++;
		}
	}
	cout << "number of rejected periods: " << cnt << endl;
	cout << "\ndone \n\n";
}
// "IT IS IMPORTANT TO GO THROUGH ALL THE SCENARIOUS IN ONE'S MIND"

void SScut_eph(Int_t detn, Int_t locut, Int_t hicut){
	unsigned int i;
	int j=0;  // i = epp counter, j = rpp counter
	for(i=0;i<epp.epp.size();i++){
		if(epp.runno[i]==0) break; // NB needed to ofset the counter so could use this way of spotting if had set epp value
		if(epp.epp[i] < locut || epp.epp[i] > hicut)
			for(j=0;j<gblperiod && (epp.runno[i]-1+j) < nruns;j++)
				rs[epp.runno[i]-1+j].flag[detn-1] = 1;
	}
}

void SScut_ephthresholds(Int_t detn, Int_t hicut, Int_t new_lothresh) {
	unsigned int i;
	Int_t locut = 0;
	int j=0;  // i = epp counter, j = rpp counter
	for(i=0;i<epp.epp.size();i++){
		if(epp.runno[i]==0) break; // NB needed to ofset the counter so could use this way of spotting if had set epp value
		if(epp.epp[i] < locut || epp.epp[i] > hicut) { // if the counts per hour for this run is not good...
			for(j=0;j<gblperiod && (epp.runno[i]-1+j) < nruns;j++){
				rs[epp.runno[i]+j-1].lo_thresh[detn-1] = new_lothresh; 
				cout << "=======> setting run " << epp.runno[i]+j << " to "<< new_lothresh << endl;
			}
		}
	}
}


void SScut_runsbytime(Int_t detn, Double_t itime, Double_t ftime, Int_t timescale = 1){
	int i;
	Int_t irun = 0, frun = 0;
	evtree->GetEvent(0); 
	Double_t time_start=rs[ev.runno].starttime;
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){  // search for start
		if( (ev.time-time_start) >= itime*timescale ){ 
			irun = ev.runno;
			break;
		}
	}
	for(;i<nentries ;i++,evtree->GetEvent(i)){  // search for end
		if( (ev.time-time_start) >  ftime*timescale ){ 
			frun = ev.runno;
			break;
		}
	}
	if(i==nentries) frun = ev.runno;
	if(frun!=0){ 
		for(i=irun;i<=frun;i++) rs[i].flag[detn-1] = 1;
		cout << "Flagged runs " << irun+1 << " to " << frun+1 << endl;
	}
	else cout << "No cuts made. (frun == 0)\n";
}

void SScut_hiRes_TimeAndEnergy(Int_t detn, Double_t itime, Double_t ftime, Int_t timescale = 1, Double_t Emin =-1, Double_t Emax =0){
	int i;
	evtree->GetEvent(0); 
	Double_t time_start=rs[ev.runno].starttime;
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){  // search for start
		if( (ev.time-time_start) > itime*timescale && (ev.time-time_start) <  ftime*timescale && (Emin==-1 || (ev.energy[detn-1] > Emin && ev.energy[detn-1]<Emax)))
			flag[i] = 1;
	}
}
	

void SScut_daynight(int hour_min, int hour_max, int day_min=0, int day_max=0){	
	struct tm tim;
	time_t unixtime;
	int i, j;
	Int_t oldrunno; 
	
	for(i=0,evtree->GetEvent(0);i<nentries;i++,evtree->GetEvent(i)){
		unixtime = time_t(ev.time);
		tim = *(localtime(&unixtime));
		if((tim.tm_hour >= hour_min && tim.tm_hour <= hour_max) && (day_min==0 || (tim.tm_wday >= day_min && tim.tm_wday <= day_max))){
			for(j=0;j<NDETS;j++) rs[ev.runno].flag[j] = 1;
			for(oldrunno=rs[ev.runno].starttime;oldrunno==rs[ev.runno].starttime && i<nentries;i++,evtree->GetEvent(i)); //forward to next run
		}
	}
	
//	struct tm {
///* ANSI standard fields */
//   int tm_sec;   /* 0 to 60 */
//   int tm_min;   /* 0 to 59 */
//   int tm_hour;  /* 0 to 23 */
//   int tm_mday;  /* 1 to 31 */
//   int tm_mon;   /* 0 to 11 */
//   int tm_year;  /* year - 1900 */
//   int tm_wday;  /* Sunday = 0 */
//   int tm_yday;  /* 0 to 365 */
//   int tm_isdst;
//       /* >0 if Daylight Savings Time,
//        *  0 if Standard,
//        * <0 if unknown */
///* extensions to ANSI standard */
//   char *tm_zone;  /* time zone name    */
//   long tm_gmtoff; /* offset from GMT   */
//};

}

void SScut_graphical_eVSe(Int_t detx=gbldetx, Int_t dety=gbldety, Char_t use_flags[2] = "__", Char_t use_lo_thresh[2] = "00"){
//  use_flags[2] = "__": use flags on detx or dety respectively eg "XY"

	int i;
	Double_t p1[2], p2[3], p3[2], p4[2];
	Double_t ex, ey, a1[2], a2[2], a3[3], a4[2];
	
	TCutG *Cut = (TCutG*)gPad->FindObject("CUTG");
	Cut->GetPoint(0,p1[0],p1[1]);
	Cut->GetPoint(1,p2[0],p2[1]);
	Cut->GetPoint(2,p3[0],p3[1]);
	Cut->GetPoint(3,p4[0],p4[1]);
	
	
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		// skip flagged groups
		while((gs[ev.groupno].flag || !gs[ev.groupno].det_on[detx-1] || !gs[ev.groupno].det_on[dety-1]) && i<nentries){
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
		}
		// skip flagged run: Skip flagged data if (flagegd && 0) || (!flagged && !0)
		while( (!rs[ev.runno].det_on || 
				(use_flags[0]!='_' && rs[ev.runno].flag[detx-1]==(use_flags[0]=='0')) || 
				(use_flags[1]!='_' && rs[ev.runno].flag[dety-1]==(use_flags[1]=='0')) 
				) && i<nentries){
			i += rs[ev.runno].nevents;
			evtree->GetEvent(i);
		}	
		ex=ev.energy[detx-1];
		ey=ev.energy[dety-1];
		if( p1[0]!=0 &&
			(use_lo_thresh[0]=='0' || (use_lo_thresh[0]=='1' && ev.energy[detx-1] > rs[ev.runno].lo_thresh[detx-1])) &&
			(use_lo_thresh[1]=='0' || (use_lo_thresh[1]=='1' && ev.energy[dety-1] > rs[ev.runno].lo_thresh[dety-1]))
			){
			llfn_poly1(a1,p1[0],p1[1],p2[0],p2[1]);
			llfn_poly1(a2,p2[0],p2[1],p3[0],p3[1]);
			llfn_poly1(a3,p3[0],p3[1],p4[0],p4[1]);
			llfn_poly1(a4,p4[0],p4[1],p1[0],p1[1]);
			if (ey <= a1[1]*ex+a1[0] && ex <= ey/a2[1]-a2[0]/a2[1] && ey >= a3[1]*ex+a3[0] && ex >= ey/a4[1]-a4[0]/a4[1]) flag[i]=1;
		}
		else if(p1[0]==0) cout << "Error: Cuts undefined\n";
	}	
}

void SScut_graphical_eVSt(Int_t dety=gbldety, Int_t time_scale=gblperiod, Char_t use_flags = '_'){
	int i;
	Double_t p1[2], p2[3], p3[2], p4[2];
	Double_t t0, tx, ey, a1[2], a2[2], a3[3], a4[2];
	
	TCutG *Cut = (TCutG*)gPad->FindObject("CUTG");
	Cut->GetPoint(0,p1[0],p1[1]);
	Cut->GetPoint(1,p2[0],p2[1]);
	Cut->GetPoint(2,p3[0],p3[1]);
	Cut->GetPoint(3,p4[0],p4[1]);
	
	evtree->GetEvent(0);
	t0 = rs[ev.runno].starttime;
	
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)){
		// skip flagged groups
		while((gs[ev.groupno].flag || !gs[ev.groupno].det_on[dety-1]) && i<nentries){
			i += gs[ev.groupno].nevents;
			evtree->GetEvent(i);
		}
		// skip flagged run: Skip flagged data if (flagegd && 0) || (!flagged && !0)
		while( (!rs[ev.runno].det_on || (use_flags!='_' && rs[ev.runno].flag[dety-1]==(use_flags=='0'))) && i<nentries){
			i += rs[ev.runno].nevents;
			evtree->GetEvent(i);
		}	
		ey=ev.energy[dety-1];
		tx=(ev.time-t0)/time_scale;
		if(p1[0]!=0){
			llfn_poly1(a1,p1[0],p1[1],p2[0],p2[1]);
			llfn_poly1(a2,p2[0],p2[1],p3[0],p3[1]);
			llfn_poly1(a3,p3[0],p3[1],p4[0],p4[1]);
			llfn_poly1(a4,p4[0],p4[1],p1[0],p1[1]);
			if (ey <= a1[1]*tx+a1[0] && tx <= ey/a2[1]-a2[0]/a2[1] && ey >= a3[1]*tx+a3[0] && tx >= ey/a4[1]-a4[0]/a4[1]) flag[i]=1;
		}
		else cout << "Error: Cuts undefined\n";
	}	
}


double SScount(int detn, double Emin, double Emax, int normalise=3600*24, char use_flags[4] = "000"){ 
// count number of events within a specified energy window. Does not include flagged events.

// 	int default_dets[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0};
// 	if(!dets) dets = default_dets;

	
	int i, j, j_start;
	int counts=0;
		
	
	for(i=0,j=0,j_start=0,evtree->GetEvent(j);i<nruns;++i) { //loop over runs
		
		// if run or group is flagged, skip to next run (set j to i+=event in run
		// skip run if...
		if( !gs[rs[i].groupno].det_on[detn-1] || 														// skip off detectors for whole group
			gs[rs[i].groupno].flag || 																// skip flagged groups
			!rs[i].det_on[detn-1] || 																			// skip off detectors
			(use_flags[1]!='_' && rs[i].flag[detn-1]==(use_flags[1]=='0'))	||							// skip flagged runs
			(rs[i].lo_thresh[detn-1]>Emin || gs[rs[i].groupno].hi_thresh[detn-1]<Emax)
			) {	
			
			j = j_start+rs[i].nevents;
			j_start = j;
			// maybe check that it is still in sync (ie that the previous event if of the previous run
			// idealy it would go like this: j = rs[i+1].evnumber which would take it to the first event of the i+1'th run
			continue;
		}
		
		for(evtree->GetEvent(j);j<j_start+rs[i].nevents;++j,evtree->GetEvent(j)) { 						// loop over events in run  (could have been j<rs[i+1].evnumber)
			// check if event is flagged
			if( Emin < ev.energy[detn-1] &&										// Fill if above lower run threshold
				Emax > ev.energy[detn-1] &&					// Fill if below upper group threshold (so far there is no upper run threshold as it is unlikely to be used)
				(use_flags[0]=='_' || !flag[j]==(use_flags[0]=='0')) &&									// Fill if using evflags and is not evflagged
				(use_flags[2]=='_' || !llfn_Decimal2Bit(BFlag[j],detn,NDETS)==(use_flags[2]=='0')) ) {  // Fill if using evBitFlagged and is not evBitFlagged (evBitFlag is a binary number which is decomposed into a flag for each detector)
			
				// count events
				++counts;
			}
		}
			
		j = j_start+rs[i].nevents;
		j_start = j;
	}
	
	if(normalise) return double(counts)/SSgetusedtimeR(detn,Emin,Emax)*normalise;
	else return counts;

}

double SScount(int *dets, double Emin, double Emax, char use_flags[4] = "000"){ 
// count number of events within a specified energy window. Does not include flagged events.

	int i;
	double total=0;
	
	for(i=0;dets[i];++i) total+= SScount(dets[i],Emin,Emax,0,use_flags);
	
	return total;
}
	


void SSexit_batch() {
	gROOT->SetBatch(kFALSE);
}
	
void SSclean_data(int detn) {
// This function removes noisy runs from the data and locates the noise threshold

	int thresholds[]={8000,500,300,250,200,150,120,100,90,80,70,60,50,40,30,20,0};
	int i, j;
	int undercut = 10; // keV  this is now used for the sliding energy window at low energies
	double safety_factor = 1.2; // this sets the level of the threshold above the tested region
	
	gROOT->SetBatch(kTRUE);
	
	for(i=0;thresholds[i+1];++i) {
				
		cout << "\n\ncleaning range: " << thresholds[i+1]-undercut << "keV to " << int(thresholds[i]) << "keV\n\n";
		SSplot_eph(detn,1,"000",100,thresholds[i+1]-undercut,thresholds[i]);
		
		// check that plot is not empty. If it is all of the thresholds are below thresholds[i+1]
		// remove under cut, and if it fails again or on the next i loop then end off
		while(epp.epp.size()==0) {
			if(undercut > 0 ) {
				undercut = 0;
				SSplot_eph(detn,1,"000",100,thresholds[i+1]-undercut,thresholds[i]);
			}
			else {
				cout << "=======> Setting baseline to " << int(safety_factor*thresholds[i]) << endl;
				for(j=0;j<nruns;++j) {
					// check if threshold falls within range and push up to upper +10keV if it does
					if(rs[j].lo_thresh[detn-1]>thresholds[i+1] && rs[j].lo_thresh[detn-1]<thresholds[i]) {
						rs[j].lo_thresh[detn-1] = safety_factor*thresholds[i];
					}
				}
				gROOT->SetBatch(kFALSE);
				return;
			}
		}
		
		cout << "fitting poisson...\n";
		SSfit_poisson();
		
		// reject full run if falls outside Poisson fit above 500 keV
		if(i==0) {
			SScut_eph(detn,0,fit_hicut);
			continue;
		}
		
		// refit only good data
		h_handle->GetXaxis()->SetRange(0,fit_hicut+1);
		SSfit_poisson();
// 		
		// check if fit is good
		if(TMath::Prob(gblfit->GetChisquare(),gblfit->GetNDF())>0.001) {
			// check if any runs had thresholds in the middle of the energy range and push to the top if there is
			cout << "Checking for midway thresholds... \n";
			for(j=0;j<nruns;++j) {
// 				cout << "Is " << rs[j].lo_thresh[detn-1] << " greater than " << thresholds[i+1]-undercut << "?  AND " << rs[j].lo_thresh[detn-1] << " < " << thresholds[i] << "?\n";
				if(rs[j].lo_thresh[detn-1]>(thresholds[i+1]-undercut) && rs[j].lo_thresh[detn-1]<thresholds[i]){
					rs[j].lo_thresh[detn-1] = thresholds[i];
					cout << "=======> Setting midway run " << j+1 << " to " << thresholds[i]<< endl;
				}
			}
			//fit is good so set thresholds
			cout << "Setting noisy run's thresholds to " << safety_factor*thresholds[i]<< endl;
			SScut_ephthresholds(detn,fit_hicut,int(safety_factor*thresholds[i]));		
		}
		else {
			cout << "=======> Fit not good. Setting runs to " << int(safety_factor*thresholds[i]) << endl;
			for(j=0;j<nruns;++j) {
				if(rs[j].lo_thresh[detn-1]<thresholds[i]) rs[j].lo_thresh[detn-1] = safety_factor*thresholds[i];
			}
			gROOT->SetBatch(kFALSE);
			return;
		}
	}
// 	gbl_draw = old_gbl_draw;
	gROOT->SetBatch(kFALSE);
	SSexit_batch();
}

/*
void llfn_cleancuts(){
	int i;
	for(i=0;i<NDETS;i++) gblcuts[i].clear();
} */

void SSclean_flags(){
	int i;
	for(i=0,evtree->GetEvent(i); i<nentries ;i++,evtree->GetEvent(i)) {flag[i] = 0; BFlag[i]=0;};
// 	llfn_cleancuts();
}

void SSclean_runflags(){
	// reloads original settings. Have to reload to keep original skipped run flags;
	int i, j;	
	TTree *rstree = (TTree*)ssfile->Get("rstats");
	run_stats rs_;
	rstree->SetBranchAddress("rs_branch",&rs_);
	for(i=0,rstree->GetEvent(i);i<nruns;i++,rstree->GetEvent(i)) for(j=0;j<NDETS;j++) rs[i].flag[j] = rs_.flag[j];
}

// void SSclean_runthresholds() {
// 	int i, j;
// 	h2_runthresholds = (TH2I*)ssfile->Get("run_thresholds");
// 	if(!h2_runthresholds) {
// 		cout << "Thresholds not found so calculating from threshold numbers...\n";
// 		h2_runthresholds = new TH2I("run_thresholds","run_thresholds",nruns,1,nruns+1,NDETS,1,NDETS+1);
// 		for(i=0;i<nruns;++i) for(j=0;j<NDETS;++j) h2_runthresholds->SetBinContent(i+1,j+1,SSch2e(SSth2ch.Eval(rs[i].thr[j]),j+1));
// 	}
// 	else cout << "Loaded thresholds.\n";
// }

// void SSreset_runthresholds(Int_t detn = 0) {
// 	int i, j;
// 
// 	if(detn==0) for(i=0;i<nruns;++i) for(j=0;j<NDETS;++j) h2_runthresholds->SetBinContent(i+1,j+1,SSch2e(SSth2ch.Eval(rs[i].thr[j]),j+1));
// 	else for(i=0;i<nruns;++i) h2_runthresholds->SetBinContent(i+1,detn,SSch2e(SSth2ch.Eval(rs[i].thr[detn-1]),detn));
// }

void SSgroups_setHelp() {
	cout << "This function sets group flags and can be used to plot one group against another\n"
	     << "or omit certain groups\n\n"
		 << "Examples: SSgroups_set(\"+2:3\") sellects groups 2 - 3\n"
		 << "          SSgroups_set(\"-4:4\") omits group 4\n"
		 << "          SSgroups_set(\" \")    remove all flags\n";
}
void SSgroups_set(Char_t use_groups[20]){
	// select which groups to use
	int i, j, group_range_sign=0, igroup=0, fgroup=0;
	char buffer[20];
	if(use_groups[0]!=' '){
		if(use_groups[0]=='+') group_range_sign = 1; 
		else if(use_groups[0]=='-') group_range_sign = -1;
		else {cout << "group_range_error, exiting...\n"; return;}
		for(i=1;use_groups[i]!=':';i++) buffer[i-1]=use_groups[i]; buffer[i-1]='\0'; igroup=atoi(buffer);
		for(j=0,i++;use_groups[i]!='\0';i++,j++) buffer[j]=use_groups[i]; buffer[j+1-1]='\0'; fgroup=atoi(buffer);
		cout << "i = " << igroup << " f = " << fgroup << endl;
		for(i=0;i<ngroups;i++){
			if(group_range_sign==1){ // set flags to use selection
				if(i < igroup-1 || i > fgroup-1) gs[i].flag = 1;
				else gs[i].flag = 0;
			}
			else{					// set flags to use outside selection
				if(i < igroup-1 || i > fgroup-1) gs[i].flag = 0;
				else gs[i].flag = 1;
			}
		}
	}
	else for(i=0;i<ngroups;i++) gs[i].flag=0;
// 	else{ // add these settings...
		
	
}

Int_t SSgroups_get(){
	int i, j;
	for(i=0,j=0;i<ngroups;i++){ 
		if(gs[i].flag==0){
			printf("%3d. ",i+1);
			llfn_getgroupname(i);
			printf(" lo_thresh: %3d,  hi_thresh %5d\n",int(gs[i].lo_thresh[0]),int(gs[i].hi_thresh[0]));
			j++;
		}
	}
	return j;
}

// this function does not work! 
void SSgetdatathresh(Int_t setnonzero = 0, Float_t addenergy = 5.){
// this runction finds the lowest energy greater than zero and adds 'addenergy' to estimate where the trigger theshold was.
	int i, j;
	Double_t lowestE[NDETS];
	Int_t group=-1;
	for(i=0,evtree->GetEvent(0);i<nentries;i++,evtree->GetEvent(i)){
		if(ev.groupno!=group){
			for(j=0;j<NDETS;j++){
				if(i!=0) evtree->GetEvent(i-1);
				if(ev.groupno!=0 && gs[ev.groupno].det_on[j] && (setnonzero || gs[ev.groupno].lo_thresh[j] == 0)) gs[ev.groupno].lo_thresh[j] = lowestE[j] + addenergy; // set the threshold counted in last group
				evtree->GetEvent(i);
				if(gs[ev.groupno].det_on[j] && ev.energy[j]>0) lowestE[j] = ev.energy[j]; // set lowest energy as first energy of group
			}
			group = ev.groupno;
		}
		for(j=0;j<NDETS;j++) if(ev.energy[j] > 0 && ev.energy[j] < lowestE[j]) lowestE[j] = ev.energy[j]; // readjust lowesstE if energy is less...
	}
	for(j=0;j<NDETS;j++) if(setnonzero || gs[ev.groupno].lo_thresh[j] == 0) gs[ev.groupno].lo_thresh[j] = lowestE[j] + addenergy; // set the lowest thresholds recorded by last group
}

void SSaddgroupmarkers(Float_t tscale=3600., Int_t Y=1){
	int i, j;
	evtree->GetEvent(0);
	Double_t starttime=rs[ev.runno].starttime;
	Double_t marker[ngroups], possition[ngroups];
	for(i=0,j=0;i<ngroups;i++){
		evtree->GetEvent(j);
		possition[i] = (rs[ev.runno].starttime - starttime +1)/tscale;
		marker[i] = Y;
		j+=gs[i].nevents;
	}
	TGraph *addmark = new TGraph(ngroups,possition,marker);
	addmark->SetMarkerColor(4);
	addmark->Draw("*");
}

///                                    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
void SSset_activedets(Char_t onoff[100]="0  0  1  0  0  0  0  1  0  0  1  1  0  0  1  0", Int_t group=0){
	int i;
	for(i=0;i<NDETS;i++) gs[group].det_on[i] = atoi(&onoff[i*3]);
}


/// ######################################### sec:auto functions ###########################################################

void SSauto_plotwaterfall(char nameID[50] = "_", Int_t markerstyle = 6, Int_t time_scale = 3600){
	int i;
	char buffer[50];
	
// Double_t E_hi, E_lo, th_hi, th_lo, t_hi;
// 
// TF1 *f1;
// TGaxis *A1; 
	
	
	gblpause = 1;
	for(i=1;i<=NDETS;i++){
// th_hi = 127.;
// th_lo = SSch2th.Eval(0);
// E_lo = gs[0].cal_eqn[0][i-1]*SSth2ch.Eval(0)+gs[0].cal_eqn[2][i-1];
// E_hi = gs[0].cal_eqn[0][i-1]*SSth2ch.Eval(th_hi)+gs[0].cal_eqn[2][i-1];
// f1 = new TF1("f1","x",th_lo,th_hi);
// 		
	SSplot_waterfall(i,markerstyle,time_scale);
// t_hi = g_handle->GetXaxis()->GetBinUpEdge(g_handle->GetXaxis()->GetLast());
// A1 = new TGaxis(t_hi,E_lo,t_hi,E_hi,"f1",1270,"+L");
// A1->Draw();
		gPad->WaitPrimitive(); 
		sprintf(buffer,"wf%s%d.gif",nameID,i);
		c1->Print(buffer);		
	}
	gblpause = 0;
}

void SSplot_highlightcorrelations(Int_t detx, Int_t dety, Int_t markerstyle=6, Int_t use_flags=1){
	SSplot_correlations(detx,dety,0,0,1,markerstyle);
	SSplot_correlations(detx,dety,0,use_flags,2,markerstyle);
	c1->ToggleToolBar();
}

void SSauto_ploteph(Int_t detn, Int_t period, Double_t Emin=500, Double_t Emax=4000){
	SSplot_eph(detn,period,"00",15.,Emin,Emax,1,1,0,1,0);
	SSplot_eph(detn,period,"0_",15.,Emin,Emax,1,1,1,1,0);
	SSplot_eph(detn,period,"00",15.,Emin,Emax,1,1,0,1,1);
	SSfit_poisson(-1,h_handle2);
	h_handle2->Draw();
	h_handle->Draw("same");
	
	c1->cd(2);
	g_handle2->Draw("A*");
	g_handle->Draw("*");
	SSfit_poly();
	gPad->WaitPrimitive();
}

void SSauto_plot2eph(Int_t detn, Int_t period, Double_t Emin=500, Double_t Emax=4000){
	char buffer[100];
	cout << "running plot_eph\n";
	SSplot_eph(detn,period,"00",20.,Emin,Emax,1,1,0,1,0);
	cout << "running fit_poisson\n";
	SSfit_poisson();
	cout << "running fit_poly\n";
	SSfit_poly();
	sprintf(buffer,"epp.lo.d%d.eps",detn);
	SSprint(1,buffer);
}
	
void SSauto_eph_fitandcut(Int_t detn, Int_t period, Double_t Emin, Double_t Emax, Double_t range_Xmax = 0, Int_t pause = 0){
	if(range_Xmax == 0) range_Xmax = 50. * period;
	SSplot_eph(detn,period,"000",range_Xmax,Emin,Emax,0,1,0,1,pause);
	cout << "plotted eph.\n";
	SSfit_poisson(-1,h_handle);
	cout << "fitted poisson. \n";
// 	cut_eph(detn,fit_locut,fit_hicut);
	SScut_eph(detn,0,fit_hicut);
	cout << "made cuts.\n";
	printf("Runs  Available: %3.0lf\n",SSgetuptime(detn)/3600.);
	printf("Runs  Used     : %3.0lf\n",(SSgetuptime(detn)-SSgetflaggedtime(detn))/3600.);
	printf("Percentage Used: %2.0lf%%\n",(SSgetuptime(detn)-SSgetflaggedtime(detn))/SSgetuptime(detn)*100.);
}

void SSalldet_auto_eph_fitandcut(Double_t Emin, Double_t Emax, Double_t range_Xmax = 20, Int_t pause =1, Int_t period =1){
	int i, j;
	c1 = new TCanvas("c1", "c1",14,30,700,500);
	for(i=0;i<ngroups;i++){ 
		for(j=0;j<NDETS;j++){ 
			if(gs[i].det_on[j]){ 
				SSauto_eph_fitandcut(j+1,period,Emin,Emax,range_Xmax,pause); 
				gPad->WaitPrimitive();
			}
		}
	}
}

void SSalldet_auto_ploteph(Int_t period =1){
	int i, j;
	c1 = new TCanvas("c1", "c1",14,30,700,500);
	for(i=0;i<ngroups;i++){ 
		for(j=0;j<NDETS;j++){ 
			if(gs[i].det_on[j]){ 
				SSauto_ploteph(j+1,period);
				gPad->WaitPrimitive();
			}
		}
	}
}	

void SSallgroups_fitandcut(Int_t detn, Int_t period, Double_t Emin, Double_t Emax, Double_t range_Xmax, Int_t pause = 0){
	int i;
	for(i=0;i<ngroups;i++) gs[i].flag = 1;
	for(i=0;i<ngroups;i++){
		gs[i].flag = 0;
		SSauto_eph_fitandcut(detn,period,Emin,Emax,range_Xmax,pause);
		cout << "done run " << i+1 << endl;
		gPad->WaitPrimitive();
		gs[i].flag = 1;
	}
	SSgroups_set();
}

void SSauto_clean(Int_t detn, Char_t file[64]="cuts.dat", Int_t pause =0){
// cleans the data. In the future could have input file with any list of energy ranges and periods to be used in cuts.
	int i;
// --------------------------get cut information ------------------------------------
	FILE * ifp;
	if(!(ifp = fopen(file,"r"))){
		cout << "Error opening file " << file << endl;
		return;
	}
	int ch = ' ', entries = 0;
	while (ch != EOF){
		ch=getc(ifp); 
		if(ch == '\n') entries++;
	}
	entries--;
	cout << "entries = " << entries << endl; 
	rewind(ifp); 
	ch=' ';
	while ( ch != EOF && (ch=getc(ifp)) != '\n' );  // goto 1st #
	
	Int_t *period = new Int_t[entries];
	Double_t *Emin = new Double_t[entries];
	Double_t *Emax = new Double_t[entries];
	
	for(i=0;i<entries;i++){
		fscanf(ifp,"%d",&period[i]); 
		fscanf(ifp,"%lf",&Emin[i]); 
		fscanf(ifp,"%lf",&Emax[i]);
		while ( ch != EOF && (ch=getc(ifp)) != '\n' ); // forward to end
	}
	cout << "Entries recieved: \n";
	for(i=0;i<entries;i++) cout << period[i] << " " << Emin[i] << " " << Emax[i] << endl;	
	

	if(cs) delete cs;
	cs = new CutStats[entries+1];
	cs[0].availruns=Int_t(SSgetusedtime(1,1,0)/3600.); cs[0].extracut=0; cs[0].mean=0.; cs[0].error=0.; cs[0].probX2=0;
	for(i=0;i<entries;i++){
		SSauto_eph_fitandcut(detn,period[i],Emin[i],Emax[i],100,pause);
		cs[i+1].availruns=Int_t(SSgetusedtime(detn,1,0)/3600.); 
		cs[i+1].extracut=cs[i].availruns-cs[i+1].availruns; 
		cs[i+1].mean=gblfit->GetParameter(0); 
		cs[i+1].error=gblfit->GetParError(0); 
		cs[i+1].probX2=TMath::Prob(gblfit->GetChisquare(),gblfit->GetNDF());
	}
	
	delete period;
	delete Emin;
	delete Emax;
	fclose(ifp);
	
}




/// ############################################ fitting tools #########################################################
 

void SSfit_poisson(Double_t mean, TH1D *hist1, Int_t min, Int_t max){
	int i;
	
	if(mean==-1) mean = hist1->GetMaximumBin();
// 	if(gblfit) gblfit->Delete();
	gblfit = new TF1("gblfit","[1]*TMath::PoissonI(x,[0])",min,max);
	gblfit->SetParName(0,"mean");
	gblfit->SetParName(1,"scale");
	gblfit->SetParameter(0,mean);
	gblfit->SetParameter(1,1);
	gblfit->SetLineColor(4);
	gblfit->SetLineWidth(2);
	hist1->Fit("gblfit");
	
	cout << "   2  ProbChi      " << TMath::Prob(gblfit->GetChisquare(),gblfit->GetNDF()) << "\n";
	
	Float_t sum;
	cout << "Suggested 99% lo and hi cuts: ";
	for(i=max,sum=0;sum<0.99;i--) sum+=TMath::PoissonI(i,gblfit->GetParameter(0)); fit_locut = i+=1;
	cout << fit_locut << " (" << sum*100 << "%) ";	
	for(i=0,sum=0;sum<0.99;i++) sum+=TMath::PoissonI(i,gblfit->GetParameter(0)); fit_hicut = i-=1;
	cout << fit_hicut << " (" << sum*100 << "%) ";
	for(i=fit_locut,sum=0;i<=fit_hicut;i++) sum+=TMath::PoissonI(i,gblfit->GetParameter(0));
	cout << " total = " << sum*100 << "%\n";
	
	// This is correct, I checked! Note the cut is made ABOVE fit_hicut and the integral is made from 0 to fit_hicut+1. 

	if(gblfit->GetParameter(0) <=0.1) cout << "WARNING, THE MEAN IS VERY LOW...\n";
}

// void fit_gaus(char help='H'){
// 	cout << "use mean = -1 to use greatest value in histogram\n";
// 	cout << "use sigma < 0 to use deltaE, (sigma = mean/deltaE)";
// 	cout << "use sigma < 0 to use value of max entry
// 	help=help;
// }

void SSfit_gaus(Double_t mean, Double_t sigma, Double_t scale, TH1D *hist1){ 
// 	if(gblfit) gblfit->Delete();
	
	if(mean==-1) mean = hist1->GetMaximumBin();
	if(sigma<0) sigma = -1.*mean*sigma/2.35;
	if(scale<0) scale = h_handle->GetBinContent(Int_t(mean));
	
	gblfit = new TF1("gblfit","[2]*TMath::Gaus(x,[0],[1])",-.5,500.);
	gblfit->SetParName(0,"mean");
	gblfit->SetParName(1,"sigma");
	gblfit->SetParName(2,"scale");
	gblfit->SetParameter(0,mean);
	gblfit->SetParameter(1,sigma);
	gblfit->SetParameter(2,scale);
	gblfit->SetLineColor(4);
	hist1->Fit("gblfit");
	
	Double_t p_sigma = TMath::Abs(gblfit->GetParameter(1));
	cout << "   3  sigma           " << p_sigma << endl;
	cout << "   4  mean + 3*sigma  " << gblfit->GetParameter(0)+p_sigma*3 << endl;
	cout << "   5  ProbChi      " << TMath::Prob(gblfit->GetChisquare(),gblfit->GetNDF()) << "\n";
	cout << "Integral(0->3sigma)*100/Entries = " << h_handle->Integral(1,Int_t(gblfit->GetParameter(0)+p_sigma*3)+2)*100/ h_handle->GetEntries()<< "\n\n";
	
	gblresults.clear();
	gblresults.push_back(gblfit->GetParameter(0));
	gblresults.push_back(p_sigma);
}

void SSfit_poly(){
	gblfit = new TF1("gblfit","[0]+[1]*x",0,1000);
	gblfit->SetParName(0,"Constant");
	gblfit->SetParName(1,"Slope");
	gblfit->SetParameter(0,0);
	gblfit->SetParameter(0,0.5);
	gblfit->SetLineColor(4);
	c1->cd(2);
	g_handle->Fit("gblfit");
}

void SSfit_exp(TH1D *th1_handle, Double_t slope=4., Double_t amplification = 4., Double_t offset=0){
// changed this function to inlcude an offset (23/02/08)
	if(offset) {
		gblfit = new TF1("gblfit","[2]+TMath::Exp([0]+[1]*x)",0,1000);
		gblfit->SetParName(2,"Offset");
		gblfit->SetParameter(2,offset);
	}
	else gblfit = new TF1("gblfit","TMath::Exp([0]+[1]*x)",0,1000);
	
	gblfit->SetParName(0,"Constant");
	gblfit->SetParName(1,"Slope");
	gblfit->SetParameter(0,amplification);
	gblfit->SetParameter(1,slope);
	
	gblfit->SetLineColor(4);
	gblfit->SetLineWidth(1);
	th1_handle->Fit("gblfit");
}

void llfn_poly1(Double_t *a, Double_t x1, Double_t y1, Double_t x2, Double_t y2){
// Converts two points to constants of a polynomial
	a[1]=(y2-y1)/(x2-x1);
	a[0]=y1-a[1]*x1;
}


/// ########################################### low level functions ####################################################

FILE *runlistp, *runnamep;
char runname[564];


void SSload(Char_t *filename){
	fname=filename;
	Int_t nbytes = 0;
	int i, j;

	
	// Import data -------------------------------------
	ssfile = new TFile(filename);
	if(!ssfile) {
		cout << "File: " << filename << " does not exist. \n";
		return;
	}
	evtree = (TTree*)ssfile->Get("events");
	//event evs; //struct event ev;
	evtree->SetBranchAddress("Event",&ev);
	
	evtree->GetEvent(0);
	//ev=evs;

	nentries = (Long64_t)evtree->GetEntries(); //nentries = (Long64_t)ev_tree->GetEntries();
	cout << "Allocating memory for flags...\n";
	flag = new Int_t[nentries];  //This will be a larg vector. Should I make something Long64?
	cout << "Reading in flags..." << endl;
	for (i = 0 ; i < nentries ; i++ ){
		evtree->GetEvent(i);
		flag[i] = ev.flag;
	}
	cout << ndetectors << " detectors found\n";
	
	// Import runflags ---------------------------------------------------
	TTree *rstree = (TTree*)ssfile->Get("rstats");
	run_stats rs_;
	rstree->SetBranchAddress("rs_branch",&rs_);
	nruns = (Int_t)rstree->GetEntries();
	rs = new run_stats[nruns]; 
	if(rs == NULL){
		cerr << "Out of Memory!" << endl;
		return;	// exit
	}
	cout << nruns << " runs found\n";
	cout << "Reading in run stats..." << endl;
	for (i = 0 ; i < nruns ; i++ ){
		nbytes += rstree->GetEvent(i);
		rs[i]=rs_;
		
		// [Oli] Bugfix for struct elements memory alignment problem:
		/// size_t mvCount = sizeof(run_stats) + size_t(&rs[i]) - size_t(&rs[i].lo_thresh[0]);
		/// memmove(&rs[i].lo_thresh[0], &rs[i].runlen+1, mvCount);		
	}
		
	// Import group stats -----------------------------------------
	TTree *gstree = (TTree*)ssfile->Get("gstats");
	group_stats gs_;
	gstree->SetBranchAddress("gs_branch",&gs_);
	ngroups = (Int_t)gstree->GetEntries();
	gs = new group_stats[ngroups]; 
	if(gs == NULL){
		cerr << "Out of Memory!" << endl;
		return;	// exit
	}
	cout << ngroups << " group(s) found\n";
	cout << "Reading in group stats..." << endl;
	for (i = 0 ; i < ngroups ; i++ ){
		nbytes += gstree->GetEvent(i);
		gs[i]=gs_;
		
		// [Oli] Bugfix for struct elements memory alignment problem:
/// 		size_t mvCount = sizeof(group_stats) + size_t(&gs[i]) - size_t(&gs[i].hi_thresh[0]);
/// 		memmove(&gs[i].hi_thresh[0], &gs[i].det_on[NDETS], mvCount);
		
		cout << i+1 << ". groupname: "; 
		llfn_getgroupname(i);
// 		cout << "   lo_thresh: " << gs[i].lo_thresh[0] << "  hi_thresh: "  << gs[i].hi_thresh[0] << " nevents: " << double(gs[i].nevents) << endl;
		printf("  hi_thresh %5d,  nruns: %2.0lf\n",int(gs[i].hi_thresh[0]),double(gs[i].nruns));
		
	}
	
	nbytes = int(nbytes*2.7); 
	nbytes += nentries*4;
	printf("\nData size loaded to RAM: %1.1lf MB\n",double(nbytes)/1.E6);

	//f->Close(); does not work if close
	
	// set rs.lo_thresh if not already set
	if(rs[0].lo_thresh[0]==-1) {
		cout << "Thresholds not set so calculating from threshold numbers...\n";
		for(i=0;i<nruns;++i) {
			for(j=0;j<NDETS;++j) {
					// the lower threshold can be set higher by the user so have to check if the run threshold is above the user thresholds (lo_thresh)
				if(SSch2e(SSth2ch.Eval(rs[i].thr[j]),j+1) > gs[rs[i].groupno].lo_thresh[j]) rs[i].lo_thresh[j]=SSch2e(SSth2ch.Eval(rs[i].thr[j]),j+1);
				else rs[i].lo_thresh[j]=gs[rs[i].groupno].lo_thresh[j];
			}
		}
	}
	
	cout << "creating binary flags...\n";
	BFlag = new Int_t[nentries];
	for(i=0;i<nentries;i++) BFlag[i] = 0;
	
	cout << "setting Figure options\n";
	llfn_setstyles(" ");
	
	
	cout << "done\n\n";
	//system("sleep 10"); 
}

void SSflag_noise_help() {
	cout << "This function flags all of the events that do not have an energy\n";
	cout << "in any active detector above the lo_thresh of its runan";
}

void SSflag_noise(int flagVal = int('n')) {
	int i, j;
	int cnt;
	for(i=0,evtree->GetEvent(0);i<nentries;i++,evtree->GetEvent(i)) {
		cnt = 0;
		for(j=0;j<ndetectors;j++) {
			if( gs[ev.groupno].det_on[j] &&
				!gs[ev.groupno].flag &&
				rs[ev.runno].det_on[j] &&
				!rs[ev.runno].flag[j]
				) {
				if(ev.energy[j] > rs[ev.runno].lo_thresh[j]) {
					cnt++;
					break;
				}
			}
		}
		if(cnt==0) flag[i] = flagVal;
	}
}


void SSsavechanges(char help = 'H') {
	help = help;
	cout << "This functions saves the data to disk,\nincluding modifications to the flags\n\n";
	cout << "Set 'reject_cuts' to 0 to save all data\n";
	cout << "Set 'reject_cuts' to -1 to set energies of flagged events to -3. NB, the event is then unflagged\n";
	cout << "Set 'reject_cuts' to 1 remove flagged events from data file.\n";
	cout << "Set 'reject_cuts' to # to remove flagged events of value # from data file. (#>1) \n";
}


void SSsavechanges(Int_t reject_cuts , Char_t *filename=fname) {
	int i, j;
	Char_t buffer[500];
	Char_t temp_fname[100];
	sprintf(buffer,"mv %s %s~", filename,filename);
	system(buffer);
	sprintf(temp_fname,"temp_%s",filename);
	TFile *hfile = new TFile(temp_fname,"RECREATE","datafile",9);
	
	event ev_;
	run_stats rs_;
	group_stats gs_;
	
	// added code to reject flagged events  02/04/2008
	Int_t *run_nevents = new Int_t[nruns];
	Int_t *run_firstevno = new Int_t[nruns];
	if(reject_cuts>0) for(i=0;i<nruns;i++) run_nevents[i]=0;
	//
	

		// event tree
		TTree *evtree2 = new TTree("events","Events");	
		EVTREE2buffer;
		evtree2->Branch("Event", &ev_,buffer);
		//run flag tree
		TTree *rstree2 = new TTree("rstats","Run Flags and stats");
		RSTREE2buffer;
		rstree2->Branch("rs_branch",&rs_,buffer);		
		// group stats tree
		TTree *gstree2 = new TTree("gstats","Group Information");
		GSTREE2buffer;
		gstree2->Branch("gs_branch",&gs_,buffer);

	int progress=0;
	cout << "Saving events...\n";
	cout << "Progress:  0%\n";
	for (i = 0,evtree->GetEvent(i) ; i < nentries ; i++,evtree->GetEvent(i) ){
		if(float(i)/nentries*100.>float(progress+10)) {
			progress+=10;
			cout << "          " << progress << "%\n";
		}
		
		if(flag[i] && reject_cuts==-1){
			ev_ = ev;
			for(j=0;j<NDETS;j++) ev_.energy[j] = -3.;
			ev_.flag = 0;
			evtree2->Fill();
		}
		// added code to reject flagged events  02/04/2008
		else if(flag[i] && reject_cuts==1) {
			flag[i] = 0;
			continue;
		}
		else if(flag[i]==reject_cuts && reject_cuts>1) {
			flag[i] = 0;
			continue;
		}
		//
		else {
			ev_ = ev;
			ev_.flag = flag[i];
			evtree2->Fill();
			if(reject_cuts>0) run_nevents[ev.runno]++;
		}
		
	}
	
	
	// added code to reject flagged events  02/04/2008
	run_firstevno[0] = 0;
	if(reject_cuts>0) for(i=1;i<nruns;i++) run_firstevno[i]=run_firstevno[i-1]+run_nevents[i-1];
	//
	
	cout << "Saving runs...\n";
	for (i = 0 ; i < nruns ; i++ ){
		rs_ = rs[i];
		// added code to reject flagged events  02/04/2008
		if(reject_cuts>0) {
			rs_.firstevno = run_firstevno[i];
			rs_.nevents = run_nevents[i];
		}
		//
		rstree2->Fill();
	}
	cout << "Saving groups...\n";
	for (i = 0 ; i < ngroups ; i++ ){
		gs_ = gs[i]; 
		gstree2->Fill();
	}

	
	hfile->Write();
	hfile->Close();
	ssfile->Close();
	sprintf(buffer,"mv %s %s",temp_fname,filename);
	system(buffer);
	// added code to reject flagged events  02/04/2008
	delete [] run_nevents;
	delete [] run_firstevno;
	delete [] flag;

	cout << "Changes saved. Reloading...\n\n";
	SSload(filename);
	
}
/*
void SSsaverunflags(Char_t *flagname){ // does not work
	int i;
	char buffer[500];
	run_stats rs_;
	
	TFile *hfile = new TFile(buffer,"UPDATE",fname,9);
	TTree *rstree = (TTree*)hfile->Get("rstats");
	sprintf(buffer,"runno/I:groupno/I:flag[%d]/I:ntevents/I",NDETS);
	cout << "hi1\n";
	rstree->Branch(flagname,&rs_,buffer);
	cout << "hi2\n";
	for(i=0;i<nruns;i++){
		cout << i << endl;
		rs_ = rs[i];
		rstree->Fill();
	}
	hfile->Write();
	hfile->Close();	
}
*/

char *llfn_getname(FILE *fd){		// fd is pointer to open file
	int c = -1;
	int count;
	do{
		c++;
		count = fread (&runname[c], 1, 1, fd ); // put file name pointed to by fd into buffer runname one byte at a time
	} while (runname[c] != '\n' && count == 1); // before end of file name
	runname[c] = 0;   // !!!! interesting :-)  
	//cout << runname << endl;
	if(count == 0){
		return NULL;
	}
	return runname;
}


Int_t llfn_time2runno(Double_t time){
	int i;
	Double_t itime;
	evtree->GetEvent(0);
	itime=rs[ev.runno].starttime;
	for(i=0,evtree->GetEvent(i);(ev.time-itime)<time && i<nentries;i++,evtree->GetEvent(i));
	return rs[ev.runno].starttime;
}

Double_t llfn_time2unixtime(Double_t time){
	int i;
	Double_t itime;
	evtree->GetEvent(0);
	itime=rs[ev.runno].starttime;
	for(i=0,evtree->GetEvent(i);(ev.time-itime)<time && i<nentries;i++,evtree->GetEvent(i));
	return ev.time;
}


void llfn_runno2runname(Int_t runno, char *runlist){
// takes in the run number of a run and returns the daq file name from list_...
	Int_t pos=0, fstart=0;
	
	char buffer[100], *runname;
	
	if(!(runlistp = fopen(runlist,"rb"))){
		cout << "Error opening file " << buffer << endl;
		exit(1);
	}
	while((runname = llfn_getname(runlistp))){
		if(!(runnamep = fopen(runname, "rb"))){
		cout << "Error opening data file " << runname << endl;
		cout << "Suggestion: Check path is correct" << endl;
		return;
		}
		pos=0;
		do{
			fread (buffer, 1, 1, runnamep );
		} while (buffer[0] != '/');
		do{
			fread (&buffer[pos], 1, 1, runnamep );
		} while (buffer[pos++] != '\n');
		fstart = atoi(buffer);
		fclose(runnamep);
		if(fstart == runno){
			cout << runname << endl;
			return;
		}
	}
	cout << "no match found\n";
	return;
}

/* out of date
Int_t getentries(Int_t group, Int_t return_output){
//  returns the number of non flagged events in a group
	int i, cntentries=0;
	for (i=0,evtree->GetEvent(i);i<nentries;i++,evtree->GetEvent(i)) if(!flag[i] && (group == 0 || (group-1) == ev.groupno)) cntentries++;
	if(!return_output) cout << "Number of non-flagged events: ";
	return cntentries;
}
*/

Double_t SSgetuptime(Int_t detn){ // returns the total amount of time a detector is on in seconds
	int i;
	Double_t tottime = 0;
	for(i=0;i<nruns;i++) if(rs[i].det_on[detn-1] && gs[rs[i].groupno].flag!=1 && gs[rs[i].groupno].det_on[detn-1]) tottime++;
	return tottime*3600.;
}

Double_t SSgetflaggedtime(Int_t detn){ // returns the total number of runs*seconds that a detector has been flagged during periods it is on
	Double_t flaggedtime=0;
	int i;
	for(i=0;i<nruns;i++) if(rs[i].det_on[detn-1] && rs[i].flag[detn-1] && gs[rs[i].groupno].flag!=1 && gs[rs[i].groupno].det_on[detn-1]) flaggedtime++;
	return flaggedtime*3600.;
}

//int llfn_USE_DATA(char useflags[2]="__", char check_type = 'e', int seq_no = -1;){
	
	

Double_t SSgetusedtime(Int_t detn, Int_t incl_flaggedtime, Int_t incl_deadtime){
// this counts the total time that is usable. Flagged time and deadtime subtracted by default
// det is used for flaggedtime
	int i;
	Double_t tottime;
	
	tottime = SSgetuptime(detn);
// 	if(incl_flaggedtime) tottime -= (SSgetflaggedtime(detn) + SSgetcuttime(detn));
	/*else */for(i=0;i<ngroups;i++) if(gs[i].flag==0 && gs[i].det_on[detn-1]) tottime -= 3600.*gs[i].runs_toobig; // subtract time of files skipped because they where to large;
	if(incl_deadtime) tottime -= SSgetdeadtime(detn);
	
	if(incl_flaggedtime) cout << "Warning, flaggedtime feature has been removed till it is updated\n";
	return tottime;
}

Double_t SSgetdeadtime(Int_t detn){
// get dead time of unflagged data by multiplying the numbe of events by the deadtime per event.

	int i;
	Double_t deadtime=0;
	for(i=0;i<nruns;i++) if(rs[i].det_on[detn-1] && rs[i].flag[detn-1]!=1 && gs[rs[i].groupno].flag!=1 && gs[rs[i].groupno].det_on[detn-1]) deadtime+=rs[i].neventsall*gbl_deadtimeperevent; // number of events in all runs used.
	return deadtime;
}

void SSgetdeadtimeR(char help='H') {
	help=help; // avoid compile warning
	cout << "This function calculates the total deadtime of a detector within the specified energies Emin and Emax (keV).\n\n" 
		 << "The default deadtime is 1.2E-4 seconds per event\n\n"
		 << "Energy range:\n"
		 << "     Emin  - Includes data with run lower thresholds below this energy (keV).\n"
		 << "             Setting Emin = 0 ignores this option.\n"
		 << "     Emax  - Include data with group upper thresholds aboe this energy (keV).\n"
		 << "             Setting Emax = 0 ignores this otption.\n"
		 << "Options (eg opt=\"vp\"): \n"
		 << "     a     - Include all events printed to the daq file in the estimate\n"
		 << "     v     - Include the veto events (eg HV eveto counts) in the estimate\n"
		 << "     p     - Plot the results to h_handle\n\n"
		 << "NB: Data included only if events, runs or groups are not flagged, there detector flags are on and\n"
		 << "    runs have thresholds within the required energy range.\n\n";
}
Double_t SSgetdeadtimeR(Int_t detn, Double_t Emin, Double_t Emax, Double_t deadtimeperevent, Char_t opt[5]="av") {
	int i;
	int incl_allevents=0, incl_vetoevents=0, plot=0;
	for(i=0;opt[i]!='\0';++i) {
		if(opt[i]=='a') incl_allevents=1;
		if(opt[i]=='v') incl_vetoevents=1;
		if(opt[i]=='p') plot=1;
	}
	
	Double_t dead_time=0;
	if(plot) h_handle = new TH1D("h2_deadtime","",nruns,0,nruns);
	
	for(i=0;i<nruns;++i) {
		if(	!gs[rs[i].groupno].flag &&											// Skip if group is flagged
			gs[rs[i].groupno].det_on[detn-1] &&									// Skip if detector is not on for this group
			rs[i].det_on[detn-1] &&												// Skip if detector is not on for this run
			!rs[i].flag[detn-1] &&												// Skop if detector is flagged during this period
			(Emax==0 || (Emax < gs[rs[i].groupno].hi_thresh[detn-1])) &&		// Skip if outside uppper range, if Emax=0 ignore this range
			(Emin==0 || (Emin > rs[i].lo_thresh[detn-1]))						// Skip if outide lower range
			) {
			dead_time += (rs[i].neventsall*incl_allevents +rs[i].nvetoevents*incl_vetoevents)*deadtimeperevent;
			if(plot) h_handle->Fill(i,(rs[i].neventsall*incl_allevents +rs[i].nvetoevents*incl_vetoevents)*deadtimeperevent/36.);
		}
		
	}
	if(plot) {
		h_handle->GetXaxis()->SetTitle("run");
		h_handle->GetYaxis()->SetTitle("Percentage dead time per hour");
		h_handle->Draw();
	}
	return dead_time;
}	

Double_t SSgetusedtimeR(Int_t detn, Double_t Emin, Double_t Emax, Double_t deadtimeperevent) { // returns the running time within the specified range
	int i; 
	
// 	h_handle2 = new TH1D("spec_time","spec time",specUpperRange,0.,specUpperRange);
	Double_t used_time=0;
	
	for(i=0;i<nruns;++i) {
		if(	!gs[rs[i].groupno].flag &&											// Skip if group is flagged
			gs[rs[i].groupno].det_on[detn-1] &&									// Skip if detector is not on for this group
			rs[i].det_on[detn-1] &&												// Skip if detector is not on for this run
			!rs[i].flag[detn-1] &&												// Skop if detector is flagged during this period
			(Emax==0 || (Emax < gs[rs[i].groupno].hi_thresh[detn-1])) &&		// Skip if outside uppper range, if Emax=0 ignore this range
			(Emin==0 || (Emin > rs[i].lo_thresh[detn-1]))						// Skip if outide lower range
			) {
			used_time += 3600. - (rs[i].neventsall+rs[i].nvetoevents)*deadtimeperevent;
		}
	}
	return used_time;
}
			

Int_t SSgetrunstoobig(Int_t detn){
	int i;
	Int_t runstoobig = 0;
	for(i=0;i<ngroups;i++) if(gs[i].flag!=1 && gs[i].det_on[detn-1]) runstoobig += gs[i].runs_toobig;
	return runstoobig;
}

Int_t SSgetusedgroupcnt(Int_t detn){
// returns the number of groups that a detector has been used
	int i;
	Int_t cnt=0;
	for(i=0;i<ngroups;i++) if(gs[i].flag!=1 && gs[i].det_on[detn-1]) cnt++;
	return cnt;
}

void SSdumprunningtimes(){
	int i;
	cout << "Available data before and after cuts in hours and (percent). \n";
	cout << "Det  nocuts    sizecuts     runcuts  deadtime  usedgroups\n";
	for(i=1;i<=NDETS;i++)
		printf("%3d  %6.0lf  %5.0lf(%3.0lf)  %5.0lf(%3.0lf)  %8.3lf  %10d \n",
			i,
			SSgetuptime(i)/3600.,
			SSgetusedtime(i,0,0)/3600.,SSgetusedtime(i,0,0)/SSgetuptime(i)*100.,
			SSgetusedtime(i,1,0)/3600.,SSgetusedtime(i,1,0)/SSgetuptime(i)*100.,
			SSgetdeadtime(i)/3600.,
			SSgetusedgroupcnt(i));
}

Int_t SSgetdatetimedif(int i, int f){
// this function can be used to report the total time between two callender dates
	Int_t start, finnish;
	evtree->GetEvent(i); start=rs[ev.runno].starttime;  
	evtree->GetEvent(f); finnish=rs[ev.runno].starttime;
	return finnish - start;
}
	
void llfn_getgroupname(int groupno){
	int i;
	char buffer[50];
	for(i=0;;i++){
		buffer[i]=Char_t(gs[groupno].group_name[i]);
		if(gs[groupno].group_name[i]=='\0') break;
	}
	printf("%-20s", buffer);
}

void SSgetResEqn(Int_t group, Int_t det){
	double m, merr, c, cerr;
	m = gs[group-1].res_eqn[0][det-1]*2.35;
	c = gs[group-1].res_eqn[2][det-1]*2.35;
	printf("Delta_E = %1.3lf*E + %1.3lf keV \n",m,c);
	merr=0;
	cerr=0;
}

void SSdumpdetstats(){ // print out useful information about the data
	int i, j;
	for(i=0;i<ngroups;i++){
		if(!gs[i].flag){
			cout << "Group " << i+1 << endl;
			cout << "Det active lo_thresh hi_thresh   calibration_eqn(m,me,c,ce)\n";
			for(j=0;j<NDETS;j++) 
				printf("%3d %6d %9.0lf %9.0lf   %4.3lf %4.3lf %6.1lf %4.1lf\n",
					j+1,gs[i].det_on[j],
					gs[i].lo_thresh[j],gs[i].hi_thresh[j],
					gs[i].cal_eqn[0][j],gs[i].cal_eqn[1][j],gs[i].cal_eqn[2][j],TMath::	Abs(gs[i].cal_eqn[3][j]));
		}
	}
}

void llfn_makewhite(TCanvas *c1){
	// set plot formatting - ie make really white
 	c1->SetFrameFillColor(0);
 	c1->SetFillColor(0);
 	c1->SetHighLightColor(0);
 	c1->SetFrameLineWidth(0);
}

void SSaddleggend(char title1[],char title2[],TH1D *h_handle1, TH1D *h_handle2){
	TLegend *leg = new TLegend(0.724138,0.769068,0.889368,0.891949,NULL,"brNDC");
	TLegendEntry *entry=leg->AddEntry(h_handle1,title1,"l");
	entry=leg->AddEntry(h_handle2,title2,"l");
	leg->SetBorderSize(0);
	leg->Draw();
}

void SSprint(Int_t opt, char name[64], TCanvas *c1){
	char buffer[64];
 	if(opt==0){
		sprintf(buffer,"%s.eps",name);
		c1->Print(buffer);
		sprintf(buffer,"%s.C",name);
		c1->Print(buffer);
	}
 	else if(opt==1){
		c1->Print(name);
		
		//TROOT simple("simple","nTuple");
		//sprintf(fname,"%s.root", argv[1]);
		//TFile *hfile = new TFile("c1.root","RECREATE","c1",9);
		//hfile->Write();
		//hfile->Close();
 	}
}

void llfn_setstyles(char OptStat[16]){
//     k :  kurtosis printed
//     K :  kurtosis and kurtosis error printed
//     s :  skewness printed
//     S :  skewness and skewness error printed
//     i :  integral of bins printed
//     o :  number of overflows printed
//     u :  number of underflows printed
//     r :  rms printed
//     R :  rms and rms error printed
//     m :  mean value printed
//     M :  mean value mean error values printed
//     e :  number of entries printed
//     n :  name of histogram is printed
	
	
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(OptStat);//gStyle->SetOptStat(10);//gStyle->SetOptStat(1000010); 	// set stat box options
	gStyle->SetOptFit(1111);
	
}

void SShelp(){

	cout << "list of functions:\n\n";
	cout << "load:........................load the file to be analysed\n\n";
	cout << "plotting functions.\n"
		 << "plot_spec:...................plot spectrum\n"
		 << "plot_waterfall:..............plot a event energies vs time\n"
		 << "plot_correlations:...........plot the energy of one detector against another\n"
		 << "plot_eph:....................plot the number of events per hour in histogram or with time\n"
		 << "plot_timediffs:..............plot the time between events of any detector or a specified detector\n"
		 << "plot_ttlcounts:..............plot the fraction of ttl flags for HF, FF, vibrations, mains-spikes\n";
		 
		 

}

void llfn_decimal2binary16(Int_t dec_number,  Int_t *bit16_number){
// adapted from M.J. Leslie.'s binary_op()
	Int_t count=16;                         /* Number of bits in a byte.    */
	Int_t MASK = 1<<(count-1);
	
	while(count--){
		*bit16_number = ( dec_number & MASK )?1:0;
		dec_number <<= 1;
		bit16_number++;
	}
}

Int_t llfn_Decimal2Bit(Int_t decimal, Int_t bitN, Int_t Nbits){
// geven a decimal number this function returns a 1 or 0 for the reqested bit.
// Left most bit is BitN=1 and the next to the right is 2 etc. 
	bitN = Nbits - bitN +1;
	Int_t MASK = 1<<(Nbits-1);
	while(Nbits-- && Nbits>=bitN) decimal <<= 1;
	return ( decimal & MASK )?1:0;
}

void SSplot_ttlcounts(){
	int i, j;
	Int_t ttlflag[16];
	h_handle = new TH1D("h1","h1",20,0.,5.);
	for(i=0,evtree->GetEvent(0);i<nentries;i++,evtree->GetEvent(i)){
		llfn_decimal2binary16(ev.veto,ttlflag);
		for(j=0;j<4;j++) if(ttlflag[j]==(j>1)) h_handle->Fill(j+1);
	}
   	Int_t ci;   // for color index setting
   	ci = TColor::GetColor("#b6b6f9");
   	h_handle->SetFillColor(ci);
	h_handle->Scale(1./nentries);
	h_handle->Draw();
}

void SSdoruntest(){
	// this run test was writen from the instrunctions at http://home.ubalt.edu/ntsbarsh/Business-stat/opre504.htm#rrunstest
	unsigned int i, runs=0, n1=0, n2=0;
	int onoff=0;
	Double_t a, s, z, expected_runs; // a = expected mean
	vector <Int_t> sequence;
	Double_t mean = TMath::Mean(epp.epp.size(),&epp.epp[0]);
	for(i=0;i<epp.epp.size();i++) {
		if(epp.epp[i] > mean) sequence.push_back(1);
		else if(epp.epp[i] < mean) sequence.push_back(0);
	}
	// get number of runs (consequative number of aboves of bellow the means...
	if(sequence.size()) onoff = !sequence[0];
	for(i=0;i<sequence.size();i++) {
		if(onoff!=sequence[i]) {
			runs++;
			onoff = sequence[i];
		}
		n1 += sequence[i];
	}
	n2 = sequence.size() - n1;
	// expected mean:
	a = 1. + 2.*n1*n2/(n1+n2);
	s = sqrt(2.*n1*n2*(2.*n1*n2-n1-n2)/((n1+n2)*(n1+n2)*(n1+n2-1.)));
	expected_runs = 1. + 2.0*n1*n2/(n1+n2);
	z = (runs - expected_runs)/s;
	
	cout << "a = " << a << endl;
	cout << "s = " << s << endl;
	cout << "runs = " << runs << endl;
	cout << "expected_runs = " << expected_runs << endl;
	cout << "z = " << z << endl;
	
	TF1 f1("f1","TMath::Gaus(x)/sqrt(2*TMath::Pi())",-10,10);
	
	printf("\nThere is a %1.2lf%% chance that the sequence is random.\n",(1.-f1.Integral(-1.*TMath::Abs(z),TMath::Abs(z)))*100.);
	
	cout << "\n\nNB. if z < 0 then too few runs. z > 0 too many runs.\n";
	
// 	cout << "\n\n\nIf the normalized test statistic is greater in absolute magnitude than 2.580, \n"
// 		 << "then this means, loosely speaking, that there is less than a one percent \n"
// 		 << "chance that the pattern under analysis was generated by a random process. \n"
// 		 << "The value 2.580 comes from statistical tables. If the test statistic is negative,\n"
// 		 << "the actual number of runs is smaller than the expected number of runs, and vice versa.\n";
// // 	return z;
}



void SScut_ttlflags(char selection[4]="HELP", Double_t cutsecsbefore=0., Double_t cutsecsafter=0.){
	// this function slectively flags data that has been tagged by TTL flags
	// selection = "A" for all or "0110111111111111" for any chosen 

	int i, j, k/*, deti*/;
	Int_t ttlflag[16];
	Double_t cur_time=0;
	Double_t time_end=0;
	
	
	if(selection[0]=='H'){
		cout 	<< "selection = A cuts all\n" 
				<< "          = 1111 cuts on hf, ff, vibration, mainsspike\n\n";
		return;
	}
	
// 	for(deti=0;deti<NDETS;deti++) newcuts[deti].clear();
		
	for(i=0,evtree->GetEvent(0),time_end=ev.time;i<nentries;i++,evtree->GetEvent(i)){
		llfn_decimal2binary16(ev.veto,ttlflag);
		for(j=0;j<4;j++){
			//(k>1) takes care of that the first 2 bits are 1 when off and the next 2 are 0 when off.
			if((selection[0]=='A' || selection[j]=='1') && (ttlflag[j]==(j>1))){ // NOT function applies to first to numbers only
				cur_time = ev.time;
				// go back "cutsecsbefore" seconds
				if((cur_time-cutsecsbefore)<ev.time) {
					for(;(cur_time-cutsecsbefore)<ev.time && ev.time>=time_end && i>=0;i--,evtree->GetEvent(i)); 
					i++;
					evtree->GetEvent(i);
				}
				// flag events till intitial position
				for(;ev.time<=cur_time && i<nentries;i++,evtree->GetEvent(i)) flag[i] = 1;
				i--; evtree->GetEvent(i);
				// flag events till "cutsecsafter" after last  ttlflag
				if((cur_time+cutsecsafter)>=ev.time){
					for(;(cur_time+cutsecsafter)>=ev.time && i<nentries;i++,evtree->GetEvent(i)){ 
						flag[i] = 1;
						llfn_decimal2binary16(ev.veto,ttlflag);
						// check that events in the period are not vetoed and if so continue cutting for another 'cutsecsafter'
						for(k=0;k<4;k++) if((selection[0]=='A' || selection[k]=='1') && (ttlflag[k]==(k>1))){
							cur_time = ev.time; // reset time of last ttle flag
							break;
						}
					}
					i--;
				}
				time_end = cur_time+cutsecsafter; // save time that the cut ended
				break;
			}
		}
	}
	
	

	
}



	

void SScut_ondummydet(Int_t dummyNo, Double_t cutsecsbefore=0., Double_t cutsecsafter=0.){

	int i, j;
	Double_t cur_time=0;
	
		
	for(i=0,evtree->GetEvent(0);i<nentries;i++,evtree->GetEvent(i)){
		for(j=0;j<4;j++){
			if(ev.energy[dummyNo]==-3.){ 
// 				cout << "hi\n";
				cur_time = ev.time;
				for(;(cur_time-cutsecsbefore)<ev.time && i>=0;i--,evtree->GetEvent(i)); // go back "cutsecsbefore" seconds
				for(;(cur_time+cutsecsafter)>=ev.time && i<nentries;i++,evtree->GetEvent(i)) flag[i]=1;
			}
		}	
	}
}


/// #################################################### BETA FUNCTIONS ##################################################



void SSbeta_getTriggerThresh(int det, char veto[3]="000"){
	// define spectral shape
	TF1 *e2 = new TF1("e2","TMath::Exp([1]*x+[0])",60,100);e2->Update();
	TF1 *g1 = new TF1("g1","gaus",100,325);g1->Update();
	TF1 *e1 = new TF1("e1","TMath::Exp([1]*x+[0])",400,1000);e1->Update();
	TF1 *fspectrum = new TF1("fspectrum","e2+g1+e1",0,2000);
	fspectrum->SetParameter(0,6.88884);
	fspectrum->SetParameter(1,-0.0140377);
	fspectrum->SetParameter(2,262.718);
	fspectrum->SetParameter(3,186.524);
	fspectrum->SetParameter(4,60.4514);
	fspectrum->SetParameter(5,3.28321);
	fspectrum->SetParameter(6,-0.00199871);
	fspectrum->SetLineColor(3);


	Int_t dets[24]={2,0};
	dets[0]=det;
	int i, j;
	int old_thr;
	Int_t rebin = 10;
	gbl_verbose=0;
// 	gbl_draw=1;
	
// 	SSlocate_Thresholds();
// 	SScut_plot_coincedences(0,1);
// 	

	double *old_run_thr = new double[nruns];
	for(i=0;i<nruns;++i) {
		old_run_thr[i] = rs[i].lo_thresh[det-1];
		rs[i].lo_thresh[det-1] = 0;
	}
	
	for(i=0;dets[i];++i) {
		old_thr = gs[0].lo_thresh[dets[i]-1];
		gs[0].lo_thresh[dets[i]-1] = 0;
		SSplot_spec2(dets[i],rebin,veto,1);
// 		h_handle->Scale(SSgetusedtime(dets[i],1,0)/(SSgetusedtime(dets[i],0,0)-SSgetusedtime(dets[i],1,0)));
		gs[0].lo_thresh[dets[i]-1] = old_thr;
		// (assuming more than 1 detector is running)
		// loop over histogram bins
		for(j=h_handle->GetMaximumBin();j<1000/rebin;++j) {
			// check is the threhold is too low and drop if needed
			if((h_handle->GetXaxis()->GetBinLowEdge(j) < gs[0].lo_thresh[dets[i]])) {
				if(h_handle->GetBinContent(j) < fspectrum->Eval(h_handle->GetXaxis()->GetBinLowEdge(j))) {
					cout << "Set threshold: " << h_handle->GetXaxis()->GetBinLowEdge(j) << "keV or setting: " 
						<< SSch2th.Eval(SSe2ch(h_handle->GetXaxis()->GetBinLowEdge(j),dets[i])) 
						<< ".  actual threshold: " // actual threshold is the threshold after ceiling the thr no.
						<< SSch2e(SSth2ch.Eval(TMath::Ceil(SSch2th.Eval(SSe2ch(h_handle->GetXaxis()->GetBinLowEdge(j),dets[i])))),dets[i]) 
						<< endl;
					cout << "old Threshold: " << rs[nruns-1].thr[dets[i]-1] << endl;
					break;
				}
			}
			// otherwise increase if needed
			else {
				if(h_handle->GetBinContent(j) < (fspectrum->Eval(h_handle->GetXaxis()->GetBinLowEdge(j)) + 100)) {
					cout << "Set threshold: " << h_handle->GetXaxis()->GetBinLowEdge(j) << "keV or setting: " 
						<< SSch2th.Eval(SSe2ch(h_handle->GetXaxis()->GetBinLowEdge(j),dets[i])) 
						<< ".  actual threshold: " 
						<< SSch2e(SSth2ch.Eval(TMath::Ceil(SSch2th.Eval(SSe2ch(h_handle->GetXaxis()->GetBinLowEdge(j),dets[i])))),dets[i]) 
						<< endl;
						cout << "old Threshold: " << rs[nruns-1].thr[dets[i]-1] << endl;
					break;
				}
			}
		}
	}
	if(gbl_draw) fspectrum->Draw("same");
	
	for(i=0;i<nruns;++i) rs[i].lo_thresh[det-1] = old_run_thr[i];
	delete[] old_run_thr;
}




/// ######################################## DEBUG ###########################################################

void SScheckorder(){
	int i=0;
	evtree->GetEvent(i);
	Double_t last_time=ev.time;
	for(i=1,evtree->GetEvent(i);i<nentries;i++) if(ev.time < last_time) cout << "ouch!\n";
}


