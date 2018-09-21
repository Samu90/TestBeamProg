# TestBeamProg

-Ntuples_wiki: basic Wiki for ntuples from T10 May 2018 Test Beam on crystal bars-see below this list.

°°°°°°°°°° The argument Conf should be substitued with the path of the .root file that one wants to analyze°°°°°°°°°°°

-plotWF.C(Conf,nevento) plots the 6 wavwforms acquired in one event if given (string runname,int eventnumber), plots iteratively for all the events if given (string runname, int 0); expects a primitive to pass from one event to the next one

-ploWF_graph.C(Conf) plots max amplitude distribuitions  and creates a file histoConf.root contaiing such distributions

-plotHisto.C(histoConf1,histoconf2) plots max amplitude distributions overlapped. histoConf1 and histoConf2 are .root files from plotWF_graph.

-plotWF_fit.C(Conf) plots max amplitude distributions and fits the distribuitions with a landau function

-plotWF_cut.C(Conf) same as plotWF_fit.C + cuts on data inbetween (0.8*MIP,3*MIP)

-plotWF_tamp.C(Conf) plots 2D histograms  t_r-t_MCP t_l-t_MCP t_ave-t_MCP vs max amplitude , fits on time averages of 1D histograms obtaining by summing 1D histograms at fixed max amplitude for bins of max amplitude. It uses time recorded with LED300 by default: AMPWALK correction.

-plot_lsign.C(Conf) plots time distributions for uncorrected (no ampwalk, no tdiff) timestamps

-plotWF_su.C(Conf) plots max amplitude vs t_left-t_right

-plotWF_corr.C(Conf) plots same graphs as _tamp.C, redefines timestamps using fits done in _tamp.C and plots time distributions for the corrected timestamps(AMPWALK corr)

-AmpCorrUniv.C(Conf) same as _corr.C with aoutomatic histograms range settings for all the configurations  

-plotWF_tdiff.C(Conf) plots same graphs as _tamp.C, then implements TDIFF correction by fitting the averages of 1D histograms in t_ave-t_MCP Vs t_left-t_right at fixed t_left-t_right with a line. Redefine timestamps in order to obtain t_ave-t_MCP independent from t_left-t_right.
Plots the distribution of t_ave-t_MCP for the new (AMPWALK+TDIFF) corrected timestamps compared with the only AMPWALK corrected one.
If ControlTdiff is TRUE prints comparison among uncorrected, ampwalk corrected, and ampwalk+tdiff corrected 2D histograms of t_ave-t_MCP vs t_left-t_right.
Plots the differential behavior of time resolution coming from these corrections in bins of t_left-t_right (15 ps mcp res subtracted).
Automatic histograms ranges for all the config

-plotWF_time.C(Conf) plots the uncorrected and fitted 2D histograms as in _tamp.C. Plots the diifferential behavior of 1) uncorrected left, right and ave timestamps in bins of t_left-t_right;
	       	     	 	     	 	   	      	    	     	       		     	      	 2) amp walk corrected left, right and ave and ampwalk+tdiff corrected ave(green) timestamps															 in bins of t_left-t_right;
Automatic histograms ranges for all the config;

-GraficiTdiff.C(Conf) same as _tdiff.C. Does not plot the differential behavior but plots ave timestamps distributions for different cuts in tdiff around the bin with maximum number of hits
Automatic histograms ranges for all the config;

-ResolutionUncorr.C(Conf) plots total uncorrected timestamps distributions

-TResAmp.C(Conf) plots the uncorrected and fitted 2D histograms as in _tamp.C. Plots the differential behavior of the time resolution as obtained in plotWF_tdiff in bins of amp_max/MIP peak and interpolates to exrapolate light needed to reach 30 ps;
Automatic histograms ranges for all the config

-RisVsMip.C(Conf) plots the 2D histograms for left and right timestamps vs amp_max performing different cuts on data (instead of cutting data around the landau fit MIP peak , scans on MIP peak value inbetwween 0.11 and 0.2 and cuts around those values). Plots the differential behavior of time resolution in bins of MIP peak values.
Automatic histograms ranges for all the config

-RisRelVsMip.C(Conf) same as RisVsMip but shows the differential behavior of timeres/gauss peak to check on the the relative values;
Automatic histograms ranges for all the config

.RisRelVsTdiff.C(Conf) plots the 6 2D histograms as in _tamp.C. Plots the time resolutions of left, right and average amp walk corrected timestamps. Plots the differential behavior of amp walk corrected timeres/gauss Peak in bins of t_left-t_right.
Automatic histograms ranges for all the config

-DarkBkg() creates and saves in controlplots the amp_max distributions of DCR runs for left and right SIPM compared with the noDCR run for 1.2 conf. Builts a graph comparing the MIP Peak at different dark current values.

-RisVsDcr(string "version") implements different type of corrections and prints comparison plots of time resolution vs dark current values for: uncorrected, tdiff corrected, old ampwalk+tdiff corr, new ampwalk corr, with new ampwalk based on the gaussian peak rather than means of the 1D histograms for correction fits (see plotWF_tamp.C). Saves all control plots and timedstamps resolutions in HDCRPlots. Prints also comparison among amp walk corrected sigma and RMS of disribution. Version should be "old" for noDCR optimized H4Analysis recos, "new" for the optimized one: the two versions take input files from different paths, please check and/or change them in the .C file

-RisVsDcr1.C(string "threshconf") same as RisVsDcr + compares old and new ampwalk fits function for left and right SiPM. Available for two "threshconf"- they should be the names of the directories containing, f.i. the configurations at a certain nino threshold value. Check the .C file to see how the paths are set.

-AmpSpatUnif.C(Conf) plots the 2D histgrams as in plotWF_tamp and prints the differential behavior of amp_max for both left and right SiPMs in bins of t_left-t_right.

-AmplitudeStudy.C(Conf) to be used with DCR runs. Plots the usual 2D histograms and studies DCR noise on amplitude distributions: plots the distributions around the MIP peak (first row) and far from it (noise, second row) and subtract one to the other to see if the noise would not spoil the MIP peak
Automatic histograms ranges for all the config

-AmpCorr.C(Conf) = AmpCorrUniv.C(Conf)?, no resolution plots.

-Copy script to move runs in lxplus from location in eos to ~/Files/.

-RecoDiff.C() compares time resolutions coming from different(DCR optimized and regular) recos and same post-reco analysis. Save all controlplots and plots a summary plot with comparisons. Paths to the runs analyzed with different reco can be modified in the .C file.








======T10 ntuples elementary wiki =======

This file contains information on the structure and main variables of the ntuples of the T10 May 2018 TesBeam (LYSO Crystal Bar time resolution).

The good ntuples(both for data and pedestals) are available on lxplus on /eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_T10_May2018/ntuples_v1/
The TestBeam logbook (info on configurations' different settings) is available at https://docs.google.com/spreadsheets/d/1ArGOxF1clg_I_9lgCssy9A57RwJ98hGXUgBIKj3aCq4/edit#gid=1831808220

Each ntuple contains 4 trees:
     --info
     --wf
     --h4
     --digi

==========TREES DESCRIPTION=============

The INFO tree contains:
_____index of reference for both run and event (Format run00000000event);
_____tableX, tableY: information on the position of the table;
_____config: code of configuration;
_____Vbias_bar: bias potential on the bar left and right SIPMs;
_____SIPM_current bar: current in SIPMs (in ADC counts);
_____NINOthr_bar: threshold (in ADC counts) on the raw SIPM signal to trigger NINO acquisition;

The tree has also a section dedicated to voltage and current on crystalbar matrices, not used in these ntuples.




The WF tree  contains the information on the raw waveforms collected from each detector in the system:
____index of reference for both run and event (format run0000000event0)
____WF_samples: total number of samples collected for a single event;
____WF_ch(vector): each sample taken is associated to the channel it belongs by the value of this vector at given sample (e.g. wf_ch[45]=2 means that the 45th sample belongs to the second channel
waveform). There are 6 or 8 channels in the configurations studied since now:
	   	       	    	     	         --0 - MCP pulse signal
						 --1/2- NINO amplified signal x 2																					    --3/4- Amplified signal x 2
						 -- Unkown channel (almost always noise)*;																				    --6/7- trigger delayed signal x 2**



To be clarified:
* only in some configurations;
**in some configurations there is only one trigger signal


The h4 tree contains parameter and settings necessary for the first reco step. This reco is operated with the H4Analysis software on wf data, and the values obtained are stored in digi tree:

It contains all integer numbers:
___index: same as in info;
___start_time;
___time_stamp;
___run: run number;
___spill: number of spill;
___event: is the event number also present in the index, +1;


The DIGI tree contains partially reconstructed and organized data. The variables used to perform the time resolution analysis are amp_max and time.
-------amp_max is a vector containing the  waveform  maximum recorded on different detectors;
-------time is the timestamps of the rising signals in different conditions;
Each different detector or detector+settings system is identified with an integer number: the system is flexible in terms of changes and upgrades of detector and/or electronics.


There are also: --- index: to retrieve run and event number (format as in info);
      	     	---n_channel: the number of channel recorded (6/8);
			 ---n_timetypes: the number of different time acquisition recorded and encoded;

amp_max and time are both n_timetypes long for each recorded event.

The first encoded values (MCP, NINOBAR1,NINOBAR2,AMPBAR1,AMPBAR2) identify the detectors: their index is to be used to retrieve max amplitudes: in fact the amp_max vectors are only full on their
first n_channel entries.

The time information is instead recorded in many different ways:
----------------timestamps on the MCP signal are always taken in CFD (Constant Fraction Discrimination, fixed at 50%) to take care of the detector noisy behavior;
----------------timestamps on SIPMs are collected both on the amplified waveforms with CFD and on NINO with different LED (Leading Edge Discriminations): in this way lots of different timestamps
			      can be taken into account to identify the right combination of detector and detector settings to obtain the better resolution. The NINOi+LED300 configuration is the one used in the macros in this analysis.
			      
