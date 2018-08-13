# TestBeamProg

-Wiki base delle ntuple utilizzate per l'analisi dei dati del Test Beam di MAGGIO 2018 A T10

      
													USEFUL PROGRAMS
													
-plotWF.C(Conf,nevent ) the program plots the 6/8 waveforms of a sigle events (nevent) if nevent==0 it plots the waveforms of the events in sequence (push enter to go on)  

-ploWF_graph.C(Conf) the program plot the max amplitute distribution and create a file (histoConf.root) with the histogram 

-plotWF_fit.C(Conf) the program plots the distribution of max amplitude and fits the distribution with a Landau function

-plotWF_cut.C(Conf) the program plots the distribution of max amplitude and fits the distribution with a Landau function and makes a cut on the events in (0.8*MIP,3*MIP)

-plotWF_tamp.C(Conf) the program plots the 2D histogram of t_r-t_MCP t_l-t_MCP vs max amplitude and t_ave-t_MCP vs t_l-t_MCP-t_r-t_MCP, it also makes a logaritmic fit on the mean values of time variable for each bin of max amplitude.

-ResolutionUncorr.C(Conf) the program calculates the time resolution of the cristal bar without any correction on amplitude or time.

-AmpCorrUniv.C(Conf) the program plots the same graph of _tamp.C whith the graph corrected with the amplitude walks, the time stamps of NINO with LED300 are used. In addiction the program plot the gaussian of t_ave-t_MCP distribution

-AmpSpatUnif.C(Conf) the program plots the values of max amplitude vs t_left-t_right corrected with amplitude walk.

-GraficiTdiff.C(Conf) the program plot the gaussian of t_ave-t_MCP distribution with correction on amp max and, with and withouth t_diff correction, calculating the different time resolution, it also calculate the time resolution for different value of t_left-t_right: |dt|<400ps, |dt|<200ps, |dt|<100ps, |dt|<50ps, |dt|<20ps, |dt|<10ps

-ResNinoT.C the program plots the time resolution values vs the NINO treshold, it plots also the gaussian distribution of t_ave-t_MCP for each treshold with related fit. 

-RisRelVsMip.C the program plots the relative time resolution (sigma/mean) vs the MIP peak value with only amplitude walk correction.

-RisRelVsTdiff.C(Conf) the program plots the relative time resolution vs t_left-t_right in bins

-RisVsMip.C the program plots the absolute time resolution vs the MIP peak value with only amplitude walk correction.

-TResAmp.C(Conf) the program plots the time resolution vs bin of (max amp)/(MIP Peak)  

-PlotWF_tdiff.C(Conf) the program plots data with tdiff corrections, it plots a comparison between raw data, data corrected with amp-walk correction and data corrected with amp walk+tdiff correction. The gaussians with and without tdiff correction are plotted and the resolutions are calculated. In addition it plots the values of time resolution for binning of tdiff and for each bin on tdiff the single histogram is plotted






======T10 ntuples elementary wiki =======

This file contains information on the structure and main variable of the ntuples of the T10 May 2018 TesBeam (LYSO Crystal Bar time resolution).

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
_____tableX, tableY: information on the position of the bar;
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


The DIGI tree contains partially reconstructed and organized data.The variables used to perform the time resolution analysis are amp_max and time.
-------amp_max is a vector containing the maximum recorded on different detectors;
-------time is the timestamps of the rising signals in different conditions;
Each different detector or detector+settings system is identified with an integer number: the system is flexible in terms of changes and upgrades of detector and/or electronics.

A stamp of the encoding is below:

There are also: --- index: to retrieve run and event number (format as in info);
      	     	---n_channel: the number of channel recorded (6/8);
			 ---n_timetypes: the number of different time acquisition recorded and encoded;

amp_max and time, are both n_timetypes long for each recorded event.

The first encoded values (MCP, NINOBAR1,NINOBAR2,AMPBAR1,AMPBAR2) identify the detectors: their index is to be used to retrieve max amplitudes: in fact the amp_max vectors are only full on their
first n_channel entries.

The time information is instead recorded in many different ways:
----------------timestamps on the MCP signal are always taken in CFD (Constant Fraction Discrimination, fixed at 50%) to take care of the detector noisy behaviour;
----------------timestamps on SIPMs are collected both on the amplified waveforms with CFD and on NINO with different LED (Leading Edge Discriminations): in this way lots of different timestamps
			      can be taken into account to identify the right combination of detector and detector settings to obtain the better resolution.
			      
