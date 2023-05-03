# GRAB-data-analysis
Two Photon Imaging Analysis Suite
Last edited: 12/17/2019
Created and Modified by Simon Trinh and Lily Zou
Designed in Matlab 2018b

Includes:
Session Analysis of single (green) channel with visualization (Simon_2p_analysis_Combo)
Session Analysis of single channel with no visualization (Simon_2p_analysis_Combo_no_visual)
Session Analysis of two channels (red and green) with visualization (Simon_2p_analysis_Red_Green)
Session Analysis of two channels with no visualization (Simon_2p_analysis_Red_Green_no_visual)
Session Struct (single channel) (Simon_Session_Struct)
Session Struct (two channel) (Simon_Session_Struct_Red_Green)
Multi-session analysis by color plots (Simon_Color_Plots)
Multi-session analysis (two channel) by color plots (Simon_Color_Plots_Red_Green)

(Note: folder tree still uses 200ms, script currently uses 500ms before stimulus. To change, go to timeFrame variable and change the multiplier)

Instructions:
1. Determine if you are reading single channel or two channel data
1a. If single channel data, use (Simon_2p_analysis_Combo) or corresponding no visual
1b. If double channel data, use (Simon_2p_analysis_Red_Green) or corresponding no visual

2. If computer graphics cannot handle visualization of several hundred figures, use no visual to generate load files 
   that can be transferred to a computer with higher processing power, then visualization code can be run at greater speeds.
2a. If load files already exist from a previous version (Prior to November 2019 or if you have made any changes to the load section)
    please delete existing load files.
2b. If you are unsure of processing power, be advised that a 16 GB RAM Intel(R) Xeon(R) CPU E3-1225 v5 @ 3.30 GHz will take appx.
    30-90 minutes for a single session of ~50K frames with no visualization and Intel(R) HD Graphics P530 cannot handle visualization.

3. To begin, please ensure each session contains the following:
	a. Mat file in the form (####_####_###.mat) (e.g. 1024_0726_000.mat)	
	b. Trial file in the form (####_####_###.trials) (e.g. 1024_0726_000.trials)
	c. Align file in the form (####_####_###.align) (e.g. 1024_0726_000.align)
	d. Behavior file in the form (behavior_####.mat) (e.g. behavior_0726.mat)
	e. Associated scanbox file (####_####_###.sbx) (e.g. 1024_0726_000.sbx)
	f. (Optional, but highly recommended) Rigid file in the form (####_####_###_rigid.mat) (e.g. 1024_0726_000_rigid.mat) and associated sbx file

4. Please ensure you are in the folder level above the session folders.
	a. Folder tree example: (start here, animal folder) AH1024 > 0721/0722/0724/0725/0726/0727/0731/0801 (all same level session folders) atlab
	   > (sample, go into 0726) files as above
	b. At the start of each run, please ensure you are in the animal folder in Matlab
	c. Should the program crash midway due to a missing file, please ensure you have the file or skip the session, and ensure you return to the 
	   animal folder level (you should see all the sessions in the sidebar)

5. At the end of the program, feel free to dbcont to clear the workspace, or dbquit to retain the workspace and continue. At any time you may load 
   the generated multi-session file (e.g. 0721-0722-0724-0725-0726-0727-0731-0801) in the animal folder to continue.

6. At this point, please run the Session Struct (or corresponding Red_Green) to generate a readable struct for the next part

7. Run the Color Plots (or corresponding Red_Green) file to generate multi-session analysis. You will get graphs with all sessions overlaid in different shades of the same color.



*** Update to exclude trials based on sample size
*** save the pole down times in another array for better session analysis
