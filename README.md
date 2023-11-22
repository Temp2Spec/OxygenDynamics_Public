# OxygenDynamics_Public
Optical investigation of oxygen dynamics in the murine cortex

This is a collection of routines for the detection and analysis of spatiotemporal oxygen dynamics in the mouse barrel cortex 

OxygenDynamics_Wrapper imports a csv file with paths and metadata of data folders that contain motion corrected (and denoised) oxygen imaging data, and if applicable behavioural data. It calls two analysis scripts:

(1) OxygenDynamics_Master perfroms analysis of oxygen dynamics 

(2) OxygenDynamics_Behaviour performs analysis of behavioural data (body movement, pupilometry and air puff). 

Some of the events detected with OxygenDynamics_Master analysis script (Oxygen Sinks) require manual curation

OxygenDynamics_Sinks_Curation app takes the unrefined data produced by OxygenDynamics_Master and with a simple yes/no response the user decides whether a putative oxygen sink is real or not. After this a new "refined" oxygen sinks data table is exported

Finally

OxygenDynamics_Stats imports the same csv file OxygenDynamics_Wrapper used, finds and collates all the refined tables for oxygen sinks as well as imports the oxygen surges and behavioural data. Overall checks at the distribution for various metrics of oxygen dynamics and correlates these with behaviouralk metrics. 

Distributions and averages are exported as csv to be plotted in other software (e.g. graphpad)


