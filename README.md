# OxygenDynamics
Optical investigation of oxygen dynamics in the murine cortex

This is a collection of analysis scripts for the detection and analysis of oxygen sinks and surges in the mouse barrel cortex

OxygenDynamics_Wrapper imports and csv file with paths and metadata of data folders that contain motion corrected (and denoised data) and if applicable behavioural data. 
It calls two analysis scripts:

OxygenDynamics_Master perfroms analysis of oxygen sinks and surges
OxygenDynamics_Behaviour performs analysis of behavioural data (body movement, pupilometry and air puff)
After this step oxygen sink data need to be manually curated

OxygenDynamics_Sinks_Curation app takes the unrefined data produced by OxygenDynamics_Master and with a simple yes/no response the user decides whether a putative oxygen sink is real or not.
After this a new "refined" oxygen sinks data table is exported

OxygenDynamics_Stats imports the same csv file OxygenDynamics_Wrapper used, finds and collates all the refined tables for oxygen sinks as well as imports the oxygen surges and behavioural data.
Overall checks at the distribution for various metrics of oxygen dynamics and correlates them with behaviour. 
Finally a generalized mixed-effects modelling approach is used for hypotheses testing.

Distributions and averages are exported as csv to be plotted in other software (e.g. graphpad)
