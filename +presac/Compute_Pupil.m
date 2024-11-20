function Info = Compute_Pupil(Info,Exp,H,plotsteps)
%******* 
%******  function Info = Compute_Pupil(Info,Exp,H)
%****
%*** Inputs:  Info - struct with information on analysis (unit name, etc)
%***          Exp - experiment struct with all data per session
%***                if Exp = [], then just plot results in Info
%***                otherwise only compute the results for return
%***          H - if integer, produce plots of the results, if figure
%***                  then gives a cell struct with panel handles to plot
%***
%*** Outputs: Info - it updates fields in the Info struct to include 
%***                  information about trial inclusion parameters
%***



% Define parameters for doing the PFR analysis
% We will want to define the time windows of interest



% Find Pupil diameter data in the Exp stucture

% Create a new structue stored in Info

% Plot pupil diameter trial by trial with x and y eye position (for time window of interest)






end

% Create seperate functions for plotting 