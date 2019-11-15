function Parameter = parseSedimentationInput(filename, Parameter)
%--------------------------------------------------------------------------
% parseConductivityInput - reads parameter from input file.
%
% Parameter = parseConductivityInput(filename, Parameter)
%--------------------------------------------------------------------------


%% reads input from the txt file
fileID = fopen(filename,'r');
in = textscan(fileID,'%s %f');
fclose(fileID);


%% Parameters
Parameter.nUnitCell = in{2}(1);           % number of unit cells 
Parameter.meshSpacingFraction = in{2}(2); % fraction of sqrt(alpha/pi)


%% Analyze input
fprintf('\nRunning testConductivity... \n');
fprintf('  n units cells =\t %d\n', Parameter.nUnitCell);
fprintf('  mesh fraction =\t %f\n', Parameter.meshSpacingFraction);

