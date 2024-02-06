% 

% The data we want to export to case format
% flow data is NX x NY x NZ x NTime x 3
flow = zeros(100,120,20,30,3);
% orientation of the data
tx = diag(ones(1,4));  % this is a 4x4 affine matrix defining rotation, translation, scaling to map matlab coordinates onto real coordiantes

sz = size(flow);
% what is the time of each timeframe
timesteps = (1:sz(4))*0.03;   %30 ms steps

%Start a new case file
cfile = CaseFile(timesteps);

%create a flow part
flowpart = cfile.CreatePart('flow');

% Create the coordiantes of all indicies
flowpart.CreateBlockGeometry(tx,sz(1:3));

% add a vector variables that are wanted
flowpart.SetVariableTime('flow',flow);


% add a constant varialbe (no time)
% Ensight does not like this, but it still works 
flow_peak = max(sqrt(sum(flow.^2,5)),[],4); 
flowpart.SetVariableTime('magnitude_flow',flow_peak(:)); % NOTE This is set to 1D, to indicate no time information


% to define 3D shapes, you cannot use the structured block Part defined above
% define a new part
% add connectivity, this is a simple square based on matlab index
shapePart = cfile.CreatePart('shape');

vertacies = [1         1        1; 
             1         1    sz(3); 
             1     sz(2)    sz(3);
             1     sz(2)        1;
             sz(1)     1        1;
             sz(1)     1    sz(3); 
             sz(1) sz(2)    sz(3);
             sz(1) sz(2)        1];
corners = [1 2 3 4 5 6 7 8];
% NOTE: cor connectivity it may be preferable not to use BlockGeometry but
% set you own using Part.SetGeometry, this however makes all variables
% point variables (not fields)
for iT = 1:sz(4) % for same connectivity for all time frame, 
    shapePart.SetGeometry(iT,vertacies(:,1),vertacies(:,2),vertacies(:,3));
    shapePart.SetConnectivity(iT,corners,'hexa8');
end

% add anouther part in the same way, perhaps a cine image
% 
mkdir('test');
cfile.WriteFile('test/myTestCase.case');

