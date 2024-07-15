%% Function to get file name input from user
%   inputs: 
%     identifier = tag to narrow down number of files, eg '*.nod' (string) 
%     question = question to prompt user (string)

function [fname]=getfilename(question,identifier)

% Display question
disp(['oo  ' question '  oo'])
% List all files matching identifier
files=dir(identifier);

% Check some files exist
if isempty(files)
    disp(['oo    ' identifier ' does not exist, listing all files..'])
    files=dir('*');
end

for ii=1:length(files);
    disp(['oo    File ' int2str(ii) ' :: ' files(ii).name])
end
inp=input('oo    Enter Number of file to use (default=1) >> ');

% set fname to this file, or the default if input is empty (1st file)
if(isempty(inp))
    fname=files(1).name;
    disp(['oo    Default file name used : ' fname])
else
    fname=files(inp).name;
end

end


    