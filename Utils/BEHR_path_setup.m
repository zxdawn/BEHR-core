function [  ] = BEHR_path_setup( )
%BEHR_path_setup Set up the paths to BEHR data.
%   This function will allow you to automatically generate the BEHR_paths
%   function which will be used to identify all the paths BEHR needs in
%   order to run.

% Paths to get:
%   Classes repo
%   Utils repo
%   SP mat dir
%   OMNO2 dir
%   OMPIXCOR dir
%   MYD06 dir
%   MCD43C3 dir
%   GLOBE database
%   BEHR mat dir
%   AMF tools path
%   NO2 profile path

file_server_ip = '128.32.208.13';

mydir = mfilename('fullpath');
[mydir, myname] = fileparts(mydir);

% Check that BEHR_paths.m does not already exist anywhere. If there's one
% in the right place, then ask if we want to overwrite it. If there's one
% somewhere else, the user will need to move or delete it so that we don't
% get two BEHR_paths.m files.
if exist(fullfile(mydir,'Constants','BEHR_paths.m'),'file')
    user_ans = questdlg('A BEHR_paths.m file already exists in BEHR/Utils/Constants. Do you wish to replace it?','File exists','Yes','No','No');
    if strcmpi(user_ans,'No')
        error('BEHR_setup:user_cancel','A BEHR_paths.m file already exists and you have chosen not to overwrite it.');
    end
elseif ~isempty(which('BEHR_paths.m'))
    uiwait(msgbox(sprintf('BEHR_paths.m exists on your Matlab search path, but at %s, not at BEHR/Utils/Constants. Please delete or move that version to BEHR/Utils/Constants.',which('BEHR_paths.m')),'modal'));
    error('BEHR_setup:user_cancel','A BEHR_paths.m file already exists outside of BEHR/Utils/Constants.\nDelete or move that file first, as MATLAB will be confused with multiple functions of the same name.')
end


% Get each path via the directory UI. If you ever need to add a new path,
% note that the field name will be how that path should be referred to in
% the call to BEHR_paths, i.e. to get paths.behr_mat_dir, one would call
% BEHR_paths('behr_mat_dir')
paths.classes = get_path('Classes repository','Please find the directory of the Matlab classes repository; it should contain at least the class JLLErrors. See the BEHR readme for information on how to clone that repo if necessary.');
paths.utils = get_path('Utils repository','Please find the directory of the general Matlab Utils repository (not the one inside the BEHR repo). See the BEHR readme for information on how to clone that repo if necessary.');
paths.sp_mat_dir = get_path('OMI SP .mat directory',sprintf('Please find the directory on the file server at %s containing OMI_SP_yyyymmdd.mat files. The file server should be mounted on your computer.',file_server_ip));
paths.omno2_dir = get_path('OMNO2 .he5 directory',sprintf('Please find the OMNO2 directory on the file server at %s. It should contain folders for each year. The file server should be mounted on your computer.',file_server_ip));
paths.ompixcor_dir = get_path('OMPIXCOR .he5 directory',sprintf('Please find the OMPIXCOR directory on the file server at %s. It should contain folders for each year. The file server should be mounted on your computer.',file_server_ip));
paths.myd06_dir = get_path('MYD06_L2 .hdf directory',sprintf('Please find the MYD06_L2 directory on the file server at %s. It should contain folders for each year. The file server should be mounted on your computer.',file_server_ip));
paths.mcd43c3_dir = get_path('MCD43C3 .hdf directory',sprintf('Please find the MCD43C3 directory on the file server at %s. It should contain folders for each year. The file server should be mounted on your computer.',file_server_ip));
paths.globe_dir = get_path('GLOBE directory',sprintf('Please find the GLOBE database directory on the file server at %s. It should contain files a10g through p10g and their .hdr files. The file server should be mounted on your computer.',file_server_ip));
paths.behr_mat_dir = get_path('BEHR .mat directory',sprintf('Please find the directory on the file server at %s containing OMI_BEHR_yyyymmdd.mat files. The file server should be mounted on your computer.',file_server_ip));
paths.amf_tools_dir = get_path('AMF tools directory','Please find the AMF_tools directory in the BEHR repository on your computer. It should contain the files damf.txt and nmcTmpYr.txt');
paths.no2_profile_path = get_path('NO2 profile directory',sprintf('Please find the directory on the file server at %s containing WRF-Chem output profiles. The file server should be mounted on your computer.',file_server_ip));
paths.website_staging_dir = get_path('Website staging directory',sprintf('Please find the directory on the file server at %s where published files are staged (WEBSITE/staging). The file server should be mounted on your computer.',file_server_ip));
fns = fieldnames(paths);

% Write the paths to a .m function
const_dir = fullfile(mydir,'Constants');
if ~exist(const_dir,'dir')
    mkdir(fullfile(mydir,'Constants'));
end

fid = fopen(fullfile(const_dir,'BEHR_paths.m'),'w');
% The opening function definition and comment block
fprintf(fid,'function p = BEHR_paths(pathname)\n');
fprintf(fid,'%% BEHR_paths() - paths used in the BEHR algorithm.\n%%\tAutomatically generated by %s on %s\n\n',myname,datestr(today,29));
fprintf(fid,'%%\tValid path names to request are:\n');
for a=1:numel(fns)
    fprintf(fid,'%%\t\t%s\n',fns{a});
end
fprintf(fid,'\n%%\tNote that this function should not be added to the BEHR git repo, as\n%%\tit must be specific to each person''s computer\n\n');

% A switch-case statement that will switch based on the path requested
fprintf(fid,'switch pathname\n');
for a=1:numel(fns)
    fprintf(fid,'\tcase ''%s''\n',fns{a});
    fprintf(fid,'\t\tp = ''%s'';\n',paths.(fns{a}));
end
% Error if pathname not recognized
fprintf(fid,'\totherwise\n');
fprintf(fid,'\t\terror(''BEHR_paths:unrecognized_path'',''Path %%s not recognized'',pathname)\n');
fprintf(fid,'end\n\n');

% This error will check if the requested path has changed since the paths
% file was created.
fprintf(fid,'if ~exist(p,''dir'')\n');
fprintf(fid,'\terror(''BEHR_paths:bad_path'',''The path %%s does not exist.\\nCheck that the file server is mounted, and that the path has not changed.\\nIf it has, rerun BEHR_path_setup or edit BEHR_paths manually'',p)\n');
fprintf(fid,'end\n');

fclose(fid);
fprintf('Be sure to add the folder BEHR/Utils/Constants to your MATLAB search path\n')

end

function folder = get_path(dir_type, dir_description)
uiwait(msgbox(dir_description, dir_type,'modal'));
folder = uigetdir();
if folder == 0
    error('BEHR_setup:user_cancel','User canceled setup, no directories have been saved.\nYou may need to re-run this script to finish BEHR setup.\n');
end

end