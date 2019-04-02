function [  ] = BEHR_publishing_uncertainty(varargin)
%BEHR_publishing_uncertainty Create the HDF files for BEHRUncert products
%   BEHR_publishing_uncertainty can accept a number of input parameters to alter its
%   behavior. All of these have default values that are set up so that
%   calling it without parameters will lead to standard behavior. The
%   parameters are:
%
%       'start': the first date to process, as a date number or a string
%       implicitly understood by datenum(). Default is '2005-01-01'
%
%       'end': the last date to process; same format requirements as
%       'start'. Default is today.
%
%       'output_type': one of the strings 'hdf' or 'txt', determines which
%       output format will be used. 'hdf' will produce HDF version 5 files.
%       Default is 'hdf'
%
%       'pixel_type': one of the strings 'native' or 'gridded', determines
%       whether the native pixels (i.e. the 'Data' structure) or the
%       gridded pixel (i.e. the 'OMI' structure) will be saved. Default is
%       'gridded'.
%
%       'reprocessed': a boolean (true or false). If true, this tells the
%       publishing algorithm to include fields that used in situ
%       measurements from the DISCOVER-AQ campaign as a priori profiles.
%       That is a specialized product that hasn't been updated in years.
%       Default is false.
%
%       'region': a string indicating which region to publish, must match
%       the directory structure in behr_paths.behr_mat_dir. Only used if
%       "mat_dir" is not specified. Default is 'us'.
%
%       'profile_mode': a string which a priori profiles' retrieval to use,
%       must match the directory structure in behr_paths.behr_mat_dir
%       (within each region). Only used if "mat_dir" is not specified.
%       Default is 'daily'.
%
%       'mat_dir': the directory from which to load the Matlab files with
%       BEHR output saved in the. If not given (or given as an empty
%       string) then behr_paths.BEHRUncertParaSubdir(region, para) is
%       called with the values of the "region" and "para"
%
%       'save_dir': the directory to which to save the resulting HDF or CSV
%       files. Default is the value returned by
%       behr_paths.website_staging_dir.
%
%       'organize': a boolean that indicates whether the output should go
%       directly in the save directory (false) or in a subdirectory named
%       behr_<pixel_type>-<output_type>-<behr version>, e.g.
%       "behr_native-hdf-v2-1C". Default is true.
%
%       'overwrite': controls the behavior of this function if one of the
%       files it is trying to output already exists. This is a number, a
%       negative value will cause it to ask you on at least the first file,
%       whether it continues to ask on successive files depends on your
%       response. 0 means do not overwrite, and a positive value means
%       always overwrite. Default is 0.
%
%       'DEBUG_LEVEL': a scalar number that controls the verbosity of this
%       function. 0 is minimum verbosity, higher numbers print more to the
%       screen.
%
%   This function can also be parallelized using the global variables
%   numThreads and onCluster.

E = JLLErrors;

p = inputParser;
p.addParameter('output_type', 'hdf');
p.addParameter('pixel_type', 'gridded');
p.addParameter('start', '2014-05-01');
p.addParameter('end', '2014-08-24');
p.addParameter('reprocessed', false);
p.addParameter('mat_dir', '');
p.addParameter('region', 'us');
p.addParameter('profile_mode', 'daily');
p.addParameter('save_dir', '');
p.addParameter('organize', true);
p.addParameter('overwrite', 0);
p.addParameter('DEBUG_LEVEL', 1);

p.parse(varargin{:});
pout = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SET OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Start and end date
start_date = pout.start;
end_date = pout.end;

% region
region = pout.region;

% Output type should be 'txt' or 'hdf'.  Text (csv) files are for native
% resolution only.
output_type = pout.output_type;

allowed_outtype = {'txt','hdf'};
if ~ismember(output_type,allowed_outtype)
    E.badinput('output_type must be one of %s',strjoin(allowed_outtype,', '));
end


% Set to 'native' to save the native OMI resolution pixels. Set to
% 'gridded' to save the 0.05 x 0.05 gridded data
pixel_type = pout.pixel_type;

allowed_pixtype = {'native','gridded'};
if ~ismember(pixel_type,allowed_pixtype)
    E.badinput('pixel_type must be one of %s',strjoin(allowed_pixtype,', '));
end

% Other options - reprocessed should be TRUE to include in situ fields
is_reprocessed = pout.reprocessed;
if ~isscalar(is_reprocessed) || ~islogical(is_reprocessed)
    E.badinput('REPROCESSED must be a scalar logical')
end

% Whether subdirectories should be created within the save directory,
% organized by pixel type, output type, and BEHR version.
organized_subdir = pout.organize;
if ~isscalar(organized_subdir) || ~islogical(organized_subdir)
    E.badinput('ORGANIZE must be a scalar logical')
end

% How to handle overwriting. 1 = overwrite, 0 = don't overwrite, -1 = ask.
overwrite = pout.overwrite;

DEBUG_LEVEL = pout.DEBUG_LEVEL;
if ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL)
    E.badinput('DEBUG_LEVEL must be a scalar number')
end

% File locations
mat_file_dir = pout.mat_dir;
save_dir = pout.save_dir;

% iterate paras
F = dir(behr_paths.BEHRUncertSubdir(region));
% Keep only directories that are not '.' and '..'
file_names = {F.name};
xx_keep = cellfun(@(x) ~regcmp(x, '\.{1,2}'), file_names) & [F.isdir];
xx_keep = xx_keep & ~strcmp(file_names, 'BaseCase');
paras = {F(xx_keep).name};

for i_field = 1:numel(paras)
    para = paras{i_field};
    disp (para);

    mat_file_dir = behr_paths.BEHRUncertParaSubdir(region, para);
    save_dir = fullfile(behr_paths.website_staging_dir(),'uncertainty',para);

    % Check that the directories exist like this so that a single error message
    % describes if both directories don't exist - handy for running on the
    % cluster so that you don't wait forever for the job to start, only to have
    % it fail b/c you forgot to make the output directory.
    dirs_dne = {}; 
    if ~exist(mat_file_dir,'dir')
        dirs_dne{end+1} = sprintf('mat_file_dir (%s)', mat_file_dir);
    end
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
        %dirs_dne{end+1} = sprintf('save_dir (%s)', save_dir);
    end
    if ~isempty(dirs_dne)
        E.dir_dne(dirs_dne);
    end

    BEHR_publishing_sub('start',start_date,'end',end_date,'region', region,'profile_mode',pout.profile_mode,...
    'mat_dir',mat_file_dir,'save_dir',save_dir,'overwrite',overwrite);

end
