function BEHR_main(date_start, date_end, prof_mode)
%Josh Laughner <joshlaugh5@gmail.com>
%Based on BEHR_nwus by Ashley Russell (02/09/2012)
%Takes "OMI_SP_yyyymmdd.m" files produced by read_omno2_v_aug2012.m as it's
%main input.

%****************************%
% CONSOLE OUTPUT LEVEL - 0 = none, 1 = minimal, 2 = all messages, 3 = times %
% Allows for quick control over the amount of output to the console.
% Choose a higher level to keep track of what the script is doing.
DEBUG_LEVEL = 2;
%****************************%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PARALLELIZATION OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specifies whether the script is executing on a cluster; this must be set
% (globally) in the calling script.  This allows for the execution of code
% needed on the cluster (i.e. adding necessary folders to the Matlab path,
% opening a parallel pool) without running them on the local machine.  If
% onCluster hasn't been defined yet, set it to false.
global onCluster;
if isempty(onCluster);
    if DEBUG_LEVEL > 0; fprintf('Assuming onCluster is false\n'); end
    onCluster = false;
end

% Defined the number of threads to run, this will be used to open a
% parallel pool. numThreads should be set in the calling run script,
% otherwise it will default to 1.
global numThreads;
if isempty(numThreads)
    numThreads = 1;
end

% Cleanup object will safely exit if there's a problem
if onCluster
    cleanupobj = onCleanup(@() mycleanup());
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DEPENDENCIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%Add the 'Utils' folder and all subfolders to MATLAB's search path. Within
%the Git repository for BEHR, this is the /Utils folder.
mpath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(mpath,'..','Utils')));


% Add the paths needed to run on the cluster. Modify these manually if
% needed.
if onCluster;
    addpath(genpath('~/MATLAB/Classes'));
    addpath(genpath('~/MATLAB/Utils'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FILE LOCATIONS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The location of the directories to read or save data to.  If onCluster is
% true, these will need to be set in the runscript - I figured this would
% be easier than setting them as shell environmental variables and using
% getenv - JLL 15 Jan 2015

if onCluster
    global behr_mat_dir;
    global sp_mat_dir;
    global amf_tools_path;
    global no2_profile_path;
    
    % Verify the paths integrity.
    nonexistant = {};
    
    if ~exist(behr_mat_dir,'dir')
        nonexistant{end+1} = 'behr_mat_dir';
    end
    if ~exist(sp_mat_dir,'dir')
        nonexistant{end+1} = 'sp_mat_dir';
    end
    if ~exist(amf_tools_path,'dir')
        nonexistant{end+1} = 'amf_tools_path';
    end
    if ~exist(no2_profile_path,'dir')
        nonexistant{end+1} = 'no2_profile_path';
    end
    
    if numel(nonexistant)>0
        string_spec = [repmat('\n\t%s',1,numel(nonexistant)),'\n\n'];
        msg = sprintf('The following paths are not valid: %s Please double check them in the run file',string_spec);
        error(E.callError('bad_cluster_path',sprintf(msg,nonexistant{:})));
    end
    
    
else
    %This is the directory where the final .mat file will be saved. This will
    %need to be changed to match your machine and the files' location.
    behr_mat_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - Scale by SCD - lw 13.5 - 18-22 UTC';
    
    %This is the directory where the "OMI_SP_*.mat" files are saved. This will
    %need to be changed to match your machine and the files' location.
    sp_mat_dir = BEHR_paths('sp_mat_dir');
    
    %Add the path to the AMF_tools folder which contains rNmcTmp2.m,
    %omiAmfAK2.m, integPr2.m and others.  In the Git repository for BEHR, this
    %is the 'AMF_tools' folder.
    amf_tools_path = BEHR_paths('amf_tools_dir');
    
    %This is the directory where the NO2 profiles are stored. This will
    %need to be changed to match your machine and the files' location.
    %no2_profile_path = '/Volumes/share/GROUP/SAT/BEHR/Monthly_NO2_Profiles';
    no2_profile_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_BEHR_13lonwt_1822UTC';
end

%Store paths to relevant files
addpath(amf_tools_path)
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');
%****************************%

%****************************%
%Process all files between these dates, in yyyy/mm/dd format
%****************************%
if nargin < 2
    date_start='2013/06/10';
    date_end='2013/06/30';
    fprintf('BEHR_main: Used hard-coded start and end dates\n');
end

% restart from last file produced (true) or run entire time period (false)
restart = false;
%****************************%

% Which WRF profiles to use. Can be 'hourly', 'daily', 'monthly', or
% 'hybrid'. For convergence, should really be set to 'monthly'.
%****************************%
if exist('prof_mode','var')
    allowed_modes = {'hourly','daily','monthly','hybrid'};
    if ~ismember(prof_mode,allowed_modes)
        error('BEHR_main:bad_input','prof_mode (if given) must be one of %s',strjoin(allowed_modes,', '));
    end
    wrf_avg_mode = prof_mode;
else
    wrf_avg_mode = 'monthly';
    fprintf('BEHR_main: Used hard-coded wrf_avg_mode = %s\n',wrf_avg_mode);
end
%****************************%

%These will be included in the file name
%****************************%
satellite='OMI';
retrieval='BEHR';
%****************************%

%****************************%
% Which cloud product to use to calculate the AMF: OMI or MODIS
%****************************%
cloud_amf = 'omi';
cloud_rad_amf = 'omi';
%****************************%


% Check that all directories given ARE directories
if ~exist(behr_mat_dir,'dir')
    E.filenotfound(behr_mat_dir)
elseif ~exist(sp_mat_dir,'dir')
    E.filenotfound(sp_mat_dir)
elseif ~exist(amf_tools_path,'dir')
    E.filenotfound(amf_tools_path)
elseif ~exist(no2_profile_path,'dir')
    E.filenotfound(no2_profile_path)
end


% Create a parallel pool if one doesn't exist and we are on a cluster
if onCluster    
    if isempty(gcp('nocreate'))
        parpool(numThreads);
    end    
    n_workers = Inf;
else
    n_workers = 0;
end

datenums = datenum(date_start):datenum(date_end);

%parfor (j=1:length(datenums), n_workers)
for j=1:length(datenums)
    month=datestr(datenums(j),'mm');
    if DEBUG_LEVEL > 0; disp(['Processing data for ', datestr(datenums(j))]); end
    filename = sprintf('OMI_SP_%s_%s.mat',BEHR_version,datestr(datenums(j),'yyyymmdd'));
    
    if DEBUG_LEVEL > 1; disp(['Looking for SP file ',fullfile(sp_mat_dir,filename),'...']); end %#ok<PFGV> % The concern with using global variables in a parfor is that changes aren't synchronized.  Since I'm not changing them, it doesn't matter.
    if isequal(exist(fullfile(sp_mat_dir,filename),'file'),0)
        if DEBUG_LEVEL > 0; disp('No SP file exists for given day'); end
        continue
    else
        if DEBUG_LEVEL > 1; fprintf('\t ...Found.\n'); end
        S=load(fullfile(sp_mat_dir,filename)); %JLL 17 Mar 2014: Will load the variable 'Data' into the workspace
        Data=S.Data;
        
        for d=1:length(Data);
            % Data is initialized in read_omno2_v_aug2012 with a single 0
            % in the Longitude field.  Since points outside the lat/lons of
            % interest are removed completely, we should also check if all
            % points are gone.
            if numel(Data(d).Longitude)==1 || isempty(Data(d).Longitude);
                if DEBUG_LEVEL > 1; fprintf('  Note: Data(%u) is empty\n',d); end
                continue %JLL 17 Mar 2014: Skip doing anything if there's really no information in this data
            else
                if DEBUG_LEVEL>0; fprintf('  Swath %u of %s \n',d,datestr(datenums(j))); end
                c=numel(Data(d).Longitude);
                
                %Data(d).MODISAlbedo(isnan(Data(d).MODISAlbedo)==1)=0; %JLL 17 Mar 2014: replace NaNs with fill values
                %Data(d).GLOBETerpres(isnan(Data(d).GLOBETerpres)==1)=1013.0000;
                
                %JLL 17 Mar 2014: Load some of the variables from 'Data' to
                %make referencing them less cumbersome. Also convert some
                %to column vectors to work with rNmcTmp2 and rDamf2
                lon = Data(d).Longitude;
                lat = Data(d).Latitude;
                loncorns=Data(d).Loncorn;
                latcorns=Data(d).Latcorn;
                time = Data(d).Time;
                
                sza = Data(d).SolarZenithAngle;
                vza = Data(d).ViewingZenithAngle;
                phi = Data(d).RelativeAzimuthAngle;
                
                mon = str2double(month)*ones(size(Data(d).Latitude));
                pressure = [1020 1015 1010 1005 1000 990 980 970 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200];% 100 50 20 5];
                if DEBUG_LEVEL > 1; disp('   Interpolating temperature data'); end
                [temperature, tmpSAVE] = rNmcTmp2(fileTmp, pressure, lon, lat, mon); %JLL 17 Mar 2014: Interpolates temperature values to the pressures and lat/lon coordinates desired
                
                surfPres = Data(d).GLOBETerpres;
                albedo = Data(d).MODISAlbedo;
                
                surfPres(surfPres>=1013)=1013; %JLL 17 Mar 2014: Clamp surface pressure to sea level or less.
                cldPres = Data(d).CloudPressure(:);
                cldPres(cldPres>=1013)=1013; % JLL 13 May 2016: Also clamp cloud pressure. Whenever this is >1013, the AMF becomes a NaN because the lookup table cannot handle "surface" pressure >1013
                
                if DEBUG_LEVEL > 1; disp('   Calculating clear and cloudy AMFs'); end
                dAmfClr = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres); %JLL 18 Mar 2014: Interpolate the values in dAmf to the albedo and other conditions input
                cloudalbedo=0.8*ones(size(Data(d).CloudFraction)); %JLL 18 Mar 2014: Assume that any cloud has an albedo of 0.8
                dAmfCld = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres); %JLL 18 Mar 2014: Interpolate dAmf again, this time taking the cloud top and albedo as the bottom pressure
                
                
                pTerr = surfPres;
                pCld = cldPres;
                if strcmpi(cloud_amf,'omi')
                    cldFrac = Data(d).CloudFraction;
                else
                    cldFrac = Data(d).MODISCloud;
                end
                
                cldRadFrac = Data(d).CloudRadianceFraction;
                
                if DEBUG_LEVEL > 1; disp('   Reading NO2 profiles'); end
                [no2_bins, apriori_bin_mode] = rProfile_WRF(datenums(j), wrf_avg_mode, loncorns, latcorns, time, pTerr, pressure, no2_profile_path); %JLL 18 Mar 2014: Bins the NO2 profiles to the OMI pixels; the profiles are averaged over the pixel
                no2Profile1 = no2_bins;
                no2Profile2 = no2_bins;
                
                if DEBUG_LEVEL > 1; disp('   Calculating BEHR AMF'); end
                % We need the tropospheric slant column to compare against
                S_behr = Data(d).ColumnAmountNO2Trop .* Data(d).AMFTrop;
                
                % Now we need to compute the SCD derived from the WRF
                % profile.
                S_wrf_clr = apply_aks_to_prof(no2Profile1, pressure, dAmfClr, pressure, pTerr);
                S_wrf_cld = apply_aks_to_prof(no2Profile1, pressure, dAmfCld, pressure, pCld);
                % The total WRF slant column will be the sum of clear and
                % above cloud columns, weighted by the cloud radiance
                % fraction, since the influence on the detector should
                % relate the to amount of light at the top of the
                % atmosphere due to clear and cloudy scenes.
                S_wrf = (1 - cldRadFrac) .* S_wrf_clr + cldRadFrac .* S_wrf_cld;
                
                noGhost=0; ak=1;
                % Calculate the initial AMFs based on the direct WRF
                % profiles
                [amf_init, ~, ~, ~, ~, no2_prof_interp_init] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak); 
                Data(d).BEHRAMFTropInitial = amf_init;
                Data(d).BEHRColumnAmountNO2TropInitial = S_behr ./ amf_init;
                % Now, for each pixel, we're going to compare the WRF SCD
                % and the BEHR SCD and scale the WRF profile's boundary
                % layer so that the SCDs match.
                chemBLH = zeros(size(amf_init));
                scaling_flags = uint8(zeros(size(amf_init)));
                scaling_warnings = uint8(zeros(size(amf_init)));
                for i=1:numel(amf_init)
                    % Compute the boundary layer height, assuming that the
                    % boundary layer can be defined as where the [NO2]
                    % drops to 1/e^2 of its original value. This function
                    % returns a NaN if it cannot find anything satisfying
                    % that criteria.
                    no2_slice = no2Profile1(:,i);
                    if all(isnan(no2_slice))
                        scaling_flags(i) = bitset(scaling_flags(i),6);
                        continue
                    end
                    chemBLH(i) = find_bdy_layer_height(no2_slice, pressure, 'exp2', 'altispres', true);
                    if isnan(chemBLH(i))
                        scaling_flags(i) = bitset(scaling_flags(i), 3);
                        continue
                    elseif chemBLH(i) - -log(pTerr/1013)*7.4 > 2
                        scaling_warnings(i) = bitset(scaling_warnings(i), 3);
                    end
                    
                    % Will use this as a test to determine how good an
                    % assumption it is that most of the NO2 is in the
                    % boundary layer. Criteria were derived by examining a
                    % number of WRF profiles and seeing what ranges of
                    % values for this quantity produced questionable BL
                    % heights.
                    prof_bottom = find(~isnan(no2_slice),1,'first');
                    bl_factor = abs(no2_slice(prof_bottom) - nanmedian(no2_slice(:)))/no2_slice(prof_bottom);
                    if bl_factor < 0.6
                        scaling_flags(i) = bitset(scaling_flags(i), 4);
                        continue
                    elseif bl_factor >= 0.6 && bl_factor <= 0.75
                        scaling_warnings(i) = bitset(scaling_warnings(i), 2);
                    end
                    
                    %perdiff = (S_wrf(i) - S_behr(i)) ./ S_behr(i);
                    %if abs(perdiff) < 0.1
                    %        scaling_flags(i) = bitset(scaling_flags(i), 2);
                    %else
                        % See notes from 12 Jan 2016.
                        % Like the previously attempted convergence method,
                        % we will scale the boundary layer by the ratio of
                        % the observed to modeled boundary layers, assuming
                        % that the free troposphere in the model is
                        % reasonably accurate
                        
                        % The catch is that, because we attempted to
                        % calculate the visible slant column, we need to
                        % compute the clear and cloudy free trop SCD and
                        % add them together. However, while the lower
                        % integration limit for clear sky is obvious (BL
                        % pressure) the cloudy one's lower limit will be
                        % either the cloud pressure or BL pressure, which
                        % ever is smaller (higher altitude). Generally, we
                        % expect the BL height should not exceed the cloud
                        % top, so flag if this is true.
                        S_ft_clr_i = apply_aks_to_prof(no2Profile1(:,i), pressure, dAmfClr(:,i), pressure, chemBLH(i));
                        S_ft_cld_i = apply_aks_to_prof(no2Profile1(:,i), pressure, dAmfCld(:,i), pressure, min([chemBLH(i), pCld(i)]));
                        if chemBLH(i) < pCld(i)
                            scaling_flags(i) = bitset(scaling_flags(i),5);
                        end
                        S_ft_i = (1-cldRadFrac(i)) .* S_ft_clr_i + cldRadFrac(i) .* S_ft_cld_i;
                        
                        S_wrf_bl_i = S_wrf(i) - S_ft_i;
                        S_behr_bl_i = S_behr(i) - S_ft_i;
                        
                        pp = pressure > chemBLH(i);
                        no2_slice(pp) = no2_slice(pp) .* (S_behr_bl_i ./ S_wrf_bl_i);
                        
                        no2Profile1(:,i) = no2_slice;
                    %end
                end
                % I don't know why the profiles are duplicated in
                % omiAmfAK2, must be a holdover from previous code.
                no2Profile2 = no2Profile1;
                
                [amf_final, amfVis_final, ~, ~, ~, scattering_weights, avg_kernels, no2_prof_interp, sw_plevels] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1);
                behr_vcd_final = S_behr ./ amf_final;
                
                sz = size(Data(d).Longitude);
                len_vecs = size(scattering_weights,1);  % JLL 26 May 2015 - find out how many pressure levels there are. Will often be 30, but might change.
                % Need this to properly reshape the scattering weights, AKs, pressure levels, and profiles
                
                Data(d).BEHRAMFTrop = reshape(amf_final,sz); %JLL 18 Mar 2014: Save the resulting AMF of the pixel
                Data(d).BEHRAMFTropVisOnly = reshape(amfVis_final,sz);
                Data(d).BEHRColumnAmountNO2Trop = reshape(behr_vcd_final, sz);
                Data(d).BEHRGhostFraction = reshape(ghost_fraction,sz);
                Data(d).BEHRScatteringWeights = reshape(scattering_weights, [len_vecs, sz]);
                Data(d).BEHRAvgKernels = reshape(avg_kernels, [len_vecs, sz]);
                Data(d).BEHRNO2apriori = reshape(no2_prof_interp_init, [len_vecs, sz]);
                Data(d).BEHRNO2ScaledApriori = reshape(no2_prof_interp, [len_vecs, sz]);
                Data(d).BEHRaprioriMode = apriori_bin_mode;
                Data(d).BEHRChemBLH = chemBLH;
                Data(d).BEHRProfileScalingFlags = scaling_flags;
                Data(d).BEHRProfileScalingWarnings = scaling_warnings;
                Data(d).BEHRPressureLevels = reshape(sw_plevels, [len_vecs, sz]);
            end
        end
        
        b=length(Data);
        for z=1:b;
            if isfield(Data,'BEHRAMFTrop')==0 || isempty(Data(z).BEHRAMFTrop)==1;
                continue
            else
                Data(z).BEHRColumnAmountNO2Trop=Data(z).ColumnAmountNO2Trop.*Data(z).AMFTrop./Data(z).BEHRAMFTrop;
                Data(z).BEHRColumnAmountNO2TropVisOnly=Data(z).ColumnAmountNO2Trop.*Data(z).AMFTrop./Data(z).BEHRAMFTropVisOnly;
                % make sure fill values in the original column or AMF are
                % fill values in BEHR.
                Data(z).BEHRColumnAmountNO2Trop(Data(z).ColumnAmountNO2Trop < -1e29 | Data(z).AMFTrop < -30000) = nan; 
                Data(z).BEHRColumnAmountNO2TropVisOnly(Data(z).ColumnAmountNO2Trop < -1e29 | Data(z).AMFTrop < -30000) = nan; 
                if DEBUG_LEVEL > 0; fprintf('   BEHR [NO2] stored for swath %u\n',z); end
            end
        end
        
        
        addpath('/Users/Josh/Documents/MATLAB/BEHR/Utils/m_map'); %JLL 18 Mar 2014: Adds the path to the m_map toolkit, needed for hdf_quadrangle_5km_new
        
        %*********************************%
        %JLL 19 Mar 2014: These will be used to define the quadrangles -
        %the quads will be smaller than the OMI pixel, and multiple quads
        %will take on the value for the same (closest) OMI pixel.  By
        %keeping the quads' centers the same over all retrievals you wish
        %to average, this will allow easier averaging over multiple OMI
        %swaths. This is a form of oversampling.
        %*********************************%
        lonmin = -125;  lonmax = -65;
        latmin = 25;   latmax = 50;
        resolution = 0.05; resolution2 = 0.05;
        %*********************************%
        %
        if lonmin > lonmax %Just in case I enter something backwards...
            error(E.badinput('Lonmin is greater than lonmax'))
        elseif latmin > latmax
            error(E.badinput('Latmin is greater than latmax'))
        end
        
        %*********************************%
        %JLL 19 Mar 2014: Save all relevant values produced by add2grid to
        %a new structure called 'OMI'
        %*********************************%
        
        if DEBUG_LEVEL > 0; disp('  Preparing OMI structure'); end
        s=numel(Data);
        
        % Prepare the OMI data structure which will receive the gridded
        % data - this will be passed to the gridding functions to keep the
        % field names in the right order.
        OMI=struct('Date','','Longitude', [], 'Latitude', [], 'Time', [], 'ViewingZenithAngle', [], 'SolarZenithAngle', [], 'ViewingAzimuthAngle', [], 'SolarAzimuthAngle', [],...
            'RelativeAzimuthAngle', [], 'AMFStrat', [], 'AMFTrop',[], 'CloudFraction', [], 'CloudRadianceFraction', [], 'CloudPressure', [], 'ColumnAmountNO2', [],...
            'SlantColumnAmountNO2', [], 'ColumnAmountNO2Trop', [], 'ColumnAmountNO2TropStd',[],'ColumnAmountNO2Strat',[],'TerrainHeight', [], 'TerrainPressure', [], 'TerrainReflectivity', [], 'vcdQualityFlags',{{}},...
            'MODISCloud', [], 'MODISAlbedo', [], 'GLOBETerpres', [], 'XTrackQualityFlags', {{}}, 'Row', [], 'Swath', [], 'TropopausePressure', [], 'BEHRColumnAmountNO2Trop',[],...
            'BEHRAMFTrop', [], 'BEHRColumnAmountNO2TropVisOnly', [], 'BEHRAMFTropVisOnly', [], 'Count', [], 'Area', [], 'Areaweight', [], 'MapData', struct);
        % Matlab treats structures as matrices, so we can duplicate our
        % structure to have the required number of entries just like a
        % matrix.
        OMI = repmat(OMI,1,s);
        hh=0;
        for d=1:s;
            if Data(d).ViewingZenithAngle==0;
            elseif numel(Data(d).ViewingZenithAngle)==1;
                continue
            else
                if DEBUG_LEVEL > 1; fprintf('   Gridding data for swath %u\n',d); end
                hh=hh+1;
                % JLL 23 Jul 2015: temporary change to study wind effects.
                % Return to add2grid_BEHR when done.
                OMI(hh) = add2grid_BEHR_winds(Data(d),OMI(hh),resolution,resolution2,[lonmin, lonmax],[latmin, latmax]); %JLL 20 Mar 2014: Superimpose data to a grid determined by lat & lon min/max and resolution above. Default resolution is 0.05 degree
            end
        end
        
        % Clean up any unused elements in OMI
        OMI(hh+1:end) = [];

        savename = sprintf('%s_%s_%s_%s.mat',satellite,retrieval,BEHR_version,datestr(datenums(j),'yyyymmdd'));
        if DEBUG_LEVEL > 0; disp(['   Saving data as',fullfile(behr_mat_dir,savename)]); end
        saveData(fullfile(behr_mat_dir,savename),Data,OMI)
    end
end
end

function saveData(filename,Data,OMI)
save(filename,'OMI','Data')
end

function mycleanup()
err=lasterror;
if ~isempty(err.message)
    fprintf('MATLAB exiting due to problem: %s\n', err.message);
    if ~isempty(gcp('nocreate'))
        delete(gcp)
    end 

    exit(1)
end
end
