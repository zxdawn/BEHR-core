function [ pix_profs ] = rProfile_pickWRF( date_in, lons, lats, sat_scds, sw_clr, sw_cld, surf_pres, cld_pres, cld_rad_frac, pressures, wrf_output_path )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%
%%%%% CONSTANTS %%%%%
%%%%%%%%%%%%%%%%%%%%%

% How far in kilometers a profile is allowed to be from the pixel center
% and still be chosen for that pixel
max_dist_km = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure that the input for the date can be understood by MATLAB as a
% date
try
    date_num_in = datenum(date_in);
catch err
    if strcmp(err.identifier,'MATLAB:datenum:ConvertDateString')
        E.badinput('%s could not be understood by MATLAB as a date')
    else
        rethrow(err);
    end
end

% LONS, LATS, SAT_SCDS, SURF_PRES, CLD_PRES, and CLD_RAD_FRAC should all be
% 2D arrays of the same size The two scattering weight arrays should have
% the vectors of scattering weights along the first dimension, and the next
% two dimensions should be the same size as LONS
sz = size(lons);
if ~isequal(size(lats), sz)
    E.badinput('LATS must be the same size as LONS')
end

if ~isequal(size(sat_scds), sz)
    E.badinput('SAT_SCDS must be the same size as LONS')
end

if ~isequal(size(cld_rad_frac), sz)
    E.badinput('CLD_RAD_FRAC must be the same size as LONS')
end

if ~isequal(size(surf_pres), sz)
    E.badinput('SURF_PRES must be the same size as LONS')
end

if ~isequal(size(cld_pres), sz)
    E.badinput('CLD_PRES must be the same size as LONS')
end

% Pressures should be the vector of standard pressures used in the OMI
% retrievals

if ~isvector(pressures)
    E.badinput('PRESSURES must be a vector')
elseif ~iscolumn(pressures)
    pressures = pressures';
end 

% The scattering weight arrays need to have the vertical coordinate along
% the first dimension

sz_sw_clr = size(sw_clr);
sz_sw_cld = size(sw_cld);
n_pres = numel(pressures);
if ~isequal(sz_sw_clr, [n_pres, sz]) || ~isequal(sz_sw_cld, [n_pres, sz])
    E.badinput('size(SW_CLR) and size(SW_CLD) must be equal to [numel(PRESSURES), size(LONS)]');
end

% Check that the WRF data path exists and has the expected type of file in
% it

if ~exist(wrf_output_path, 'dir')
    E.badinput('WRF_OUTPUT_PATH %s does not exist', wrf_output_path)
end

wrf_files = dir(fullfile(wrf_output_path, 'WRF_BEHR_monthly*.nc'));
if isempty(wrf_files)
    E.filenotfound(sprintf('WRF_BEHR_monthly*.nc files in %s', wrf_output_path));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

pix_profs = nan([numel(pressures), size(lons)]);
pix_model_scds = nan(size(lons));

[wrf_no2, wrf_pres, wrf_lon, wrf_lat] = load_wrf_vars(date_in, wrf_output_path);

% Rearrange WRF arrays so that the vertical coordinate is first, will make
% subsetting easy
wrf_no2 = permute(wrf_no2, [3 1 2]);
wrf_pres = permute(wrf_pres, [3 1 2]);

wrf_no2 = interp_wrf_to_pres_vec(wrf_no2, wrf_pres, pressures);

for a=1:numel(lons)
    [closest_no2, closest_lon, closest_lat] = subset_profiles_by_distance(wrf_no2, wrf_lon, wrf_lat, lons(a), lats(a), max_dist_km);
    closest_scds = calc_wrf_scds(closest_no2, pressures, sw_clr(:,a), sw_cld(:,a), surf_pres(a), cld_pres(a), cld_rad_frac(a)); 
    [~, i] = min(abs(sat_scds(a) - closest_scds));
    pix_profs(:, a) = closest_no2(:, i);
    pix_model_scds(a) = closest_scds(i);
end

pix_profs = remove_below_surface_no2(pix_profs, pressures, surf_pres);

end

function [wrf_no2, wrf_pres, wrf_lon, wrf_lat] = load_wrf_vars(date_in, wrf_output_path)

E=JLLErrors;
E.addCustomError('ncvar_not_found','The variable %s is not defined in the file %s. Likely this file was not processed with (slurm)run_wrf_output.sh, or the processing failed before writing the calculated quantites.');

% Find the file for this day or this month
year_in = year(date_in);
month_in = month(date_in);

file_pat = sprintf('WRF_BEHR_monthly_%04d-%02d-*.nc', year_in, month_in);


F = dir(fullfile(wrf_output_path,file_pat));
if numel(F) < 1
    E.filenotfound(file_pat);
elseif numel(F) > 1
    E.toomanyfiles(file_pat);
else
    wrf_info = ncinfo(fullfile(wrf_output_path,F(1).name));
end

wrf_vars = {wrf_info.Variables.Name};

% Load NO2 and check what unit it is - we'll use that to convert to
% parts-per-part later. Make the location of units as general as possible
try
    wrf_no2 = ncread(wrf_info.Filename, 'no2');
catch err
    if strcmp(err.identifier,'MATLAB:imagesci:netcdf:unknownLocation')
        E.callCustomError('ncvar_not_found','no2',F(1).name);
    else
        rethrow(err);
    end
end
ww = strcmp('no2',wrf_vars);
no2_attr = {wrf_info.Variables(ww).Attributes.Name};
aa = strcmp('units',no2_attr);
wrf_no2_units = wrf_info.Variables(ww).Attributes(aa).Value;

if isempty(wrf_no2_units)
    if DEBUG_LEVEL > 1; fprintf('\tWRF NO2 unit not identified. Assuming ppm\n'); end
    wrf_no2_units = 'ppm';
end

% Convert to be an unscaled mixing ratio (parts-per-part)
wrf_no2 = convert_units(wrf_no2, wrf_no2_units, 'ppp');

% Load the remaining variables
try
    varname = 'pres';
    wrf_pres = ncread(wrf_info.Filename, varname);
    varname = 'XLONG';
    wrf_lon = ncread(wrf_info.Filename, varname);
    varname = 'XLAT';
    wrf_lat = ncread(wrf_info.Filename, varname);
catch err
    if strcmp(err.identifier,'MATLAB:imagesci:netcdf:unknownLocation')
        E.callCustomError('ncvar_not_found',varname,F(1).name);
    else
        rethrow(err);
    end
end
end

function [closest_no2, closest_lon, closest_lat] = subset_profiles_by_distance(wrf_no2, wrf_lon, wrf_lat, pix_lon, pix_lat, max_dist)

E = JLLErrors;

if ndims(wrf_no2) > 3
    E.notimplemented('4D WRF NO2 or PRES array')
end

% Fairly slow: 0.2 sec/pixel
% prof_dist = nan(size(wrf_lon));
% for a=1:numel(wrf_lon)
%     prof_dist(a) = m_lldist([wrf_lon(a), pix_lon], [wrf_lat(a), pix_lat]);
% end

% Much much faster, 0.0002 sec/pixel to just take the absolute difference
% in degrees, however this means fewer profiles will be available at higher
% latitudes. To account for this, we scale the longitude difference by the
% cosine of the average latitude. It's a very ad hoc correction, but the
% exact distance isn't terribly important as long as we get a similar
% number of pixels. See notes from 13 Jan 2017 for reasoning.
del_lon = wrf_lon - pix_lon;
mean_lat = (wrf_lat + pix_lat) ./ 2;
prof_dist = sqrt((del_lon ./ cosd(mean_lat)).^2 + (wrf_lat - pix_lat).^2);

xx = prof_dist <= max_dist/111;

if sum(xx) < 1
    %E.callError('rProfile:subset_by_dist','No profiles found within %f km of the pixel at lon %f, lat %f', max_dist, pix_lon, pix_lat);
    closest_no2 = nan(size(wrf_no2,1),1);
    closest_lon = nan;
    closest_lat = nan;
else
    closest_no2 = wrf_no2(:,xx);
    closest_lon = wrf_lon(xx);
    closest_lat = wrf_lat(xx);
end
end

function wrf_no2_interp = interp_wrf_to_pres_vec(wrf_no2, wrf_pres, pressures)
% Interpolate all the NO2 profiles to the input pressures, then average
% them together. Extrapolate so that later we can be sure to have one
% bin below the surface pressure for omiAmfAK2 and integPr2.
% Interpolate in log-log space to account for the exponential
% dependence of pressure on altitude and the often exponential decrease
% of concentration with altitude.
sz = size(wrf_no2);
wrf_no2_interp = nan([numel(pressures), sz(2:end)]);

for a=1:prod(sz(2:end))
    interp_no2 = interp1(log(wrf_pres(:,a)), log(wrf_no2(:,a)), log(pressures), 'linear', 'extrap');
    interp_no2 = exp(interp_no2);
    wrf_no2_interp(:,a) = interp_no2;
end
end

function wrf_scds = calc_wrf_scds(wrf_no2, pressures, sw_clr, sw_cld, surf_pres, cld_pres, cld_rad_frac)
S_wrf_clr = nan(1, size(wrf_no2,2));
S_wrf_cld = nan(1, size(wrf_no2,2));
for a=1:size(wrf_no2,2)
    S_wrf_clr(a) = apply_aks_to_prof(wrf_no2(:,a), pressures, sw_clr, pressures, surf_pres);
    S_wrf_cld(a) = apply_aks_to_prof(wrf_no2(:,a), pressures, sw_cld, pressures, cld_pres);
end
wrf_scds = (1 - cld_rad_frac) .* S_wrf_clr + cld_rad_frac .* S_wrf_cld;
end

function profs = remove_below_surface_no2(profs, pressures, surf_pres)
for a=1:numel(surf_pres)
    last_below_surf = find(pressures > surf_pres(a),1,'last')-1;
    profs(1:last_below_surf, a) = nan;
end
end
