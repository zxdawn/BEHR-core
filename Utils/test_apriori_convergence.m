function [ converged_profs, true_prof, start_prof, amfs_converged, amfs_true, amfs_start ] = test_apriori_convergence( true_prof, true_pres, true_lon, true_lat, start_prof, start_pres, init_wrf_path, sza, vza, raa, alb, surfpres, cldradfrac, cldpres, date_in )
%TEST_APRIORI_CONVERGENCE Ideal case of scaling a monthly profile to match daily
%   CONVERGED_PROFS = TEST_APRIORI_CONVERGENCE( TRUE_PROF, TRUE_PRES,
%   START_PROF, START_PRES, SZA, VZA, RAA, ALB, SURFP, CLDRADFRAC, CLDPRES)
%   - applies the scaling approach to START_PROF to try to make them look
%   like TRUE_PROF. TRUE_PRES and START_PRES are the corresponding pressure
%   coordinates. SZA, VZA, RAA, ALB, SURFP are the inputs for the clear sky
%   AMF, CLDRADFRAC and CLDPRES determine the cloudy AMF. The profiles and
%   pressures should be given as 3D arrays with the vertical coordinate as
%   the third dimension (so they can be read directly from WRF files). The
%   other inputs should all be scalars. DATE_IN can specify a date in
%   yyyy-mm-dd or datenumber format. This also assumes that the true and
%   starting profiles are on the same grid.

E = JLLErrors;
if ndims(true_prof) ~= 3
    E.badinput('TRUE_PROF must be a 3d array')
end
sz = size(true_prof);
if ~isequal(size(true_pres), sz)
    E.badinput('TRUE_PRES must be the same size as TRUE_PROF (a 3d array)')
end
if ~isequal(size(true_lon), sz(1:2))
    E.badinput('TRUE_LON must be the same size as the first two dimensions of TRUE_PROF')
end
if ~isequal(size(true_lat), sz(1:2))
    E.badinput('TRUE_LAT must be the same size as the first two dimensions of TRUE_PROF')
end

if ~ischar(init_wrf_path)
    E.badinput('INIT_WRF_PATH must be a string')
elseif ~exist(init_wrf_path, 'dir')
    E.badinput('INIT_WRF_PATH must be a directory')
else
    F = dir(fullfile(init_wrf_path, 'WRF_BEHR_*.nc'));
    if isempty(F)
        E.badinput('No WRF_BEHR files found in %s', init_wrf_path);
    end
end

if ~isscalar(sza)
    E.badinput('SZA must be a scalar')
end
if ~isscalar(vza)
    E.badinput('VZA must be a scalar')
end
if ~isscalar(raa)
    E.badinput('RAA must be a scalar')
end
if ~isscalar(alb)
    E.badinput('ALB must be a scalar')
end
if ~isscalar(surfpres)
    E.badinput('SURFP must be a scalar')
end
if ~isscalar(cldradfrac)
    E.badinput('CLDRADFRAC must be a scalar')
end
if ~isscalar(cldpres)
    E.badinput('CLDPRES must be a scalar')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

bplevs = behr_pres_levels;
bplevs = bplevs(:);

% First, interpolate the two profile sets to the standard pressure levels.
% Also reorders the dimensions so that vertical is first
true_prof = interp_no2_to_pres(true_prof, true_pres);
start_prof = interp_no2_to_pres(start_prof, start_pres);

% Calculate the AMF scattering weights for the inputs
mon = month(date_in);
mon_mat = repmat(mon, size(true_prof,2), size(true_prof,3));
sza_mat = repmat(sza, size(true_prof,2), size(true_prof,3));
vza_mat = repmat(vza, size(true_prof,2), size(true_prof,3));
raa_mat = repmat(raa, size(true_prof,2), size(true_prof,3));
alb_mat = repmat(alb, size(true_prof,2), size(true_prof,3));
surfpres_mat = repmat(surfpres, size(true_prof,2),size(true_prof,3));
cldpres_mat = repmat(cldpres, size(true_prof,2),size(true_prof,3));
cldradfrac_mat = repmat(cldradfrac, size(true_prof,2),size(true_prof,3));

[dAmfClr, dAmfCld, temperature] = compute_behr_sc_weights(true_lon, true_lat, mon_mat, sza_mat, vza_mat, raa_mat, alb_mat, surfpres_mat, cldpres_mat);

alpha = 1 - 0.003 * (temperature - 220);   % temperature correction factor vector
alpha_i=max(alpha,0.1);
alpha = min(alpha_i,10);
 
% Don't think this made it into the actual scaling code! Be sure to update.
dAmfClr = dAmfClr .* alpha;
dAmfCld = dAmfCld .* alpha;

% This should be handled by integPr2 in apply_aks_to_prof
% dAmfClr(bplevs > surfpres) = 1e-30;
% dAmfCld(bplevs > cldpres) = 1e-30;


% Calculate the SCDs for the true and starting profiles
S_wrf_clr_true = apply_aks_to_prof(true_prof, bplevs, dAmfClr, bplevs, surfpres_mat);
S_wrf_cld_true = apply_aks_to_prof(true_prof, bplevs, dAmfCld, bplevs, cldpres_mat);
S_wrf_true = (1 - cldradfrac) .* S_wrf_clr_true + cldradfrac .* S_wrf_cld_true;


converged_profs = rProfile_pickWRF(date_in, true_lon, true_lat, S_wrf_true, dAmfClr, dAmfCld, surfpres_mat, cldpres_mat, cldradfrac_mat, bplevs, init_wrf_path);


amfs_true = omiAmfAK2(surfpres_mat, cldpres_mat, cldradfrac_mat, cldradfrac_mat, bplevs, dAmfClr, dAmfCld, temperature, true_prof);
amfs_start = omiAmfAK2(surfpres_mat, cldpres_mat, cldradfrac_mat, cldradfrac_mat, bplevs, dAmfClr, dAmfCld, temperature, start_prof);
amfs_converged = omiAmfAK2(surfpres_mat, cldpres_mat, cldradfrac_mat, cldradfrac_mat, bplevs, dAmfClr, dAmfCld, temperature, converged_profs);


end

function no2 = interp_no2_to_pres(no2_in, pres_in)
% Interpolate a given NO2 array to the standard pressure levels, also make
% the vertical coordinate the first dimension of the array

plevs = behr_pres_levels;
no2 = nan(numel(plevs), size(no2_in,1), size(no2_in,2));

for a=1:size(no2,2)
    for b=1:size(no2,3)
        tmp_no2 = squeeze(no2_in(a,b,:));
        tmp_pres = squeeze(pres_in(a,b,:));
        interp_no2 = interp1(log(tmp_pres), log(tmp_no2), log(plevs), 'linear', 'extrap');
        no2(:,a,b) = exp(interp_no2);
    end
end


end


