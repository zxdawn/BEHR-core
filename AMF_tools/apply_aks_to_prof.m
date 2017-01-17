function [ column_out ] = apply_aks_to_prof( prof_no2, prof_pres, aks, aks_pres, surf_p )
%APPLY_AKS_TO_WRF_PROF Applies averaging kernel (or scattering weights) to a WRF NO2 profile
%   V = APPLY_AKS_TO_PROF( PROF_NO2, PROF_PRES, AKS, AKS_PRES, SURF_P,
%   TROP_P ) applies the averaging kernel or scattering weight vector AKS
%   to the profile PROF_NO2. PROF_NO2 should be in unscaled mixing ratio.
%   (If AKs is the averaging kernel, then the output will be a vertical
%   column "retrived" from the modeled profile; if it is scattering
%   weights, then the output will be the modeled slant column corresponding
%   to the modeled vertical profile.) PROF_PRES should be the pressure
%   coordinate corresponding to the profiles, and AKS_PRES the pressure
%   coordinate corresponding to the AKS vector.
%
%   This function assumes that each profile and AK vector is a slice along
%   the first dimension of those variables, i.e. PROF_NO2(:,1) and AKS(:,1)
%   are the first profile and AK vector. If arrays are provided, then all
%   must be the same size. If a vector is provided for either of the
%   pressure inputs, it will be used for all profile and AK vectors.
%
%   This function will interpolate the profile to the AK pressure levels
%   then multiply the AKS vector element-wise by the interpolated profile.
%   It will then integrate from SURF_P to the top of the interpolated
%   vector using integPr2.

E = JLLErrors;

% How similar the pressures have to be to skip interpolating
pres_diff_threshold = 0.1;

% Input checking
prof_pres_isvec = isvector(prof_pres);
aks_pres_isvec = isvector(aks_pres);

if ~isequal(size(prof_no2), size(aks))
    E.badinput('PROF_NO2 and AKS must be the same size array');
end
if prof_pres_isvec
    prof_pres = prof_pres(:); % force to be column vector
    if numel(prof_pres) ~= size(prof_no2,1)
        E.badinput('If a vector PROF_PRES must be the same length as the first dimension of PROF_NO2');
    end
else
    if ~isequal(size(prof_no2), size(prof_pres))
        E.badinput('If an array, PROF_PRES must be the same size as PROF_NO2');
    end
end
if aks_pres_isvec
    aks_pres = aks_pres(:); % force to be column vector
    if numel(aks_pres) ~= size(aks,1)
        E.badinput('If a vector AKS_PRES must be the same length as the first dimension of AKS');
    end
else
    if ~isequal(size(aks), size(aks_pres))
        E.badinput('If an array, AKS_PRES must be the same size as AKS');
    end
end

sz = size(prof_no2);
% This kludge is needed to compare the sizes of the profile and surf_p
% arrays properly, since Matlab's size() always returns at least a two
% element vector, sz always needs to have at least 3 elements.
if numel(sz) == 2; sz = [sz, 1]; end
if ~isequal(size(surf_p), sz(2:end)) && ~isscalar(surf_p) && prod(sz(2:end)) > 1
    E.badinput('SIZE(SURF_P) must equal the size of the second and higher dimensions of PROF_NO2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate the model profiles to the AK pressures. We do it this way
% rather than AKS to profile pressure because AKs can change rapidly with
% altitude in a way not capture by the discrete levels given.
%
% Interpolate in log-log space because both concentration and pressure tend
% to have an exponential relationship to altitude.

sz = size(prof_no2);
column_out = nan(sz(2:end));

for a=1:prod(sz(2:end))
    if prof_pres_isvec
        a_pp = 1;
    else
        a_pp = a;
    end
    if aks_pres_isvec
        a_ap = 1;
    else
        a_ap = a;
    end
    
    if numel(prof_pres(:, a_pp)) == numel(aks_pres(:, a_ap)) && all(abs(prof_pres(:, a_pp) - aks_pres(:, a_ap)) < pres_diff_threshold)
        % The pressures are similar; we don't need to interpolate (they've
        % probably already been interpolated)
    else
        tmp_prof = exp(interp1(log(prof_pres(:,a_pp)), log(prof_no2(:,a)), log(aks_pres(:,a_ap))));
        prof_no2(:,a) = tmp_prof;
    end
    
    column_out(a) = integPr2(prof_no2(:,a) .* aks(:,a), aks_pres(:,a_ap), surf_p(a));
end



end

