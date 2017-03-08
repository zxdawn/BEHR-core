function [ dAmfClr, dAmfCld, temperature ] = compute_behr_sc_weights( lon, lat, temperature_month, sza, vza, raa, albedo, surfPres, cldPres )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

sz = size(lon);
if ~ismatrix(lon)
    E.badinput('LON must be a matrix or vector')
elseif ~isequal(sz, size(lat))
    E.badinput('LAT must be a matrix or vector the same size as LON');
elseif ~isequal(sz, size(temperature_month))
    E.badinput('TEMPERATURE_MONTH must be a matrix or vector the same size as LON');
elseif ~isequal(sz, size(sza))
    E.badinput('SZA must be a matrix or vector the same size as LON');
elseif ~isequal(sz, size(vza))
    E.badinput('VZA must be a matrix or vector the same size as LON');
elseif ~isequal(sz, size(raa))
    E.badinput('RAA must be a matrix or vector the same size as LON');
elseif ~isequal(sz, size(albedo))
    E.badinput('ALBEDO must be a matrix or vector the same size as LON');
elseif ~isequal(sz, size(surfPres))
    E.badinput('SURFPRES must be a matrix or vector the same size as LON');
elseif ~isequal(sz, size(cldPres))
    E.badinput('CLDPRES must be a matrix or vector the same size as LON');
end

% Find the lookup tables for temperature and pressure
fileTmp = fullfile(behr_repo_dir,'AMF_tools','nmcTmpYr.txt');
fileDamf = fullfile(behr_repo_dir, 'AMF_tools', 'damf.txt');
pressure = behr_pres_levels();

if ~exist(fileTmp,'file')
    E.filenotfound('Could not find the temperature profile file at %s', fileTmp);
elseif ~exist(fileDamf,'file')
    E.filenotfound('Could not find the temperature profile file at %s', fileDamf)
end 

temperature = rNmcTmp2(fileTmp, pressure, lon, lat, temperature_month); %JLL 17 Mar 2014: Interpolates temperature values to the pressures and lat/lon coordinates desired

% Surface pressures greater than (i.e. below) sea level produce NaNs in the
% scattering weights, because the loop up table does not include points
% with surface pressure beyond 1013 hPa. Clamping the surface pressures to
% 1013 ensures that reasonable scattering weights are chosen.

if any(surfPres(:) > 1013)
    warning('Surface pressures greater than 1013 hPa found, clamping to 1013')
end
surfPres(surfPres > 1013) = 1013; 

if any(cldPres(:) > 1013)
    warning('Cloud pressures greater than 1013 hPa found, clamping to 1013')
end
cldPres(cldPres > 1013) = 1013;

% Get the scattering weights
dAmfClr = rDamf2(fileDamf, pressure, sza, vza, raa, albedo, surfPres); 
cloudalbedo=0.8*ones(size(lon)); % Assume that any cloud has an albedo of 0.8
dAmfCld = rDamf2(fileDamf, pressure, sza, vza, raa, cloudalbedo, cldPres);

end

