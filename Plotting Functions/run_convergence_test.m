function [ profiles, amfs, amf_settings ] = run_convergence_test( date_in )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

sza = 0;
vza = 0;
raa = 0;
alb = 0.04;
surfp = 1013;

cldradfrac = 0.0;
cldp = 600;

lon = -95;
lat = 37.5;

wrf_path_daily = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_BEHR/hourly-14.0-lonwt-1822UTC';
wrf_path_monthly = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_BEHR/monthly-13.5-lonwt-1822UTC';

% Load no2 and pressure arrays, convert no2 to mixing ratio (unscaled) from
% ppmv.
file_pattern = 'WRF_BEHR_%s_%s.nc';
daily_file = fullfile(wrf_path_daily, sprintf(file_pattern, 'hourly', datestr(date_in, 'yyyy-mm-dd')));
monthly_file = fullfile(wrf_path_monthly, sprintf(file_pattern, 'monthly', datestr(eomdate(date_in), 'yyyy-mm-dd')));

no2_daily = double(ncread(daily_file, 'no2'))/1e6;
pres_daily = double(ncread(daily_file, 'pres'));
lon_daily = double(ncread(daily_file,'XLONG'));
lat_daily = double(ncread(daily_file,'XLAT'));

no2_daily = no2_daily(:,:,:,1);
pres_daily = pres_daily(:,:,:,1);
lon_daily = lon_daily(:,:,1);
lat_daily = lat_daily(:,:,1);

no2_monthly = double(ncread(monthly_file, 'no2'))/1e6;
pres_monthly = double(ncread(monthly_file, 'pres'));

[profiles.converged, profiles.daily, profiles.monthly, amfs.converged, amfs.daily, amfs.monthly] =... 
    test_apriori_convergence(no2_daily, pres_daily, lon_daily, lat_daily, no2_monthly, pres_monthly, wrf_path_monthly, sza, vza, raa, alb, surfp, cldradfrac, cldp, date_in);

amf_settings.sza = sza;
amf_settings.vza = vza;
amf_settings.raa = raa;
amf_settings.alb = alb;
amf_settings.surfPres = surfp;
amf_settings.cldradfrac = cldradfrac;
amf_settings.cldPres = cldp;
amf_settings.lon = lon;
amf_settings.lat = lat;
amf_settings.date = datestr(date_in,'yyyy-mm-dd');

% Sample plots
% Find the location with the biggest increase and decrease in the bottom
% bin
del = no2_daily - no2_monthly;
[~,xx_pos] = max(reshape(del(:,:,1),[],1));
[~,xx_neg] = min(reshape(del(:,:,1),[],1));
[~,xx_sm] = min(reshape(abs(del(:,:,1)),[],1));
[~,xx_postop] = max(reshape(del(:,:,end),[],1));
[~,xx_negtop] = min(reshape(del(:,:,end),[],1));
[~,xx_bigdelA] = max(abs(amfs.converged(:) - amfs.daily(:)));

convergence_prof_plot(profiles.daily, profiles.monthly, profiles.converged, amfs.daily, amfs.monthly, amfs.converged, xx_pos, 'Largest positive surface difference');
convergence_prof_plot(profiles.daily, profiles.monthly, profiles.converged, amfs.daily, amfs.monthly, amfs.converged, xx_neg, 'Largest negative surface difference');
convergence_prof_plot(profiles.daily, profiles.monthly, profiles.converged, amfs.daily, amfs.monthly, amfs.converged, xx_sm, 'Smallest absolute surface difference');
convergence_prof_plot(profiles.daily, profiles.monthly, profiles.converged, amfs.daily, amfs.monthly, amfs.converged, xx_postop, 'Largest positive top difference');
convergence_prof_plot(profiles.daily, profiles.monthly, profiles.converged, amfs.daily, amfs.monthly, amfs.converged, xx_negtop, 'Largest negative top difference');
convergence_prof_plot(profiles.daily, profiles.monthly, profiles.converged, amfs.daily, amfs.monthly, amfs.converged, xx_bigdelA, 'Largest absolute amf difference');

end


