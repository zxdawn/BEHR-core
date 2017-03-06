function [ converged_profs ] = test_apriori_convergence( true_prof, true_pres, start_prof, start_pres, sza, vza, raa, alb, surfpres, cldradfrac, cldpres, date_in )
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
if ~isequal(size(start_prof), sz)
    E.badinput('START_PROF must be the same size as TRUE_PROF (a 3d array)')
end
if ~isequal(size(start_pres), sz)
    E.badinput('START_PRES must be the same size as TRUE_PROF (a 3d array)')
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


end



