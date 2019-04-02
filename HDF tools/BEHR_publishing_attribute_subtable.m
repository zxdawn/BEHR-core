function [ attr_table ] = BEHR_publishing_attribute_table( varargin )
%BEHR_PUBLISHING_ATTRIBUTE_TABLE Returns table of HDF attributes for BEHR fields
%   Detailed explanation goes here

E = JLLErrors;

as_struct = ismember('struct', varargin);

allowed_subsets = {'all', 'sp', 'behr', 'behr-insitu', 'pub', 'pub-insitu'};
xx = ismember(varargin, allowed_subsets);
if sum(xx) > 1
    E.badinput('Only one of %s may be an input', strjoin(allowed_subsets, ', '));
elseif sum(xx) < 1
    subset = 'all';
else
    subset = varargin{xx};
end 


% This cell array will have the variable name, unit, range, fill, product
% (SP or BEHR) and description in that order.
longfill = single(-1.267650600228229401496703205376e30);
shortfill = single(-32767);
pixcorfill = single(-1.00000001504747e+30);
behrfill = single(behr_fill_val());
behrflagfill = bitset(0, 32);
nofill = NaN;

% Variable name, unit, range, fill value, product, description

attr_table = {'PercentChangeNO2', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeNO2Vis', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeLNO2', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeLNO2_pickering', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';... 
    'PercentChangeLNO2Vis', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeLNOx', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeLNOx_pickering', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeLNOxVis', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeAMF', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeAMFVis', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeAMFLNO2', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeAMFLNO2_pickering', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeAMFLNO2Vis', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeAMFLNOx', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeAMFLNOx_pickering', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'PercentChangeAMFLNOxVis', 'unitless', [-Inf, Inf], longfill, 'BEHR', 'PercentChangeNO2';...
    'Swath', 'unitless', [0, Inf], nofill, 'SP', 'Orbit number since OMI launch';...
                };
            
%attr_table = add_psm_weight_fields(attr_table);
attr_table = choose_subset(attr_table, subset);
            
if numel(unique(attr_table(:,1))) < size(attr_table, 1)
    E.callError('attr_def', 'One or more attributes is multiply defined in the attributes table');
end
            
if as_struct
    fields = {'unit', 'range', 'fillvalue', 'product', 'description'};
    S = struct;
    for a = 1:size(attr_table,1)
        S.(attr_table{a,1}) = cell2struct(attr_table(a,2:end), fields, 2);
    end
    attr_table = S;
end

%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%
    function attr_table = add_psm_weight_fields(attr_table)
        weights_vars = BEHR_publishing_gridded_fields.psm_weight_vars;
        table_lines = cell(numel(weights_vars), size(attr_table, 2));
        for i=1:numel(weights_vars)
            if ismember(weights_vars{i}, BEHR_publishing_gridded_fields.psm_weight_vars('insitu'))
                product = 'BEHR-InSitu';
            else
                product = 'BEHR';
            end

            table_lines(i,:) = [weights_vars(i), {'unitless', [0, Inf], behrfill, product, sprintf('Weight field for the %s field', BEHR_publishing_gridded_fields.all_psm_vars{i})}];
        end
        attr_table = cat(1, attr_table, table_lines);
    end

    
    function attr_table = choose_subset(attr_table, subset)
        if strcmpi(subset, 'all')
            return
        elseif strcmpi(subset, 'pub')
            products = {'SP', 'BEHR', 'PIXCOR'};
        elseif strcmpi(subset, 'pub-insitu')
            products = {'SP', 'BEHR', 'PIXCOR', 'BEHR-InSitu'};
        else
            products = {subset};
        end

        xx = false(size(attr_table,1));

        for a=1:size(attr_table,1)
            xx(a) = any(strcmpi(attr_table{a,5}, products));
        end

        attr_table = attr_table(xx,:);
    end

end


