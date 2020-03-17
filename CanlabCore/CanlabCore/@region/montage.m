function o2 = montage(obj, varargin)
% This function displays a region object on a standard slice montage
%
% :Usage:
% ::
%
%    montage(obj, [optional inputs])
%
% - takes all optional inputs to canlab_results_fmridisplay and addblobs method
%
% :Input:
%
%   **obj:**
%        a region object
%
% :Optional Inputs:  - see help canlab_results_fmridisplay
%
%   **o2***
%        An existing fmridisplay object, with no keyword strings
%
%   **'solid' or 'color'**
%        Followed by a single [r g b] triplet to define color for all regions. Passed
%        through to addblobs.
%
%   **'colors'**
%        Followed by a cell array of [r g b] triplets, one color for each
%        region. Note: If 'symmetric' option is on, will reorganize color
%        to region assignments.
%
%   **symmetric** [default]
%       Mirror left/right blobs with same colors
%       See match_colors_left_right
%
%   **nosymmetric**
%       Standard color map, no L/R color-matching for symmetry 
%
%   **'noblobs':**
%        do not display blobs
%
%   **'nooutline':**
%        do not display blob outlines
%
%   **'addmontages':**
%        when entering existing fmridisplay obj, add new montages
%
%   **'noremove':**
%        do not remove current blobs when adding new ones
%
%   **'outlinecolor:**
%        followed by new outline color
%
%   **'splitcolor':**
%        followed by 4-cell new split colormap colors (help fmridisplay or edit code for defaults as example)
%
%   **'montagetype':**
%        'full' for full montages of axial and sagg slices.
%
%        'full hcp' for full montage, but with surfaces and volumes from
%        HCP data
%
%        'compact' [default] for single-figure parasagittal and axials slices.
%
%        'compact2': like 'compact', but fewer axial slices.
%
%        'multirow': followed by number of rows
%           e.g., o2 = canlab_results_fmridisplay([], 'multirow', 2);
%
%        'regioncenters': Show one slice per region, centered on region center
%               Note: 'nozoom' option omits zooming-in on regions
%
%   **'noverbose':**
%        suppress verbose output, good for scripts/publish to html, etc.
%
%   **'overlay':**
%        specify anatomical image for montage (not surfaces), followed by
%        image name
%        e.g., o2 = canlab_results_fmridisplay([], 'overlay', 'icbm152_2009_symmetric_for_underlay.img')';
%
%         The default brain for overlays is based on Keuken et al. 2014
%         For legacy SPM8 single subject, enter as arguments:
%         'overlay', which('SPM8_colin27T1_seg.img')
%
% Other inputs to addblobs (fmridisplay method) are allowed, e.g., 'cmaprange', [-2 2], 'trans'
%
% See help fmridisplay
% e.g., 'color', [1 0 0]
%
%
% Examples:
%
% o2 = canlab_results_fmridisplay([], 'noverbose');
% o2 = montage(r, o2);      % symmetric colors left/right
% o2 = removeblobs(o2);
% o2 = montage(r, o2, 'map');
% 
% create_figure('slices'); axis off
% o2 = canlab_results_fmridisplay([], 'multirow', 2);
% brighten(.6)
% hcp152t1 = which('HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii');
% r = region(fmri_data(hcp152t1), 'unique_mask_values');
% o2 = montage(r, o2, 'wh_montages', 3:4);
% o2 = montage(r(1:20), o2, 'wh_montages', 1:2, 'color', [1 .5 0]);
%
% Use a different montage type in canlab_results_fmridisplay:
% o2 = montage(r, 'compact2', 'nosymmetric');
%
% Plot one region blob per slice on a series of montages.  
% o2 = montage(r, 'regioncenters', 'nosymmetric');

% edited: Tor, 1/2018, unique color option/default

% Defaults and inputs
% -----------------------------------------------------------------------

methodtype = 'symmetric';
colortype = 'unique';
one_blob_per_slice = false; 
dozoom = true;
colors = scn_standard_colors(length(obj));

if any(strcmp(varargin, 'map')), methodtype = 'map'; end
if any(strcmp(varargin, 'nosymmetric')), methodtype = 'map'; end
if any(strcmp(varargin, 'old')), methodtype = 'old'; end
if any(strcmp(varargin, 'solid')), colortype = 'solid'; end
if any(strcmp(varargin, 'color')), colortype = 'solid'; end
if any(strcmp(varargin, 'colors'))
    wh = find(strcmp(varargin, 'colors'));
    colors = varargin{wh + 1}; varargin{wh} = []; varargin{wh + 1} = [];  
end
if any(strcmp(varargin, 'regioncenters')), one_blob_per_slice = true; end
if any(strcmp(varargin, 'nozoom')), dozoom = false; end

% Initialize display if needed
% -----------------------------------------------------------------------

funhan = @(x) isa(x, 'fmridisplay');
whfmridisplay = cellfun(funhan, varargin);

o2 = [];

if any(whfmridisplay)
    wh = find(whfmridisplay);
    o2 = varargin{wh(1)};
end

if ~exist('o2', 'var') || ~isa(o2, 'fmridisplay')
    create_figure('fmridisplay'); axis off

        % some inputs control which slices are shown: standard slices, or specific locations for these blobs
        % regioncenters is passed through, so it determines the type of montage to create: 
        o2 = canlab_results_fmridisplay(obj, 'noverbose', varargin{:}, 'noblobs', 'nooutline');

end


% Define colors for 'unique' option
% And render blobs
% -------------------------------------------------------------------------
switch colortype
    
    case 'unique' % different color for each blob (some reuse)
        
        switch methodtype
            case 'symmetric'
                colors = match_colors_left_right(obj);
        end
        
        if one_blob_per_slice
            % One blob on each slice
            
            for i = 1:length(obj)
                
                o2 = addblobs(o2, obj(i), 'color', colors{i}, 'noverbose', varargin{:}, 'wh_montages', i);
                drawnow
                
            end
            
            if dozoom
                zoom_in_on_regions(o2, obj, 'axial');  % hard-coded for now, could change orientation...must make flexible in canlab_results_fmridisplay
            end
            
        else
            % All visible blobs on each slice
            
            % Redefine blob groups based on color categories, to reduce number of regions, combining those with same color too many function calls are bad news!! overloads handle graphics system?
            [region_groups, colors] = redefine_colors_based_on_groups(obj, colors);
            
            for i = 1:length(region_groups)
                
                o2 = addblobs(o2, region_groups{i}, 'color', colors{i}, 'noverbose', varargin{:});
                drawnow
                
            end
            
        end
        
    case 'solid' % solid color, user-entered
        
        if one_blob_per_slice
            
            for i = 1:length(obj)
                
                o2 = addblobs(o2, obj(i), varargin{:}, 'wh_montages', i);
                drawnow
                
            end
            
        else
            
            % Just render blobs
            o2 = addblobs(o2, obj, varargin{:});
        end
        
    case 'old'
        
        montage_clusters([], obj, varargin{:})
        
end


end % function


function [region_groups, color_groups] = redefine_colors_based_on_groups(r, colors)

% regroup to reduce number of regions, combining those with same color
% too many function calls are bad news!! overloads handle graphics system?

catcolors = cat(1, colors{:});
uniquecolors = unique(catcolors, 'rows');
ncolors = size(uniquecolors, 1);

region_groups = {};
color_groups = {};

for i = 1:ncolors
    
    mycolor = uniquecolors(i, :);
    cmatchfun = @(x) ~any(x - mycolor);
    
    wh_regions = cellfun(cmatchfun, colors);
    
    region_groups{i} = r(wh_regions);
    
    color_groups{i} = mycolor;
    
    % o2 = addblobs(o2, r(wh_regions), 'maxcolor', mycolor, 'mincolor', mycolor, 'wh_montages', 1:2);
    
end

end
