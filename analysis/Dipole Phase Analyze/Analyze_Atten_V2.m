% ======================= Analyze_Atten_Rewrite.m =======================
% Non‑OOP rewrite of your script. Keep your ComputePeakIntensity.m on path.
% What this script does:
%   1) Load input/output TIFFs, subtract background
%   2) Compute peak intensity maps I_in, I_out (W/cm^2)
%   3) Find centroids (weighted / geometric)
%   4) Align output to input; compute alpha(x,y) = -(1/dz) ln(Iout/Iin)
%   5) Figures:
%       - I maps with centers (pixels)
%       - alpha(x,y) (pixels) + radial profile
%       - SEPARATE calibrated figures in micrometers:
%           * I_in(x,y)
%           * I_out(x,y)
%           * alpha(x,y)
%       - Optional interactive view (I_in + alpha)
%
% USAGE: run this file as a script.

close all; clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'theory'));

%% -------- Constants (SI) --------
% (Only those actually used below; others removed for clarity)
c = 2.99792458e8;            % m/s

%% -------- User parameters --------
pixel_um = 0.88;                       % µm
pixel_m  = pixel_um * 1e-6;
tau_fs   = 40;                          % fs
tau_s    = tau_fs * 1e-15;

E1_mJ = 4*0.8;    E1_J = E1_mJ * 1e-3;   % input.tiff
E2_mJ = 4*0.8*0.8;    E2_J = E2_mJ * 1e-3;   % output.tiff

q         = 101;                        % (kept for context)
lambda_nm = 808;                       % nm
omega_d   = 2*pi*c/(lambda_nm*1e-9);   %#ok<NASGU>

% Files
name_bg  = 'background_He.tiff';
name_in1 = 'input_He.tiff';
name_in2 = 'output_He.tiff';

% Centroid mode and (optional) manual override for output center
centroid_mode = 'geometric';           % 'weighted' or 'geometric'
override_c2   = [];            % [] to disable, or [x,y]

% Separation along z for attenuation
dz_m = 8e-3;                           % 8 mm

%% -------- Load images once --------
bg  = im2double(imread(name_bg));
in1 = im2double(imread(name_in1));
in2 = im2double(imread(name_in2));

% Background subtraction (clip at 0)
corr1 = max(in1 - bg, 0);
corr2 = max(in2 - bg, 0);

%% -------- Peak intensity maps (W/cm^2) --------
[~, I1_Wcm2] = ComputePeakIntensity(corr1, pixel_m, E1_J, tau_s);
[~, I2_Wcm2] = ComputePeakIntensity(corr2, pixel_m, E2_J, tau_s);

Ipk1 = max(I1_Wcm2(:)); Ipk2 = max(I2_Wcm2(:));
fprintf('Peak intensity (input.tiff,  %.3f mJ): %.3g W/cm^2\n', E1_mJ, Ipk1);
fprintf('Peak intensity (output.tiff, %.3f mJ): %.3g W/cm^2\n', E2_mJ, Ipk2);

%% -------- Find centroids --------
[c1_xy, s1] = locateSpot(I1_Wcm2, centroid_mode);
[c2_xy, s2] = locateSpot(I2_Wcm2, centroid_mode);
if ~isempty(override_c2) && numel(override_c2)==2 && all(isfinite(override_c2))
    c2_xy = override_c2;
end

% Quick visualization with overlays (pixels)
figure('Name','I maps with centers (pixels)','NumberTitle','off');
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
ax1 = nexttile; imagesc(ax1, I1_Wcm2); axis(ax1,'image'); set(ax1,'YDir','normal');
colormap(ax1, hot); cb=colorbar(ax1); cb.Label.String='Peak Intensity (W/cm^2)';
title(ax1, sprintf('%s  (%.3f mJ) — %s centroid', name_in1, E1_mJ, centroid_mode));
xlabel(ax1,'X pixel'); ylabel(ax1,'Y pixel'); drawCrossAt(ax1, c1_xy, 12, 1.8, 'w'); drawEllipse(ax1, s1, 'c', 1.4);
ax2 = nexttile; imagesc(ax2, I2_Wcm2); axis(ax2,'image'); set(ax2,'YDir','normal');
colormap(ax2, hot); cb=colorbar(ax2); cb.Label.String='Peak Intensity (W/cm^2)';
title(ax2, sprintf('%s (%.3f mJ) — %s centroid', name_in2, E2_mJ, centroid_mode));
xlabel(ax2,'X pixel'); ylabel(ax2,'Y pixel'); drawCrossAt(ax2, c2_xy, 12, 1.8, 'w'); drawEllipse(ax2, s2, 'c', 1.4);
sgtitle(tl, 'Peak Intensity Maps with Spot Centers');

%% -------- Align output to input; compute alpha(x,y) --------
shift_xy   = c1_xy - c2_xy;                       % [dx, dy] in pixels
I2_aligned = imtranslate(I2_Wcm2, shift_xy, 'linear','OutputView','same','FillValues',0);

frac = 0.02;                                      % overlap mask
mask = (I1_Wcm2 > frac*max(I1_Wcm2(:))) & (I2_aligned > frac*max(I2_aligned(:)));

eps0 = 1e-12;                                     % guard
Att_lin   = I2_aligned ./ max(I1_Wcm2, eps0);
alpha_map = -(1/dz_m) * log(max(Att_lin, eps0));
alpha_map(~mask) = NaN;

% --- alpha(x,y) in pixels
figure('Name','alpha(x,y) (pixels)','NumberTitle','off');
ax = axes; imagesc(ax, alpha_map); axis(ax,'image'); set(ax,'YDir','normal');
try, if exist('turbo','file'), colormap(ax,turbo); else, colormap(ax,jet); end; catch, colormap(ax,jet); end
cb = colorbar(ax); cb.Label.String = '\alpha (1/m)';
title(ax, sprintf('\\alpha(x,y), %s → %s (dz=%.3f m)', name_in1, name_in2, dz_m));
xlabel(ax,'X pixel'); ylabel(ax,'Y pixel');

% --- Radial profile of alpha
nbins = 200; [rad_um, alpha_rad] = radial_profile_scalar(alpha_map, c1_xy, pixel_um, nbins, mask);
figure('Name','alpha radial profile','NumberTitle','off');
plot(rad_um, alpha_rad, 'LineWidth', 1.8); grid on;
xlabel('Radius (\mum)'); ylabel('\alpha (1/m)'); title('Radial profile of \alpha (annular mean)');

%% -------- Calibrated images in micrometers (separate figures) --------
[H,W] = size(I1_Wcm2); x_um = (0:W-1) * pixel_um; y_um = (0:H-1) * pixel_um;

% Shared intensity limits for fair comparison
vmax = max([I1_Wcm2(:); I2_aligned(:)],[],'omitnan'); if ~isfinite(vmax)||vmax<=0, vmax = 1; end
clim_I = [0, vmax];

% Robust alpha limits
alpha_vec = alpha_map(isfinite(alpha_map));
if ~isempty(alpha_vec)
    lo = prctile(alpha_vec,1); hi = prctile(alpha_vec,99);
    if ~(isfinite(lo)&&isfinite(hi)&&lo<hi), lo=min(alpha_vec); hi=max(alpha_vec); end
    clim_A = [lo, hi];
else
    clim_A = [0,1];
end

% I_in(x,y)
figure('Name','I_{in}(x,y) in \mum','NumberTitle','off');
ax = axes; imagesc(ax, x_um, y_um, I1_Wcm2); axis(ax,'image'); set(ax,'YDir','normal');
colormap(ax,'hot'); caxis(ax,clim_I); cb=colorbar(ax); cb.Label.String='Peak Intensity (W/cm^2)';
xlabel(ax,'x (\mum)'); ylabel(ax,'y (\mum)'); title(ax,'Input peak intensity I_{in}(x,y)');
try, impixelinfo; end

% I_out(x,y)
figure('Name','I_{out}(x,y) in \mum','NumberTitle','off');
ax = axes; imagesc(ax, x_um, y_um, I2_aligned); axis(ax,'image'); set(ax,'YDir','normal');
colormap(ax,'hot'); caxis(ax,clim_I); cb=colorbar(ax); cb.Label.String='Peak Intensity (W/cm^2)';
xlabel(ax,'x (\mum)'); ylabel(ax,'y (\mum)'); title(ax,'Output peak intensity I_{out}(x,y)');
try, impixelinfo; end

% alpha(x,y)
figure('Name','alpha(x,y) in \mum','NumberTitle','off');
ax = axes; imagesc(ax, x_um, y_um, alpha_map); axis(ax,'image'); set(ax,'YDir','normal');
try, if exist('turbo','file'), colormap(ax,turbo); else, colormap(ax,jet); end; catch, colormap(ax,jet); end
caxis(ax,clim_A); cb=colorbar(ax); cb.Label.String='\alpha (1/m)';
xlabel(ax,'x (\mum)'); ylabel(ax,'y (\mum)'); title(ax,'\alpha(x,y)');
try, impixelinfo; end

%% -------- Optional interactive figure (I_in + alpha) --------
%  (Click to drop permanent datatips; hover for values.)
makeInteractive_I_alpha(I1_Wcm2, alpha_map, c1_xy, pixel_um);


% ========================== Local functions ==========================
function makeInteractive_I_alpha(ref_map, alpha_map, ref_center, pixel_um)
    f = figure('Name','Interactive I(x,y) & \alpha(x,y)','NumberTitle','off');
    tl = tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');

    % --- Left: I_in
    ax1 = nexttile(tl); hIm1 = imagesc(ax1, ref_map);
    axis(ax1,'image'); set(ax1,'YDir','normal'); colormap(ax1,'hot');
    cb1 = colorbar(ax1); cb1.Label.String = 'I_{in} (W/cm^2)';
    title(ax1,'Input peak intensity I_{in}(x,y)'); xlabel(ax1,'X pixel'); ylabel(ax1,'Y pixel');

    % --- Right: alpha
    ax2 = nexttile(tl); hIm2 = imagesc(ax2, alpha_map);
    axis(ax2,'image'); set(ax2,'YDir','normal');
    try, if exist('turbo','file'), colormap(ax2,turbo); else, colormap(ax2,jet); end; catch, colormap(ax2,jet); end
    cb2 = colorbar(ax2); cb2.Label.String = '\alpha (1/m)';
    title(ax2,'\alpha(x,y)'); xlabel(ax2,'X pixel'); ylabel(ax2,'Y pixel');

    % Pixel readouts
    try, impixelinfo(hIm1); impixelinfo(hIm2); end

    % Stash data for datacursor callback
    setappdata(f,'Iin_map',ref_map);
    setappdata(f,'alpha_map',alpha_map);
    setappdata(f,'center_xy',ref_center);
    setappdata(f,'pixel_um',pixel_um);

    % IMPORTANT: set click handler BEFORE enabling datacursormode
    set(f,'WindowButtonDownFcn',@(~,~) localClickToDatatip());

    % Data cursor
    dcm = datacursormode(f);
    set(dcm,'Enable','on','SnapToDataVertex','on','UpdateFcn',@(obj,evt) localDatatip(evt,f));
end

function txt = localDatatip(evt, fig)
    Iin_map   = getappdata(fig,'Iin_map');
    alpha_map = getappdata(fig,'alpha_map');
    center_xy = getappdata(fig,'center_xy');
    pixel_um  = getappdata(fig,'pixel_um');
    pos = evt.Position; j = round(pos(1)); i = round(pos(2));
    [H,W] = size(Iin_map);
    if i<1 || i>H || j<1 || j>W, txt = {'(outside image)'}; return; end
    Iin   = Iin_map(i,j); alpha = alpha_map(i,j);
    r_um  = hypot(j-center_xy(1), i-center_xy(2)) * pixel_um;
    x_um  = (j-1) * pixel_um; y_um = (i-1) * pixel_um;
    txt = {sprintf('Pixel: (x=%d, y=%d)', j, i), ...
           sprintf('Coord: (%.2f \xB5m, %.2f \xB5m)', x_um, y_um), ...
           sprintf('r from center: %.2f \xB5m', r_um), ...
           sprintf('I_{in}: %.4g W/cm^2', Iin), ...
           sprintf('\\alpha: %.4g 1/m', alpha)};
end

function localClickToDatatip()
    ax = gca; cp = get(ax,'CurrentPoint'); x = cp(1,1); y = cp(1,2);
    xl = xlim(ax); yl = ylim(ax);
    if x>=xl(1) && x<=xl(2) && y>=yl(1) && y<=yl(2)
        try, datatip(ax,x,y); end
    end
end

function [r_um, prof] = radial_profile_scalar(img, center_xy, pixel_um, nbins, mask)
    if nargin < 5 || isempty(mask), mask = ~isnan(img); end
    [H,W] = size(img); [X,Y] = meshgrid(1:W, 1:H);
    R_um = hypot(X - center_xy(1), Y - center_xy(2)) * pixel_um;
    Rmax_um = max(R_um(mask)); edges = linspace(0, Rmax_um, nbins+1); r_um = 0.5*(edges(1:end-1)+edges(2:end));
    prof = nan(1, nbins);
    for k = 1:nbins
        in_bin = mask & R_um >= edges(k) & R_um < edges(k+1) & isfinite(img);
        if any(in_bin(:)), prof(k) = mean(img(in_bin),'omitnan'); end
    end
end

function [center_xy, stats] = locateSpot(I, centroid_mode)
    I = double(I); I_s = imgaussfilt(I, 1.0); In = I_s / max(I_s(:));
    T  = graythresh(In); BW = In > max(0.35, 0.8*T); BW = bwareafilt(BW, 1);
    if ~any(BW(:)), p = prctile(In(:), 99.5); BW = In >= p; BW = bwareafilt(BW, 1); end
    switch lower(string(centroid_mode))
        case "weighted"
            S = regionprops(BW, In, 'WeightedCentroid','Centroid','Area','BoundingBox','Orientation','MajorAxisLength','MinorAxisLength','PixelIdxList');
            if isempty(S), [yy, xx] = find(In == max(In(:)), 1, 'first'); center_xy = [xx, yy];
            else, center_xy = S(1).WeightedCentroid; end; stats = S(1);
        case "geometric"
            S = regionprops(BW, 'Centroid','Area','BoundingBox','Orientation','MajorAxisLength','MinorAxisLength','PixelIdxList');
            if isempty(S), [yy, xx] = find(In == max(In(:)), 1, 'first'); center_xy = [xx, yy];
            else, center_xy = S(1).Centroid; end; stats = S(1);
    end
    stats.Mask = BW; %#ok<STRNU>
end

function drawCrossAt(ax, center_xy, L, lw, color)
    if nargin < 4, lw = 1.6; end
    if nargin < 5, color = 'w'; end
    x = center_xy(1); y = center_xy(2); hold(ax,'on');
    line(ax, [x-L, x+L], [y, y], 'Color', color, 'LineWidth', lw);
    line(ax, [x, x], [y-L, y+L], 'Color', color, 'LineWidth', lw);
    plot(ax, x, y, '+', 'Color', color, 'MarkerSize', 10, 'LineWidth', lw);
end

function drawEllipse(ax, stats, color, lw)
    if nargin < 3, color = 'c'; end
    if nargin < 4, lw = 1.2; end
    if ~isfield(stats,'MajorAxisLength') || stats.MajorAxisLength==0, return; end
    a = stats.MajorAxisLength/2; b = stats.MinorAxisLength/2; th = deg2rad(-stats.Orientation);
    t = linspace(0, 2*pi, 200); R = [cos(th) -sin(th); sin(th) cos(th)]; xy = R * [a*cos(t); b*sin(t)];
    x  = xy(1,:) + stats.Centroid(1); y = xy(2,:) + stats.Centroid(2);
    hold(ax,'on'); plot(ax, x, y, '-', 'Color', color, 'LineWidth', lw);
end
% ===================== end of Analyze_Atten_Rewrite.m =====================
