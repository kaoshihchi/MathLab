close all; clear;
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'theory'));
%% -------- Constants (SI) --------
global q_e m_e hbar
q_e  = 1.602176634e-19;      % C
m_e  = 9.1093837015e-31;     % kg
hbar = 1.054571817e-34;      % J*s
mu_0 = 4*pi*1e-7;            % H/m
c    = 2.99792458e8;         % m/s

%% -------- User parameters --------
pixel_um = 0.88;                       % µm
pixel_m  = pixel_um * 1e-6;
tau_fs   = 40;                         % fs
tau_s    = tau_fs * 1e-15;

E1_mJ = 8.3;      E1_J = E1_mJ * 1e-3; % input.tiff
E2_mJ = 7.3;    E2_J = E2_mJ * 1e-3; % output.tiff

q         = 45;                        % harmonic order
lambda_nm = 808;                       % nm
omega_d   = 2*pi*c/(lambda_nm*1e-9);
I_p_eV    = 15.7596;                   % eV (Ar)
I_p       = I_p_eV * q_e;              % J

%% -------- Load images once --------
name_bg = 'background.tiff'; 
name_in1 = 'input.tiff'; 
name_in2 = 'output_0.64psi.tiff'; 

bg  = im2double(imread(name_bg));
in1 = im2double(imread(name_in1));
in2 = im2double(imread(name_in2));

% Background subtraction (clip at 0)
corr1 = max(in1 - bg, 0);
corr2 = max(in2 - bg, 0);

%% -------- Peak intensity maps (W/cm^2) --------
% If you already have these helpers, keep them. Otherwise:
% [~, I_Wcm2] = ComputePeakIntensity(corr, pixel_m, E_J, tau_s);
[~, I1_Wcm2] = ComputePeakIntensity(corr1, pixel_m, E1_J, tau_s);
[~, I2_Wcm2] = ComputePeakIntensity(corr2, pixel_m, E2_J, tau_s);

% Display intensity maps
figure; tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(I1_Wcm2); axis image; colormap hot; cb=colorbar;
cb.Label.String='Peak Intensity (W/cm^2)';
title(sprintf('input.tiff  (%.3f mJ)', E1_mJ)); xlabel('X pixel'); ylabel('Y pixel');

nexttile; imagesc(I2_Wcm2); axis image; colormap hot; cb=colorbar;
cb.Label.String='Peak Intensity (W/cm^2)';
title(sprintf('output.tiff (%.3f mJ)', E2_mJ)); xlabel('X pixel'); ylabel('Y pixel');

Ipk1 = max(I1_Wcm2(:)); Ipk2 = max(I2_Wcm2(:));
fprintf('Peak intensity (input.tiff,  %.3f mJ): %.3g W/cm^2\n', E1_mJ, Ipk1);
fprintf('Peak intensity (output.tiff, %.3f mJ): %.3g W/cm^2\n', E2_mJ, Ipk2);

%% -------- Choose centroid mode --------
% Options: 'weighted' or 'geometric'
centroid_mode = 'weighted';   % <- change to 'weighted' when desired

% -------- Use it on your two intensity maps --------
[c1_xy, s1] = locateSpot(I1_Wcm2, centroid_mode);
[c2_xy, s2] = locateSpot(I2_Wcm2, centroid_mode);
% c2_xy = [615, 705]; 

x1_um = (c1_xy(1) - 1) * pixel_um;  y1_um = (c1_xy(2) - 1) * pixel_um;
x2_um = (c2_xy(1) - 1) * pixel_um;  y2_um = (c2_xy(2) - 1) * pixel_um;

% Re-plot with overlays
figure; tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile; imagesc(I1_Wcm2); axis image; colormap(ax1, hot); cb=colorbar(ax1);
cb.Label.String='Peak Intensity (W/cm^2)';
title(ax1, sprintf('input.tiff  (%.3f mJ) — %s centroid', E1_mJ, centroid_mode));
xlabel(ax1,'X pixel'); ylabel(ax1,'Y pixel'); set(ax1,'YDir','normal');
drawCrossAt(ax1, c1_xy, 12, 1.8, 'w'); drawEllipse(ax1, s1, 'c', 1.4);

ax2 = nexttile; imagesc(I2_Wcm2); axis image; colormap(ax2, hot); cb=colorbar(ax2);
cb.Label.String='Peak Intensity (W/cm^2)';
title(ax2, sprintf('output.tiff (%.3f mJ) — %s centroid', E2_mJ, centroid_mode));
xlabel(ax2,'X pixel'); ylabel(ax2,'Y pixel'); set(ax2,'YDir','normal');
drawCrossAt(ax2, c2_xy, 12, 1.8, 'w'); drawEllipse(ax2, s2, 'c', 1.4);

sgtitle(tl, 'Peak Intensity Maps with Spot Centers (+ optional ellipse)');

fprintf('\n--- Spot centers (%s) ---\n', centroid_mode);
fprintf('input.tiff  center: (x=%.2f px, y=%.2f px)  -> (%.2f µm, %.2f µm)\n', ...
    c1_xy(1), c1_xy(2), x1_um, y1_um);
fprintf('output.tiff center: (x=%.2f px, y=%.2f px)  -> (%.2f µm, %.2f µm)\n\n', ...
    c2_xy(1), c2_xy(2), x2_um, y2_um);

%% -------- Compute alpha(x,y) = -(1/dz)*ln(Iout/Iin) --------
dz_m = 8e-3;       % distance between planes (8 mm)

% Pick input as reference and output as the one to move
ref_map    = I1_Wcm2;  ref_center  = c1_xy;   ref_title  = name_in1;
move_map   = I2_Wcm2;  move_center = c2_xy;   move_title = name_in2;

% Align output onto input by centroid shift
shift_xy   = ref_center - move_center;     % [dx, dy] in pixels
I2_aligned = imtranslate(move_map, shift_xy, ...
                 'linear','OutputView','same','FillValues',0);

% Build mask where both maps have signal
frac = 0.02;   % threshold as fraction of each max
mask = (ref_map > frac*max(ref_map(:))) & (I2_aligned > frac*max(I2_aligned(:)));

% Compute attenuation factor and alpha
eps0 = 1e-12;
Att_lin   = I2_aligned ./ max(ref_map, eps0);
alpha_map = -(1/dz_m) * log(max(Att_lin, eps0));
alpha_map(~mask) = NaN;

% -------- Plot alpha(x,y) --------
figure;
imagesc(alpha_map); axis image; set(gca,'YDir','normal');
colormap(jet); cb = colorbar; cb.Label.String = '\alpha (1/m)';
title(sprintf('\\alpha(x,y), %s → %s (dz=%.3f m)', ref_title, move_title, dz_m));
hold on; 
% drawCrossAt(gca, ref_center, 12, 1.8, 'w'); 
hold off;
xlabel('X pixel'); ylabel('Y pixel');

% -------- Radial profile of alpha --------
nbins = 200;
[rad_um, alpha_rad] = radial_profile_scalar(alpha_map, ref_center, pixel_um, nbins, mask);

figure;
plot(rad_um, alpha_rad, 'LineWidth', 1.8); grid on;
xlabel('Radius (µm)'); ylabel('\alpha (1/m)');
title('Radial profile of \alpha (annular mean)');


%%
makeInteractive_I_alpha(ref_map, alpha_map, ref_center, pixel_um);


%% ===== Interactive I(x,y) with data cursor showing I and alpha =====
function makeInteractive_I_alpha(ref_map, alpha_map, ref_center, pixel_um)
    % ref_map: input intensity (W/cm^2), aligned
    % alpha_map: attenuation alpha (1/m), NaN where invalid
    % ref_center: [x,y] centroid (pixels)
    % pixel_um: pixel size (µm)

    f = figure('Name','Interactive I(x,y) & \alpha(x,y)','NumberTitle','off');
    tl = tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');

    % --- Left: I_in
    ax1 = nexttile(tl);
    hIm1 = imagesc(ax1, ref_map);  % capture IMAGE HANDLE
    axis(ax1,'image'); set(ax1,'YDir','normal');
    colormap(ax1,'hot'); cb1 = colorbar(ax1); cb1.Label.String = 'I_{in} (W/cm^2)';
    title(ax1,'Input peak intensity I_{in}(x,y)');
    xlabel(ax1,'X pixel'); ylabel(ax1,'Y pixel');
    % drawCrossAt(ax1, ref_center, 12, 1.6, 'w');

    % --- Right: alpha
    ax2 = nexttile(tl);
    hIm2 = imagesc(ax2, alpha_map); % capture IMAGE HANDLE
    axis(ax2,'image'); set(ax2,'YDir','normal');
    if exist('turbo','file'), colormap(ax2,turbo); else, colormap(ax2,jet); end
    cb2 = colorbar(ax2); cb2.Label.String = '\alpha (1/m)';
    title(ax2,'\alpha(x,y)');
    xlabel(ax2,'X pixel'); ylabel(ax2,'Y pixel');
    % drawCrossAt(ax2, ref_center, 12, 1.6, 'w');

    % --- Pixel readout widgets (need IMAGE handles, not axes)
    try
        impixelinfo(hIm1);
        impixelinfo(hIm2);
    catch
        % Older MATLAB: try the figure/handle form
        try
            impixelinfo(f, hIm1);
            impixelinfo(f, hIm2);
        catch
            % If unavailable, silently continue
        end
    end

    % Stash data for datacursor callback
    setappdata(f,'Iin_map',ref_map);
    setappdata(f,'alpha_map',alpha_map);
    setappdata(f,'center_xy',ref_center);
    setappdata(f,'pixel_um',pixel_um);

    % Custom data cursor for BOTH panels
    dcm = datacursormode(f);
    set(dcm,'Enable','on','SnapToDataVertex','on',...
        'UpdateFcn',@(obj,evt) localDatatip(evt,f));

    % Optional: left-click to place permanent datatip
    set(f,'WindowButtonDownFcn',@(~,~) localClickToDatatip(f));
end

%% --- Data tip text builder ---
function txt = localDatatip(evt, fig)
    Iin_map   = getappdata(fig,'Iin_map');
    alpha_map = getappdata(fig,'alpha_map');
    center_xy = getappdata(fig,'center_xy');
    pixel_um  = getappdata(fig,'pixel_um');

    pos = evt.Position;  % [x y] in pixels (image coords)
    j = round(pos(1));   % column (x)
    i = round(pos(2));   % row    (y)

    [H,W] = size(Iin_map);
    if i<1 || i>H || j<1 || j>W
        txt = {'(outside image)'}; return;
    end

    Iin   = Iin_map(i,j);
    alpha = alpha_map(i,j);
    r_um  = hypot(j - center_xy(1), i - center_xy(2)) * pixel_um;
    x_um  = (j - 1) * pixel_um;
    y_um  = (i - 1) * pixel_um;

    txt = {
        sprintf('Pixel: (x=%d, y=%d)', j, i)
        sprintf('Coord: (%.2f µm, %.2f µm)', x_um, y_um)
        sprintf('r from center: %.2f µm', r_um)
        sprintf('I_{in}: %.4g W/cm^2', Iin)
        sprintf('\\alpha: %.4g 1/m', alpha)
    };
end

%% --- Click anywhere to drop a permanent datatip (optional) ---
function localClickToDatatip(fig)
    ax = gca;
    cp = get(ax,'CurrentPoint');
    x = cp(1,1); y = cp(1,2);
    xl = xlim(ax); yl = ylim(ax);
    if x>=xl(1) && x<=xl(2) && y>=yl(1) && y<=yl(2)
        try, datatip(ax,x,y); end
    end
end



function [r_um, prof] = radial_profile_scalar(img, center_xy, pixel_um, nbins, mask)
    if nargin < 5 || isempty(mask), mask = ~isnan(img); end
    [H,W] = size(img);
    [X,Y] = meshgrid(1:W, 1:H);
    R_um = hypot(X - center_xy(1), Y - center_xy(2)) * pixel_um;

    Rmax_um = max(R_um(mask));
    edges = linspace(0, Rmax_um, nbins+1);
    r_um = 0.5*(edges(1:end-1)+edges(2:end));

    prof = nan(1, nbins);
    for k = 1:nbins
        in_bin = mask & R_um >= edges(k) & R_um < edges(k+1) & isfinite(img);
        if any(in_bin(:))
            prof(k) = mean(img(in_bin),'omitnan');
        end
    end
end

%% -------- Locate focal spot (single particle) with mode switch --------
function [center_xy, stats] = locateSpot(I, centroid_mode)
    % I: intensity map (double, nonnegative)
    I = double(I);
    if ~isfinite(max(I(:))) || max(I(:)) <= 0
        error('locateSpot: invalid intensity map');
    end

    % 1) Light smoothing to suppress noise; normalize
    I_s = imgaussfilt(I, 1.0);          % sigma ~ 1 pixel
    In  = I_s / max(I_s(:));            % [0,1]

    % 2) Threshold -> keep largest blob
    T  = graythresh(In);                 % Otsu
    BW = In > max(0.35, 0.8*T);          % robust for various contrasts
    BW = bwareafilt(BW, 1);

    % Fallback if thresholding found nothing
    if ~any(BW(:))
        warning('Thresholding found no blob; falling back to top 0.5%% pixels.');
        p  = prctile(In(:), 99.5);
        BW = In >= p;
        BW = bwareafilt(BW, 1);
    end

    % 3) Compute center based on mode
    switch lower(string(centroid_mode))
        case "weighted"
            % intensity-weighted centroid (sub-pixel)
            S = regionprops(BW, In, ...
                'WeightedCentroid','Centroid','Area','BoundingBox', ...
                'Orientation','MajorAxisLength','MinorAxisLength','PixelIdxList');
            if isempty(S)
                [yy, xx] = find(In == max(In(:)), 1, 'first');
                center_xy = [xx, yy];
                stats = struct('WeightedCentroid', center_xy, 'Centroid', center_xy, ...
                    'Area', 1, 'BoundingBox', [xx yy 1 1], ...
                    'Orientation', 0, 'MajorAxisLength', 0, 'MinorAxisLength', 0);
            else
                center_xy = S(1).WeightedCentroid;   % [x,y]
                stats     = S(1);
            end

        case "geometric"
            % pure geometry centroid of the binary mask
            S = regionprops(BW, ...
                'Centroid','Area','BoundingBox','Orientation', ...
                'MajorAxisLength','MinorAxisLength','PixelIdxList');
            if isempty(S)
                [yy, xx] = find(In == max(In(:)), 1, 'first');
                center_xy = [xx, yy];
                stats = struct('Centroid', center_xy, ...
                    'Area', 1, 'BoundingBox', [xx yy 1 1], ...
                    'Orientation', 0, 'MajorAxisLength', 0, 'MinorAxisLength', 0);
            else
                center_xy = S(1).Centroid;           % [x,y]
                stats     = S(1);
            end

        otherwise
            error('centroid_mode must be ''weighted'' or ''geometric''.');
    end

    % Also return the binary mask via stats for optional visualization
    stats.Mask = BW;
end
%% -------- Simple cross & (optional) ellipse overlay helpers --------
function drawCrossAt(ax, center_xy, L, lw, color)
    if nargin < 4, lw = 1.6; end
    if nargin < 5, color = 'w'; end
    x = center_xy(1);  y = center_xy(2);
    hold(ax, 'on');
    line(ax, [x-L, x+L], [y, y], 'Color', color, 'LineWidth', lw);
    line(ax, [x, x], [y-L, y+L], 'Color', color, 'LineWidth', lw);
    plot(ax, x, y, '+', 'Color', color, 'MarkerSize', 10, 'LineWidth', lw);
end

function drawEllipse(ax, stats, color, lw)
    % Visualize fitted ellipse from regionprops (major/minor + orientation)
    if nargin < 3, color = 'c'; end
    if nargin < 4, lw = 1.2; end
    if ~isfield(stats,'MajorAxisLength') || stats.MajorAxisLength==0
        return;
    end
    a = stats.MajorAxisLength/2;
    b = stats.MinorAxisLength/2;
    th = deg2rad(-stats.Orientation); % regionprops angle: CCW from x-axis, display needs minus
    t = linspace(0, 2*pi, 200);
    R = [cos(th) -sin(th); sin(th) cos(th)];
    xy = R * [a*cos(t); b*sin(t)];
    x  = xy(1,:) + stats.Centroid(1);
    y  = xy(2,:) + stats.Centroid(2);
    hold(ax,'on'); plot(ax, x, y, '-', 'Color', color, 'LineWidth', lw);
end
