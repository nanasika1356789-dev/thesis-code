clear; clc; close all;

%% ===== User inputs =====
videoFile = 're-entrant 50% consolidated.MOV';    
%% ===== Read video =====
v = VideoReader(videoFile);
fps = v.FrameRate;
N   = max(1, floor(v.Duration*fps));

I1 = readFrame(v);
G1 = rgb2gray(I1);
G1 = adapthisteq(G1,'ClipLimit',0.02);    

figure(1); imshow(G1); title('Click LEFT-group points (top→bottom), then press Enter');
[xL,yL] = ginput; hold on; plot(xL,yL,'go','MarkerSize',8,'LineWidth',1.2);
title('Click RIGHT-group points (top→bottom), then press Enter');
[xR,yR] = ginput; plot(xR,yR,'ro','MarkerSize',8,'LineWidth',1.2); drawnow;

if numel(xL)~=numel(xR)
    error('The number of points should be same.');
end

L0 = sortrows([xL yL], 2);
R0 = sortrows([xR yR], 2);
pts0 = [L0; R0];
m    = size(L0,1);                

%% ===== (Optional) mm calibration =====
answer = questdlg('Do you want to calibrate pixels to mm?', ...
                  'Calibration', 'Yes','No','Yes');
if strcmp(answer,'Yes')
    figure(2); imshow(G1); title('Click two points with a KNOWN real distance (mm)');
    [xC,yC] = ginput(2);
    pxDist  = hypot(diff(xC), diff(yC));
    prompt  = {'Enter known distance (mm):'};
    dlgttl  = 'Calibration distance';
    defAns  = {'10'};
    resp    = inputdlg(prompt, dlgttl, [1 40], defAns);
    known_mm = str2double(resp{1});
    mm_per_px = known_mm / pxDist;
else
    mm_per_px = NaN; 
end

%% ===== KLT tracker =====
tracker = vision.PointTracker('NumPyramidLevels',4, ...
                              'BlockSize',[31 31], ...
                              'MaxBidirectionalError',3);
initialize(tracker, pts0, G1);

nPts  = size(pts0,1);
trajX = nan(N, nPts);
trajY = nan(N, nPts);
trajX(1,:) = pts0(:,1)';  trajY(1,:) = pts0(:,2)';

v.CurrentTime = 0;
k = 1;
while hasFrame(v) && k<=N
    if k==1
        I = I1; G = G1;
    else
        I = readFrame(v);
        G = rgb2gray(I);
        G = adapthisteq(G,'ClipLimit',0.02);
    end

    if k>1
        [pts, found] = step(tracker, G);
        prev = [trajX(k-1,:)' trajY(k-1,:)'];
        pts(~found,:) = prev(~found,:);  
        trajX(k,:) = pts(:,1)'; trajY(k,:) = pts(:,2)';
    end

    if mod(k,10)==0
        imshow(I); hold on; plot(trajX(k,1:end), trajY(k,1:end),'yx'); drawnow;
    end
    k = k+1;
end

trajX = trajX(1:k-1,:); trajY = trajY(1:k-1,:);
t = (0:size(trajX,1)-1)'/fps;  

%% ===== Multi-point averaging: lateral width =====
XL = trajX(:,1:m);  YL = trajY(:,1:m);        
XR = trajX(:,m+1:end); YR = trajY(:,m+1:end); 

rowDist_px = hypot(XR - XL, YR - YL);         
W_px       = mean(rowDist_px, 2, 'omitnan');  

if ~isnan(mm_per_px)
    W_mm = W_px * mm_per_px;
    W0   = mean(W_mm(1:min(5,end)),'omitnan');
    eps_x = (W_mm - W0) / W0;                 
else
    W0   = mean(W_px(1:min(5,end)),'omitnan');
    eps_x = (W_px - W0) / W0;                  
end

%% ===== Load axial strain from Excel/CSV =====
[fn,pn] = uigetfile({'*.xlsx;*.xls;*.csv','Excel/CSV Files'}, ...
                     'Select file with axial strain (time + epsilon_y)');
if isequal(fn,0)
    error('Wrong file');
end
T = readtable(fullfile(pn,fn));

varNamesLower = lower(string(T.Properties.VariableNames));

itime = find(contains(varNamesLower, "time") | varNamesLower=="t", 1);

istr  = find(contains(varNamesLower, "epsilon_y") | ...
             contains(varNamesLower, "eps_y")     | ...
             contains(varNamesLower, "strain_y")  | ...
             contains(varNamesLower, "axialstrain"), 1);
if isempty(istr)
    istr = find(contains(varNamesLower, "epsilon") | ...
                contains(varNamesLower, "strain"), 1);
end

if isempty(itime) || isempty(istr)
    [idx,ok] = listdlg('PromptString','Select TIME column:', ...
                       'ListString', T.Properties.VariableNames, ...
                       'SelectionMode','single');
    if ~ok, error('Did not choose time colomn'); end
    itime = idx;

    [idx,ok] = listdlg('PromptString','Select AXIAL STRAIN column:', ...
                       'ListString', T.Properties.VariableNames, ...
                       'SelectionMode','single');
    if ~ok, error('Did not choose time colomn'); end
    istr = idx;
end

time_ext = T{:,itime};
eps_y_ext = T{:,istr};

good = isfinite(time_ext) & isfinite(eps_y_ext);
time_ext = time_ext(good); eps_y_ext = eps_y_ext(good);
if ~isempty(time_ext)
    time_ext = time_ext - time_ext(1);
end

eps_y = interp1(time_ext, eps_y_ext, t, 'linear', 'extrap');

eps_y(abs(eps_y) < 1e-8) = NaN;

%% ===== Poisson's ratio =====
nu = eps_x ./ eps_y;

%% ===== Plots =====
figure('Name','Results','Color','w');

subplot(3,1,1);
if ~isnan(mm_per_px)
    plot(t, W_mm, 'LineWidth',1.3);
    ylabel('Width (mm)');
else
    plot(t, W_px, 'LineWidth',1.3);
    ylabel('Width (px)');
end
grid on; title('Average lateral width');

subplot(3,1,2);
plot(t, eps_x, 'LineWidth',1.3);
ylabel('\epsilon_x'); grid on;

subplot(3,1,3);
plot(eps_y, nu, 'o-','LineWidth',1.2);
xlabel('Axial strain \epsilon_y'); ylabel('Poisson''s ratio \nu');
grid on; title('Poisson ratio vs axial strain');

%% ===== Export CSV =====
if ~isnan(mm_per_px)
    out = table(t, W_mm, eps_x, eps_y, nu, ...
        'VariableNames', {'time_s','width_mm','epsilon_x','epsilon_y','nu'});
else
    out = table(t, W_px, eps_x, eps_y, nu, ...
        'VariableNames', {'time_s','width_px','epsilon_x','epsilon_y','nu'});
end
writetable(out, 'poisson_extY_re_entrant_0%.csv');
disp('Saved: poisson_from_video_with_extY.csv');
if ~isnan(mm_per_px)
    fprintf('Calibration: %.6f mm/px\n', mm_per_px);
end
