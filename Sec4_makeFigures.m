% Get data from R somehow...
tv = load('delTV.txt');
loc = load('delLoc.txt');
cElev = load('delElev.txt');
cGrad = load('delGrad.txt');
border = load('delBorder.txt');
obs = load('delObs.txt');
oLoc = load('delOLoc.txt');
thetaMean = load('delMean.txt');
modSD  = load('delSD.txt');
modSD2 = load('delSD2.txt');
covMat = load('delCovMat.txt');
corMat = load('delCorMat.txt');
sdMat  = load('delSDMat.txt');
covMatMAP = load('delMAPCovMat.txt');
corMatMAP = load('delMAPCorMat.txt');
sdMatMAP  = load('delMAPSDMat.txt');
avRange = load('delAvRange.txt');
avStdDev = load('delAvStdDev.txt');
predOM = load('delPredOnMesh.txt');


% Get mask to remove sea
mask = makeSeaMask(loc, tv, border);

% Plot Prediction
hmm = figure;
predOM(mask == 0) = nan;
h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(loc(:,1)), 1), predOM);
set(h, 'edgecolor', 'none')
set(gca, 'DataAspectRatio', [1, 1, 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14)
ylabel('Northing (km)', 'FontSize', 14)
colorbar;
caxis([0.35, 3.7])
view(0, 90);
hold on;
xlim([-170, 730]);
ylim([6350, 7550]);
shading interp;
plot(border(:,1), border(:,2), 'k');
print('-dpng', '-r300', 'results/Precipitation/predPrec.png')
close(hmm)

% Plot elevation covariate
hmm = figure;
cElev(mask == 0) = nan;
h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(loc(:,1)), 1), cElev);
set(h, 'edgecolor', 'none')
set(gca, 'DataAspectRatio', [1, 1, 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14)
ylabel('Northing (km)', 'FontSize', 14)
colorbar;
view(0, 90);
hold on;
xlim([-170, 730]);
ylim([6350, 7550]);
shading interp;
plot(border(:,1), border(:,2), 'k');
print('-dpng', '-r300', 'results/Precipitation/covElev.png')
close(hmm)

% Plot gradient covariate
hmm = figure;
cGrad(mask == 0) = nan;
h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(loc(:,1)), 1), cGrad);
set(h, 'edgecolor', 'none')
set(gca, 'DataAspectRatio', [1, 1, 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14)
ylabel('Northing (km)', 'FontSize', 14)
colorbar;
view(0, 90);
hold on;
xlim([-170, 730]);
ylim([6350, 7550]);
shading interp
plot(border(:,1), border(:,2), 'k');
print('-dpng', '-r300', 'results/Precipitation/covGrad.png')
close(hmm)

% Plot observations
hmm = figure;
scatter(oLoc(:,1), oLoc(:,2), 20, obs, 'filled');
hold on;
plot(border(:,1), border(:,2), 'k');
xlim([-150, 550]);
ylim([6400, 7300]);
set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14);
ylabel('Northing (km)', 'FontSize', 14);
colorbar;
print('-dpng', '-r300', 'results/Precipitation/pNorwayObservation.png');
close(hmm)

% Plot mean effect in range
hmm = figure;
h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(loc(:,1)), 1), thetaMean(2)+thetaMean(5)*cElev+thetaMean(6)*cGrad);
set(h, 'edgecolor', 'none');
set(gca, 'DataAspectRatio', [1, 1, 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14)
ylabel('Northing (km)', 'FontSize', 14)
colorbar;
view(0, 90);
hold on;
xlim([-150, 550]);
ylim([6400, 7300]);
shading interp
plot(border(:,1), border(:,2), 'k');
print('-dpng', '-r300', 'results/Precipitation/effectRange.png')
close(hmm)

% Plot mean effect in std.dev.
hmm = figure;
h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(loc(:,1)), 1), thetaMean(3)+thetaMean(8)*cElev+thetaMean(9)*cGrad);
set(h, 'edgecolor', 'none');
set(gca, 'DataAspectRatio', [1, 1, 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14)
ylabel('Northing (km)', 'FontSize', 14)
colorbar;
view(0, 90);
hold on;
xlim([-150, 550]);
ylim([6400, 7300]);
shading interp
plot(border(:,1), border(:,2), 'k');
print('-dpng', '-r300', 'results/Precipitation/effectStdDev.png')
close(hmm)

% Check separation of variance and correlation
hmm = figure;
h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(modSD), 1), modSD);
set(h, 'edgecolor', 'none');
set(gca, 'DataAspectRatio', [1, 1, 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14)
ylabel('Northing (km)', 'FontSize', 14)
colorbar;
view(0, 90);
hold on;
xlim([-500, 800]);
ylim([6100, 7550]);
shading interp
plot(border(:,1), border(:,2), 'k');
print('-dpng', '-r300', 'results/Precipitation/StdDev1.png')
close(hmm)

% Check separation of variance and correlation
hmm = figure;
h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(modSD2),1), modSD2);
set(h, 'edgecolor', 'none');
set(gca, 'DataAspectRatio', [1, 1, 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14)
ylabel('Northing (km)', 'FontSize', 14)
colorbar;
view(0, 90);
hold on;
xlim([-500, 800]);
ylim([6100, 7550]);
shading interp
plot(border(:,1), border(:,2), 'k');
print('-dpng', '-r300',  'results/Precipitation/StdDev2.png')
close(hmm)


    % Expected standard deviation
    hmm = figure;
    sdMat(mask == 0) = nan;

    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(sdMat),1), sdMat);
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 14);
    xlabel('Easting (km)', 'FontSize', 14)
    ylabel('Northing (km)', 'FontSize', 14)
    colorbar;
    view(0, 90);
    hold on;
  xlim([-170, 730]);
ylim([6350, 7550]);

shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300', 'results/Precipitation/expSD.png')
    close(hmm)
    
    % Expected covariance
    num = 1
    hmm = figure;
    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(covMat(:,num)),1), covMat(:,num));
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 14);
    xlabel('Easting (km)', 'FontSize', 14)
    ylabel('Northing (km)', 'FontSize', 14)
    colorbar;
    view(0, 90);
    hold on;
    xlim([-500, 800]);
    ylim([6100, 7550]);
    shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300', 'results/Precipitation/expCov.png')
    close(hmm)
    
    % Expected correlation
    num = 1
    hmm = figure;
    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(corMat(:,num)),1), corMat(:,num));
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 14);
    xlabel('Easting (km)', 'FontSize', 14)
    ylabel('Northing (km)', 'FontSize', 14)
    colorbar;
    view(0, 90);
    hold on;
    xlim([-500, 800]);
    ylim([6100, 7550]);
    shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300', 'results/Precipitation/expCor1.png')
    close(hmm)
    
    % Expected correlation
    num = 3
    hmm = figure;
    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(corMat(:,num)),1), corMat(:,num));
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 14);
    xlabel('Easting (km)', 'FontSize', 14)
    ylabel('Northing (km)', 'FontSize', 14)
    colorbar;
    view(0, 90);
    hold on;
    xlim([-500, 800]);
    ylim([6100, 7550]);
    shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300',  'results/Precipitation/expCor1.png')
    close(hmm)
    
    % Plot contour plot
    ll = [0.905; 0.567; 0.356; 0.223];
    loc1 = loc(365, 1:2);
    loc2 = loc(120, 1:2);
    hmm = figure;
    tricontour(tv, loc(:,1), loc(:,2), corMat(:,1), ll);
    hold on;
    tricontour(tv, loc(:,1), loc(:,2), corMat(:,3), ll);
    Fig1Ax1 = get(1, 'Children');
    Fig1Ax1Line1 = get(Fig1Ax1, 'Children');
    set(Fig1Ax1Line1, 'LineWidth', 2);
    plot(loc1(1), loc1(2), 'xr', 'LineWidth', 2);
    plot(loc2(1), loc2(2), 'xr', 'LineWidth', 2);
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 30);
    xlabel('Easting (km)', 'FontSize', 30)
    ylabel('Northing (km)', 'FontSize', 30)
    colorbar;
    caxis([0, 1]);
    view(0, 90);
    hold on;
    xlim([-500, 800]);
    ylim([6100, 7550]);
    plot(border(:,1), border(:,2), 'k');
    print('-depsc2', 'results/Precipitation/corCurves.eps')
    close(hmm)
 
    
    % MAP standard deviation
    hmm = figure;
    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(sdMat),1), sdMatMAP);
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 14);
    xlabel('Easting (km)', 'FontSize', 14)
    ylabel('Northing (km)', 'FontSize', 14)
    colorbar;
    view(0, 90);
    hold on;
    xlim([-500, 800]);
    ylim([6100, 7550]);
    shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300', 'results/Precipitation/expSD.png')
    close(hmm)
    
    % MAP covariance
    num = 2
    hmm = figure;
    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(covMat(:,num)),1), covMatMAP(:,num));
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 14);
    xlabel('Easting (km)', 'FontSize', 14)
    ylabel('Northing (km)', 'FontSize', 14)
    colorbar;
    view(0, 90);
    hold on;
    xlim([-500, 800]);
    ylim([6100, 7550]);
    shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300', 'results/Precipitation/expCov.png')
    close(hmm)
    
    % MAP correlation
    num = 1
    hmm = figure;
    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(corMat(:,num)),1), corMatMAP(:,num));
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 30);
    xlabel('Easting (km)', 'FontSize', 30)
    ylabel('Northing (km)', 'FontSize', 30)
    colorbar;
    view(0, 90);
    hold on;
    xlim([-500, 800]);
    ylim([6100, 7550]);
    shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300',  'results/Precipitation/expCor.png')
    close(hmm)
    
    % Average range
    hmm = figure;
        avRange(mask == 0) = nan;

    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(avRange),1), avRange);
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 14);
    xlabel('Easting (km)', 'FontSize', 14)
    ylabel('Northing (km)', 'FontSize', 14)
    colorbar;
    view(0, 90);
    hold on;
xlim([-170, 730]);
ylim([6350, 7550]);

    shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300', 'results/Precipitation/avRange.png')
    close(hmm)
    
    % Average std.dev.
    hmm = figure;
    avStdDev(mask == 0) = nan;

    h = trisurf(tv, loc(:,1), loc(:,2), zeros(length(avStdDev),1), avStdDev);
    set(h, 'edgecolor', 'none');
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    set(gca, 'FontSize', 14);
    xlabel('Easting (km)', 'FontSize', 14)
    ylabel('Northing (km)', 'FontSize', 14)
    colorbar;
    view(0, 90);
    hold on;
  xlim([-170, 730]);
ylim([6350, 7550]);

    shading interp
    plot(border(:,1), border(:,2), 'k');
    print('-dpng', '-r300', 'results/Precipitation/avStdDev.png')
    close(hmm)
    
    
% Plot mesh
hmm = figure;
h = trimesh(tv, loc(:,1), loc(:,2));
set(h, 'Color', 'black');
set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'FontSize', 14);
xlabel('Easting (km)', 'FontSize', 14)
ylabel('Northing (km)', 'FontSize', 14)
hold on;
plot(border(:,1), border(:,2), 'k')
xlim([-500, 800]);
ylim([6100, 7550]);
print('-depsc2', 'results/Precipitation/mesh.eps')
close(hmm);