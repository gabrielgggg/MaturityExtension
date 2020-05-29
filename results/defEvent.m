clear; clc;
load results.mat;
load simulation.mat;

window = 5;
simSz = size(bLixPath, 1);

ccases = zeros(size(bLixPath));
for ix = (window+1):simSz-(window+1)
    if dPath(ix) == 1 && dPath(ix-1) == 0
        ccases(ix) = 1;
    end
end

totalCases = sum(ccases);

xAxis = (-window):window;
evSz = size(xAxis, 2);

selected = zeros([simSz, evSz]);
for ix = 1:evSz
    shft = xAxis(ix);
    selected((window+1):simSz-(window+1), ix) = (ccases((window+1-shft):(simSz-(window+1)-shft)) == 1)';
end

dev = nan([evSz, 1]);
exev = nan([evSz, 1]);
gdpev = nan([evSz, 1]);
cev = nan([evSz, 1]);
bSev = nan([evSz, 1]);
bLev = nan([evSz, 1]);
spSev = nan([evSz, 1]);
spLev = nan([evSz, 1]);
durPrev = nan([evSz, 1]);
durIssev = nan([evSz, 1]);
nxgdpev = nan([evSz, 1]);

for ix = 1:evSz
    dev(ix) = mean(dPath(selected(:, ix) == 1) == 1);
    exev(ix) = mean(dPath(selected(:, ix) == 1) == 2);
    
    gdpev(ix) = mean(gdpPath(selected(:, ix) == 1));
    cev(ix) = mean(cPath(selected(:, ix) == 1));
    
    bSev(ix) = mean(bSpath(selected(:, ix) == 1));
    bLev(ix) = mean(bLpath(selected(:, ix) == 1));
    
    if xAxis(ix) < 0
        spSev(ix) = mean(spSsim(selected(:, ix) == 1));
        spLev(ix) = mean(spLsim(selected(:, ix) == 1));
    end

    durPrev(ix) = mean(durPrPath(selected(:, ix) == 1));
    durIssev(ix) = mean(durIssPath(selected(:, ix) == 1));
    nxgdpev(ix) = mean(nxGdp(selected(:, ix) == 1));
end

figure(1); set(gcf,'units','normalized','position',[0.25 0.25 0.5 0.5]);
plot(xAxis, dev*100, 'b-o', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]); title('default & exclusion');
hold on; plot(xAxis, exev*100, 'r-x', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]);
xlabel('Years'); set(gca,'FontSize',18); orient(gcf, 'landscape');

figure(2);  set(gcf,'units','normalized','position',[0.25 0.25 0.5 0.5]);
plot(xAxis, gdpev*100, 'b-o', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]); title('GDP');
xlabel('Years'); set(gca,'FontSize',18); orient(gcf, 'landscape');

figure(3);  set(gcf,'units','normalized','position',[0.25 0.25 0.5 0.5]);
plot(xAxis, nxgdpev*100, 'b-o', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]); title('NX/GDP');
xlabel('Years'); set(gca,'FontSize',18); orient(gcf, 'landscape');

figure(4);  set(gcf,'units','normalized','position',[0.25 0.25 0.5 0.5]);
yyaxis left; plot(xAxis, bSev*100, 'b-o', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]); title('b_S and b_L');
yyaxis right; plot(xAxis, bLev*100, 'r-x', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]);
xlabel('Years'); set(gca,'FontSize',18); orient(gcf, 'landscape');

figure(5);  set(gcf,'units','normalized','position',[0.25 0.25 0.5 0.5]);
plot(xAxis, spSev*100, 'b-o', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]); title('Spread S and L (r)');
hold on;
plot(xAxis, spLev*100, 'r-x', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]);
xlabel('Years'); set(gca,'FontSize',18); orient(gcf, 'landscape');

figure(6);  set(gcf,'units','normalized','position',[0.25 0.25 0.5 0.5]);
yyaxis left; plot(xAxis, durIssev, 'b-o', 'LineWidth', 2); xlim([xAxis(1), xAxis(end)]); title('Duration');
yyaxis right; plot(xAxis, durPrev, 'r-x', 'LineWidth', 2); 
xlabel('Years'); set(gca,'FontSize',18); orient(gcf, 'landscape');

