load results.mat;

% yPath = dlmread('yPath.txt');
% dPath = dlmread('dPath.txt');
% 
% gdpPath = dlmread('gdpPath.txt');
% cPath = dlmread('cPath.txt');
% nxGdp = (gdpPath - cPath) ./ gdpPath;
% 
% qSpath = dlmread('qSpath.txt');
% qLpath = dlmread('qLpath.txt');
% 
% bSpath = dlmread('bSpath.txt');
% bLpath = dlmread('bLpath.txt');
% 
% bSsim = bS(bSpath);
% bLsim = bL(bLpath);


delimiter = '\t';
dataArray = dlmread("simulation.tab", delimiter, 2, 0);

ix = 1;
timePath = dataArray(:, ix); ix = ix + 1;
yIxPath = dataArray(:, ix); ix = ix + 1;
yPath = dataArray(:, ix); ix = ix + 1;
bSixPath = dataArray(:, ix); ix = ix + 1;
bSpath = dataArray(:, ix); ix = ix + 1;
bLixPath = dataArray(:, ix); ix = ix + 1;
bLpath = dataArray(:, ix); ix = ix + 1;
dPath = dataArray(:, ix); ix = ix + 1;
qSpath = dataArray(:, ix); ix = ix + 1;
qLpath = dataArray(:, ix); ix = ix + 1;
gdpPath = dataArray(:, ix); ix = ix + 1;
cPath = dataArray(:, ix); % ix = ix + 1;

clear delimiter dataArray ix;

nxGdp = (gdpPath - cPath) ./ gdpPath;
durPrPath = dSS + (dLL - dSS) * bLpath ./ (bSpath + bLpath);

bSstock = zeros(size(bSpath));
bLstock = zeros(size(bLpath));
bSstock(2:end) = bSpath(1:end-1);
bLstock(2:end) = bLpath(1:end-1);

ellL = bLpath - (1.0 - deltaL) * bLstock;
ellL(ellL < 0.0) = 0.0;
durIssPath = dSS + (dLL - dSS) * ellL ./ (bSpath + ellL);

spSsim = kS * (1 ./ qSpath - 1);
spLsim = kL * (1 ./ qLpath - 1);

simSz = size(bSpath, 1);

validWindow = 3;
valid = false(size(bLpath));
for tt = (validWindow+1):simSz
    valid(tt) = (sum(dPath((tt-validWindow):tt)) == 0);
end

validDur = valid & (bSpath + bLpath > 0.0);
validDurI = validDur & (bSpath + ellL > 0.0);

ssppS = zeros(size(y));
ssppL = zeros(size(y));
mm = zeros(size(y));
ssttdd = zeros(size(y));
p90 = zeros(size(y));
p10 = zeros(size(y));
mmStock = zeros(size(y));
ssttddStock = zeros(size(y));
p90stock = zeros(size(y));
p10stock = zeros(size(y));
for yIx = 1:ySz
    trgt = valid & yIxPath == yIx & (bSpath + bLpath) > 0.0;

    ssppS(yIx) = nanmean(spSsim(trgt));
    ssppL(yIx) = nanmean(spLsim(trgt));

    mmStock(yIx) = nanmean(durPrPath(trgt));
    ssttddStock(yIx) = nanstd(durPrPath(trgt));
    p90stock(yIx) = quantile(durPrPath(trgt), 0.9);
    p10stock(yIx) = quantile(durPrPath(trgt), 0.1);
    
    mm(yIx) = nanmean(durIssPath(trgt));
    ssttdd(yIx) = nanstd(durIssPath(trgt));
    p90(yIx) = quantile(durIssPath(trgt), 0.9);
    p10(yIx) = quantile(durIssPath(trgt), 0.1);
end
if usejava('desktop')
  figure(5);
  plot(y, mm, 'bo-', y, p90,'b--', y, p10, 'b--');
  hold on;
  plot(y, mmStock, 'ro-', y, p90stock,'r--', y, p10stock, 'r--');
  xlabel('y');
  ylabel('duration');
  
  figure(25);
  plot(y, ssppS, 'bo-', y, ssppL, 'ro--');
  xlabel('y');
  ylabel('spread');
end

labels = true;

poss = 0;
acted = 0;
for ix = 2:size(dPath, 1)
    if dPath(ix-1) == 0
        poss = poss + 1;
    end
    if dPath(ix) == 1 && dPath(ix-1) == 0
        acted = acted + 1;
    end
end
defProbSim = acted / poss * 100;

labels = true;

if (labels), fprintf('%25s ', 'Def Probability (%): '); end
fprintf('%10.5f \n', defProbSim);

if (labels), fprintf('%25s ', 'Time in def (%): '); end
fprintf('%10.5f \n', mean(dPath == 1) * 100);

if (labels), fprintf('%25s ', 'Time in arrears (%): '); end
fprintf('%10.5f \n', mean(dPath == 2) * 100);

if (labels), fprintf('%25s ', 'Sp S (%): '); end
fprintf('%10.5f \n', mean(spSsim(valid)) * 100);
if (labels), fprintf('%25s ', 'Std Sp S (%): '); end
fprintf('%10.5f \n', std(spSsim(valid)) * 100);

if (labels), fprintf('%25s ', 'Sp L (%): '); end
fprintf('%10.5f \n', mean(spLsim(valid)) * 100);
if (labels), fprintf('%25s ', 'Std Sp L (%): '); end
fprintf('%10.5f \n', std(spLsim(valid)) * 100);

if (labels), fprintf('%25s ', 'Debt/GDP (%): '); end
fprintf('%10.5f \n', nanmean((bSstock(valid) + bLstock(valid)) ./ gdpPath(valid)) * 100);

if (labels), fprintf('%25s ', 'S Debt/GDP (%): '); end
fprintf('%10.5f \n', nanmean(bSstock(valid) ./ gdpPath(valid)) * 100);

if (labels), fprintf('%25s ', 'L Debt/GDP (%): '); end
fprintf('%10.5f \n', nanmean(bLstock(valid) ./ gdpPath(valid)) * 100);

if (labels), fprintf('%25s ', 'S / (S+L) (%): '); end
fprintf('%10.5f \n', nanmean(bSstock(valid) ./ (bSstock(valid) + bLstock(valid))) * 100);

if (labels), fprintf('%25s ', 'Corr durPr, GDP (%): '); end
fprintf('%10.5f \n', corr(gdpPath(validDur), durPrPath(validDur)) * 100);

if (labels), fprintf('%25s ', 'Corr durPr, Sp S (%): '); end
fprintf('%10.5f \n', corr(spSsim(validDur), durPrPath(validDur)) * 100);

if (labels), fprintf('%25s ', 'Corr durPr, Sp L (%): '); end
fprintf('%10.5f \n', corr(spLsim(validDur), durPrPath(validDur)) * 100);

if (labels), fprintf('%25s ', 'Corr DurI, GDP (%): '); end
fprintf('%10.5f \n', corr(gdpPath(validDurI), durIssPath(validDurI)) * 100);

if (labels), fprintf('%25s ', 'Corr DurI, Sp S (%): '); end
fprintf('%10.5f \n', corr(spSsim(validDurI), durIssPath(validDurI)) * 100);

if (labels), fprintf('%25s ', 'Corr DurI, Sp L (%): '); end
fprintf('%10.5f \n', corr(spLsim(validDurI), durIssPath(validDurI)) * 100);


if (labels), fprintf('%25s ', 'Std C / Std GDP (%): '); end
fprintf('%10.5f \n', std(cPath(valid)) / std(gdpPath(valid)) * 100);
if (labels), fprintf('%25s ', 'Corr NX/GDP, GDP (%): '); end
fprintf('%10.5f \n', corr(gdpPath(valid), nxGdp(valid)) * 100);

ttmp = squeeze(chPr(:, :, :, end));
if (labels), fprintf('%25s ', 'Max min chPr (%): '); end
fprintf('%10.5f \n', max(ttmp(:)) * 100);

ttmp = squeeze(gammaPr(:, end));
if (labels), fprintf('%25s ', 'Max min gammaPr (%): '); end
fprintf('%10.5f \n', max(ttmp(:)) * 100);

clear ttmp;


