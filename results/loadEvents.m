load results.mat;

evMat = dlmread('events.tab');

retState = evMat(:, 2);
defState = evMat(:, 7);
noEv = size(retState, 1);

yRetEv = y(retState);
yDefEv = y(defState);
% yRetEv = yRetEv - max(0.0, ldb0 .* yRetEv + ldb1 .* yRetEv.^2);

bSev = bS(evMat(:, 3));
bLev = bL(evMat(:, 4));

gammaSev = bS(evMat(:, 5));
gammaLev = bL(evMat(:, 6));

mrkt_gamma_qS = zeros(size(gammaSev));
mrkt_gamma_qL = zeros(size(gammaSev));
for ix = 1:noEv
    mrkt_gamma_qS(ix) = qSa(evMat(ix, 2), evMat(ix, 5), evMat(ix, 6));
    mrkt_gamma_qL(ix) = qLa(evMat(ix, 2), evMat(ix, 5), evMat(ix, 6));
end

durDef = dSS + (dLL - dSS) *  bLev ./ (bSev + bLev);
durRet = dSS + (dLL - dSS) *  gammaLev ./ (gammaSev + gammaLev);

DeltaCred = EdCre(retState);

labels = true;

if (labels), fprintf('%25s ', 'Dur Before: '); end
fprintf('%10.5f \n', nanmean(durDef));
if (labels), fprintf('%25s ', 'Dur After: '); end
fprintf('%10.5f \n', nanmean(durRet));

if (labels), fprintf('%25s ', 'B/Y Before (%): '); end
fprintf('%10.5f \n', mean((bSev + bLev) ./ yDefEv) * 100);
if (labels), fprintf('%25s ', 'B/Y After (%): '); end
fprintf('%10.5f \n', mean((gammaSev + gammaLev) ./ yRetEv) * 100);

if (labels), fprintf('%25s ', 'Agg Haircut (%): '); end
fprintf('%10.5f \n', (1.0 - nanmean( (1 + rf) * (gammaSev + gammaLev) ./ (bSev + bLev) )) * 100);
if (labels), fprintf('%25s ', 'Agg Hair, Mrkt (%): '); end
fprintf('%10.5f \n', (1.0 - nanmean( DeltaCred ./ (mrkt_gamma_qS .* bSev + mrkt_gamma_qL .* bLev) )) * 100);

if (labels), fprintf('%25s ', 'S Haircut (%): '); end
fprintf('%10.5f \n', (1.0 - muS * nanmean( (gammaSev + gammaLev) ./ (muS * bSev + bLev) )) * 100);
if (labels), fprintf('%25s ', 'S Hair Mrkt (%): '); end
fprintf('%10.5f \n', ...
    (1.0 - muS * nanmean( DeltaCred ./ (muS * bSev + bLev) ./ mrkt_gamma_qS ) ) * 100);

if (labels), fprintf('%25s ', 'L Haircut (%): '); end
fprintf('%10.5f \n', (1.0 - nanmean( (1 + rf) * (gammaSev + gammaLev) ./ (muS * bSev + bLev) )) * 100);
if (labels), fprintf('%25s ', 'L Hair Mrkt (%): '); end
fprintf('%10.5f \n', ...
    (1.0 - nanmean( DeltaCred ./ (muS * bSev + bLev) ./ mrkt_gamma_qL )) * 100);


mm = zeros(size(y));
ssttdd = zeros(size(y));
p90 = zeros(size(y));
p10 = zeros(size(y));
for yIx = 1:ySz
    trgt = retState == yIx;
    mm(yIx) = nanmean(durRet(trgt));
    ssttdd(yIx) = std(durRet(trgt));
    p90(yIx) = quantile(durRet(trgt), 0.9);
    p10(yIx) = quantile(durRet(trgt), 0.1);
end
if usejava('desktop')
  figure(555)
  plot(y, mm, 'bo-', y, p90,'b--', y, p10, 'b--');
  xlabel('y');
  ylabel('duration');
  title('\gamma');
end
