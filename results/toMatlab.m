clear; clc;

params = dlmread('parameters.txt');
ix = 1;
ySz = params(ix); ix = ix + 1;
bSSz = params(ix); ix = ix + 1;
bLSz = params(ix); ix = ix + 1;
chPrSz = params(ix); ix = ix + 1;
crra = params(ix); ix = ix + 1;
beta = params(ix); ix = ix + 1;
rf = params(ix); ix = ix + 1;
eta = params(ix); ix = ix + 1;
etaA = params(ix); ix = ix + 1;
deltaS = params(ix); ix = ix + 1;
deltaL = params(ix); ix = ix + 1;
alpha = params(ix); ix = ix + 1;
muS = params(ix); ix = ix + 1;
ldb0 = params(ix); ix = ix + 1;
ldb1 = params(ix); ix = ix + 1;
rho = params(ix); ix = ix + 1;
rhoDef = params(ix); ix = ix + 1;
rhoNash = params(ix); ix = ix + 1;
adj = params(ix); ix = ix + 1;
bSshare = params(ix);
clear params ix;

kS = deltaS + rf;
kL = deltaL + rf;
dSS = (1 + rf) / (deltaS + rf);
dLL = (1 + rf) / (deltaL + rf);

y = dlmread('../shocks/y.tab');

exoTran = reshape(dlmread('../shocks/PI.tab'), [ySz, ySz]);

bS = dlmread('bS.txt');
bL = dlmread('bL.txt');

Vaut = dlmread('Vaut.txt');

defPr = reshape(dlmread('defPr.txt'), [ySz, bSSz, bLSz]);
V = reshape(dlmread('V.txt'), [ySz, bSSz, bLSz]);
Va = reshape(dlmread('Va.txt'), [ySz, bSSz, bLSz]);
Vd = dlmread('Vd.txt');
qS = reshape(dlmread('qS.txt'), [ySz, bSSz, bLSz]);
qL = reshape(dlmread('qL.txt'), [ySz, bSSz, bLSz]);

qSa = reshape(dlmread('qSa.txt'), [ySz, bSSz, bLSz]);
qLa = reshape(dlmread('qLa.txt'), [ySz, bSSz, bLSz]);

xiS = reshape(dlmread('xiS.txt'), [ySz, bSSz, bLSz]);
xiL = reshape(dlmread('xiL.txt'), [ySz, bSSz, bLSz]);

chPr = reshape(dlmread('chPr.txt'), [ySz, bSSz, bLSz, chPrSz]);
chS = reshape(dlmread('chS.txt'), [ySz, bSSz, bLSz, chPrSz]);
chL = reshape(dlmread('chL.txt'), [ySz, bSSz, bLSz, chPrSz]);

chPrA = reshape(dlmread('chPrA.txt'), [ySz, bSSz, bLSz, chPrSz]);
chSa = reshape(dlmread('chSa.txt'), [ySz, bSSz, bLSz, chPrSz]);
chLa = reshape(dlmread('chLa.txt'), [ySz, bSSz, bLSz, chPrSz]);

gammaPr = reshape(dlmread('gammaPr.txt'), [ySz, chPrSz]);
gammaS = reshape(dlmread('gammaS.txt'), [ySz, chPrSz]);
gammaL = reshape(dlmread('gammaL.txt'), [ySz, chPrSz]);

EdCre = dlmread('EdCre.txt');

save results.mat;
% delete gamma*.txt
% delete ch*.txt
% delete xiS.txt
% delete xiL.txt
% delete qS.txt
% delete qL.txt
% delete qSa.txt
% delete qLa.txt
% delete V.txt
% delete Vd.txt
% delete Vaut.txt
% delete defPr.txt
