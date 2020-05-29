yIx = 11;
bSix = 5;
bLix = 25;

sface = zeros([bSSz, bLSz]);
for ix = 1:chPrSz
    sface(chS(yIx, bSix, bLix, ix), chL(yIx, bSix, bLix, ix)) = ...
        chPr(yIx, bSix, bLix, ix) ./ (1-defPr(yIx, bSix, bLix));
end
figure(5);
contourf(bL, bS, sface, 3); colorbar();
hold on;
scatter(bS(bSix), bL(bLix), 'or');
xlabel('b_S');
ylabel('b_L');

gface = zeros([bSSz, bLSz]);
for ix = 1:chPrSz
    gface(gammaS(yIx, ix), gammaL(yIx, ix)) = ...
        gammaPr(yIx, ix);
end
figure(10);
contourf(bL, bS, gface, 3); colorbar();
xlabel('\gamma_S');
ylabel('\gamma_L');