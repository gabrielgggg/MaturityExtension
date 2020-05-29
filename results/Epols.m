clear; clc;
load results.mat;

EbSpr = zeros([ySz bSSz bLSz]);
EbLpr = zeros([ySz bSSz bLSz]);
EdurPr = zeros([ySz bSSz bLSz]);

for yIx = 1:ySz
    for bSix = 1:bSSz
        for bLix = 1:bLSz
            if defPr(yIx, bSix, bLix) < 0.95
                probs = squeeze(chPr(yIx, bSix, bLix, :)) ./ (1-defPr(yIx, bSix, bLix));
                Ses = squeeze(bS(chS(yIx, bSix, bLix, :)));
                Les = squeeze(bL(chL(yIx, bSix, bLix, :)));
                lRatio = Les ./ (Les + Ses);
                lRatio(isnan(lRatio)) = 0.0;
                
                EbSpr(yIx, bSix, bLix) = dot( Ses, probs );
                EbLpr(yIx, bSix, bLix) = dot( Les, probs );
                EdurPr(yIx, bSix, bLix) = dSS + (dLL - dSS) * dot( lRatio , probs );
            else
                EbSpr(yIx, bSix, bLix) = NaN;
                EbLpr(yIx, bSix, bLix) = NaN;
                EdurPr(yIx, bSix, bLix) = NaN;
            end
        end
    end
end

EbSaPr = zeros([ySz bSSz bLSz]);
EbLaPr = zeros([ySz bSSz bLSz]);
EdurApr = zeros([ySz bSSz bLSz]);

for yIx = 1:ySz
    for bSix = 1:bSSz
        for bLix = 1:bLSz
            probs = squeeze(chPrA(yIx, bSix, bLix, :));
            Ses = squeeze(bS(chSa(yIx, bSix, bLix, :)));
            Les = squeeze(bL(chLa(yIx, bSix, bLix, :)));
            lRatio = Les ./ (Les + Ses);
            lRatio(isnan(lRatio)) = 0.0;

            EbSaPr(yIx, bSix, bLix) = dot( Ses, probs );
            EbLaPr(yIx, bSix, bLix) = dot( Les, probs );
            EdurApr(yIx, bSix, bLix) = dSS + (dLL - dSS) * dot( lRatio , probs );
        end
    end
end

EgammaS = zeros([ySz 1]);
EgammaL = zeros([ySz 1]);
EgammaPr = zeros([ySz 1]);

for yIx = 1:ySz
    probs = squeeze(gammaPr(yIx, :));
    Ses = squeeze(bS(gammaS(yIx, :)));
    Les = squeeze(bL(gammaL(yIx, :)));
    lRatio = Les ./ (Les + Ses);
    lRatio(isnan(lRatio)) = 0.0;

    EgammaS(yIx) = dot( Ses, probs );
    EgammaL(yIx) = dot( Les, probs );
    EgammaPr(yIx) = dSS + (dLL - dSS) * dot( lRatio , probs );
end
