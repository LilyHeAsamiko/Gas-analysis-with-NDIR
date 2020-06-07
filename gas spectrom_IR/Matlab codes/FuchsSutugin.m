function [ out ] = FuchsSutugin( dp, lambda )
% [ out ] = FuchsSutugin( dp, lambda ) Calculates Fuchs-Sutugin
% correction factor for particle diameter dp and free path length
% lambda
%
% INPUT
% dp - particle diameter (m)
% lambda - free path length (m)
%
% OUTPUT
% out - Fuchs-Sutugin correction factor (-)

Kn = 2*lambda./dp;
out = (1 + Kn) ./ (1 + 1.71*Kn + 1.33*Kn.^2);


end

