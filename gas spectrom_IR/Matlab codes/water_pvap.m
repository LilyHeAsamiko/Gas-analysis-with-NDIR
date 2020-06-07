function [ pvap ] = water_pvap( T )
% [ pvap ] = water_pvap( T ) calculates partial pressure of water vapour 
% at temperature T (K), result in Pascals.
%
% INPUT
% T - Temperature (K)
%
% OUTPUT
% pvap - partial pressure of water vapour (Pa)

pvap = exp(77.34491296-7235.424651./T-8.2*log(T)+0.0057113*T);
end

