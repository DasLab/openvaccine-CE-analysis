function [rmse,fit] = expfit_varyamp(p,times,data);
% [rmse,fit] = expfit(p,times,data);

tau = p(1);
amp = p(2);

fit = amp * exp( - times/tau );

rmse = 0.0;
if ~exist( 'data','var') return; end;

rmse = sqrt( sum(( fit-data ).^2 ));
%rmse = sqrt( sum(( log(fit/data)/log(2) ).^2 ));

