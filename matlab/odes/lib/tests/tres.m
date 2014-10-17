clear all;

addpath('../');

spts = linspace(-1,1,10);

% choos a quadratic:
q = spts.^2;

% pick some sub interval
local_spts = linspace( -.25, 0.1, 5 );

[EvalCoeffs Integrate_Coeffs] = New_Res_Coeffs( spts, local_spts );

qvals = (EvalCoeffs * q')';
max( abs( qvals - local_spts(2:end).^2 ) )

qint = ( local_spts.^3 - local_spts(1)^3 ) / 3;
qint = qint(2:end);

myqint = (Integrate_Coeffs * q')';

max( abs( qint - myqint ) )
