function [b,bse,bt,r2,F] = regress_jc(y,X,display)
% This is the MATLAB function regress modified by John Cochrane to 
% display more useful output
% 
%
%  REGRESS_JC Multiple linear regression using least squares.
%   b = regress_jc(y,X,display) 
%   returns the vector of regression coefficients, b,
%   in the linear model  y = Xb, (X is an nxp matrix, y is the nx1  vector of observations). 
%   NOTE if you want a constant in the regression, the first column of X needs to be 
%   a column of ones. 
%   Display = 1 prints the regression. Display = 0 suppresses printout 
%   Display = 0 is useful when you're runnning a lot of regressions and you want to collect
%   the results in a table. 
%   If you omit the "display" option it will default to one and display results
% 
%    [B,BSE,BT,R2,F] = REGRESS_JC(y,X,display) 
%    BSE: standard erorrs
%    BT: t statistics
%    R2: Rsquared
%    F: F statistic
%
%   The X matrix should include a column of ones so that the model
%   contains a constant term.  The F and p values are computed under
%   the assumption that the model contains a constant term, and they
%   are not correct for models without a constant.  The R-square
%   value is the ratio of the regression sum of squares to the
%   total sum of squares.
% 
%   Changes from matlab function: 
%   returns standard errors and t statistics, not confidence
%   intervals for betas. 
%   removes silly "confidence errors for the residuals" 
%   adds option to display results
%   removes p values so it will run without statistics toolbox. 

%   References:
%      [1] Samprit Chatterjee and Ali S. Hadi, "Influential Observations,
%      High Leverage Points, and Outliers in Linear Regression",
%      Statistical Science 1986 Vol. 1 No. 3 pp. 379-416. 
%      [2] N. Draper and H. Smith, "Applied Regression Analysis, Second
%      Edition", Wiley, 1981.

%   B.A. Jones 3-04-93
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.13 $  $Date: 2002/01/17 21:31:53 $


if nargin < 3
    display = 1;
end

if  nargin < 2,              
    error('REGRESS requires at least two input arguments, y and X.');      
end 

alpha = 0.05; % all tests are 5% now. 

% Check that matrix (X) and left hand side (y) have compatible dimensions
[n,p] = size(X);
[n1,collhs] = size(y);
if n ~= n1, 
    error('The number of rows in Y must equal the number of rows in X.'); 
end 

if collhs ~= 1, 
    error('Y must be a vector, not a matrix'); 
end

% Remove missing values, if any
wasnan = (isnan(y) | any(isnan(X),2));
if (any(wasnan))
   y(wasnan) = [];
   X(wasnan,:) = [];
   n = length(y);
end

% Find the least squares solution.
[Q, R]=qr(X,0);
b = R\(Q'*y);

% Find a confidence interval for each component of x
% Draper and Smith, equation 2.6.15, page 94

if (size(R,1)>=size(R,2))
   RI = R\eye(p);
   xdiag=sqrt(sum((RI .* RI)',1))';
   T = X*RI;
else
   xdiag = NaN*ones(size(b));
   T = NaN*ones(size(y));
end
nu = max(0,n-p);                % Residual degrees of freedom
yhat = X*b;                     % Predicted responses at each data point.
r = y-yhat;                     % Residuals.
if nu ~= 0
   rmse = norm(r)/sqrt(nu);        % Root mean square error.
else
   rmse = Inf;
end
s2 = rmse^2;                    % Estimator of error variance.
if (nargout>=2)|display; 
   %tval = tinv((1-alpha/2),nu);  % Commented out so it will run without
   %statistics toolbox. If you want to compute t statistics, this is how to
   %do it. 
   bse = xdiag*rmse;             % jc addition
   bt = b./bse;                  % jc addition
   %bint = [b-tval*xdiag*rmse, b+tval*xdiag*rmse];
end

% Calculate R-squared.
if (nargout>3)|display,
   RSS = norm(yhat-mean(y))^2;  % Regression sum of squares.
   TSS = norm(y-mean(y))^2;     % Total sum of squares.
   r2 = RSS/TSS;                % R-square statistic.
   if (p>1)
      F = (RSS/(p-1))/s2;       % F statistic for regression
   else
      F = NaN;
   end
%   prob = 1 - fcdf(F,p-1,nu);   % Significance probability for regression
                                 % commented out so it will run without
                                 % statistics toolbox. 
   % All that requires a constant.  Do we have one?
   if (~any(all(X==1)))
      % Apparently not, but look for an implied constant.
      b0 = R\(Q'*ones(n,1));
      if (sum(abs(1-X*b0))>n*sqrt(eps))
         warning(sprintf(['R-square is not well defined unless X has' ...
                       ' a column of ones.\nType "help regress" for' ...
                       ' more information.']));
      end
   end
end

% Restore NaN so inputs and outputs conform
if (nargout>2 & any(wasnan))
   tmp = ones(size(wasnan));
   tmp(:) = NaN;
   tmp(~wasnan) = r;
   r = tmp;
end
if (nargout>3 & any(wasnan))
   tmp = ones(length(wasnan),2);
   tmp(:) = NaN;
   tmp(~wasnan,:) = rint;
   rint = tmp;
end

if display; 
    fprintf('Coeffs   '); fprintf(' %8.3g ',b'); fprintf('\n'); 
    fprintf('Std errs '); fprintf(' %8.3g ',bse'); fprintf('\n'); 
    fprintf('t stats  '); fprintf(' %8.3g ',bt'); fprintf('\n');  
    fprintf('R2: %8.3g F stat: %8.3g  \n',  [r2 F]); 
end; 
