%*****IEEE CEC 2022 Competition on Dynamic Optimization Problems Generated by Generalized Moving Peaks Benchmark******
%Author: Danial Yazdani
%Last Edited: December 06, 2021
%
% ------------
% Reference:
% ------------
%  
%  D. Yazdani et al.,
%            "Benchmarking Continuous Dynamic Optimization: Survey and Generalized Test Suite,"
%            IEEE Transactions on Cybernetics (2020).
% 
%  D. Yazdani et al.,
%            "Generalized Moving Peaks Benchmark," arXiv:2106.06174, (2021).
% 
%  T. Blackwell and J. Branke,
%            "Multiswarms, exclusion, and anti-convergence in dynamic environments"
%            IEEE Transactions on Evolutionary Computation (2006).
% ------------
% Notification:
% ------------
% This code solves Generalized Moving Peaks Benchmark (GMPB) by mQSO.
% It is assumed that the environmental changes are VISIBLE, therefore,
% mQSO is informed about changes (i.e., mQSO does not need to detect
% environmental changes). Also note that mQSO does not access to a prior knowledge
% about the shift severity value. The shift severity is learned in this code.
%
% 
% -------
% Inputs:
% -------
%
%    The Participants can set peak number, change frequency, dimension, 
%    and shift severity in lines 59-62 of "main.m" according to the
%    competition instractions available in the following link:
%                                
%                 https://www.danialyazdani.com/CEC-2022
% 
%
% ------------
% Output:
% ------------
% 
% Offline error
% 
% --------
% License:
% --------
% This program is to be used under the terms of the GNU General Public License
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Danial Yazdani
% e-mail: danial.yazdani AT gmail dot com
%         danial.yazdani AT yahoo dot com
% Copyright notice: (c) 2021 Danial Yazdani
%*********************************************************************************************************************
function [result,Problem] = fitness(X,Problem)
[SolutionNumber,~] = size(X);
result = NaN(SolutionNumber,1);
for jj=1 : SolutionNumber
    if Problem.FE >= Problem.MaxEvals || Problem.RecentChange == 1
        return;
    end
    x = X(jj,:)';
    f=NaN(1,Problem.PeakNumber);
    for k=1 : Problem.PeakNumber
        a = Transform((x - Problem.PeaksPosition(k,:,Problem.Environmentcounter)')'*Problem.RotationMatrix{Problem.Environmentcounter}(:,:,k)',Problem.tau(Problem.Environmentcounter,k),Problem.eta(k,:,Problem.Environmentcounter));
        b = Transform(Problem.RotationMatrix{Problem.Environmentcounter}(:,:,k) * (x - Problem.PeaksPosition(k,:,Problem.Environmentcounter)'),Problem.tau(Problem.Environmentcounter,k),Problem.eta(k,:,Problem.Environmentcounter));
        f(k) = Problem.PeaksHeight(Problem.Environmentcounter,k) - sqrt( a * diag(Problem.PeaksWidth(k,:,Problem.Environmentcounter).^2) * b);
    end
    result(jj) = max(f);
    Problem.FE = Problem.FE + 1;
    SolutionError = Problem.OptimumValue(Problem.Environmentcounter) - result(jj);
    if rem(Problem.FE , Problem.ChangeFrequency)~=1
        if Problem.CurrentError(Problem.FE-1)<SolutionError
            Problem.CurrentError(Problem.FE) = Problem.CurrentError(Problem.FE-1);
        else
            Problem.CurrentError(Problem.FE) = SolutionError;
        end
    else
        Problem.CurrentError(Problem.FE) =  SolutionError;
    end
    if ~rem(Problem.FE , Problem.ChangeFrequency) && Problem.FE < Problem.MaxEvals %% 每ChangeFrequency代改变一次环境
        Problem.Environmentcounter = Problem.Environmentcounter+1;
        Problem.RecentChange = 1;
    end
end
end

function Y = Transform(X,tau,eta)
Y = X;
tmp = (X > 0);
Y(tmp) = log(X(tmp));
Y(tmp) = exp(Y(tmp) + tau*(sin(eta(1).*Y(tmp)) + sin(eta(2).*Y(tmp))));
tmp = (X < 0);
Y(tmp) = log(-X(tmp));
Y(tmp) = -exp(Y(tmp) + tau*(sin(eta(3).*Y(tmp)) + sin(eta(4).*Y(tmp))));
end