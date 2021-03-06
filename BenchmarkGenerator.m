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
function [Problem] = BenchmarkGenerator(PeakNumber,ChangeFrequency,Dimension,ShiftSeverity,EnvironmentNumber)
Problem                     = [];
Problem.FE                  = 0;
Problem.PeakNumber          = PeakNumber;
Problem.ChangeFrequency     = ChangeFrequency;
Problem.Dimension           = Dimension;
Problem.ShiftSeverity       = ShiftSeverity;
Problem.EnvironmentNumber   = EnvironmentNumber;
Problem.Environmentcounter  = 1;
Problem.RecentChange        = 0;
Problem.MaxEvals            = Problem.ChangeFrequency * Problem.EnvironmentNumber;
Problem.Ebbc                = NaN(1,Problem.EnvironmentNumber);
Problem.CurrentError        = NaN(1,Problem.MaxEvals);
Problem.CurrentPerformance  = NaN(1,Problem.MaxEvals);
Problem.MinCoordinate       = -100;
Problem.MaxCoordinate       = 100;
Problem.MinHeight           = 30;
Problem.MaxHeight           = 70;
Problem.MinWidth            = 1;
Problem.MaxWidth            = 12;
Problem.MinAngle            = -pi;
Problem.MaxAngle            = pi;
Problem.MaxTau              = 1;
Problem.MinTau              = -1;
Problem.MaxEta              = 20;
Problem.MinEta              = -20;
Problem.HeightSeverity      = 7;
Problem.OptimumValue        = NaN(Problem.EnvironmentNumber,1);
Problem.OptimumID           = NaN(Problem.EnvironmentNumber,1);
Problem.PeaksHeight         = NaN(Problem.EnvironmentNumber,Problem.PeakNumber);
Problem.PeaksPosition       = NaN(Problem.PeakNumber,Problem.Dimension,Problem.EnvironmentNumber);
Problem.PeaksWidth          = NaN(Problem.PeakNumber,Problem.Dimension,Problem.EnvironmentNumber);
Problem.PeaksPosition(:,:,1)= Problem.MinCoordinate + (Problem.MaxCoordinate-Problem.MinCoordinate)*rand(Problem.PeakNumber,Problem.Dimension);
Problem.PeaksHeight(1,:)    = Problem.MinHeight + (Problem.MaxHeight-Problem.MinHeight)*rand(Problem.PeakNumber,1);
Problem.PeaksWidth(:,:,1)   = Problem.MinWidth + (Problem.MaxWidth-Problem.MinWidth)*rand(Problem.PeakNumber,Problem.Dimension);
[Problem.OptimumValue(1), Problem.OptimumID(1)]     = max(Problem.PeaksHeight(1,:));
Problem.AngleSeverity       = pi/9;
Problem.TauSeverity         = 0.2;
Problem.EtaSeverity         = 2;
Problem.EllipticalPeaks     = 1;
Problem.InitialRotationMatrix = NaN(Problem.Dimension,Problem.Dimension,Problem.PeakNumber);
for ii=1 : Problem.PeakNumber
    [Problem.InitialRotationMatrix(:,:,ii) , ~] = qr(rand(Problem.Dimension));
end
Problem.WidthSeverity       = 1;
Problem.PeaksAngle          = Problem.MinAngle + (Problem.MaxAngle-Problem.MinAngle)*rand(Problem.PeakNumber,1);
Problem.tau                 = Problem.MinTau + (Problem.MaxTau-Problem.MinTau)*rand(Problem.PeakNumber,1);
Problem.RotationMatrix      = cell(1,Problem.EnvironmentNumber);%NaN(Problem.PeakNumber,Problem.Dimension,Problem.EnvironmentNumber);
Problem.RotationMatrix{1}   = Problem.InitialRotationMatrix;
Problem.PeaksAngle          = NaN(Problem.EnvironmentNumber,Problem.PeakNumber);
Problem.tau                 = NaN(Problem.EnvironmentNumber,Problem.PeakNumber);
Problem.eta                 = NaN(Problem.PeakNumber,4,Problem.EnvironmentNumber);
Problem.PeaksAngle(1,:)     = Problem.MinAngle + (Problem.MaxAngle-Problem.MinAngle)*rand(Problem.PeakNumber,1);
Problem.tau(1,:)            = Problem.MinTau + (Problem.MaxTau-Problem.MinTau)*rand(Problem.PeakNumber,1);
Problem.eta(:,:,1)          = Problem.MinEta + (Problem.MaxEta-Problem.MinEta)*rand(Problem.PeakNumber,4);
for ii=2 : Problem.EnvironmentNumber%Generating all environments
    ShiftOffset = randn(Problem.PeakNumber,Problem.Dimension);
    Shift          = (ShiftOffset ./ pdist2(ShiftOffset,zeros(1,Problem.Dimension))).* Problem.ShiftSeverity;
    PeaksPosition  = Problem.PeaksPosition(:,:,ii-1) + Shift;
    PeaksWidth  = Problem.PeaksWidth(:,:,ii-1) + (randn(Problem.PeakNumber,Problem.Dimension).* Problem.WidthSeverity);
    PeaksHeight = Problem.PeaksHeight(ii-1,:) + (Problem.HeightSeverity*randn(1,Problem.PeakNumber));
    PeaksAngle  = Problem.PeaksAngle(ii-1,:) + (Problem.AngleSeverity.*randn(1,Problem.PeakNumber));
    PeaksTau    = Problem.tau(ii-1,:) + (Problem.TauSeverity.*randn(1,Problem.PeakNumber));
    PeaksEta    = Problem.eta(:,:,ii-1) + (randn(Problem.PeakNumber,4).* Problem.EtaSeverity);
    tmp = PeaksAngle > Problem.MaxAngle;
    PeaksAngle(tmp) = (2*Problem.MaxAngle)- PeaksAngle(tmp);
    tmp = PeaksAngle < Problem.MinAngle;
    PeaksAngle(tmp) = (2*Problem.MinAngle)- PeaksAngle(tmp);
    tmp = PeaksTau > Problem.MaxTau;
    PeaksTau(tmp) = (2*Problem.MaxTau)- PeaksTau(tmp);
    tmp = PeaksTau < Problem.MinTau;
    PeaksTau(tmp) = (2*Problem.MinTau)- PeaksTau(tmp);
    tmp = PeaksEta > Problem.MaxEta;
    PeaksEta(tmp) = (2*Problem.MaxEta)- PeaksEta(tmp);
    tmp = PeaksEta < Problem.MinEta;
    PeaksEta(tmp) = (2*Problem.MinEta)- PeaksEta(tmp);
    tmp = PeaksPosition > Problem.MaxCoordinate;
    PeaksPosition(tmp) = (2*Problem.MaxCoordinate)- PeaksPosition(tmp);
    tmp = PeaksPosition < Problem.MinCoordinate;
    PeaksPosition(tmp) = (2*Problem.MinCoordinate)- PeaksPosition(tmp);
    tmp = PeaksHeight > Problem.MaxHeight;
    PeaksHeight(tmp) = (2*Problem.MaxHeight)- PeaksHeight(tmp);
    tmp = PeaksHeight < Problem.MinHeight;
    PeaksHeight(tmp) = (2*Problem.MinHeight)- PeaksHeight(tmp);
    tmp = PeaksWidth > Problem.MaxWidth;
    PeaksWidth(tmp) = (2*Problem.MaxWidth)- PeaksWidth(tmp);
    tmp = PeaksWidth < Problem.MinWidth;
    PeaksWidth(tmp) = (2*Problem.MinWidth)- PeaksWidth(tmp);
    Problem.PeaksPosition(:,:,ii)= PeaksPosition;
    Problem.PeaksHeight(ii,:)    = PeaksHeight;
    Problem.PeaksWidth(:,:,ii)   = PeaksWidth;
    Problem.PeaksAngle(ii,:)     = PeaksAngle;
    Problem.tau(ii,:)            = PeaksTau;
    Problem.eta(:,:,ii)          = PeaksEta;
    for jj=1 : Problem.PeakNumber
        Problem.RotationMatrix{ii}(:,:,jj) = Problem.InitialRotationMatrix(:,:,jj) * Rotation(Problem.PeaksAngle(ii,jj),Problem.Dimension);
    end
    Problem.Environmentcounter = Problem.Environmentcounter + 1;
    [Problem.OptimumValue(ii), Problem.OptimumID(ii)] = max(PeaksHeight);
end
Problem.Environmentcounter = 1;
end

function output = Rotation(teta,Dimension)
counter = 0;
PageNumber = Dimension * ((Dimension-1)/2);
X = NaN(Dimension,Dimension,PageNumber);
for ii=1 : Dimension
    for jj=(ii+1) : Dimension
        if ii~=jj
            TmpMatrix = eye(Dimension);
            TmpMatrix(ii,ii) = cos(teta);
            TmpMatrix(jj,jj) = cos(teta);
            TmpMatrix(ii,jj) = sin(teta);
            TmpMatrix(jj,ii) = -sin(teta);
            counter = counter + 1;
            X(:,:,counter) = TmpMatrix;
        end
    end
end
output = eye(Dimension);
tmp = randperm(PageNumber);
for ii=1 : PageNumber
    output = output * X(:,:,tmp(ii));
end
end