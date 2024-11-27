function [fcDistance] = fc_distance(fc1,fc2) 

%[fcDistance] = fc_distance(fc1,fc2)
%Inputs: fc1 and fc2 are empirical and simulated FCs of size N x N.
%Output: fcd= Euclidean distance between the two FCs. A scaler value. 

%[fcDistance] = fc_distance(fc1,fc2,lesionAreas)
% computes the euclidean distance between given two FC matrices. If lesionAreas is not empty
% then those areas will not be used while calculating the distance
% fc1(lesionAreas,:) = 0;
% fc1(:,lesionAreas) = 0;
% fc2(lesionAreas,:) = 0;
% fc2(:,lesionAreas) = 0;
% fcDistance = sqrt(mean(mean((fc1 - fc2).^2)));
% fcDistance = sum(sum(abs(fc1 - fc2)));


fcDistance = sqrt(sum(sum((fc1 - fc2).^2)));
