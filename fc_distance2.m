function [fcDistance,fcc] = fc_distance2(fc1,fc2,lesionAreas) %fc_distance(fc1,fc2)
%[fcDistance] = fc_distance(fc1,fc2,lesionAreas)
% computes the euclidean distance between given two FC matrices. If lesionAreas is not empty
% then those areas will not be used while calculating the distance
fc1(lesionAreas,:) = 0;
fc1(:,lesionAreas) = 0;
fc2(lesionAreas,:) = 0;
fc2(:,lesionAreas) = 0;
% fcDistance = sqrt(mean(mean((fc1 - fc2).^2)));
fcDistance = sqrt(sum(sum((fc1 - fc2).^2)));
fccc = corrcoef(fc1, fc2);
fcc=fccc(1,2);
% fcDistance = sum(sum(abs(fc1 - fc2)));
end