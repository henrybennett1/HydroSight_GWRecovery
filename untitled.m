close all; clear; clc;
load("downloads\cleaned\46854_boreData.mat");

temp1 = boreDataWL(1:168,:);
temp2 = boreDataWL(169:336,:);
temp3 = boreDataWL(337:504,:);
plot(temp1(:,5));
hold on;
plot(temp2(:,5));
plot(temp3(:,5));
hold off;

years = boreDataWL(:,1);

breakpoint = 1;
j = 2;
for i = 2:length(years)
    if years(i) < years(i-1)
        breakpoint(j) = i;
        j = j + 1;
    end
end

for k = 1:(length(breakpoint)-1)
    temp = boreDataWL(breakpoint(k):breakpoint(k+1),:);
end