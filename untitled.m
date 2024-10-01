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
breakpoint = [breakpoint, length(years)];

tempStruct = struct;

for k = 1:(length(breakpoint)-1)
    % tempStruct(k).year = boreDataWL(breakpoint(k):breakpoint(k+1),1);
    % tempStruct(k).month = boreDataWL(breakpoint(k):breakpoint(k+1),2);
    % tempStruct(k).day = boreDataWL(breakpoint(k):breakpoint(k+1),3);
    % tempStruct(k).idk = boreDataWL(breakpoint(k):breakpoint(k+1),4);
    tempStruct(k).value = boreDataWL(breakpoint(k):breakpoint(k+1),5);
end