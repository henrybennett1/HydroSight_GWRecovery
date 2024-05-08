%V1
%Emission probabilities = 1 if no obs value
%Emission probabilities are the observed likelihoods
data = ;
alpha1 = zeros(size(data,1),1);            %Potentially setting size of Alpha

alpha1(1) = obj.parameters.probabilities.initial_1 * emission.probs1(1);       %DEFINING ALPHA FOR FIRST LOOP BASED ON INITIAL PROBABILITIES
scalefactor = sum(alpha1);              %Alternatives include max(alpha1) and mean(alpha1)                                                
lscale1 = log(scalefactor);
alpha1(1) = alpha1(1) / scalefactor;        %alpha(1) is reset to 1

for i = 2:size(data,1)
    alpha1(i) = alpha1 .* obj.parameters.probabilities.trans_1 * emission.probs1(i);        %Dont know where emission probabilities are stored/found?  emission.probs1 is a placeholder for syntax
    scalefactor = sum(alpha1);
    lscale1 = lscale1 + log(scalefactor);
    alpha1(i) = alpha1(i) / scalefactor; 
end

%Log Likelihood is returned as lscale. The origin of this formula is
%Hydrostate, but it also comes from Page 49 in Zucchini textbook

%V2 (Same as V1 with changed number)
