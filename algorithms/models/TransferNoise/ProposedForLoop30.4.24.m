%V1
%Emission probabilities are set to 1 if no obs data
data = ;
alpha1 = zeros(size(data,1),1);            %Potentially setting size of Alpha
alpha1(1) = obj.parameters.probabilities.initial_1 * emission.probs1(1);       %DEFINING ALPHA FOR FIRST LOOP BASED ON INITIAL PROBABILITIES
sumalpha1 = sum(alpha1);                                                   
lscale1 = log(sumalpha1);
alpha1(1) = alpha1(1) / sumalpha1;        %alpha(1) is reset to 1

for i = 2:size(data,1)
    alpha1(i) = alpha1 .* obj.parameters.probabilities.trans_1 * emission.probs1(i);        %Dont know where emission probabilities are stored/found?  emission.probs1 is a placeholder for syntax
    sumalpha1 = sum(alpha1);
    lscale1 = lscale1 + log(sumalpha1);
    alpha1(i) = alpha1(i) / sum(alpha1); 
end

%V2 (Same as V1 with changed number)
