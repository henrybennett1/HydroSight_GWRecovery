emis = [0.2,0.2,0.2,0.2,0.2,0.2;0.1,0.1,0.1,0.1,0.1,0.1];
tprobs = [0.2,0.8;0.1,0.9];
[seq,states] = hmmgenerate(100,tprobs,emis);
states2 = hmmviterbi(seq,tprobs,emis);

accuracy = sum(states == states2);

plot(1:length(states),states,'-')