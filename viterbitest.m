%Viterbi Decoding Test
%Viterbi takes a dataset which has error, and then 

%Trellis diagram is the transition matrix? or at least it functions similarly
codein = ;                          % 
trellis = poly2trellis();           % 
tbdepth = 2;                        % Generally set to 2 or 3, but check with Tim
opmode = 'trunc';                   % 'cont' preserves values between instances, and 'term' assumes that the state ends in the same state it starts, all opmodes assume it starts in state 0, whichever that is
dectype = 'soft';                   % 'hard' only works if inputs are 0 or 1, and 'unquant' works when inputs are + or -
nsdec = 3;                          % Integer between 1 and 13 used for 'soft' dectype to improve error decoding recovery?

decodeout = vitdec(codein,trellis,tbdepth,opmode,dectype);