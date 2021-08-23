function ICc = IC12(fp)
%  this function assigns initial conditions for the ODE solver. 

% 1	N
% 2	Sm0
% 3	Eu0
% 4	Em0
% 5	Iu0
% 6	Im0
% 7	Au0
% 8	Am0
% 9	Ih0
% 10	Iicu0
% 11	R
% 12	D


% refer to Equation 1 in the manuscript. 
% Su =   N    -  Sm       Eu      Em      Iu      Im      Au      Am      Ih      Iicu    R       D                                                                        
Su =  fp(1) - (fp(2) +  fp(3)+  fp(4)+  fp(5)+  fp(6)+  fp(7)+  fp(8)+  fp(9)+  fp(10)+ fp(11)+ fp(12));

% Initial conditions. 
% IC = (su  sm      Eu      em       Iu      Im        Au      Am      Ih      Iicu    R       D      )
ICc = [Su   fp(2)    fp(3)    fp(4)    fp(5)    fp(6)    fp(7)    fp(8)    fp(9)    fp(10)    fp(11)    fp(12)  ];

end



