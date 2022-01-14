OPTIONS MPRINT MLOGIC;
%INCLUDE "C:\Users\atw14\Documents\LCA_Covariates_3Step_v10.sas";
LIBNAME sasf "C:\Users\atw14\Documents\";

PROC LCA DATA=sasf.simple_simulated_example   
         OUTPOST=Noninclusive_Post  ;
    NCLASS 3;
    ITEMS I1 I2 I3 I4 I5;
    CATEGORIES 2 2 2 2 2;  
    SEED 11111;
    NSTARTS 5;
    ID id X1 X2;
    RHO prior=1; 
RUN;

%LCA_Covariates_3Step(postprobs = Noninclusive_Post,
    id = id, 
    Covariates = X1 X2);
