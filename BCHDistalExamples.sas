
/*------------------------------------------------------------------*
| The documentation and code below is supplied by The Methodology Center, Penn State University.             
*-------------------------------------------------------------------*/
																				  



%INCLUDE "C:\SAS\project files\LCA_distal_BCH_v110.sas";  /* The %INCLUDE statement is necessary for SAS Macros. The file path must match where the files exist ON YOUR MACHINE. */

/* LIBNAME sasdata "C:\SAS\project files\"; */ /* Identify the directory where your datafiles exist in this file path. */



 
 

    
    
PROC LCA DATA = sasdata.SimData_Binary OUTPARAM = Binary_param OUTPOST = Binary_post;
  ID id;
  NCLASS 5; 
  ITEMS item001-item008; 
  CATEGORIES 2 2 2 2 2 2 2 2;  
  SEED 12345;    
  RHO PRIOR = 1;  
  NSTARTS 20;
  MAXITER 5000;          
  CRITERION 0.000001; 
RUN;

    
    
    
%LCA_Distal_BCH(input_data = sasdata.SimData_Binary,  
                param = Binary_param,
                post = Binary_post, 
                id = id,
                distal = z, 
                metric = binary ); 

                      
    

PROC LCA DATA = sasdata.SimData_conti OUTPARAM = conti_param OUTPOST = conti_post ;
    ID id;
    NCLASS 5; 
    ITEMS item001-item008; 
    CATEGORIES 2 2 2 2 2 2 2 2;  
    SEED 12345;    
    RHO PRIOR = 1; 
    NSTARTS 20;
MAXITER 5000;          
CRITERION 0.000001; 
RUN;




%LCA_Distal_BCH(input_data = sasdata.SimData_conti, 
            param = conti_param,
            post = conti_post,
            id = id,
            distal = z,
            metric = Continuous );




PROC LCA DATA = sasdata.SimData_Count OUTPARAM = Count_param OUTPOST = Count_post;
    ID id;
    NCLASS 5; 
    ITEMS item001-item008; 
    CATEGORIES 2 2 2 2 2 2 2 2; 
    SEED 12345;    
    RHO PRIOR = 1; 
    NSTARTS 20;
MAXITER 5000;  
CRITERION 0.000001; 
RUN;



%LCA_Distal_BCH(input_data = sasdata.SimData_Count, 
            param = Count_param,
			post = Count_post,
			id=id,
			distal = z,  
            metric = Count );





PROC LCA DATA = sasdata.SimData_Categ OUTPARAM = Categ_param OUTPOST = Categ_post;
    ID id;
    NCLASS 5; 
    ITEMS item001-item008; 
    CATEGORIES 2 2 2 2 2 2 2 2; 
    SEED 12345;    
    RHO PRIOR = 1; 
    NSTARTS 20;
    MAXITER 5000;          
    CRITERION 0.000001; 
    RUN;

    
    
    
%LCA_Distal_BCH(input_data = sasdata.SimData_Categ,  
        param = Categ_param,
        post = Categ_post, 
        id = id,
        distal = z, 
        metric = categorical );


                  
