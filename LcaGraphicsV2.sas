/*---------------------------------------------------------------------*
 |  The documentation and code below is supplied by the Methodology 
 |  Center, at The Pennsylvania State University.
 *---------------------------------------------------------------------
 |  SAS Macros for Latent Class Analysis Graphics, Version 2.
 |  For use as supplements to PROC LCA.
 |  (c) 2015 The Pennsylvania State University.
 *---------------------------------------------------------------------*/


%MACRO IdentificationPlot(SeedsDataset=);
/*-------------------------------------------------------------------- 
 | MACRO NAME  : %IdentificationPlot
 | SHORT DESC  : Plot fitted log-likelihoods for different starting 
 |               values for an LCA model, based on output data from 
 |               PROC LCA.  This can be used to help assess how well
 |               a model is identified (how plausible it is that the
 |               global maximum of the likelihood as a function of
 |               the data is well-defined and has been found.) 
 *-------------------------------------------------------------------- 
 | CREATED BY  : John DZIAK                                  (12/20/2010)
 *-------------------------------------------------------------------- 
 | 1. This macro requires as input the name of a dataset created by
 |    PROC LCA using the OUTSEEDS option.  For instance, it can be called
 |    like:     
 |     %IdentificationPlot(SeedsDataset=seeds1);
 *------------------------------------------------------------------
 | Copyright 2015 The Pennsylvania State University.
 |
 | This program is free software; you can redistribute it and/or
 | modify it under the terms of the GNU General Public License as
 | published by the Free Software Foundation; either version 2 of
 | the License, or (at your option) any later version.
 |
 | This program is distributed in the hope that it will be useful,
 | but WITHOUT ANY WARRANTY; without even the implied warranty of
 | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 | General Public License for more details.
 *------------------------------------------------------------------*/
	DATA LcaMacroTemp1;
		SET &SeedsDataset;
		LogLik = ROUND(log_likelihood,.01); 
		CALL SYMPUT("LowestLogLik",Loglik);
		CALL SYMPUT("NumSeeds",_N_);
	RUN; 
	PROC SORT DATA=LcaMacroTemp1; BY DESCENDING Loglik; RUN;
	DATA LcaMacroTemp2; SET LcaMacroTemp1; RUN;
	DATA LcaMacroTemp2;
		SET LcaMacroTemp2; 
		CALL SYMPUT("HighestLogLik",Loglik);
	RUN; 
	DATA LcaMacroTemp3;
		LowestLoglik = SYMGET("LowestLoglik");
		HighestLoglik = SYMGET("HighestLoglik");
		MinForGraph = FLOOR(LowestLoglik-.1);
		MaxForGraph = CEIL(HighestLogLik+.1);
		CALL SYMPUT("MinForGraph",MinForGraph);
		CALL SYMPUT("MaxForGraph",MaxForGraph);
	RUN;  
	PROC CHART DATA=LcaMacroTemp1;
		TITLE "Frequency distribution of log-likelihoods for multiple starting values";
		HBAR LogLik / DISCRETE ;
	RUN;
	DATA LcaMacroTemp1; SET LcaMacroTemp1;
		LogLik = PUT(LogLik,8.2);
	RUN;
	PROC GCHART DATA=LcaMacroTemp1;
		TITLE " ";
		NOTE "Frequency distribution of log-likelihoods for multiple starting values";
		HBAR LogLik / DISCRETE ;
	RUN; QUIT;
%MEND IdentificationPlot;



%MACRO OddsRatioPlot(ParamDataset =,
                    StdErrDataset =);
/*-------------------------------------------------------------------- 
 | MACRO NAME  : %OddsRatioPlot
 | SHORT DESC  : Plot the confidence intervals for odds ratios of
 |               class membership as a function of the covariates in
 |               an LCA model. 
 *-------------------------------------------------------------------- 
 | CREATED BY  : John DZIAK                                  (05/17/2010)
 *-------------------------------------------------------------------- 
 | This macro requires as input the names of two datasets created by
 | PROC LCA using the OUTPARAM and OUTSTDERR options.  For instance, 
 | it can be called like:     
 |     %OddsRatioPlot(ParamDataset=pars1,StdErrDataset =sds1);
 *------------------------------------------------------------------
 | 1. Thanks to Dr. Mariya Shiyko for suggesting including this plot.
 | 2. This macro is not for use with estimates created using the BINARY option.
 *------------------------------------------------------------------
 | Copyright 2010 The Pennsylvania State University.
 |
 | This program is free software; you can redistribute it and/or
 | modify it under the terms of the GNU General Public License as
 | published by the Free Software Foundation; either version 2 of
 | the License, or (at your option) any later version.
 |
 | This program is distributed in the hope that it will be useful,
 | but WITHOUT ANY WARRANTY; without even the implied warranty of
 | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 | General Public License for more details.
 *------------------------------------------------------------------*/
  /* Find names and number of variables ; */
	%LET NumVariables = 0;
    DATA LcaMacroTemp1;
        SET &ParamDataset;
        WHERE PARAM="BETA" AND VARIABLE ~="" AND GROUP=1; 
        CALL SYMPUT("VariableID",_N_);
        CALL SYMPUT(CAT("VariableName",_N_),Variable);
        CALL SYMPUT("NumVariables",_N_);
    RUN;
	%IF %EVAL(&NumVariables = 0) %THEN %DO;
		%PUT OddsRatioPlot only works for models with covariates.;
		PROC IML;
				PRINT "OddsRatioPlot only works for models with covariates.";
		RUN;
	%END;
	%ELSE %DO;
	    /* Find number of groups ; */
	    DATA LcaMacroTemp2;
	        SET &ParamDataset;
	        WHERE PARAM="BETA" AND VARIABLE = "&VariableName1";
	        CALL SYMPUT("NumGroups",_N_);
	    RUN;  
	    %DO GroupID = 1 %TO %EVAL(&NumGroups); 
		/* Read in point estimates */ 
        DATA LcaMacroTemp1;
            SET &ParamDataset;
            WHERE PARAM="BETA" AND VARIABLE ~="" AND GROUP=&GroupID; 
            CALL SYMPUT("VariableID",_N_);
            CALL SYMPUT(CAT("VariableName",_N_),Variable);
            CALL SYMPUT("NumVariables",_N_);
        RUN;
		PROC TRANSPOSE DATA=LcaMacroTemp1 OUT=LcaMacroPointEstimates;
            ID Variable;
        RUN;   
		DATA LcaMacroPointEstimates;
			SET LcaMacroPointEstimates;
            WHERE _NAME_ ~= "GROUP" & _NAME_ ~= "RESPCAT" & &VariableName1 ~= . ; 
            RowID = _N_; 
		RUN;
		/* Read in standard errors */
        DATA LcaMacroTemp2;
            SET &StdErrDataset;
            WHERE PARAM="BETA" AND VARIABLE ~="" AND GROUP=&GroupID;
        RUN;
        PROC TRANSPOSE DATA=LcaMacroTemp2 OUT=LcaMacroStandardErrors;
            ID Variable;
        RUN;			
		DATA LcaMacroStandardErrors;
			SET LcaMacroStandardErrors;
            WHERE _NAME_ ~= "GROUP" & _NAME_ ~= "RESPCAT" & &VariableName1 ~= . ; 
            RowID = _N_; 
		RUN;
		* Set up the file with information for the plot;
		DATA LcaMacroPointEstimatesAll; STOP; RUN;
		DATA LcaMacroStandardErrorsAll; STOP; RUN;
		DATA _NULL_; CALL SYMPUT("StdErrsAvailable",1); RUN;
		%DO VariableID = 1 %TO %EVAL(&NumVariables);
			DATA _NULL_;
				temp = CAT("VariableName",&VariableID);
				CALL SYMPUT("ThisVariable",SYMGET(TRIM(temp)));
			RUN; 
			DATA temp&VariableID;  
				SET LcaMacroPointEstimates;
                Estimate = &ThisVariable;
				Label = CAT(CAT("Class",SUBSTR(_Label_,LENGTH(_Label_),1)), ",&ThisVariable");
			RUN;    
            DATA LcaMacroPointEstimatesAll; 
				SET LcaMacroPointEstimatesAll temp&VariableID;
			RUN;     
			DATA temp&VariableID;
				SET LcaMacroStandardErrors;
                StdError = &ThisVariable; 
				Label = CAT(CAT("Class",SUBSTR(_Label_,LENGTH(_Label_),1)), ",&ThisVariable");
            RUN; 
            DATA LcaMacroStandardErrorsAll; 
				SET LcaMacroStandardErrorsAll temp&VariableID;
			RUN; 
			DATA LcaMacroStandardErrorsAll; 
				SET LcaMacroStandardErrorsAll;  
				IF  StdError = . THEN CALL SYMPUT("StdErrsAvailable",0);
			RUN; 
			PROC IML;
				USE LcaMacroStandardErrorsAll; 
				READ ALL INTO X;
				CLOSE LcaMacroStandardErrorsAll; 
				Y = NROW(X);
				IF (Y<1) THEN CALL SYMPUT("StdErrsAvailable","0");
			QUIT; 
		%END;
        DATA LcaMacroDoTheWork;
            MERGE LcaMacroPointEstimatesAll LcaMacroStandardErrorsAll;  
            BetaMid = Estimate;
			IF (StdError > 1000000) THEN StdError = 1000000;
			IF (StdError = .) THEN StdError = 1000000;
            BetaHigh = Estimate + 1.96*StdError;
            BetaLow = Estimate - 1.96*StdError;
            RatioMid = EXP(MIN(BetaMid,20));
			RatioHigh = EXP(MIN(BetaHigh,20));
            RatioLow = EXP(MIN(BetaLow,20)); 
        RUN; 
		DATA LcaMacroRatioLow;    SET LcaMacroDoTheWork;    Ratio = RatioLow;    RUN;
        DATA LcaMacroRatioMid;    SET LcaMacroDoTheWork;    Ratio = RatioMid;    RUN;
        DATA LcaMacroRatioHigh;    SET LcaMacroDoTheWork;    Ratio = RatioHigh;    RUN;
        DATA LcaMacroRatioAll; SET LcaMacroRatioLow LcaMacroRatioMid LcaMacroRatioHigh; RUN;
		* Do the plot; 
        AXIS1 LABEL=("") OFFSET=(10,10); 
        SYMBOL1 I=HILOCB C=RED VALUE=NONE WIDTH=3 MODE=INCLUDE; 
        SYMBOL2 C=GREEN VALUE=DOT INTERPOL=NONE WIDTH=10 MODE=INCLUDE; 
        AXIS2 LABEL=("Odds Ratio") ORDER=(.01 .1 1 10 100) LOGBASE=10 LOGSTYLE=EXPAND;
		%IF %EVAL(&StdErrsAvailable = 1) %THEN %DO;
	        PROC GPLOT DATA=LcaMacroRatioAll;
	            TITLE " ";   
	            NOTE m=(10pct,95pct) "95% Confidence Intervals for Odds Ratios";
				%IF %EVAL(&numGroups > 1) %THEN %DO;
					NOTE m=(10pct,90pct) "(For Group &GroupID)";
				%END;
				PLOT Ratio*Label=1 RatioMid*Label=2 / OVERLAY HAXIS=AXIS1 VAXIS=AXIS2 VREF=1 VMINOR=9 CVREF=BLUE ;
	        RUN; 
			QUIT;
		%END;
		%ELSE %DO;
			%PUT The plot could not be shown because standard errors could not be computed.;
			PROC IML;
				PRINT "The plot could not be shown because standard errors could not be computed.;";
			RUN;
			QUIT;
		%END;
		%END;
	%END;
%MEND OddsRatioPlot;

 
%MACRO ItemResponsePlot(ParamDataset=);
/*-------------------------------------------------------------------- 
 | MACRO NAME  : %ItemResponsePlot
 | SHORT DESC  : Plot the response probabilities for each class and item
 |               in an LCA model. 
 *-------------------------------------------------------------------- 
 | CREATED BY  : John DZIAK                                  (05/17/2010)
 *-------------------------------------------------------------------- 
 | This macro requires as input the names of a datasets created by
 | PROC LCA using the OUTPARAM option.  For instance, 
 | it can be called like:     
 |     %ItemResponsePlot(ParamDataset=pars1); 
 *------------------------------------------------------------------*
 | Copyright 2010 The Pennsylvania State University.
 |
 | This program is free software; you can redistribute it and/or
 | modify it under the terms of the GNU General Public License as
 | published by the Free Software Foundation; either version 2 of
 | the License, or (at your option) any later version.
 |
 | This program is distributed in the hope that it will be useful,
 | but WITHOUT ANY WARRANTY; without even the implied warranty of
 | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 | General Public License for more details.
 *------------------------------------------------------------------*/
    PROC IML; 
        USE &ParamDataset;
            READ ALL INTO Params;
        CLOSE &ParamDataset;
        Params[1,2] = -999;
        A = Params[<>,1];
        B = Params[<>,2];
        CALL SYMPUT("NumGroups",CAT(" ",A));
        CALL SYMPUT("NumResponseCategories",CAT(" ",B));
    QUIT; 
    * Check for measurement invariance (because PROC LCA output datasets
    don't explicitly indicate this and we don't want to inconvenience
    the user by asking);
    * So pretend to be doing a data analysis on the output dataset.
    * Is any variability in the rho estimate for category 1 accounted for by group?;
	ODS SELECT NONE;
	PROC GLM DATA=&ParamDataset;
        WHERE PARAM="RHO" & RESPCAT=1;
        CLASS GROUP VARIABLE;
        MODEL ESTLC1 = VARIABLE GROUP GROUP*VARIABLE;
  	  ODS OUTPUT ModelANOVA=LcaMacroTempAnova1;
    RUN; 
    ODS OUTPUT CLOSE;
	ODS SELECT ALL; 
    DATA LcaMacroTempAnova1; 
        SET LcaMacroTempAnova1;
        WHERE HypothesisType = 1 & Source = "GROUP";
    RUN;
    DATA LcaMacroTempAnova2;
        SET LcaMacroTempAnova1;
        OneGroupOnly = 1;
        IF SS > 0.0001 THEN OneGroupOnly = 0;
        CALL SYMPUT('OneGroupOnly',OneGroupOnly);
    RUN;   
	%PUT _ALL_;
    %IF %EVAL(&OneGroupOnly)>0 %THEN %LET NumGroups = 1;
     /* If there is measurement invariance, only plot the first group,
        since the others will be identical. */ 
    %DO RespCatNum=1 %TO &NumResponseCategories;
        %DO GroupNum=1 %TO &NumGroups;
            %DO RespCatNum = 1 %TO %EVAL(&NumResponseCategories);
                %LET VarNamesString=;
                DATA LcaMacroTemp1;
                    SET &ParamDataset;
                    WHERE ((Param="RHO") & (Group=&GroupNum) & RespCat=&RespCatNum); 
                    CALL SYMPUT(CAT('LcaVar',_N_),TRIM(Variable)); 
                    CALL SYMPUT('NumVars',_N_);   
                    CALL SYMPUT("VarNamesString",CAT(CAT(SYMGET("VarNamesString")," '"),CAT(TRIM(Variable),"' "))); 
                RUN;
                PROC TRANSPOSE DATA=&ParamDataset OUT=LcaMacroTemp2;RUN; 
                %LET EstNamesString=;
                %LET EstNamesStringB=;
                DATA LcaMacroTemp2;
                    SET LcaMacroTemp2;
                    WHERE ((_NAME_ ~= "GROUP") & (_NAME_ ~= "RESPCAT")); 
                    CALL SYMPUT('NumClasses',_N_);
                    CALL SYMPUT("EstNamesString",CAT(SYMGET("EstNamesString"),CAT(TRIM(_Name_)," "))); 
                    CALL SYMPUT("EstNamesStringB",CAT(SYMGET("EstNamesStringB"),CAT(TRIM(_Name_),"*Variable "))); 
                RUN; 
                AXIS1 ORDER=(&VarNamesString) LABEL=('');
                AXIS2 ORDER=(0.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0) LABEL=('');
                SYMBOL1 CI=BLUE CV=BLACK INTERPOL=JOIN WIDTH=1 VALUE='1';
                SYMBOL2 CI=GREEN CV=BLACK INTERPOL=JOIN WIDTH=1 VALUE='2';
                SYMBOL3 CI=RED  CV=BLACK INTERPOL=JOIN WIDTH=1 VALUE='3';
                SYMBOL4 CI=BIP CV=BLACK INTERPOL=JOIN WIDTH=1  VALUE='4';
                SYMBOL5 CI=CHARCOAL CV=BLACK INTERPOL=JOIN WIDTH=1 VALUE='5';
                SYMBOL6 CI=ORANGE CV=BLACK INTERPOL=JOIN WIDTH=1 VALUE='6';
                SYMBOL7 CI=BIPPK  CV=BLACK INTERPOL=JOIN WIDTH=1  VALUE='7';
                SYMBOL8 CI=BROWN CV=BLACK INTERPOL=JOIN WIDTH=1  VALUE='8';
                SYMBOL9 CI=MAGENTA CV=BLACK INTERPOL=JOIN WIDTH=1  VALUE='9';
				PROC GPLOT DATA=LcaMacroTemp1;
                    TITLE " ";   
                    PLOT &EstNamesStringB  / OVERLAY HAXIS=AXIS1 VAXIS=AXIS2 CVREF=LTGRAY;
                    NOTE m=(10pct,95pct)  "Probability of response=&RespCatNum as a function of class and item"; 
                RUN; QUIT; 
            %END;
        %END;
    %END;
%MEND ItemResponsePlot; 
