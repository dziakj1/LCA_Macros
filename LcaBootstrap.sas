%MACRO LcaBootstrap(    null_outest=,
                        alt_outest=,
                        null_outparam=,
                        alt_outparam=,
                        n=,
                        num_bootstrap=99,
                        num_starts_for_null=20,
                        num_starts_for_alt=20,
                        cores=1
);
/*  Version 4            April 8, 2016     
    By John Dziak, The Methodology Center, The Pennsylvania State University
    (c) 2016 The Pennsylvania State University
 Dr. YoungKyoung Min and Dr. Shu Xu provided helpful feedback for 
 testing and improving this code.
 This SAS macro does a parametric bootstrap test for
 latent class analysis using PROC LCA. 
 This program is distributed in the hope that it will be useful,
 but without any warranty; without even the implied warranty of
 merchantability or fitness for a particular purpose.  
Input:
    -- null_outest and alt_outest are the OUTEST= datasets 
        from PROC LCA run on the null and alternative models 
        respectively.
    -- null_outparam and alt_outparam are the OUTPARAM= dataset 
        from PROC LCA run on the null and alternative models.
    -- n is the original sample size used for the LCA analyses.
    -- num_bootstrap is the number of bootstrap 
        replications, which should be at least 99 (although 
        999 is preferable).  
    -- num_starts_for_null is the number of starting values
        to fit under the null hypothesis in each bootstrap
        replication to help find the global maximum 
        of the likelihood under this hypothesis.  It should be at least
        10.  It must be at least 2.
    -- num_starts_for_alt is the number of starting values
        to fit under the alternative hypothesis in each bootstrap
        replication to help find the global maximum 
        of the likelihood under this hypothesis.  It should be 
        least 10, and no less than num_starts_for_null, in order
        to reduce the risk of getting a nonsensical negative log-
        likelihood ratio statistic due to suboptimal fitting of the
        alternative hypothesis.  It must be at least 2.
    -- cores is the number of processor cores to use if more 
        than one is available (defaults to 1)
Output: 
    Three datasets are created:
    -- CalculatedLrtForBootstrap contains the -2 log likelihood 
        test statistic on the original provided dataset.  
    -- LcaBootstraps contains the -2 log likelihood test 
        statistic for each bootstrap-simulated dataset.
    -- BootstrapResults contains the p-value for the test, which
        is based on how many of the numbers in LcaBootstraps 
        are greater than the numbers in CalculatedLrtForBootstrap.
        It also contains information on the average proportion of 
        starting values agreeing with the best fit value in the 
        bootstrap replications (a measure of how well-identified
        the bootstrap datasets were).
Limitations: Currently the model cannot have covariates,
    groups or customized restrictions.  The items must be
    dichotomous.  There is little error checking 
    for inappropriate input so the user is on the honor system.
    Prints a lot of repetitive stuff to the log.  The alternative 
      hypothesis must involve exactly one class more than the null.
Notes:
    A "parametric bootstrap" implies repeatedly simulating 
    new datasets under the estimated null hypothesis (k-class) 
    distribution and then testing the null hypothesis for 
    each simulated dataset to approximate the null hypothesis 
    distribution of the test statistic.  The result is a test 
    statistic in the CalculatedLRTforBootstrap dataset, a list 
    of the bootstrapped LRTs in the LcaBootstraps dataset, 
    and a p-value in the BootstrapPValue dataset.
    The macro needs to be fed the results from the fits under
    the null and alternative hypotheses.  The null hypothesis
    results in null_outparam are used both to generate the 
    bootstrapped datasets and to provide starting values 
    for fitting the null hypothesis model to the bootstrapped 
    datasets.  The alternative hypothesis results 
    in alt_outparam are needed to provide starting 
    values for fitting the alternative hypothesis
    model to the bootstrapped datasets.  The log-likelihoods
    from the null and alternative hypothesis results (in
    null_outest and alt_outest) are also used to calculate 
    the negative-two log-likelihood ratio statistic for the 
    original observed dataset, which is compared to the
    empirical distribution of the negative-two log-likelihood
    ratio statistics for the bootstrapped datasets in order
    to get a p-value.
    Note:  
    While performing the computations, this macro creates/
    overwrites and deletes some datasets in the Work library,
    including: 
       a_bootstrap_dataset, boot_outest, boot_outest_extra, 
       boot_outest_main, boot_seeds, boot_summary, boot_temp1, 
       bootstraplca_macro_info, bootstraplikelihoodsh0, 
       bootstraplikelihoodsh1, estim_gamma_dataset, 
       loglikelihood_extra, null_gamma_dataset, null_rho_dataset
*/
    /* We first include a nested macro, namely a version of SimulateLcaDataset. */
    %MACRO SimulateLcaDataset(  true_gamma_dataset=,
                                true_rho_dataset=, 
                                output_dataset_name=, 
                                total_n= );  * Dec 3, 2010 version;
            PROC IML;  
                * Process input...    ;
                USE &true_gamma_dataset;
                    READ ALL INTO true_gamma;
                CLOSE &true_gamma_dataset;
                USE &true_rho_dataset;
                    READ ALL INTO true_rho_yes;
                CLOSE &true_rho_dataset;
                total_n = &total_n;
                true_nclass = NCOL(true_rho_yes); 
                n_items = NROW(true_rho_yes); 
                ***************************************************** ;
                items = ('Item' + CHAR((1:n_items),3)); 
                                                * assumes < 1000 items;
                CALL CHANGE(items," ","0",0);
                true_class = J(&total_n,1,0);  
                DO this_person = 1 TO &total_n ;
                    temp = 0;
                    CALL RANDGEN(temp, 'TABLE', true_gamma);
                    true_class[this_person] = temp;
                END;
                responses = J(&total_n, n_items, 0);
                DO this_person = 1 TO &total_n ;
                    DO this_item = 1 TO n_items;
                        temp = 0;
                        this_persons_class = true_class[this_person];
                        CALL RANDGEN(temp, 'BERNOULLI', 
                                true_rho_yes[this_item,this_persons_class]); 
                             responses[this_person,this_item] = 2 - temp;     
                    END;
                END;
                ResponsePatterns = J(2**n_items, n_items, 0); 
                DO this_item = 1 TO n_items;  
                    temp = REPEAT(1, 2**(n_items-this_item), 1) // REPEAT(2, 2**(n_items-this_item), 1);
                    temp = REPEAT(temp, 2**(this_item-1),1);
                    ResponsePatterns[,this_item] = temp;
                END; 
                Count = J((2**n_items),1,0);
                DO pattern_index = 1 TO (2**n_items);
                    Matches = J(&total_n, 1, 1);
                    DO this_item = 1 TO n_items;
                        Matches =     Matches #
                                    (responses[,this_item]=ResponsePatterns[pattern_index,this_item]); 
                        Count[pattern_index] = Matches[+];
                    END;  
                END;  
                 responses_matrix =  ResponsePatterns || Count; 
                cols = items || 'Count'; 
                CREATE &output_dataset_name FROM responses_matrix [ COLNAME = cols ];
                    APPEND  FROM responses_matrix;
                        * The name of the append command can be misleading. ;
                        * It does not append to the data from the previous ;
                        * replication, which is replaced.  It just appends ;
                        * to the new empty dataset, i.e., writes the data. ;
                CLOSE &output_dataset_name;
            QUIT;
            DATA &output_dataset_name;
                SET &output_dataset_name;
                WHERE Count > 0; * omit non-occurring patterns;
            RUN;  
        %MEND;
    OPTIONS NONOTES;
    %LET LcaBootstrapDebug = 0;
    %LET LcaBootstrapIdentificationCrit = .01;
    * INITIAL PROCESSING:                                            ;
    PROC IML;
        CALL SYMPUT("bootstrap_macro_okay","y");
        USE &null_outparam; 
            READ ALL INTO null_outparam;
        CLOSE &null_outparam;  
        nr = NROW(null_outparam);
        nc = NCOL(null_outparam); 
        nItems = (nr-1)/2;
        nclass_null = nc-2;
        USE &alt_outparam; 
            READ ALL INTO alt_outparam;
        CLOSE &alt_outparam; 
        IF (nItems ^= (nr-1)/2) THEN DO;
            PRINT "ERROR IN MACRO: The input datasets are incompatible.";
            PRINT "Null #items:"  nItems;
            PRINT "Alt. #items:"  (nr-1)/2;
            PRINT "The number of items should be the same.";
            CALL SYMPUT("bootstrap_macro_okay","n");
        END;
        IF (nclass_null > 10) THEN DO;
            PRINT "ERROR IN MACRO: There are too many classes.";
            CALL SYMPUT("bootstrap_macro_okay","n");
        END;
        IF (nItems <2) THEN DO;
            PRINT "ERROR IN MACRO: There is not enough data.";
            CALL SYMPUT("bootstrap_macro_okay","n");
        END;
        IF (nItems > 20) THEN DO;
            PRINT "ERROR IN MACRO: There are too many items.";
            CALL SYMPUT("bootstrap_macro_okay","n");
        END;
        nr = NROW(alt_outparam);
        nc = NCOL(alt_outparam); 
        nclass_alt = nc-2;
        IF (nclass_alt ^= nclass_null + 1) THEN DO;
            PRINT "ERROR IN MACRO: The input datasets are incompatible.";
            CALL SYMPUT("bootstrap_macro_okay","n");
            IF (nclass_null > 0) THEN DO;
                PRINT "Null #classes:"  nclass_null;
            END;
            IF (nclass_null <= 0) THEN DO;
                PRINT("The null hypothesis parameters dataset seems to be empty or not from PROC LCA.");
            END;
            IF (nclass_alt > 0) THEN DO;
                PRINT "Alt. #classes:"  nclass_alt;
                PRINT "This should be one greater than the null #classes.";
            END;
            IF (nclass_alt <= 1) THEN DO;
                PRINT("The null hypothesis parameters dataset seems to be empty, or not from PROC LCA.");
            END;
            CALL SYMPUT("bootstrap_macro_okay","n");
        END;
        IF (%EVAL(&num_starts_for_null)<10) THEN DO;
            PRINT "ERROR IN MACRO: A higher value for num_starts_for_null is needed.";
            CALL SYMPUT("bootstrap_macro_okay","n");            
        END;
        IF (%EVAL(&num_starts_for_alt)<10) THEN DO;
            PRINT "ERROR IN MACRO: A higher value for num_starts_for_alt is needed.";
            CALL SYMPUT("bootstrap_macro_okay","n");
        END;
        IF (%EVAL(&num_starts_for_alt)<%EVAL(&num_starts_for_null)) THEN DO;
            PRINT "ERROR IN MACRO: A higher value for num_starts_for_alt is needed.";
            CALL SYMPUT("bootstrap_macro_okay","n");
        END;
        num_starts_for_null_minus_one = %EVAL(&num_starts_for_null-1);
         CREATE bootstraplca_macro_info VAR {nclass_null 
                                             nclass_alt 
                                             nItems
                                             num_starts_for_null_minus_one };
            APPEND;
         CLOSE bootstraplca_macro_info ;
    QUIT; 
    DATA bootstraplca_macro_info; SET bootstraplca_macro_info;
        CALL SYMPUT("nclass_null",PUT(nclass_null,4.));
        CALL SYMPUT("nclass_alt",PUT(nclass_alt,4.));
        CALL SYMPUT("num_starts_for_null_minus_one",PUT(num_starts_for_null_minus_one,6.));
		CALL SYMPUT("nItems",PUT(nItems,6.));
         IF nItems = 2 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 ");
            CALL SYMPUT("categories","2 2 ");
        END;
         IF nItems = 3 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 ");
            CALL SYMPUT("categories","2 2 2 ");
        END;
         IF nItems = 4 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 ");
            CALL SYMPUT("categories","2 2 2 2 ");
        END;
         IF nItems = 5 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 ");
            CALL SYMPUT("categories","2 2 2 2 2 ");
        END;
         IF nItems = 6 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 ");
        END;
         IF nItems = 7 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 ");
        END;
         IF nItems = 8 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 9 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 10 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 11 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 12 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 13 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 Item013");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 14 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 Item013 Item014");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 15 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 Item013 Item014 Item015");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 16 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 Item013 Item014 Item015 Item016 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 17 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 Item013 Item014 Item015 Item016 Item017 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 18 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 Item013 Item014 Item015 Item016 Item017 Item018 ");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;
         IF nItems = 19 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 Item013 Item014 Item015 Item016 Item017 Item018 Item019");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;
        IF nItems = 20 THEN DO; 
            CALL SYMPUT("items","Item001 Item002 Item003 Item004 Item005 Item006 Item007 Item008 Item009 Item010 Item011 Item012 Item013 Item014 Item015 Item016 Item017 Item018 Item019 Item020");
            CALL SYMPUT("categories","2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ");
        END;

    RUN; 
        %IF &bootstrap_macro_okay=n %THEN %ABORT; 
    * STEP ONE:                                              ;
    * First we organize the estimate datasets in a form      ;     
    * that is easy for the SimulateLcaDataset macro to read. ;      
    DATA estim_gamma_dataset;
        SET &null_outparam;
        WHERE PARAM="GAMMA";
        DROP PARAM GROUP VARIABLE RESPCAT;
    RUN;
    PROC TRANSPOSE DATA=estim_gamma_dataset 
                    OUT=estim_gamma_dataset;
    RUN;
    DATA null_gamma_dataset;
        SET estim_gamma_dataset;
        KEEP COL1;
    RUN;
    DATA null_rho_dataset;
        SET &null_outparam;
        WHERE PARAM="RHO" & RESPCAT=1;
        DROP PARAM GROUP VARIABLE RESPCAT;
    RUN;
    * STEP TWO:                                                 ; 
    * Initialize datasets to hold the lists of log-likelihoods  ;
    * which will be calculated for each bootstrap dataset.      ;
    DATA BootstrapLikelihoodsH0; boot = .; RUN;
    DATA BootstrapLikelihoodsH1; boot = .; RUN;

    * STEP THREE:                                               ;
    * Run the main loop in which we do the repeated resampling  ;
    * according to the parametric bootstrap idea                ;
    %DO bootIndex = 1 %TO &num_bootstrap;
		%PUT Doing bootstrap number &bootIndex;
        * We simulate a replication dataset under the null      ;
        * hypothesis assumptions.  ;
        %SimulateLcaDataset(  
                true_gamma_dataset=null_gamma_dataset,
                true_rho_dataset=null_rho_dataset, 
                output_dataset_name=a_bootstrap_dataset, 
                total_n=&n );
        *     -------    The Null Hypothesis Part -------        ;
        * Now we analyze the simulated dataset under H0 using the original fit as starting values.;
        PROC LCA DATA=a_bootstrap_dataset START=&null_outparam OUTEST=boot_outest_main NOPRINT;
           NCLASS &nclass_null ;
           ITEMS &items ;
           FREQ count;
           CATEGORIES &categories ;    
           ESTIMATION EM; 
           CORES &cores ;   
        RUN;     
        * Now we analyze the simulated dataset under H0 using random starting values.;
        PROC LCA DATA=a_bootstrap_dataset OUTEST=boot_outest_extra OUTSEEDS=boot_seeds NOPRINT;
           NCLASS &nclass_null ;
           ITEMS &items ;
           FREQ count ;
           CATEGORIES &categories ;    
           ESTIMATION EM; 
           CORES &cores;   
           SEED 123456;
           NSTARTS &num_starts_for_null_minus_one;
        RUN;     
        DATA loglikelihood_extra;
            SET boot_outest_extra;
            loglik_extra = log_likelihood;
            KEEP loglik_extra;
        RUN;
        DATA boot_temp1; 
            SET boot_outest_main ;
            seed = -1;
            KEEP log_likelihood seed;
        RUN;
        PROC SORT DATA=boot_seeds; BY seed; RUN;
        DATA boot_outest;
            MERGE boot_seeds boot_temp1; 
            BY seed;
        RUN;  
        PROC SORT DATA=boot_outest; BY DESCENDING log_likelihood; RUN;
        PROC IML;             
            USE boot_outest;
                READ ALL VAR { log_likelihood } ;
            CLOSE boot_outest; 
            best_log_likelihood = log_likelihood[1];  
            diff = log_likelihood - best_log_likelihood; 
            startvals_agree = (SUM(ABS(diff)< &LcaBootstrapIdentificationCrit)-1) / (NROW(log_likelihood)-1); 
            IF (NROW(log_likelihood)=1) THEN startvals_agree = .; 
            boot = &bootIndex ;
            CREATE boot_summary  VAR { boot best_log_likelihood startvals_agree};
                APPEND;
            CLOSE boot_summary; 
        RUN; 
        DATA BootstrapLikelihoodsH0; 
            * Record the H0 likelihood for this bootstrap iteration;
            SET BootstrapLikelihoodsH0 boot_summary ; 
        RUN;
        *   -------  The Alternative Hypothesis Part -------  ;
        * Now we analyze the simulated dataset under H1.; 
        PROC LCA DATA=a_bootstrap_dataset OUTEST=boot_outest_extra OUTSEEDS=boot_seeds NOPRINT;
           NCLASS &nclass_alt ;
           ITEMS &items ;
           FREQ count;
           CATEGORIES &categories ;    
           ESTIMATION EM; 
           CORES &cores;   
           SEED 123456;
           NSTARTS &num_starts_for_alt;
        RUN;     
        DATA loglikelihood_extra;
            SET boot_outest_extra;
            loglik_extra = log_likelihood;
            KEEP loglik_extra;
        RUN;
        DATA boot_temp1; 
            SET boot_outest_main ;
            seed = -1;
            KEEP log_likelihood seed;
        RUN;
        PROC SORT DATA=boot_seeds; BY seed; RUN;
        DATA boot_outest;
            MERGE boot_seeds boot_temp1; 
            BY seed;
        RUN;  
        PROC SORT DATA=boot_outest; BY DESCENDING log_likelihood; RUN;
        PROC IML;             
            USE boot_outest;
                READ ALL VAR { log_likelihood } ;
            CLOSE boot_outest; 
            best_log_likelihood = log_likelihood[1];  
            diff = log_likelihood - best_log_likelihood; 
            startvals_agree = (SUM(ABS(diff) < &LcaBootstrapIdentificationCrit)-1) / (NROW(log_likelihood)-1); 
            IF (NROW(log_likelihood)=1) THEN startvals_agree = .; 
            boot = &bootIndex;
            CREATE boot_summary  VAR { boot best_log_likelihood startvals_agree};
                APPEND;
            CLOSE boot_summary; 
        RUN; 
        DATA BootstrapLikelihoodsH1; 
            * Record the H1 likelihood for this bootstrap iteration;
            SET BootstrapLikelihoodsH1 boot_summary ; 
        RUN;
        *   -------  End of the main loop -------  ;
    %END;   
    * STEP FOUR:                                               ;
    * Now we calculate the test statistic for each of the      ; 
    * bootstrap replications we created.                       ;
    DATA BootstrapLikelihoodsH0;
        SET BootstrapLikelihoodsH0;
        loglikH0 = best_log_likelihood; 
        startvals_agree_H0 = startvals_agree;
        KEEP boot loglikH0 startvals_agree_H0;
    RUN;
    DATA BootstrapLikelihoodsH1;
        SET BootstrapLikelihoodsH1;
        loglikH1 = best_log_likelihood; 
        startvals_agree_H1 = startvals_agree;
        KEEP boot loglikH1 startvals_agree_H1;
    RUN; 
    DATA BootstrapLikelihoodsH0; SET BootstrapLikelihoodsH0 ; WHERE boot > 0; RUN;
    DATA BootstrapLikelihoodsH1; SET BootstrapLikelihoodsH1 ; WHERE boot > 0; RUN;
    DATA LcaBootstraps;
        MERGE BootstrapLikelihoodsH0 BootstrapLikelihoodsH1 ;
        BY boot;
    RUN;
    DATA LcaBootstraps; SET LcaBootstraps;
        loglikRatioTest = -2 * loglikH0 + 2 * loglikH1;
        * By multiplying by 2 we imitate the usual                 ;
        * G-squared statistic.                                      ;
    RUN;

* STEP FIVE:                                                     ;
*    Now we get the p-value, which is just a slightly adjusted   ;
*    version of the proportion of bootstrap replications having  ;
*    test statistics higher than the actually obtained            ;
*    test statistic.                                             ;
    OPTIONS NOTES;
DATA BootstrapResult; RUN;
PROC IML;
        USE &null_outest;
            READ ALL VAR {log_likelihood} INTO data_loglikH0;
        CLOSE &null_outest;
        USE &alt_outest;
            READ ALL VAR {log_likelihood} INTO data_loglikH1;
        CLOSE &alt_outest;
        data_loglikRatioTest = -2 * data_loglikH0 + 2 * data_loglikH1;
        CREATE CalculatedLRTforBootstrap FROM data_loglikRatioTest;
            APPEND FROM data_loglikRatioTest;
        CLOSE CalculatedLRTforBootstrap;
        USE LcaBootstraps;
            READ ALL VAR {loglikRatioTest} INTO boot_loglikRatioTests;
            READ ALL VAR {startvals_agree_H0} INTO boot_startvals_agree_H0;
            READ ALL VAR {startvals_agree_H1} INTO boot_startvals_agree_H1;
        CLOSE LcaBootstraps;
        * Now data_loglikRatioTest is a single number giving the
        * test statistic for the data we actually got.  
        * boot_loglikRatioTests is a vector of test statistics simulated
        * from a distribution in which the null hypothesis is true.
        * We reject the null hypothesis if and only if 
        * data_loglikRatioTest is larger than almost all of the 
        * numbers in boot_loglikRatioTests.; 
        IF ANY(boot_loglikRatioTests=.) THEN DO;
			PRINT("ERROR IN MACRO: The test statistic could not be calculated ");
			PRINT("in some bootstrap replication.");
			CALL SYMPUT("bootstrap_macro_okay","n");
        END; 
		IF ANY(boot_loglikRatioTests<0) THEN DO;
			PRINT("ERROR IN MACRO: Alternative model was estimated poorly in at least one bootstrap");
			PRINT("replication.  A larger value for num_starts_for_alt is recommended.");
			CALL SYMPUT("bootstrap_macro_okay","n");
        END; 
        count = boot_loglikRatioTests > data_loglikRatioTest;
        count2 = boot_loglikRatioTests <= data_loglikRatioTest; 
        %IF &bootstrap_macro_okay=n %THEN %DO;            
             PRINT("ERROR IN MACRO -- p-value cannot be calculated");
             %ABORT; 
        %END;
        IF (count[+]+count2[+] = &num_bootstrap) THEN DO;
            pvalue = (count[+] + 1) / (&num_bootstrap + 1);
            answer_from_bootstrap = "The bootstrap p-value was" + CHAR(pvalue);
            startvals_agree_H0 = boot_startvals_agree_H0[+] / &num_bootstrap ;
            startvals_agree_H1 = boot_startvals_agree_H1[+] / &num_bootstrap ;
            CREATE BootstrapResult VAR {pvalue startvals_agree_H0 startvals_agree_H1};
                APPEND; 
            CLOSE BootstrapResult;
            PRINT answer_from_bootstrap;
        END; ELSE DO;
            PRINT("ERROR IN MACRO -- p-value cannot be calculated");
            CALL SYMPUT("bootstrap_macro_okay","n");
        END;
     QUIT; 
* STEP SIX:                                               ;
*     Clean up by removing the temporary datasets that  ;
*      had been created.                                 ;
     %IF %EVAL(&LcaBootstrapdebug = 0) %THEN %DO;
         PROC DATASETS  NOLIST; 
             DELETE  a_bootstrap_dataset
                    bootstraplikelihoodsh0
                    bootstraplikelihoodsh1
                    bootstraplca_macro_info
                    boot_outest
                    boot_outest_extra
                    boot_outest_main
                    boot_seeds
                    boot_summary
                    boot_temp1
                    estim_gamma_dataset
                    loglikelihood_extra
                    null_gamma_dataset
                    null_rho_dataset;
        QUIT;      
    %END;
%MEND;
