%MACRO SimulateLcaDataset(  true_gamma_dataset=,
							true_rho_dataset=, 
							output_dataset_name=, 
							total_n= ); /*
/*  Version 3.0          December 3, 2010     
	By John Dziak, The Methodology Center, The Pennsylvania State University
	(c) 2010 The Pennsylvania State University  
	Thanks to YoungKyoung Min for her very valuable input.
	This SAS macro creates simulated categorical data suitable
	for testing out latent class analysis and related methods. 
	Input:  
	 -- true_gamma_dataset is the name of a dataset
		containing a single column of k numbers whose
		sum is one.  These will be treated as the 
		population proportions of k latent classes.
	 -- true_rho_dataset is the name of a numeric dataset
		with m rows and k columns.  rho(m,k) should give 
		the population probability of the "1" response 
		on item m for a member of class k.  Heuristically
		we think of "1" as meaning "yes" (vs. 2 for "no").   
	 -- output_dataset_name is the name which will be given
		to the simulated dataset.
	 -- total_n is the number of subjects in the simulated
		dataset (not the number in each class).
Limitations: Currently the model cannot have
	groups or customized restrictions.  The items must be
	dichotomous.  There is practically no error checking 
	for inappropriate input so the user is on the honor system.
*/
	PROC IML;  
		* Process input...	;
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
				Matches = 	Matches #
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
