def concat_matfiles(exp_code,exp_ids,pathname,subject_names,exp_num=0):
    # Create a DataFrame from mat files by looping over subjects (subject_names)
    file_names = os.listdir(pathname)
    
    df = DataFrame()
    for exp_id in exp_ids:
        temp2 = DataFrame()
        for subject_name in subject_names:

            # Create a tag to find files based on a subject name
            if exp_num == 0:
                subject_tag = subject_name + str(exp_code) + '_' + str(exp_id) + '_'
            else:
                subject_tag = subject_name + str(exp_code) + '_' + str(exp_code) + '_' + str(exp_id) + '_'
            subject_files = fnmatch.filter(file_names, subject_tag + '*.mat') # Filter the files

            # Loop over filtered files
            temp1 = DataFrame()
            for subject_file in subject_files:

                # Convert a MATLAB matrix to a DataFrame
                temp0 = DataFrame(loadmat(pathname+subject_file)['datamatrix'])
                temp1 = pd.concat([temp1,temp0],ignore_index=True) # Concatenate DataFrames

            # Create a column to assign subject names
            temp1['Name'] = DataFrame([subject_name]*len(temp1))
            temp2 = pd.concat([temp2,temp1],ignore_index=True) # Concatenate

            # Create a column to assign experiment index
            if exp_id == 1:
                feat_name = 'OriA'
            elif exp_id == 2:
                feat_name = 'ColA'
            elif exp_id == 3:
                feat_name = 'OriC'
            elif exp_id == 4:
                feat_name = 'ColC'
            elif exp_id == 5:
                feat_name = 'OriD2'
            elif exp_id == 6:
                feat_name = 'ColD2'
            elif exp_code == 6 and exp_id == 1:
                feat_name = 'OriD1'
            elif exp_code == 6 and exp_id == 2:
                feat_name = 'ColD1'
            else:
                feat_name = 'B'                
            temp2['Feature'] = DataFrame([feat_name]*len(temp2))
                
        # Concatenate all subject data to create a single DataFrame            
        df = pd.concat([df,temp2],ignore_index=True)
    return df