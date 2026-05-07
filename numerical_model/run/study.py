# -*- coding: utf-8 -*-

#reading the study
#this is where the users has to define the study to use, the runs to compute and enter a few information
with open(preprod_folder+'study.csv', 'r') as f:
    reader = csv.reader(f,  delimiter=';')
    for row in reader:
        if row[0]=='to_run':
            to_run=row[1]
        elif row[0]=='study':
            study=row[1]
        else:
            len_param_list_prev=len(row)
            for l in list(range(1,len(row[1:])+1)):
                if len_param_list_prev==len(row) and row[l]=='':
                    len_param_list_prev=l
            exec(row[0]+'=row[1:len_param_list_prev]')
