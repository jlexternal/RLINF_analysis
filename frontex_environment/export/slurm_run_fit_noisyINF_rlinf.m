taskID = getenv('SLURM_ARRAY_TASK_ID');
iter_num = str2num(taskID);
run_fit_noisyINF_rlinf(iter_num);
