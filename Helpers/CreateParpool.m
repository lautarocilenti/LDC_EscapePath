function [] = CreateParpool()
%CREATEPARPOOL 
    if CheckIfCluster()
        if isempty(gcp('nocreate'))
            pc = parcluster('local');
            nproc = str2num(getenv('SLURM_CPUS_PER_TASK'));
            pc.NumWorkers = nproc;
            job_folder = fullfile(getenv('TMPDIR'),getenv('SLURM_JOB_ID'));
            mkdir(job_folder);
            pc.JobStorageLocation = job_folder;
            parpool(pc,nproc)
        end

    end
end

function [result] = CheckIfCluster()
    if (contains(cd,'lcilenti@umd.edu'))
        result = 1;
    else
        result = 0;
    end
end