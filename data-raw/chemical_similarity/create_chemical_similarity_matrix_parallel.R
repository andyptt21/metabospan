options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library(rcdk)
library(tidyverse)
library(parallel)

detectBatchCPUs <- function() { 
    ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
    if (is.na(ncores)) { 
        ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
    } 
    if (is.na(ncores)) { 
        return(4) # for helix
    } 
    return(ncores) 
}

ncpus <- detectBatchCPUs() 
options(mc.cores = ncpus) 

## Do NOT parallelize get.fingerprints

ramp_structures<-readRDS("metabHasStructure.Rds") %>%
    rename(smiles=iso_smiles)

##Adding information

start=Sys.time()
atom_containers<-sapply(ramp_structures$smiles,function(x){
    out<-parse.smiles(x)
    print(match(x,ramp_structures$smiles))
    ##Sys.sleep(.5)
    return(out)
})
print(Sys.time()-start)

print("atom containers generated")
##saveRDS(atom_containers,"atom_containers_interactive.Rds")

##atom_containers=readRDS("atom_containers_interactive.Rds")
counter<-0
fingerprints = lapply(atom_containers, function(mol){
    counter<<-counter+1
    print(counter)
    ret_val = NULL
    if(length(mol) > 0){
        tryCatch(
        {ret_val = get.fingerprint(mol, type = "pubchem")},
        error=function(cond) {
            message(cond)
            # Choose a return value in case of error
            return(NA)
        }
        )
    }
    return(ret_val)
})

##saveRDS(fingerprints, file = "fingerprints_interactive.rds")
print("fingerprints generated")

ramp_structures$fingerprint = fingerprints
saveRDS(ramp_structures, file = "all_ramp_structures.rds")

ramp_structures<-readRDS("all_ramp_structures.rds")

compute_tanimoto_sim = function(dat,threshold=0.9){

    #Read in file and get fingerprint bit vectors.
    fingerprint_vecs = lapply(dat$fingerprint, function(fprint){as.numeric(unlist(strsplit(as.character(fprint),"")))})
    
    #Open output files.
    ## namefile = file(n, "w")
    ## weightfile = file(w, "w")
    
    #For each molecule, compute its similarity to all other molecules
    #and write to a file.
    weights_ids = mclapply(seq(1,length(fingerprint_vecs)), function(mol){
    cat(paste0(mol,"... "))
        #Compute all similarities to this molecule.
    sims = sapply(seq(1,length(fingerprint_vecs)),
                  function(mol2){
            sim = 0
            ##Simliarity of a molecule to itself is set to 0 so it will be filtered out.
            if (dat$ramp[[mol]] != dat$ramp[[mol2]]){
               sim = tanimoto(fingerprint_vecs[[mol]], fingerprint_vecs[[mol2]]) 
            }
            return(sim)
        })
            
        #Retain only similarities greater than the threshold.
        weights_keep = unlist(sims[which(sims >= threshold)])
        ids_keep = dat$ramp[which(sims >= threshold)]
        
        #Return ids and weights.
        return(data.frame(as.data.frame(weights_keep),
                          as.data.frame(ids_keep,
                                        stringsAsFactors = FALSE)
                          ))
    },mc.cores=ncpus
    )
    
    names(weights_ids)<-dat$ramp

    ## out_matrix<-weights_ids[[1]][,1]
    
    ## for(i in 2:length(weights_ids)){
    ##     out_matrix<-cbind(out_matrix,weights_ids[[i]][,1])
    ## }

    ## colnames(out_matrix)=rownames(out_matrix)=weights_ids[[1]][,2]
    out_matrix<-bind_rows(weights_ids,
                          .id = "ref_id")
    return(out_matrix)
}

tanimoto = function(x, y){
    intersection = length(which((x * y) == 1))
    union1 = length(which(x == 1)) + length(which(y == 1)) - intersection
    return(intersection / union1)
}

start=Sys.time()
chemical_similarity_matrix<-compute_tanimoto_sim(ramp_structures,0.9)
print(Sys.time()-start)


saveRDS(chemical_similarity_matrix,"RaMP_chemical_similarity_matrix.Rds")
