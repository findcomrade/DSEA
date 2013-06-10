collect_drugs <- function(start, end){
  drugs <- c()
  viability <- c()
  genders <- c()
  diagnosis <- c()
  
  samples <- rep(sample_names, end-start+1)
  for(it in start:end){
    drugs <- c(drugs, rep(drug_names[it], 151))
    viability <- c(viability, cell_viability[,it])  
    diagnosis <- c(diagnosis, as.character(pd$Diagnosis))
    genders <- c(genders, as.character(pd$Gender))
  }  
  return (data.frame(sample=samples, drug=drugs, viability=viability, diagnosis=diagnosis, gender=genders))
}

samples_by_drug <- function(dname, viability_matrix){
  samples <- row.names(viability_matrix)
  drugs <- rep(dname, length(samples))
  viability <- viability_matrix[,dname]
  diagnosis <- c()
  for (sm in 1:length(samples)){
    diagnosis <- c(diagnosis, as.character(pd$Diagnosis[pd$Specimen.ID == as.character(samples[sm])]))
  }
  # genders <- as.character(pd$Gender)

  return (data.frame(sample=samples, drug=drugs, viability=viability, diagnosis=diagnosis))
}

filter_cell_viability <- function(diagnosis){
  filtered <- matrix(ncol=length(drug_names))
  for (sample in 1:nrow(pd)){
    if (pd$Diagnosis[sample] == diagnosis) {
      filtered <- rbind( filtered, as.matrix(cell_viability[as.character(pd$Specimen.ID[sample]),]) ) 
    }
  }
  filtered <- filtered[2:dim(filtered)[1],]
  return (filtered)
}

samples_by_diagnosis <- function(diagnosis){
  
  return(1)
}