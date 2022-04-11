I should compile this into a pipeline/single script but for now here's instructions on how to build sysdata.Rda, which contains all the precomputed data structures necessary for running MetaboSPAN.

1. Run _extract\_ramp\_data.r_. This generates metabHasPathway and metabHasStructure Rds files which are stored in their respective subdirectories 
2. Run scripts in subdirectories (_pathway\_similarity_ and _chemical\_similarity_). It's best to submit these jobs as they take most of the day to run. Ontologies are currently deprecated
3. **Optional:** _create\_ramp\_knowledge\_graph_ is used soley to make the fused pathway+chemical knowledge graph. It's currently not in use by the package._
4. Run _create\_sysData.R_ to convert intermediate files into adjacency matrices and compile into sysData.Rda
