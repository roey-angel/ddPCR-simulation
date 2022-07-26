# Simulating the bias in ddPCR caused by DNA fagmentation and multiple gene copies

## Usage:
Model_ddPCR_underestimation_V2.sh will run the simulaiton while Run_model_ddPCR_underestimation_V2.sh is a wrapper for running the simulation n times.
When done, run Summarise_ddPCR_underestimation.sh to create a summary table.

## Output header description

| Column               | Description                                                             | 
|----------------------|-------------------------------------------------------------------------|
|# Genomes             |Number of analysed genomes                                               |${nGENOMES}
|# Frags per genome    |Number of fragments generated per genome                                 |${nFRAG}
|Frags w. SSU          |Number of total fragments with one or more SSU genes in them             |${nHITs}
|P frags w. SSU        |Fraction of Frags w. SSU / # Genomes * # Frags per genome                |${pHITs}
|Frags w. PCR          |Number of total fragments with one or more PCR hits in them              |${nPCRHits}
|P frags w. PCR        |Fraction of total fragments with one or more PCR hits in them            |${pPCRHits}
|Frags w. multi SSU    |Number of total fragments with multiple SSU genes in them                |${nMultHits}
|P frags w. multi SSU  |Fraction of total fragments with multiple SSU genes in them              |${pMultHits}
|Frags w. multi PCR    |Number of total fragments with multiple PCR hits in them                 |${nMultPCRHits}
|P frags w. multi PCR  |Fraction of total fragments with multiple PCR hits in them               |${pMultPCRHits} 
|P SSU frags with PCR  |Fraction of PCR hits in total number of fragments with SSU               |${pFragWPCR}
|P multi PCR           |Fraction of multiple PCR hits in total number of fragments with a PCR hit|${pFragWMultPCR}
|Corrected ddPCR value |Number of true SSU genes in all fragments from all genomes               |${cddPCR}
