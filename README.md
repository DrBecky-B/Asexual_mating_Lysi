Data and code for 'Nevertheless sex persisted: facultative sex is common but costly in a parthenogenetic wasp' (working title).
This read me provides details about the headers in each of the files in this repository:

File name: Lysi_asex_mate-vFeb-24.R: this is the R code used to run the analyses and make the figures in the manuscript

File name: G1.mating.rates.txt: this is the .txt file which contains data about how many female wasps mated and remained virgin in the first generation (G0). Date: date the mating observations took place. 
Line: asexual line (female wasp). Mated: how many female wasps accepted a mating for that line. Virgin: how many female wasps were not offered a mate (remained virgin) for that line. 
Refuse: how many females were offered a male but kicked him off without mating. No.attempt: how many females were offered a male but the male did not attempt to mate. 

File name: G1-raw-check.txt: this is the .txt file that contains all the fitness data for the first generation of females (G0). Number: female ID number. 
Set up date: date that the female emerged and was set up in her mating treatment. Block: block number (which day from 1-3 the female emerged on and was set up on). 
Line: which line the female came from. Mated/Virgin: whether the female was mated or virgin. N nymphs: how many aphid nymphs were on the leaf disc when the wasp was transferred.
Survival (hours): how many hours (up to a maximum of 72) did the female survive for on the leaf disc (after 72 hours wasps were transferred with aphids to the plant and survival could not be checked).
Failed to parasitise: did the female fail to parasitise any aphids? 1 = failed, 0 = successful. Failed to produce offspring: did the female fail to produce any adult offspring? 1 = failed, 0 = successful.
Mating failure: for mated sexual females only (all others NA) - did mated sexual females only produce male offspring
Mummies: number of mummies counted. Wasps: how many adult offspring counted. Unemerged: estimate - how many mummies are unemerged? (calculated from mummies-wasps). 
Males: number of male offspring. Females: number of female offspring. Exclude (lost or killed mother): should data from this wasp be removed because the mother was lost or killed during handling. 
SA: for genotyped broods only - did the brood contain any sexual alleles? Alleles: for genotyped broods only - which alleles were recovered for this brood (asexual only 184.184; sexual only 192.192 or asexual and sexual 184.192)

File name: G2-raw-check.txt: this is the .txt file that contains all the fitness data for the second generation of females (G1). Number G1: female ID number of mother. Number G2: female ID number for this female.
Set up date: date that the female emerged and was given hosts. Block 1: block number of mother (which day from 1-3 the female emerged on and was set up on). Block 2: block number of this female (which day from 1-3 the female emerged on and was set up on).
Line: which line the female came from. Mated/Virgin: whether the female was mated or virgin. N nymphs: how many aphid nymphs were on the leaf disc when the wasp was transferred.
Survival (hours): how many hours (up to a maximum of 72) did the female survive for on the leaf disc (after 72 hours wasps were transferred with aphids to the plant and survival could not be checked).
Failed to parasitise: did the female fail to parasitise any aphids? 1 = failed, 0 = successful. Failed to produce offspring: did the female fail to produce any adult offspring? 1 = failed, 0 = successful.
Mating failure: for mated sexual females only (all others NA) - did mated sexual females only produce male offspring
Mummies: number of mummies counted. Wasps: how many adult offspring counted. Unemerged: estimate - how many mummies are unemerged? (calculated from mummies-wasps). 
Exclude (no G2): should this row be removed because the G0 mother didn't produce any offspring. Exclude (lost or killed mother): should data from this wasp be removed because the mother was lost or killed during handling. 
SA: for genotyped broods only - did the brood contain any sexual alleles? Alleles: for genotyped broods only - which alleles were recovered for this brood (asexual only 184.184; sexual only 192.192 or asexual and sexual 184.192)
