
## CURRENT

1. Determine SNP in VCF
2. Feed VCF data into the newick tree.

## QOL

 * Mutect2 needs a proper --tmp-dir set or it will throw warnings

 * Cache the hash of query
 * Should all directory lookups follow the same hierarchy?
     * Reference Hierarchy: Company database > User homebrew database > Current working database
     * Proposed hierarchy: [Install Directory] > [Working Directory] > [Home/User Directory] > [Absolute Path]
     * Should they choose only one of these and use as absolute reference only as a backup?
      * Should absolute pathing be the first attempt?
 * One program takes all references or just one at a time?
 * "Download" client separate from CanSNPer

## STRUCTURE

* modules/
*   CanSNPer2.py
*   DirectoryLibrary.py
*   DownloadReferences.py
*   Wrappers.py
*   Aligners.py
*   Mappers.py
*   SNPCallers.py
*   ParseXMFA2.py
*   LogKeeper.py
*   ErrorFixes.py
*   DatabaseConnection.py
* defaultFlags.toml
* Databases/
*   francisella_tularensis.db
* References/
* 
* 

## DUMP

 * Minimap2
    * @RG PL: [CAPILLARY, DNBSEQ (MGI/BGI), ELEMENT, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT (Oxford Nanopore), PACBIO (Pacific Bio-sciences), SINGULAR, SOLID, and ULTIMA]
