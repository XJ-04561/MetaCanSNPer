
## CURRENT

1. SNP naming scheme changed to 

## QOL

 * Mutect2 needs a proper --tmp-dir set or it will throw warnings

 * Should all directory lookups follow the same hierarchy?
     * Reference Hierarchy: Company database > User homebrew database > Current working database
     * Proposed hierarchy: [Install Directory] > [Working Directory] > [Home/User Directory] > [Absolute Path]
     * Should they choose only one of these and use as absolute reference only as a backup?
      * Should absolute pathing be the first attempt?
 * One program takes all queries or just one at a time?

## DUMP

 * Minimap2
    * @RG PL: [CAPILLARY, DNBSEQ (MGI/BGI), ELEMENT, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT (Oxford Nanopore), PACBIO (Pacific Bio-sciences), SINGULAR, SOLID, and ULTIMA]
