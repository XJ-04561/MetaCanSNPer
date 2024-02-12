
[CURRENT]

1. Implement minimap2 and GATK4

[QOL]

 * Cache the hash of query
 * Should all directory lookups follow the same hierarchy?
     * Reference Hierarchy: Company database > User homebrew database > Current working database
     * Proposed hierarchy: [Install Directory] > [Working Directory] > [Home/User Directory] > [Absolute Path]
     * Should they choose only one of these and use as absolute reference only as a backup?
      * Should absolute pathing be the first attempt?
 * One program takes all references or just one at a time?
 * "Download" client separate from CanSNPer