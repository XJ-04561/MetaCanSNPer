

 * Should all directory lookups follow the same hierarchy?
     * Proposed hierarchy: [Install Directory] > [Working Directory] > [Home/User Directory] > [Absolute Path]
     * Should they choose only one of these and use as absolute reference only as a backup?
      * Should absolute pathing be the first attempt?

 * Make sure that SNP data or pre-mapped/aligned data can be given.
 * One program takes all references or just one at a time?
 * Make combined method for aligning and mapping
 * How should references be handled as a hierarchy of actions?
    * Example:
        1. Check whether it exists already on the disk.
        2. Download from internet-based servers.
            a. Remove or keep the downloaded references based on argument flag.
