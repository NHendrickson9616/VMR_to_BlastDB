#!/usr/bin/env awk
#
# summarize mismatch-by-genus table in one line
#
BEGIN {
    # mismatch-by-gneus columns
    KEY=1;
    GENUS=2;
    ERR=3;
    CT=4
    
    # pre-populate possible error list
    # must match code in "extract_results.awk"
    split("no_hits,mismatch,mismatch_seqsame,match_segdiff,match_noQseg,match_noSseg",err_list,",")
    for(i in err_list) {
	key = err_list[i]
	#print "SETUP ["key"]"
	cts[key]=0
    }
}
# each line has an error class and a count
{
    sum=sum+$CT;
    cts[$ERR]=cts[$ERR]+$CT
}
END{
    # one line summary
    printf "total mismatches:"sum",";
    for(key in cts) { 
	printf key":"cts[key]","
    };
    print ""
}
