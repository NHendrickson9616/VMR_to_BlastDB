#
# parse hit lines grep'ed from BLAST (fmt=7)
#
# output looks like:
# Aquamavirus	MH760796	Aquamavirus_A	Aquamavirus_A	73.049	4267	1001	129	62	4253	59	4251	0.0	1371
# Aquareovirus	AF450318	# BLAST processed 1 queries
# Aquareovirus	AF450319	Aquareovirus_A#Seg3	Aquareovirus_E#L3	81.570	1172	214	2	12	1182	1486	2656	0.0	966
# ...
BEGIN {
    FS="\t"
    OFS="\t"
    GENUS=1
    ACCESSION=2
    QUERY=3
    SUBJECT=4
    IDENTITY=5
    ALIGN_LEN=6
    MISMATCHES=7
    GAPS=8
    Q_START=9
    Q_END=10
    S_START=11
    S_END=12
    EVALUE=13
    BITSCORE=14
}
{ 
    if(NF==3) { 
        print "no_hits\t"$GENUS"\t"$ACCESSION;
        next;
    };
    if($QUERY==$SUBJECT) {
	# everything (species#seg) matches
        print "match",$0
    } else {
	# split by SEGMENT (#name)
	split($QUERY,q_arr,"#")
	split($SUBJECT,s_arr,"#")
	if( q_arr[1] == s_arr[1] ) {
	    # species names match
	    if( q_arr[2] == "" || s_arr[2] == "" ) {
		# species names are the same, but one does not HAVE a segment name
		print "match_noseg",$0
	    } else {
		# species names are the same, but segment names aren't
		print "match_badseg",$0
	    }
	} else {
	    # nothing matches
	    print "mismatch",$0
	}
    }
}

