#
# parse hit lines grep'ed from BLAST (fmt=7)
#
# output looks like:
# Aquamavirus	MH760796	MH760796#Aquamavirus_A	EU142040.?#Aquamavirus_A	73.049	4267	1001	129	62	4253	59	4251	0.0	1371
# Aquareovirus	AF450318	# BLAST processed 1 queries
# Aquareovirus	AF450319	AF450319#Aquareovirus_A#Seg3	HM989932.?#Aquareovirus_E#L3	81.570	1172	214	2	12	1182	1486	2656	0.0	966
# ...
BEGIN {
    FS="\t"
    OFS="\t"
    
    # BLAST OUTPUT fmt=7 vs VMR_E database
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

    # fields in ICTV's fasta headers
    SEQ_ACC=1
    SEQ_SPECIES=2
    SEQ_SEG=3
}
{ 
    if(NF==3) { 
        print "no_hits\t"$GENUS"\t"$ACCESSION;
        next;
    };
    if($QUERY==$SUBJECT) {
	# everything (access#species#seg) matches
        print "match",$0
    } else {
	# extract species_name and segment_name (split by #)
	split($QUERY,q_arr,"#")
	split($SUBJECT,s_arr,"#")
	print "# QUERY=" q_arr[SEQ_SPECIES] " [" q_arr[SEQ_SEG] "]; SUBJECT=" s_arr[SEQ_SPECIES] " [" s_arr[SEQ_SEG]"]"
	if( q_arr[SEQ_SPECIES] == s_arr[SEQ_SPECIES] ) {
	    # species names match
	    if( q_arr[SEQ_SEG] == "" && s_arr[SEQ_SEG] == "" ) {
		# species names are the same, neither has a segment
		print "match",$0
	    } else if( q_arr[SEQ_SEG] == s_arr[SEQ_SEG] ) {
		# species names are the same, and segment names a the same
		print "match_seg",$0
	    } else if( q_arr[SEQ_SEG] == "" && s_arr[SEQ_SEG] != "" ) {
		# species names are the same, but only QUERY does NOT have a segment name
		print "match_noQseg",$0
	    } else if( q_arr[SEQ_SEG] != "" && s_arr[SEQ_SEG] == "" ) {
		# species names are the same, but only SUBJECT does NOT have a segment name
		print "match_noSseg",$0
	    } else {
		# species names are the same, but segment names aren't
		print "match_segdiff",$0
	    }
	} else if( q_arr[SEQ_SEG] == "" && s_arr[SEQ_SEG] == "" && s_arr[SEQ_SEG] != "" ) {
	    # ONLY segment names match, species name does not match
	    print "mismatch_segsame",$0
	} else {
	    # nothing matches
	    print "mismatch",$0
	}
    }
}

