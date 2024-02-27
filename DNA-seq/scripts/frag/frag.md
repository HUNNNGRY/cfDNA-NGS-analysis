
# get bam frag
get-FragmentSizeBam.sh {sample}

# get ratio
get-FragmentSizeDepth.sh {sample}

## join to matrix
get-FragmentSizeMat.R

## (test) 
frag_size.R
