```
# sbr2smc.aligned.bam and it's corresponding .bam file will be generated
asts -q subreads.bam -t smc.bam -p sbr2smc.aligned
```

# chagelog

## 0.9.1

* reporter bug fix

## 0.9.0

* reporter
* gsmm2(0.20 -> 0.22)
* identity_gap_compressed instead of identity_without_long_indel(10)

## 0.8.0

* using hpc
* default kmer/wins and fallback kmer/wins


## 0.7.1

* length filter when add sbr to SubreadsAndSmc
* log file

## 0.4.0

* dw,cr,ar and qual from sbr.bam will dump the result bam
* worked for smc reads that length < 200 && >200

## 0.3.0

* add ch tag to bam record

## 0.2.1

* worked

## 0.1.1

## 0.1.0
