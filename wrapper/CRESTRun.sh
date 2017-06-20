#!/bin/sh

EXT_CLIP=/srv/public/kjk/noname/bin/crest/extractSClip.pl
CREST_BIN=/srv/public/kjk/noname/bin/crest/CREST.pl
BLAT_PORT=$1
BAM_FILE=$2
GENOME_FA=$3
GENOME_2BIT=$4
OUT_DIR=$5

# move
cd /srv/public/kjk/noname/bin/crest

echo "[start] `date`"
start_time=`date +%s`

# clip ext.
$EXT_CLIP -i $BAM_FILE --ref_genome $GENOME_FA -o $OUT_DIR

echo "[end-clip] `date`"

# CREST
$CREST_BIN -f $OUT_DIR/out.sorted.bam.cover -d $BAM_FILE --ref_genome $GENOME_FA  -t $GENOME_2BIT --cap3 /srv/public/kjk/noname/bin/CAP3/cap3 --blatclient /srv/public/kjk/noname/bin/blat/bin/gfClient --blatserver localhost --blatport $BLAT_PORT --blat /srv/public/kjk/noname/bin/blat/bin/blat -o $OUT_DIR --min_sclip_reads 1

echo "[end] `date`"
end_time=`date +%s`
runtime=$((end_time-start_time))
echo $runtime

exit 0
