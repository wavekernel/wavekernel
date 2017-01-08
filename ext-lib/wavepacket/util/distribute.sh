DATE=`date +'%Y%m%d'`
ARCHIVE_NAME=WavePacket_${DATE}.zip
git archive HEAD -o ${ARCHIVE_NAME}
