DATE=`date +'%Y%m%d'`
ARCHIVE_NAME=ELSES${DATE}.tar.gz
git archive --format=tar.gz --prefix=ELSES${DATE}/ HEAD -o ${ARCHIVE_NAME}
scp ${ARCHIVE_NAME} k009914@kashiwa.issp.u-tokyo.ac.jp:~/
