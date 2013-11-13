DATE=`date +'%Y%m%d'`
ARCHIVE_NAME=eigen_test_${DATE}.tar.gz
git archive --format=tar.gz HEAD -o ${ARCHIVE_NAME} eigen_test/
