#!/bin/sh

echo "==> Building sax"
(cd sax ; make)
echo "==> Building wxml"
(cd wxml ; make)
echo "==> Building cml"
(cd cml ; make)
echo "==> Building xpath"
(cd xpath ; make)
#
# *** Comment out strings and dom if you do not
#     have a f95 compiler
#     (and complain to your vendor!!)
#
echo "==> Building strings"
(cd strings ; make)
echo "==> Building dom"
(cd dom ; make)
echo "-----------------------------------"
echo "==> Building SAX examples"
(cd Examples/sax ; sh build.sh)
echo "==> Building XPATH examples"
(cd Examples/xpath ; make)
echo "==> Building WXML examples"
(cd Examples/wxml ; make)
echo "==> Building CML examples"
(cd Examples/cml ; make)
#
# *** Comment out the dom examples if
#     you do not have a f95 compiler
#     (and complain to your vendor!!)
#
echo "==> Building DOM examples"
(cd Examples/dom ; make)
#
echo "==> Building SAX Tutorial"
(cd Tutorial/sax ; make)
echo "==> Building XPATH Tutorial"
(cd Tutorial/xpath ; make)


