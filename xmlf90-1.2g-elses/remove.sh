#!/bin/sh

echo "==> Cleaning strings"
(cd strings ; make clean)
echo "==> Cleaning sax"
(cd sax ; make clean)
echo "==> Cleaning wxml"
(cd wxml ; make clean)
echo "==> Cleaning cml"
(cd cml ; make clean)
echo "==> Cleaning xpath"
(cd xpath ; make clean)
echo "==> Cleaning dom"
(cd dom ; make clean)
echo "----------------"
echo "==> Cleaning SAX examples"
(cd Examples/sax ; sh remove.sh)
echo "==> Cleaning DOM examples"
(cd Examples/dom ; make clean)
echo "==> Cleaning XPATH examples"
(cd Examples/xpath ; make clean)
echo "==> Cleaning WXML examples"
(cd Examples/wxml ; make clean)
echo "==> Cleaning CML examples"
(cd Examples/cml ; make clean)
echo "==> Cleaning SAX Tutorial"
(cd Tutorial/sax ; make clean)
echo "==> Cleaning XPATH Tutorial"
(cd Tutorial/xpath ; make clean)
echo "==> Cleaning $FLIB_ROOT/modules:"
rm -f $FLIB_ROOT/modules/*
echo "==> Cleaning $FLIB_ROOT/lib/libflib.a:"
rm -f $FLIB_ROOT/lib/libflib.a


