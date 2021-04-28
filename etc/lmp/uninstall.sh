#!/bin/bash


loc=`pwd`

echo ""
echo "Removing the following folders and contents:"
echo "${loc}/build"
echo "${loc}/exe"
echo ""

rm -rf ./build ./exe

echo "Uninstall complete."
echo ""
