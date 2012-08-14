#!/bin/bash
LOCALVER=`git rev-list HEAD | wc -l`
if [ $LOCALVER \> 1 ] ; then
    VER=`git rev-list origin/master | wc -l`
    VER_DIFF=$(($LOCALVER-$VER))
    if [ $VER_DIFF != 0 ] ; then
        VER="$VER+$VER_DIFF"
    fi
    if git status | grep -q "modified:" ; then
        VER="${VER}M"
    fi
	VER=`echo ${VER} | sed 's/\ //'`
	echo "Revision=${VER}"
    echo "Creating fluidReg.hpp."
    OUT="#ifndef VERSION_H_\n #define VERSION_H_ \n\n #define REVISION \"${VER}\" \n\n #endif" 
    echo -e $OUT > ./src/fluidReg.hpp
fi
