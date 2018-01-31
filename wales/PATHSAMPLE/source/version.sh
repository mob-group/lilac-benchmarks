#!/bin/bash
# SVN VERSION SCRIPT
# To be run on compilation - returns the current version of the code.
# This script should be called from within the Makefile, and the results passed to the preprocessor
# For example:
# DEFS+=-DSVNVERSION="`./version.sh`"

tmpfile=`mktemp`

if [ $# -gt 1 ]; then
  outfile=$1
else
  outfile=SVNREV
fi

# First, we try to run svnversion to get the version number
svnversion . | sed 's/.*://' | sed 's/M//' > $tmpfile
# svnversion returns 'exported' if you have it installed but run it in a non svn directory.
# If svnversion is not installed at all, the file will be empty. In either case, we want to
# use the version number recorded previously. 
if [ "`cat $tmpfile`" != "exported" ] && [ "`cat $tmpfile`" != "" ]; then 
# If we are working with svn however, we want to overwrite the version number with the current one.
    cp $tmpfile $outfile 
else
    #If we're not using svn then we might be using git for the repository.
    #Try getting the git hash key
    #git rev-parse --short HEAD > $tmpfile
    git describe --always --long --dirty > $tmpfile
    #if this fails it will set $? to something other than 0
    #An additional check is that if it fails it will contain the word "fatal"
    stat="$?"
    fatal=`grep -i "fatal" $tmpfile`
    if [ "$stat" -eq 0 -a -z "$fatal" ]; then
      #git hash worked, so save it in file $outfile
      cp $tmpfile $outfile 
    fi
fi
# Remove the temporary file
rm $tmpfile
# Return the version from the VERSION file
cat $outfile
