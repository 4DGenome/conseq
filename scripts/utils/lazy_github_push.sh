#!/bin/bash


#==================================================================================================
# Created on: 2017-04-26
# Usage: ./lazy_github_push.sh
# Author: javier.quilez@crg.eu
# Goal: combines git add, commit and push
#==================================================================================================



# git add
if [[ -z "$2" ]]; then
	git add *
else
	git add $2
fi

# git commit
git commit -m "$1"

# git push
git push