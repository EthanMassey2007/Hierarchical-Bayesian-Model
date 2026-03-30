#!/bin/bash

DATE="$1"
MESSAGE="$2"

GIT_AUTHOR_DATE="$DATE" \
GIT_COMMITTER_DATE="$DATE" \
git commit -m "$MESSAGE"
