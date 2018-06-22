#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#
# Unstage several files which might confuse git and which don't neccessarily need to be added to the repo everytime.
files="$(basename $PWD).jar buildversion.txt src/de/dkfz/b080/co/*Plugin.java"

git reset HEAD $files
git checkout -- $files

