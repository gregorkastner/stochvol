# Source: David Winterbottom, https://codeinthehole.com/tips/tips-for-using-a-git-pre-commit-hook/
# Accessed on 2020-07-10
# Licence for this file: CC BY-NC-SA 4.0
# Changes have been made to adapt the code to our project.
# Please note that this file is not part of the CRAN package stochvol, only part of the development repository hosted on GitHub.

# To install this file as a local pre-commit hook in the repo execute
# chmod ug+x pre-commit && ln -s `pwd -P`/pre-commit.sh .git/hooks/pre-commit

# Stash staged files that aren't included in this commit
STASH_NAME="pre-commit-$(date +%s)"
git stash save -q --keep-index $STASH_NAME

# Test prospective commit
R -e 'quit(save = "no", status = length(devtools::check(vignettes = FALSE, cran = FALSE)$errors) > 0)'
RESULT=$?

# Restore stashed files in the workspace
git stash pop

# Return an appropriate value: 0 => everything's fine, other value => error
[ $RESULT -ne 0 ] && exit 1
exit 0

