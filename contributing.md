# Contributing

Follow these steps for adding a feature or fixing a bug. This will update the
master branch on github (origin master)

1. pull the master branch from github
2. checkout the devel branch.  If no devel branch exists, create iterature
3. do all work on the devel branch
## updating the master branch
4. run check on devel to make sure everything is working
5. checkout the master branch.  
6. pull changes from github
7. checkout the devel branch again
8. rebase the 


```
git checkout master
git checkout devel
## do work on devel, making frequent commits

# look in git log to find the commit hashes
git cherry-pick oldest-new-commit..newest-new-commit 
git svn info # should just give repository information
git svn rebase # pulls down changes from svn repository, may be conflicts
git svn dcommit # "pushes" changes to svn
```
