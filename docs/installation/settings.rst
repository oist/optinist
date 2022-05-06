# GitHub SSH access settings
**You only need to do the following once.**

**CAUTION for Windows Users** : If you don't have git installed, you need to install  [git-for-windows](https://git-scm.com/download/win) or, git command in WSL2.

On Terminal or Git Bash or WSL2 terminal,
```
git config --global user.name "user.name"
git config --global user.email "user@oist.jp"
git config --global core.quotepath false
```
On your shell,
```
git config --global user.name "user.name"
git config --global user.email "user@oist.jp"
git config --global core.quotepath false
```
make your ssh key
```
mkdir ~/.ssh
cd ~/.ssh
ssh-keygen -t rsa
#Enter
#Input passphrase (empty for no passphrase):
#Input passphrase (again)
```
Open the generated public key (rsa.pub) with a text editor and copy all the contents.
Go to GitHub and follow the steps below to register your public key.

1. Log in to GitHub and select `Settings` from the menu on the top right
2. Select `SSH and GPG keys`
3. Press `New SSH Key`
4. Enter the `Title` (optional) and `Key` (paste the copied content) and press `Add SSH` key.

This completes the SSH connection settings!
