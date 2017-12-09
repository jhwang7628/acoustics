sudo apt-get install ssh -y
cd /home/jui-hsien/code/acoustics
git remote set-url origin git@bitbucket.org:jw969/acoustics.git
ssh-keygen
echo -e "\n" 
cat ~/.ssh/id_rsa.pub
echo -e "\n"
echo "copy the key to bitbucket and hit [ENTER]"
read 
git pull
cd build_release 
tmux new -s make-fdtd -d "make fdtd-acoustic-simulator-viewer -j32"
