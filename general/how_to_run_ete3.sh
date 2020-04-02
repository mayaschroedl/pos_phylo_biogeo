########
# How to run the brnln_outliers.py script from https://github.com/mossmatters/phyloscripts/blob/master/brlenoutliers/brlen_outliers.py
# or any other ete3 script on a Server
########

# Install Minconda  (you can ignore this step if you already have Anaconda/Miniconda)
wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/anaconda_ete/
export PATH=~/anaconda_ete/bin:$PATH;

#To install this package with conda run:
conda install -c etetoolkit ete3

#sometimes ete3 doesn't install pyqt correctly --> do it manually
conda install pyqt


#if run on server without display:
#install virtual display like this: https://kovyrin.net/2007/10/01/how-to-run-gui-programs-on-a-server-without-any-monitor/
apt-get install xvfb
Xvfb -shmem -screen 0 1280x1024x24
DISPLAY=:0 xdpyinfo
xvfb:2:respawn:/usr/bin/Xvfb :0 -ac -screen 0 2048x1536x24


#activate display
display=$(shuf -i 100-200 -n 1)
export DISPLAY=:${display}
Xvfb :${display} -screen 0 1024x768x16 > /dev/null 2>&1 &
echo "export DISPLAY=:${display}" > ~/.xvfb


##run script
python brlen_outliers.py tree