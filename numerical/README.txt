Installation instructions
-------------------------

1) Install StochKit

cd /path/to/this/directory
wget http://downloads.sourceforge.net/project/stochkit/StochKit2/StochKit2.0.11/StochKit2.0.11.tgz
tar xvfz StochKit2.0.11.tgz
cd StochKit2.0.11
sh install.sh

sudo mkdir /usr/local/lib
sudo cp libs/boost_1_53_0/stage/lib/libboost_* /usr/local/lib/


2) CMEpy : Package to compute the CME.

cd /path/to/this/directory
git clone git://github.com/fcostin/cmepy.git
cd cmepy
python setup.py install
