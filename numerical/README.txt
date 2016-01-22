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


Example of commands
-------------------

Compare number of robots per species:

python compare_species_size.py \
  --alpha 1 \
  --beta 1 \
  --nrobots 10 \
  --nspecies 3 \
  --cme \
  --duration=10 \
  --save_epsilons data_1_1.bin

Compare species behavior by varying alpha:

python compare_species_behavior.py \
  --nrobots 10 \
  --nspecies 3 \
  --alpha 1 \
  --beta 1 \
  --duration=10 \
  --save_epsilons=data_alpha_1_1.bin \
  --cme \
  --sweep_type alpha

Compare species behavior by varying beta:

python compare_species_behavior.py \
  --nrobots 10 \
  --nspecies 3 \
  --alpha 1 \
  --beta 1 \
  --duration=10 \
  --save_epsilons=data_beta_1_1.bin \
  --cme \
  --sweep_type beta
