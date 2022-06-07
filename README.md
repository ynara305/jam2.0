# JAM2: a Monte-Carlo event generator for high-energy nuclear collisions.

This code simulates high energy nuclear collisions based on the transport theoretical models. 
Models include hadronic cascade model, relativistic quantum molecular dynamics (RQMD.RMF, RQMDs, RQMDv), and hydrodynamics.
Hybrid simulation of hydrodynamics + cascade model is available.

[](
https://drive.google.com/file/d/1vN1UdJFpUwivv0fhIMB4NSSYierVmM2b/pythia-8.244mod.tar.gz
)

First download modified pythia8.244 (pythia8.244mod2.tar.gz) from
https://drive.google.com/file/d/10ZKyqbM0tE-Ca6EUCRCgvHniImAchK0v/view?usp=sharing

Then Install Pythia8 to your preferred directory:

./configure --prefix=$HOME/lib/pythia8

make install

Then go to JAM source directory and do

autoreconf -i

./configure Pythia8=$HOME/lib/pythia8

make

