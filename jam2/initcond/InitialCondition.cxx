#include <jam2/initcond/InitialCondition.h>

namespace jam2 {

int InitialCondition::findAZ(string nucl)
{
    int mass = atoi(nucl.c_str());

    string element[109]={"0",
       "H",  "He", "Li" ,"Be", "B",  "C",  "N",  "O",
       "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",
       "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr",
       "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
       "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
       "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
       "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",
       "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
       "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
       "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
       "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
       "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",
       "Bk", "Cf", "Es", "Fm", "Md", "No", "Lw", "Ku",
       "105","106","107","108"};

    string pcode[33]={"e-",  "e+",  "nu_e",  "nu_e~",
      "mu-",  "mu+",  "nu_mu",  "nu_mu~",  "tau-",
      "tau+",  "nu_tau",  "nu_tau~",  "pi+",  "pi-",
      "n0",  "n0bar",  "p+",  "pbar-",  "Gamma",
      "Lambda0",  "Sigma-",  "Sigma0",  "Sigma+",  "Xi-",
      "Xi0",  "Omega-",  "pi0", "Reggeon", "Pomeron",
      "K-",  "K+", "K0", "K~0"};

     int kcode[33]={11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,
       211,-211,2112,-2112,2212,-2212,22,3122,3112,3212,3222,
       3312,3322,3334,111,28,29,-321,321,311,-311};

     int pdgID=-1;

    int iz=0, ip=-1;
    if(mass>1) {
        string nuc=nucl;
        while(atoi(nuc.c_str())>0) {
            nuc.erase(nuc.begin());
        }
          for(int i=1;i<=108;i++) {
            if(nuc==element[i]) {iz=i; pdgID=0; break;}
          }

          /*
          for(int i=1;i<=108;i++) {
            int j=nucl.find(element[i]);
            if(j>0) { iz=i; break;}
          }
          */

	pdgID=1000000000 + mass*10 + iz*10000;

	//cout << " A= " << mass << " iz= " << iz << " id= " << pdgID <<endl;

    } else {
        mass=1;
        for(int i=0;i<33;i++) {
            //int j=nucl.find(pcode[i]);
            //if(j>=0) { ip=i; break;}
            if(nucl==pcode[i]) {ip=i; pdgID=kcode[i]; break;}
        }
        int i1= pcode[ip].find("+");
        int i2=pcode[ip].find("-");
	        int i3=pcode[ip].find("0");
        if(i1 > 0) iz=1;
        if(i2 > 0) iz=-1;
        if(i3 > 0) iz=0;
    }


        //atomic = mass;
        //charge = iz;

      //cout << "A= " << mass << " z= " << iz << endl;
      //if(ip>=0) cout << "beam = " << pcode[ip] << endl;

    return pdgID;
}

} // end namespace jam2
