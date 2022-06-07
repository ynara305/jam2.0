#include <jam2/xsection/HadronContent.h>
#include <cmath>

namespace jam2 {

using namespace std;

//***********************************************************************

//void HadronContent::attflv(int kf,int& ifla,int& iflb)
void HadronContent::findFlavor(int kf,int& ifla,int& iflb)
{
//...Purpse: to give spin and quarkflavour to the ends of
//...the excited strings.
//...For mesons, the order of the end flavors is randomly given;
//...For baryons, where a quark-diquark combination,
//...the diquark is always assigned to iflb.
//...Use SU(6) weight proton=1/3d(uu)1 + 1/6u(ud)1 + 1/2u(ud)0
//...                 nurtron=1/3u(dd)1 + 1/6d(ud)1 + 1/2d(ud)0

      const double vfr13=0.167;
      const double vfr14=0.333;
      const double vfr15=0.5;

//...for proton
//...vfr13 probability of finding a diquark ud with spin 1 1/6
//...vfr14 probability of finding a diquark uu with spin 1 1/3
//...vfr15 probability of finding a diquark ud with spin 0 1/2
//...54321
//...    x quark
//... xx0x di-quark
//...   1x leptons
//...   xx gauge and Higgs bosons
//...  xxx meson
//...x0xxx meson
//... xxxx baryon
//...xxxxx baryon

    int j = abs(kf);
    int  jc1=j         % 10;      // spin
    int  jc2=(j/10)    % 10;
    int  jc3=(j/100)   % 10;
    int  jc4=(j/1000)  % 10;
    //int  jc5=(j/10000) % 10;

    ifla=0;
    iflb=0;
    if(j < 100) return;

//...Di-quark
    if(jc2 == 0) return;

    int id1a=jc4*100+jc3*10+jc2;
    int  id1= kf > 0 ?  id1a : -id1a;


//..Mesons
    if(jc4 == 0) {
//...Identify the quark and antiquark in mesons:
/*
	int j100=jc3;
        int j10=jc2;
        int isgn = pow(-1, max(j100, j10));
        if(kf < 0) isgn = -isgn ;
        if(isgn > 0) j10 = -j10;
        if(isgn < 0) j100 = -j100;
*/

	int j100, j10;
	int sgn = kf > 0 ? 1 : -1;
	if(jc3%2 == 0) {
	    j100=jc3*sgn;
	    j10=-jc2*sgn;
	} else {
	    j100=jc2*sgn;
	    j10=-jc3*sgn;
	}

        if(id1 == 11 || id1 == 22) {   // pi0 eta rho0 omega...
	    double ranx=rndm->flat();
	    if(ranx < 0.5) {
		j100=1;
		j10=-1;
	     } else {
		j100=2;
		j10=-2;
	     }
	}

        if(j == 331) {  // eta'  (1/6 1/6 2/3)
          double ranx=rndm->flat();
          if(ranx < 0.1667) {
             j100=1;
             j10=-1;
	  } else if(ranx < 0.3333) {
             j100=2;
             j10=-2;
	  } else {
             j100=3;
             j10=-3;
	  }
	}

        //if(spin < 0.5) {
        if(j10 < 0) {
          ifla=j100;
          iflb=j10;
	} else {
         ifla=j10;
         iflb=j100;
	}

//...Baryons
    } else if(jc4 != 0) {

        int j1000=jc4;
        int j100=jc3;
        int j10=jc2;
        if(kf < 0) {
          j1000=  -j1000;
          j100 = -j100;
          j10 = -j10;
	}
        double spin=rndm->flat();

//... spin 1/2 baryons
        if(jc1 == 2) {

        if(spin < vfr13) {
          ifla=j1000;
          iflb=ifrkfc(j100,j10,0,1.0);
	} else if(spin < vfr13+vfr14) {
          ifla=j10;
          iflb=ifrkfc(j1000,j100,0,1.0);
	} else if(spin < vfr13+vfr14+vfr15) {
          ifla=j100;
          double s=0.0;
          if(j1000 == j10) s=1.0;
          iflb=ifrkfc(j1000,j10,0,s);
	}

//...Certain Lambda-like hadrons have two lightest quarks in spin-0:
        if(abs(j100) < abs(j10)) {
        if(spin < vfr13) {
          ifla=j1000;
          double s=0.0;
          if(j100 == j10) s=1.0;
          iflb=ifrkfc(j100,j10,0,s);
	} else if(spin < vfr13+vfr14) {
          ifla=j10;
          iflb=ifrkfc(j1000,j100,0,1.0);
	} else if(spin < vfr13+vfr14+vfr15) {
          ifla=j100;
          iflb=ifrkfc(j1000,j10,0,1.0);
	}
	}

//...spin 3/2 baryons
	} else {
      if(j1000 == j100 && j100 == j10) {
         ifla=j1000;
         iflb=ifrkfc(j100,j10,0,1.0);
      } else if(j1000 == j100 && j100 != j10) {
         if(spin < 0.3333) {
           ifla=j10;
           iflb=ifrkfc(j1000,j100,0,1.0);
	 } else {
           ifla=j100;
           iflb=ifrkfc(j1000,j10,0,1.0);
	 }
      } else if(j1000 == j10 && j10 != j100) {
         if(spin < 0.3333) {
           ifla=j100;
           iflb=ifrkfc(j1000,j10,0,1.0);
	 } else {
           ifla=j10;
           iflb=ifrkfc(j1000,j100,0,1.0);
	 }
      } else {
         if(spin < 0.3333) {
           ifla=j1000;
           iflb=ifrkfc(j100,j10,0,1.0);
	 } else if(spin < 0.6667) {
           ifla=j100;
           iflb=ifrkfc(j1000,j10,0,1.0);
	 } else {
           ifla=j10;
           iflb=ifrkfc(j1000,j100,0,1.0);
	 }
      }
	}

	// in case of anti-baryon, qq-bar + qbar.
        if(kf < 0) {
	    int tmp = ifla;
	    ifla=iflb;
	    iflb=tmp;
	}

    } else {
        cout << "(attflv:) Unrecognized particle code KF" <<endl;
	exit(1);
    }

}

//***********************************************************************

int HadronContent::ifrkfc(int ia,int ib,int ic,double s)
{

//...Purpose: to return the kf code for flavor having ia ib ic.
//...The kf code for a 2- or 3-quark system of spin s composed by 
//...flavor ia, ib, ic: (the system must be qq or qqq, not qqbar, etc).
//...it corresponds to a diquark system if ic=0.........................

    int ia0 = max( abs(ia), max(abs(ib),abs(ic)));
    int ic0 = min( abs(ia), min(abs(ib),abs(ic)));
    int ib0 = abs(ia+ib+ic)-ia0-ic0;
    int ifrkfc = 1000*ia0 + 100*ib0 + 10*ic0 + int(2.0*(s+0.2))+ 1;
    if(ia != abs(ia) || ib != abs(ib)) ifrkfc = -ifrkfc;

    return ifrkfc;

}

//****************************************************************

double HadronContent::jamemjet(int kfl10,int kfl20)
{
    int kfl1=kfl10;
    int kfl2=kfl20;
    if(abs(kfl1) > 10) {
        int kftmp=kfl1;
        kfl1=kfl2;
        kfl2=kftmp;
    }

    const double parj32=1.0;
    double m1 = particleData->constituentMass(kfl1);
    double m2 = particleData->constituentMass(kfl2);

    // meson.
    if(abs(kfl1) <= 10 && abs(kfl2) <= 10) {
	
      //cout << " mth= " <<1.0*parj32+0.001 + m1 + m2 << " m1= "<< m1 << " m2= "<<m2 <<endl;
      //cin.get();

        double fac=0.0;
        if(abs(kfl1)+abs(kfl2)  > 6) fac=0.3;
        //cout << "mjet= "<< max(1.0, 1.0*parj32+0.001 + m1 + m2 + fac) <<endl;
	//cin.get();
        return max(1.0, 1.0*parj32+0.001 + m1 + m2 + fac);
        //return max(2.0, 1.0*parj32+0.001 + m1 + m2);

    } else {

        int kfla=(abs(kfl2)/1000) % 10;
        int kflb=(abs(kfl2)/100)  % 10;
        int ns=0;
        if(kfla == 3) ns++;
        if(kflb == 3) ns++;
        if(abs(kfl1) == 3) ns++;

	double parc51=2.0;
        return max(parc51+0.001 + ns*0.35,1.0*parj32 + m1 + m2);

    }

}

void HadronContent::findFlavorMeson(int kf,int* ifla,int* iflb, double* prob,int& n)
{
//...    x quark
//... xx0x di-quark
//...   1x leptons
//...   xx gauge and Higgs bosons
//...  xxx meson
//...x0xxx meson
//... xxxx baryon
//...xxxxx baryon

    int j = abs(kf);
    //int  jc1=j         % 10;      // spin
    int  jc2=(j/10)    % 10;
    int  jc3=(j/100)   % 10;
    int  jc4=(j/1000)  % 10;
    //int  jc5=(j/10000) % 10;

    if(jc4 != 0 || j < 100) {
	cout << " findFlavorMeson jc4= " << jc4
	     << " kf= " << kf << endl;
	exit(1);
    }

    int id1a=jc3*10+jc2;
    int id1= kf > 0 ?  id1a : -id1a;

//...Identify the quark and antiquark in mesons:
    int j100, j10;
    int sgn = kf > 0 ? 1 : -1;
    if(jc3%2 == 0) {
	j100=jc3*sgn;
	j10=-jc2*sgn;
    } else {
	j100=jc2*sgn;
        j10=-jc3*sgn;
    }

    if(j10 < 0) {
	ifla[0]=j100;
        iflb[0]=j10;
    } else {
        ifla[0]=j10;
        iflb[0]=j100;
    }
    n=1;

    if(id1 == 11 || id1 == 22) {   // pi0 eta rho0 omega...
	ifla[0]=-1;
	iflb[0]=1;
	ifla[1]=-2;
	iflb[1]=2;
	prob[0]=0.5;
	n=2;
    } else if (j == 331) {  // eta'  (1/6 1/6 2/3)
	ifla[0]=-1;
	iflb[0]=1;
	ifla[1]=-2;
	iflb[1]=2;
	ifla[2]=-3;
	iflb[2]=3;
	prob[0]=0.1667;
	prob[1]=0.3333;
	n=3;
    } 

}

// Combine two flavours (including diquarks) to produce a hadron.
// The weighting of the combination may fail, giving output 0.

//int StringFlav::combine(FlavContainer& flav1, FlavContainer& flav2) {
int HadronContent::combine(int flav1, int flav2) {

  // Recognize largest and smallest flavour.
  int id1Abs = abs(flav1);
  int id2Abs = abs(flav2);
  int idMax  = max(id1Abs, id2Abs);
  int idMin  = min(id1Abs, id2Abs);

  // Construct a meson.
  if (idMax < 9 || idMin > 1000) {

    // Pick spin state and preliminary code.
    int flav = (idMax < 3) ? 0 : idMax - 2;
    //double rndmSpin = mesonRateSum[flav] * rndm->flat();
    //int spin = -1;
    //do rndmSpin -= mesonRate[flav][++spin];
    //while (rndmSpin > 0.);
    //int idMeson = 100 * idMax + 10 * idMin + mesonMultipletCode[spin];
    int idMeson = 100 * idMax + 10 * idMin + 1;

    // For nondiagonal mesons distinguish particle/antiparticle.
    if (idMax != idMin) {
      int sign = (idMax%2 == 0) ? 1 : -1;
      if ( (idMax == id1Abs && flav1 < 0)
        || (idMax == id2Abs && flav2 < 0) ) sign = -sign;
      idMeson *= sign;

    // For light diagonal mesons include uubar - ddbar - ssbar mixing.
    } else if (flav < 2) {
      //double rMix = rndm->flat();
      //if      (rMix < mesonMix1[flav][spin]) idMeson = 110;
      //else if (rMix < mesonMix2[flav][spin]) idMeson = 220;
      //else                                   idMeson = 330;
      //idMeson += mesonMultipletCode[spin];
      idMeson += 1;

      // Additional suppression of eta and eta' may give failure.
      //const double etaSup = 0.6;
      //const double etaPrimeSup = 0.12;
      //if (idMeson == 221 && etaSup < rndm->flat()) return 0;
      //if (idMeson == 331 && etaPrimeSup < rndm->flat()) return 0;
    }

    // Finished for mesons.
    return idMeson;
  }

  // SU(6) factors for baryon production may give failure.
  int idQQ1 = idMax / 1000;
  int idQQ2 = (idMax / 100) % 10;
  int spinQQ = idMax % 10;
  int spinFlav = spinQQ - 1;
  if (spinFlav == 2 && idQQ1 != idQQ2) spinFlav = 4;
  if (idMin != idQQ1 && idMin != idQQ2) spinFlav++;
  //if (baryonCGSum[spinFlav] < rndm->flat() * baryonCGMax[spinFlav])
  //  return 0;

  // Order quarks to form baryon. Pick spin.
  int idOrd1 = max( idMin, max( idQQ1, idQQ2) );
  int idOrd3 = min( idMin, min( idQQ1, idQQ2) );
  int idOrd2 = idMin + idQQ1 + idQQ2 - idOrd1 - idOrd3;
  //int spinBar = (baryonCGSum[spinFlav] * rndm->flat()
  //  < baryonCGOct[spinFlav]) ? 2 : 4;
  int spinBar=2;

  // Distinguish Lambda- and Sigma-like.
  bool LambdaLike = false;
  if (spinBar == 2 && idOrd1 > idOrd2 && idOrd2 > idOrd3) {
    LambdaLike = (spinQQ == 1);
    if (idOrd1 != idMin && spinQQ == 1) LambdaLike = (rndm->flat() < 0.25);
    else if (idOrd1 != idMin)           LambdaLike = (rndm->flat() < 0.75);
  }

  // Form baryon code and return with sign.
  int idBaryon = (LambdaLike)
    ? 1000 * idOrd1 + 100 * idOrd3 + 10 * idOrd2 + spinBar
    : 1000 * idOrd1 + 100 * idOrd2 + 10 * idOrd3 + spinBar;
   return (flav1 > 0) ? idBaryon : -idBaryon;

}

} // end namespace jam2
