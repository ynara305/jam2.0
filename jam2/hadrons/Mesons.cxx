#include <jam2/hadrons/Mesons.h>

namespace jam2 {

using namespace Pythia8;

void Mesons::addMesons(ParticleData* particleData, ParticleTable* m, 
	ParticleTable* l0, ParticleTable* l10, ParticleTable* l1p,
	ParticleTable* str0, ParticleTable* strp)
{
    addLightIsoSpin0Meson(particleData,m,l0);
    addLightIsoSpin1Meson(particleData,m,l10,l1p);
    addStrangeMeson(particleData,m,str0,strp);
}

void Mesons::addLightIsoSpin0Meson(ParticleData* particleData,
	ParticleTable* m, ParticleTable* light0)
{
    //int id=JamHadronData::id_light0;

    const int nkf=13;
    int kf[nkf]={221, 223, 225, 331, 333,335,10221,10223,10331,10333,20223,
	20333, 9010221 };

    for(int i=0;i<nkf;i++) {
	ParticleDataEntry* p=particleData->findParticle(kf[i]);
	if(p) {
	    m->add(p); light0->add(p);
	}
    }

    ParticleDataEntry* p;

    particleData->addParticle(9000221,"f_0(500)",1, 0, 0,
          0.5500    ,0.5000    ,0.2900    , 1.000    ,0.0);
    p = particleData->particleDataEntryPtr(9000221);
    p->addChannel( 1,0.3333333333,3, 211, -211);
    p->addChannel( 1,0.3333333333,3, 111,  111);
    p->addChannel( 1,0.3333333333,3, -211, 211);
    m->add(p);
    light0->add(p);

    //p = new JParticle("eta(1295)","",   20221,id, 0, 0, 0, 0, 0, 0, 0, 0,
    //       1.295    ,0.5300E-01,0.5000E-01, 0.000    ,1, 0);
    particleData->addParticle(100221,"eta(1295)",1, 0, 0,
		1.294, 0.055, 1.24, 1.4,0.0);
    p = particleData->particleDataEntryPtr(100221);
    p->addChannel( 1,0.1666600000,3, 9000111, 111);
    p->addChannel( 1,0.1666600000,3, 9000211, -211);
    p->addChannel( 1,0.1666600000,3, -9000211, 211);
    p->addChannel( 1,0.5000000000,3, 221, 9000221);
    m->add(p);
    light0->add(p);


    //p = new JParticle("omega(1420)","",   30223,id, 2, 0, 0, 0, 0, 0, 0, 0,
    //       1.419    ,0.1740    ,0.2000    , 0.000    ,1, 0);
    particleData->addParticle(100223,"omega(1420)",  3, 0, 0,
           1.419    ,0.1740    ,0.9200    , 1.500,0.0);
    p = particleData->particleDataEntryPtr(100223);
    p->addChannel( 1,0.3333300000,3, 213, -211);
    p->addChannel( 1,0.3333300000,3, -213, 211);
    p->addChannel( 1,0.3333300000,3, 113, 111);
    m->add(p);
    light0->add(p);

    //p = new JParticle("phi(1680)","",   30333,id, 2, 0, 0, 0, 0, 0, 0, 0,
      //     1.680    ,0.1500    ,0.2000    , 0.000    ,1, 0);
    particleData->addParticle(100333,"phi(1680)",3, 0, 0,
          1.680    ,0.1500    ,1.400    , 1.800, 0.0);
    p = particleData->particleDataEntryPtr(100333);
    p->addChannel( 1,0.5000000000,3, 321, -323);
    p->addChannel( 1,0.5000000000,3, -321, 323);
    m->add(p);
    light0->add(p);

    particleData->addParticle(9000223,"f_1(1510)",3, 0, 0,
           1.512    ,0.3500E-01, 1.50    , 2.000    ,0.0);
    p = particleData->particleDataEntryPtr(9000223);
    p->addChannel( 1,0.5000000000,3, 323, -321);
    p->addChannel( 1,0.5000000000,3, -323, 321);
    m->add(p);
    light0->add(p);

    particleData->addParticle(30223,"omega(1650)",3, 0, 0,
           1.67    ,0.315    ,0.50    , 2.000    ,0.0);
    p = particleData->particleDataEntryPtr(30223);
    p->addChannel( 1,0.3333300000,3, 213, -211);
    p->addChannel( 1,0.3333300000,3, -213, 211);
    p->addChannel( 1,0.3333300000,3, 113, 111);
    m->add(p);
    light0->add(p);

}

void Mesons::addLightIsoSpin1Meson(ParticleData* particleData,
	ParticleTable* m, ParticleTable* light10,ParticleTable* light1p)
{
    //int id=JamHadronData::id_pi;
    //            pi0  pi+  rho0 rho+ a_20 a_2+ a_0(1450)0 a_0(1450)+
    const int kf[16]={111, 211, 113, 213, 115, 215, 10111,10211,
	//a_1(1260)0   b_1(1235)0       rho(1450)      a_0(980)
	20113, 20213,  10113, 10213,   100113,100213, 9000111, 9000211};

    for(int i=0;i<16;i++) {
	ParticleDataEntry* p=particleData->findParticle(kf[i]);
	if(p) {
	  m->add(p);
          int kf2=(kf[i]/100)%10;
          int kf1=(kf[i]/10)%10;
	  if(kf1==1 && kf2==1) 
	    light10->add(p);
	  else
	    light1p->add(p);
	}
    }

    ParticleDataEntry* p;

   particleData->addParticle(100111,"pi(1300)0", 1, 0, 0,
           1.300    ,0.4000    ,1.200    , 0.000 );
    p = particleData->particleDataEntryPtr(100111);
    p->addChannel( 1,0.2500000000,4, 213, -211); // l=1
    p->addChannel( 1,0.2500000000,4, -213, 211); // l=1
    p->addChannel( 1,0.5000000000,3, 9000221, 111); //  l=0
    m->add(p);
    light10->add(p);

    particleData->addParticle(100211,"pi(1300)+","pi(1300)-",1, 3, 0,
           1.300    ,0.4000    ,1.2000    , 0.000);
    p = particleData->particleDataEntryPtr(100211);
    p->addChannel( 1,0.2500000000,4, 213, 111);   // l=1
    p->addChannel( 1,0.2500000000,4, 113, 211);   // l=1
    p->addChannel( 1,0.5000000000,3, 9000221,211); //l=0
    m->add(p);
    light1p->add(p);


    particleData->addParticle(10115,"pi_2(1670)0", 5, 0, 0,
           1.670    ,0.2400    ,1.6000    , 2.000);
    p = particleData->particleDataEntryPtr(10115);
    p->addChannel( 1,0.5620000000,5, 225, 111);
    p->addChannel( 1,0.1550000000,4, 213, -211);
    p->addChannel( 1,0.1550000000,4, -213, 211);
    p->addChannel( 1,0.0800000000,5, 9000221, 111);
    p->addChannel( 1,0.0210000000,4, 321, -323);
    p->addChannel( 1,0.0210000000,4, -321, 323);
    p->addChannel( 1,0.0060000000,6, 221, 111);
    m->add(p);
    light10->add(p);

    particleData->addParticle(10215,"pi_2(1670)+","pi_2(1670)-",5,3,0,
           1.670    ,0.2400    ,1.600    , 2.000);
    p = particleData->particleDataEntryPtr(10215);
    p->addChannel( 1,0.5620000000,5, 225, 211);
    p->addChannel( 1,0.1550000000,4, 213, 111);
    p->addChannel( 1,0.1550000000,4, 113, 211);
    p->addChannel( 1,0.0800000000,5, 9000221, 211);
    p->addChannel( 1,0.0210000000,4, 321, -313);
    p->addChannel( 1,0.0210000000,4, -311, 323);
    p->addChannel( 1,0.0060000000,6, 221, 211);
    m->add(p);
    light1p->add(p);

    particleData->addParticle(30113,"rho(1700)0",3, 0, 0,
           1.720    ,0.250    ,1.60    , 2.000);
    p = particleData->particleDataEntryPtr(30113);
    p->addChannel( 1,1.0000000000,3, 113, -211, 211);
    m->add(p);
    light10->add(p);

    particleData->addParticle(30213,"rho(1700)+","rho(1700)-",3, 3, 0,
           1.720    ,0.250    ,1.600    , 2.000);
    p = particleData->particleDataEntryPtr(30213);
    p->addChannel( 1,1.0000000000,3, 213, 111, 111);
    m->add(p);
    light1p->add(p);

}


void Mesons::addStrangeMeson(ParticleData* particleData,
	ParticleTable* m, ParticleTable* str0,ParticleTable* strp)
{
    //int id=JamHadronData::id_str;

    int kf[16]={130,310, 311, 321, 313,323, 10313,10323, 10321,10311,
	315,325, 20313, 20323, 30313,30323};

    // skip 130 and 310
    for(int i=2;i<16;i++) {
	ParticleDataEntry* p=particleData->findParticle(kf[i]);

	if(kf[i]==10313 || kf[i]==10323) {
	    p->setMMin(1.2);
	}
	if(p) {
	  m->add(p); 
          int kf1=(kf[i]/10)%10;
	  if(kf1==1) 
	    str0->add(p);
	  else
	    strp->add(p);
	}
    }

    ParticleDataEntry* p;
    /*
    particleData->addParticle(30313,"K*(1680)0","K*(1680)bar0",3,0, 0,
           1.717    ,0.3220    ,0.0000    , 0.000 );
    p = particleData->particleDataEntryPtr(30313);
    p->addChannel( 1,0,1,  0.2580    ,4, 321, -211);
    p->addChannel( 1,0,1,  0.1290    ,4, 311, 111);
    p->addChannel( 1,0,0,  0.2093    ,3, 321, -213);
    p->addChannel( 1,0,0,  0.1047    ,3, 311, 113);
    p->addChannel( 1,0,0,  0.1993    ,3, 323, -211);
    p->addChannel( 1,0,0,  0.0997    ,3, 313, 111);
    m->add(p);

    particleData->addParticle(30323,"K*(1680)+","K*(1680)-",3, 3, 0,
           1.717    ,0.3220    ,0.0000    , 0.000);
    p = particleData->particleDataEntryPtr(30323);
    p->addChannel( 1,  0.2580    ,4, 311, 211);
    p->addChannel( 1,  0.1290    ,4, 321, 111);
    p->addChannel( 1,  0.2093    ,3, 311, 213);
    p->addChannel( 1,  0.1047    ,3, 321, 113);
    p->addChannel( 1,  0.1993    ,3, 313, 211);
    p->addChannel( 1,  0.0997    ,3, 323, 111);
    m->add(p);
    */

//---------------------------------------------------------------------------
    particleData->addParticle(10315,"K_2(1770)0","K_2(1770)bar0",5,0,0,
           1.773    ,0.1860    ,1.6000    , 2.200);
    p = particleData->particleDataEntryPtr(10315);
    p->addChannel( 1, 0.5333,3, 325, -211);
    p->addChannel( 1, 0.2667,3, 315, 111);
    p->addChannel( 1, 0.0333,4, 323, -211);
    p->addChannel( 1, 0.0167,4, 313, 111);
    p->addChannel( 1, 0.0500,3, 311, 225);
    p->addChannel( 1, 0.0500,4, 311, 333);
    p->addChannel( 1, 0.0500,4, 311, 223);
    m->add(p);
    str0->add(p);

    particleData->addParticle(10325,"K_2(1770)+","K_2(1770)-",5, 3, 0,
           1.773    ,0.1860    ,1.6000    , 2.200);
    p = particleData->particleDataEntryPtr(10325);
    p->addChannel( 1, 0.5333,3, 315, 211);
    p->addChannel( 1, 0.2667,3, 325, 111);
    p->addChannel( 1, 0.0333,4, 313, 211);
    p->addChannel( 1, 0.0167,4, 323, 111);
    p->addChannel( 1, 0.0500,3, 321, 225);
    p->addChannel( 1, 0.0500,4, 321, 333);
    p->addChannel( 1, 0.0500,4, 321, 223);
    m->add(p);
    strp->add(p);
//---------------------------------------------------------------------------

    particleData->addParticle(317,"K_3(1780)0","K_3(1780)bar0",7,0, 0,
           1.776    ,0.1590    ,1.6000    , 2.200);
    p = particleData->particleDataEntryPtr(317);
    p->addChannel( 1,0.3000,5, 321, -213);
    p->addChannel( 1,0.1500,5, 311, 113);
    p->addChannel( 1,0.1820,5, 323, -211);
    p->addChannel( 1,0.0910,5, 313, 111);
    p->addChannel( 1,0.1287,6, 321, -211);
    p->addChannel( 1,0.0643,6, 311, 111);
    p->addChannel( 1,0.0800,6, 311, 221);
    p->addChannel( 1,0.0027,4, 325, -211);
    p->addChannel( 1,0.0013,4, 315, 111);
    m->add(p);
    str0->add(p);

    particleData->addParticle(327,"K_3(1780)+","K_3(1780)-",7, 3, 0,
           1.776    ,0.1590    ,1.6000    , 2.200);
    p = particleData->particleDataEntryPtr(327);
    p->addChannel( 1,0.3000,5, 311, 213);
    p->addChannel( 1,0.1500,5, 321, 113);
    p->addChannel( 1,0.1820,5, 313, 211);
    p->addChannel( 1,0.0910,5, 323, 111);
    p->addChannel( 1,0.1287,6, 311, 211);
    p->addChannel( 1,0.0643,6, 321, 111);
    p->addChannel( 1,0.0800,6, 321, 221);
    p->addChannel( 1,0.0027,4, 315, 211);
    p->addChannel( 1,0.0013,4, 325, 111);
    m->add(p);
    strp->add(p);
//---------------------------------------------------------------------------

    particleData->addParticle(20315,"K_2(1820)0","K_2(1820)bar0",5,0,0,
           1.816    ,0.2760    ,1.600    , 2.200 );
    p = particleData->particleDataEntryPtr(20315);
    p->addChannel( 1, 0.5333,3, 325, -211);
    p->addChannel( 1, 0.2667,3, 315, 111);
    p->addChannel( 1, 0.0333,4, 323, -211);
    p->addChannel( 1, 0.0167,4, 313, 111);
    p->addChannel( 1, 0.0500,3, 311, 225);
    p->addChannel( 1, 0.0500,4, 311, 333);
    p->addChannel( 1, 0.0500,4, 311, 223);
    m->add(p);
    str0->add(p);

    particleData->addParticle(20325,"K_2(1820)+","K_2(1820)-",5, 3, 0,
           1.816    ,0.2760    ,1.600    , 2.200);
    p = particleData->particleDataEntryPtr(20325);
    p->addChannel( 1, 0.5333,3, 315, 211);
    p->addChannel( 1, 0.2667,3, 325, 111);
    p->addChannel( 1, 0.0333,4, 313, 211);
    p->addChannel( 1, 0.0167,4, 323, 111);
    p->addChannel( 1, 0.0500,3, 321, 225);
    p->addChannel( 1, 0.0500,4, 321, 333);
    p->addChannel( 1, 0.0500,4, 321, 223);
    m->add(p);
    strp->add(p);

}

}
