#include <jam2/hadrons/Baryons.h>
#include <string>
#include <jam2/hadrons/BaryonTable.h>

namespace jam2 {

using namespace std;
using Pythia8::ParticleData;
using Pythia8::ParticleDataEntry;

static const double Mnucl=0.9383, Mpion=0.138, eKinMin=0.004;
static const double Mlambda = 1.11568, Msigma=1.19744;
static const double mMinD  = Mnucl + Mpion + eKinMin;
static const double mMinNs = Mnucl + 2*Mpion + eKinMin;
static const double mMinDs = Mnucl + 3*Mpion + eKinMin;
//static const double mMinNs = Mnucl + Mpion + eKinMin;
//static const double mMinDs = Mnucl + Mpion + eKinMin;


///////////////////////////////////////////////////////////////////////////

/*
void Baryons::setNucleons(ParticleData* table)
{
    //int id = id_nucl;
    table->addParticle(2112,"nbar0",2112,2,0,0,0.93957,0.0,0.0,0.0,2.655e+14);
    table->addParticle(2212,"pbar-",2212,2,0,0,0.93827);
}

void Baryons::setDeltas(ParticleData* table)
{
    //int id=id_delt;
    table->addParticle(
	   "Delta-","Deltabar+",1114,id, 3,-3, 3, 0, 0, 0, 0, 0,
           1.234    ,0.1200    ,0.1400    , 0.000    ,1, 0);
    table->add(p);
    p->addChannel( 1,0,1,  1.0000    , 2112, -211);

    p = new JParticle("Delta0","Deltabar0",2114,id, 3, 0, 3, 0, 0, 0, 0, 0,
           1.233    ,0.1200    ,0.1400    , 0.000    ,1, 0);
    table->add(p);
    p->addChannel( 1,0,1,  0.3310    , 2212, -211);
    p->addChannel( 1,0,1,  0.6630    , 2112, 111);
    p->addChannel( 1,0,0,  0.0060    , 2112, "gamma");

    p = new JParticle("Delta+","Deltabar-",    2214,id, 3, 3, 3, 0, 0, 0, 0, 0,
           1.232    ,0.1200    ,0.1400    , 0.000    ,1, 0);
    table->add(p);
    p->addChannel( 1,0,1,  0.6630    , 2212, 111);
    p->addChannel( 1,0,1,  0.3310    , 2112, 211);
    p->addChannel( 1,0,0,  0.0060    , 2212, "gamma");

    p = new JParticle("Delta++","Deltabar--",    2224,id, 3, 6, 3, 0, 0, 0, 0, 0,
           1.231    ,0.1200    ,0.1400    , 0.000    ,1, 0);
    table->add(p);
    p->addChannel( 1,0,1,  1.0000    , 2212, 211);

    for(int i=0; i < delta_resonance::num; i++) {
	addDeltaResonance(i,table);
    }
}

void Baryons::setLambdas(ParticleData* table)
{
    int id=ParticleTable::id_lambda;
    JParticle* p = new JParticle(
	    "Lambda0","Lambdabar0",    3122,id, 1, 0, 3,-1, 0, 0, 0, 0,
             1.116    , 0.000    , 0.000    ,0.1093E-42,1, 0);
    table->add(p);
    p->addChannel( 1,0,0,  0.6390    , 2212, -211);
    p->addChannel( 1,0,0,  0.3580    , 2112, 111);
    p->addChannel( 1,0,0,  0.0020    , 2112, "gamma");
    p->addChannel( 1,0,0,  0.0010    , "nu_ebar", "e-", 2212);

}


void Baryons::setSigmas(ParticleTable* table)
{
    int id=ParticleTable::id_sigma;
    JParticle* p1 = new JParticle("Sigma-","Sigmabarr+",3112,id,
	    1,-3,3,-1,0,0,0,0,1.19744,0.0,0.0,44.34,1);
    JParticle* p2 = new JParticle("Sigma0","Sigmabar0",3212,id,
	    1,0,3,-1,0,0,0,0,1.19255,0.0,0.0, 0.0   ,1);
    JParticle* p3 = new JParticle("Sigma+","Sigmabar-",3222,id,
	    1,3,3,-1,0,0,0,0,1.18937,0.0,0.0,23.96,1);

    table->add(p1);
    table->add(p2);
    table->add(p3);
    p1->addChannel( 1,0,0,0.9990000000, 2112, -211);
    p1->addChannel( 1,0,0,0.0010000000, "nu_ebar", "e-", 2112);
    p2->addChannel( 1,0,0,1.0000000000, 3122, "gamma");
    p3->addChannel( 1,0,0,0.5160000000, 2212, 111);
    p3->addChannel( 1,0,0,0.4830000000, 2112, 211);
    p3->addChannel( 1,0,0,0.0010000000, 2212, "gamma");

}

void Baryons::setXis(ParticleTable* table)
{
    int id=ParticleTable::id_xi;
    JParticle* p1 = new JParticle("Xi0","Xibar0",3322,id,
	    1,0,3,-2,0,0,0,0,1.3149,0.0,0.0,87.1,1);
    JParticle* p2 = new JParticle("Xi-","Xibar-",3312,id,
	    1,-3,3,-2,0,0,0,0,1.3213,0.0,0.0,49.1,1);
    table->add(p1);
    table->add(p2);

    p1->addChannel( 1,0,0,0.9954, 3122, 111);
    p1->addChannel( 1,0,0,0.0011, 3122, "gamma");
    p1->addChannel( 1,0,0,0.0035, 3212, "gamma");

    p2->addChannel( 1,0,0,0.9988000000, 3122, -211);
    p2->addChannel( 1,0,0,0.0001000000, 3112, "gamma");
    p2->addChannel( 1,0,0,0.0006000000, "nu_ebar", "e-", 3122);
    p2->addChannel( 1,0,0,0.0004000000, "nu_mubar", "mu-", 3122);
    p2->addChannel( 1,0,0,0.0001000000, "nu_ebar", "e-", 3212);
}
*/

///////////////////////////////////////////////////////////////////////////

void Baryons::addNuclResonance(int inuc,ParticleData* table,
	ParticleTable* nucl,
	ParticleTable* nstar, ParticleTable* pstar)
{
    using namespace nucleon_resonance;
    //int id=id_nucls;
    ParticleDataEntry *p = table->findParticle(pdg1[inuc]);
    if(p) {
	ParticleDataEntry *p2 = table->findParticle(pdg2[inuc]);
	nucl->add(p); nucl->add(p2);
	nstar->add(p);
	pstar->add(p2);
	return;
    }

//     (1/2,1/2) -> (1/2,1/2)+(1,0)
//...2) n+ --> 1/3|p+,pi0>  + 2/3|n0,pi+>
//...3) n+ --> 1/2|d++,pi-> + 1/3|d+,pi0> + 1/6|d0,pi+>
//......n+ --> 1/3|sigma0,k+> + 2/3|sigma+,k0>

//...2) n0 --> 1/3|n0,pi0>  + 2/3|p+,pi->
//...3) n0 --> 1/2|d-,pi+>  + 1/3|d0,pi0> + 1/6|d+,pi->
//......n0 --> 1/3|sigma0,k0> + 2/3|sigma-,k+>

    //double mMin=mMinNs;
    double mMin=mMinD;
    double mMax=2.5;

    // N(1440) 
    if(inuc==0) mMin=mMinD;

    string nam = name[inuc] + '0';
    string anti = name[inuc] + "bar0"; 
    table->addParticle(pdg1[inuc],nam,anti,spin[inuc],0, 0,
	    mass[inuc],width[inuc],mMin, mMax);

    string nam2  = name[inuc] + '+';
    string anti2 = name[inuc] + "bar-"; 
    table->addParticle(pdg2[inuc],nam2,anti2,spin[inuc],3, 0,
		mass[inuc],width[inuc], mMin,mMax);

    ParticleDataEntry* pa1 = table->particleDataEntryPtr(pdg1[inuc]);
    ParticleDataEntry* pa2 = table->particleDataEntryPtr(pdg2[inuc]);
    nucl->add(pa1); nucl->add(pa2);
    nstar->add(pa1); pstar->add(pa2);

    // 1) n0 -> n   + sigma
    if(branch[inuc][0] > 0.0) {
        double br = branch[inuc][0];
        int l = angmom[inuc][0]+3;
	pa1->addChannel(1,br,l,2112,9000221);
	pa2->addChannel(1,br,l,2212,9000221);
    }

    // N -> N + pi
    if(branch[inuc][1] > 0) {
	double brat1 = branch[inuc][1]*1./3.;     // 2) n0 -> n   + pi0
	double brat2 = branch[inuc][1]*2./3.;     // 2) n0 -> p   + pi-
        int l = angmom[inuc][1]+3;
       	pa1->addChannel(1,brat1,l,2112,111);
       	pa1->addChannel(1,brat2,l,2212,-211);
       	pa2->addChannel(1,brat1,l,2212,111);
       	pa2->addChannel(1,brat2,l,2112,211);
    }

    // N -> Delta + pi
    if(branch[inuc][2] > 0) {
	double brat1 = branch[inuc][2]*1./2.;     // 3) n0 -> d-  + pi+
	double brat2 = branch[inuc][2]*1./3.;     // 3) n0 -> d0  + pi0
	double brat3 = branch[inuc][2]*1./6.;     // 3) n0 -> d+  + pi-
        int l = angmom[inuc][2]+3;
       	pa1->addChannel(1,brat1,l,1114,211);
       	pa1->addChannel(1,brat2,l,2114,111);
       	pa1->addChannel(1,brat3,l,2214,-211);
       	pa2->addChannel(1,brat1,l,2224,-211);
       	pa2->addChannel(1,brat2,l,2214,111);
       	pa2->addChannel(1,brat3,l,2114,211);
    }

    // N -> N + rho
    if(branch[inuc][3] > 0) {
	double brat1 = branch[inuc][3]*1./3.;     // 4) n0 -> n   + rho0
	double brat2 = branch[inuc][3]*2./3.;     // 4) n0 -> p   + rho-
        int l = angmom[inuc][3]+3;
       	pa1->addChannel(1,brat1,l,2112,113);
       	pa1->addChannel(1,brat2,l,2212,-213);
       	pa2->addChannel(1,brat1,l,2212,113);
       	pa2->addChannel(1,brat2,l,2112,213);
    }

    // N -> N + eta
    if(branch[inuc][4] > 0) {
        int l = angmom[inuc][4]+3;
	double brat1 = branch[inuc][4];           // 5) n0 -> n   + eta
       	pa1->addChannel(1,brat1,l,2112,221);
       	pa2->addChannel(1,brat1,l,2212,221);
    }

    // N -> lambda   + kaon
    if(branch[inuc][5] > 0) {
        int l = angmom[inuc][5]+3;
       	pa1->addChannel(1,branch[inuc][5],l,3122,311);
       	pa2->addChannel(1,branch[inuc][5],l,3122,321);
    }

    if(branch[inuc][6] > 0) {
	double brat1 = branch[inuc][6]*1./3.;    // 6) n0 -> sigma0   + kaon0
	double brat2 = branch[inuc][6]*2./3.;    // 6) n0 -> sigma-   + kaon+
        int l = angmom[inuc][6]+3;
       	pa1->addChannel(1,brat1,l,3212,311);
       	pa1->addChannel(1,brat2,l,3112,321);
       	pa2->addChannel(1,brat1,l,3212,321);
       	pa2->addChannel(1,brat2,l,3222,311);
    }

     // 6) n0 -> n + omega
    if(branch[inuc][7] > 0) {
        int l = angmom[inuc][7]+3;
       	pa1->addChannel(1,branch[inuc][7],l,2112,223);
       	pa2->addChannel(1,branch[inuc][7],l,2212,223);
    }

    // 7) n0 -> n + eta'
    if(branch[inuc][8] > 0) {
        int l = angmom[inuc][8]+3;
	pa1->addChannel(1,branch[inuc][8],l,2112,331);
	pa2->addChannel(1,branch[inuc][8],l,2212,331);
    }


}

void Baryons::addDeltaResonance(int ix,ParticleData* table,ParticleTable* del,
	ParticleTable* dmstar, ParticleTable* d0star, ParticleTable* dpstar,
	ParticleTable* dppstar)
{
    using namespace delta_resonance;
    //int id=ParticleTable::id_delts;
    ParticleDataEntry *p = table->findParticle(pdg1[ix]);
    if(p) {
	del->add(p);
	ParticleDataEntry *p2=table->findParticle(pdg2[ix]);
	ParticleDataEntry *p3=table->findParticle(pdg3[ix]);
	ParticleDataEntry *p4=table->findParticle(pdg4[ix]);
	del->add(p2); del->add(p3); del->add(p4);
	dmstar->add(p);
	d0star->add(p2);
	dpstar->add(p3);
	dppstar->add(p4);
	//cout << "Delta resonance pa = " << p->id() <<endl;
	//exit(1);
	return;
    }

    //double mMin=mMinDs;
    double mMin=mMinD;
    double mMax=2.5;

    string nam1 = name[ix] + '-';
    string anti1 = name[ix] + "bar+"; 
    table->addParticle(pdg1[ix],nam1,anti1,spin[ix],-3,0,
  	    mass[ix],width[ix],mMin,mMax);

    string nam2 = name[ix] + '0';
    string anti2 = name[ix] + "bar0"; 
    table->addParticle(pdg2[ix],nam2,anti2,spin[ix],0,0,
  	    mass[ix],width[ix],mMin,mMax);

    string nam3  = name[ix] + '+';
    string anti3 = name[ix] + "bar-"; 
    table->addParticle(pdg3[ix],nam3,anti3,spin[ix],3,0,
  	    mass[ix],width[ix],mMin,mMax);

    string nam4  = name[ix] + "++";
    string anti4 = name[ix] + "bar--"; 
    table->addParticle(pdg4[ix],nam4,anti4,spin[ix],6,0,
  	    mass[ix],width[ix],mMin,mMax);

    Pythia8::ParticleDataEntry* pa1 = table->particleDataEntryPtr(pdg1[ix]);
    Pythia8::ParticleDataEntry* pa2 = table->particleDataEntryPtr(pdg2[ix]);
    Pythia8::ParticleDataEntry* pa3 = table->particleDataEntryPtr(pdg3[ix]);
    Pythia8::ParticleDataEntry* pa4 = table->particleDataEntryPtr(pdg4[ix]);

    del->add(pa1);
    del->add(pa2);
    del->add(pa3);
    del->add(pa4);

    dmstar->add(pa1);
    d0star->add(pa2);
    dpstar->add(pa3);
    dppstar->add(pa4);

    // 1)  D* -> n* + pi
    // d-*  --> |n,pi->
    // d0*  --> 2/3|n,pi0> + 1/3|p,pi->
    // d+*  --> 2/3|p,pi0> + 1/3|n,pi+>
    // d++* --> |p,pi+>
    double br = branch[ix][0];
    //int ag = angmom[ix][0];
    if(br > 0.0) {
        int l = angmom[ix][0]+3;
	pa1->addChannel(1,br,        l,12112,-211);
	pa2->addChannel(1,br*2.0/3.0,l,12112, 111);
	pa2->addChannel(1,br*1.0/3.0,l,12212,-211);
	pa3->addChannel(1,br*2.0/3.0,l,12212, 111);
	pa3->addChannel(1,br*1.0/3.0,l,12112, 211);
	pa4->addChannel(1,br,        l,12212, 211);
    }

    // 2)  D* -> n + pi
    br = branch[ix][1];
    if(br > 0) {
        int l = angmom[ix][1]+3;
	pa1->addChannel(1,br,         l, 2112,-211);
	pa2->addChannel(1,br*2.0/3.0, l, 2112,111);
	pa2->addChannel(1,br*1.0/3.0, l, 2212,-211);
	pa3->addChannel(1,br*2.0/3.0, l, 2212,111);
	pa3->addChannel(1,br*1.0/3.0, l, 2112,211);
	pa4->addChannel(1,br,         l, 2212,211);
    }
    // 3)  D* -> d + pi
    // d-*  --> 3/5 |d-,pi0> + 2/5|d0,pi->
    // d0*  --> 1/15|d0,pi0>  + 8/15|d+,pi-> + 2/5|d-,pi+>
    // d+*  --> 1/15|d+,pi0>  + 8/15|d0,pi+> + 2/5|d++,pi->
    // d++* --> 3/5 |d++,pi0> + 2/5|d+,pi+>
    br = branch[ix][2];
    if(br > 0) {
        int l = angmom[ix][2]+3;
	pa1->addChannel(1,br*3.0/5.0,  l, 1114,111);
	pa1->addChannel(1,br*2.0/5.0,  l, 2114,-211);
	pa2->addChannel(1,br*1.0/15.0, l, 2114,111);
	pa2->addChannel(1,br*8.0/15.0, l, 2214,-211);
	pa2->addChannel(1,br*2.0/5.0,  l, 1114,211);
	pa3->addChannel(1,br*1.0/15.0, l, 2214,111);
	pa3->addChannel(1,br*8.0/15.0, l, 2114,211);
	pa3->addChannel(1,br*2.0/5.0,  l, 2224,-211);
	pa4->addChannel(1,br*3.0/5.0,  l, 2224,111);
	pa4->addChannel(1,br*2.0/5.0,  l, 2214,211);
    }
    // 4)  D* -> n + rho
    br = branch[ix][3];
    if(br > 0) {
        int l = angmom[ix][3]+3;
	pa1->addChannel(1,br,         l, 2112,-213);
	pa2->addChannel(1,br*2.0/3.0, l, 2112,113);
	pa2->addChannel(1,br*1.0/3.0, l, 2212,-213);
	pa3->addChannel(1,br*2.0/3.0, l, 2212,113);
	pa3->addChannel(1,br*1.0/3.0, l, 2112,213);
	pa4->addChannel(1,br,         l, 2212,213);
    }
    // 5)  D* -> sigma + kaon
    // d-*  --> |k0,sigma->
    // d0*  --> 2/3|k0,sigma0> + 1/3|k+,sigma->
    // d+*  --> 2/3|k+,sigma0> + 1/3|k0,sigma+>
    // d++* --> |k+,sigma+>
    br = branch[ix][4];
    if(br > 0) {
        int l = angmom[ix][4]+3;
	pa1->addChannel(1,br,         l, 3112,311);
	pa2->addChannel(1,br*2.0/3.0, l, 3212,311);
	pa2->addChannel(1,br*1.0/3.0, l, 3112,321);
	pa3->addChannel(1,br*2.0/3.0, l, 3212,321);
	pa3->addChannel(1,br*1.0/3.0, l, 3222,311);
	pa4->addChannel(1,br,         l, 3222,321);
    }

};

//...S=-1, I=0 (uds)
void Baryons::addLambdaResonance(int ix,ParticleData* table,ParticleTable* lam)
{
    using namespace lambda_resonance;
    //int id=ParticleTable::id_lambdas;

    //if(table->findParticle(pdg[ix])==NULL) return;
    ParticleDataEntry *p = table->findParticle(pdg[ix]);
    if(p) {
	lam->add(p);
	if(p->id()==3124) p->setMMin(1.44); // Lambda(1520)
	if(p->id()==13122) p->setMMin(1.34); // Lambda(1405)
         //particleData->mMin(13122,1.34); // Lambda(1405) 
	return;
    }
    //if(pdg[ix]==13122) return; // Lambda(1405)
    //if(pdg[ix]== 3124) return; // Lambda(1520)
    //if(pdg[ix]==23122) return; // Lambda(1600)
    //if(pdg[ix]==33122) return; // Lambda(1670)

    double minM=mass[ix]-width[ix];
    double maxM=mass[ix]+width[ix];
    //double minM=Msigma + Mpion + eKinMin;
    //double maxM=2.0;

    string nam = name[ix] + '0';
    string anti = name[ix] + "bar0"; 
    table->addParticle(pdg[ix],nam,anti,spin[ix],0,0,
	    mass[ix],width[ix],minM,maxM);
    Pythia8::ParticleDataEntry* pa = table->particleDataEntryPtr(pdg[ix]);
    lam->add(pa);

    //cout << "id= "<< pa->id() << " mmin= "<< pa->mMin()<<endl;

    // (1) Lambda* -> N + Kbar
    //....|Lambda*> -> 1/2|n,Kb0>  + 1/2|p,Kb->
    if(branch[ix][0] > 0.0) {
	double  brat1=branch[ix][0]*0.5; // n+kbar0
        double  brat2=branch[ix][0]*0.5; // p+kbar-
        int l = angmom[ix][0]+3;
       	pa->addChannel(1,brat1,l, 2112,-311);
       	pa->addChannel(1,brat2,l, 2212,-321);
    }
    // (2) Lambda* -> Sigma + pi
    //....|Lambda*> -> 1/3|S0,pi0> + 1/3|S+,pi-> + 1/3|S-,pi+>
    if(branch[ix][1] > 0.0) {
	double  brat1=branch[ix][1]*1.0/3.0; // Sigma0 + pi0
	double  brat2=branch[ix][1]*1.0/3.0; // Sigma+ + pi-
	double  brat3=branch[ix][1]*1.0/3.0; // Sigma- + pi+
        int l = angmom[ix][1]+3;
       	pa->addChannel(1,brat1,l, 3212,111);
       	pa->addChannel(1,brat2,l, 3222,-211);
       	pa->addChannel(1,brat3,l, 3112,211);
    }
    // (3) Lambda* -> Lambda + eta
    if(branch[ix][2] > 0.0) {
	double  brat1=branch[ix][2];
        int l = angmom[ix][2]+3;
       	pa->addChannel(1,brat1,l, 3122,221);
    }
    // (4) Lambda* -> Xi + kaon
    //....|Lambda*> -> 1/2|X0,K0> + 1/2|X-,K+>
    if(branch[ix][3] > 0.0) {
	double  brat1=branch[ix][3]*0.5;  // Lambda -> Xi0 + K0
	double  brat2=branch[ix][3]*0.5;  // Lambda -> Xi- + K+
        int l = angmom[ix][3]+3;
       	pa->addChannel(1,brat1,l, 3322,311);
       	pa->addChannel(1,brat2,l, 3312,321);
    }
    // (5) Lambda* -> Sigma* + pi
    if(branch[ix][4] > 0.0) {
	double  brat1=branch[ix][4]*1.0/3.0; // Sigma0 + pi0
	double  brat2=branch[ix][4]*1.0/3.0; // Sigma+ + pi-
	double  brat3=branch[ix][4]*1.0/3.0; // Sigma- + pi+
        int l = angmom[ix][4]+3;
       	pa->addChannel(1,brat1,l, 3214,  111);
       	pa->addChannel(1,brat2,l, 3224, -211);
       	pa->addChannel(1,brat3,l, 3114,  211);
    }
    // (6) Lambda* -> Lambda + omega
    if(branch[ix][5] > 0.0) {
	double  brat1=branch[ix][5];
        int l = angmom[ix][5]+3;
       	pa->addChannel(1,brat1,l, 3122,223);
    }
    // (7) Lambda* -> N + Kbar*
    if(branch[ix][6] > 0.0) {
	double  brat1=branch[ix][6]*0.5; // n+kbar0
        double  brat2=branch[ix][6]*0.5; // p+kbar-
        int l = angmom[ix][6]+3;
       	pa->addChannel(1,brat1,l, 2112,  -313);
       	pa->addChannel(1,brat2,l, 2212,  -323);
    }

};


void Baryons::addSigmaResonance(int ix,ParticleData* table,ParticleTable* sig,
	ParticleTable* smstar, ParticleTable* s0star, ParticleTable* spstar)
{
    using namespace sigma_resonance;
    ParticleDataEntry *p = table->findParticle(pdg1[ix]);
    if(p) {
	sig->add(p);
        ParticleDataEntry *p2=table->findParticle(pdg2[ix]);
        ParticleDataEntry *p3=table->findParticle(pdg3[ix]);
	sig->add(p2); sig->add(p3);
	smstar->add(p); s0star->add(p2); spstar->add(p3);
	return;
    }
    //int id=ParticleTable::id_sigmas;

    //double mMin=1.35, mMax=2.0;
    //double mMin=mass[ix]-width[ix];
    double mMin=1.35;
    double mMax=mass[ix]+width[ix];
    //double mMin=Mlambda + Mpion + eKinMin;
    //double mMax=2.0;

    static const double c[3][19]={
         {1.0,0.0,  1.0, 0.5,0.5,  1.0, 1.0,0.0, 1.0,0.0,
          1.0, 1.0, 0.75,0.25, 1.0,0.0, 1.0,0.0, 1.0},
         {0.5,0.5,  1.0, 0.5,0.5,  1.0, 0.5,0.5, 0.5,0.5,
          1.0, 1.0, 0.5,0.5,  0.5,0.5, 0.5,0.5, 1.0},
	  {1.0,0.0,  1.0, 0.5,0.5,  1.0, 1.0,0.0, 1.0,0.0,
           1.0, 1.0, 0.75,0.25, 1.0,0.0, 1.0,0.0, 1.0} };

    string nam1  = name[ix] + '-';
    string anti1 = name[ix] + "bar+"; 
    table->addParticle(pdg1[ix],nam1,anti1,spin[ix],-3,0,
		mass[ix],width[ix],mMin,mMax);

    string nam2 = name[ix] + '0';
    string anti2 = name[ix] + "bar0"; 
    table->addParticle(pdg2[ix],nam2,anti2,spin[ix],0,0,
	    mass[ix],width[ix],mMin,mMax);

    string nam3  = name[ix] + '+';
    string anti3 = name[ix] + "bar-"; 
    table->addParticle(pdg3[ix],nam3,anti3,spin[ix],3,0,
		mass[ix],width[ix],mMin,mMax);

    ParticleDataEntry* pa1 = table->particleDataEntryPtr(pdg1[ix]);
    ParticleDataEntry* pa2 = table->particleDataEntryPtr(pdg2[ix]);
    ParticleDataEntry* pa3 = table->particleDataEntryPtr(pdg3[ix]);
    sig->add(pa1);
    sig->add(pa2);
    sig->add(pa3);
    smstar->add(pa1); s0star->add(pa2); spstar->add(pa3);

    // Sigma* -> N + Kbar
    //....|S-*> -> |n,Kb->
    //....|S0*> -> 1/2|n,Kb0>   + 1/2|p,Kb->
    //....|S+*> -> |p,Kb0>
    if(branch[ix][0] > 0.0) {
      //brat(1)=rdec(1)*c(ich,1)
      //brat(2)=rdec(1)*c(ich,2)
        int l = angmom[ix][0]+3;
       	pa1->addChannel(1,branch[ix][0]*c[0][0],l, 2112, -321);
       	pa2->addChannel(1,branch[ix][0]*c[1][0],l, 2112, -311);
       	pa2->addChannel(1,branch[ix][0]*c[1][1],l, 2212, -321);
       	pa3->addChannel(1,branch[ix][0]*c[2][0],l, 2212, -311);
    }
    // Sigma* -> Lambda + pion
    if(branch[ix][1] > 0.0) {
      //brat(3)=rdec(2)*c(ich,3)
        int l = angmom[ix][1]+3;
       	pa1->addChannel(1,branch[ix][1]*c[0][2],l, 3122,-211);
       	pa2->addChannel(1,branch[ix][1]*c[1][2],l, 3122,111);
       	pa3->addChannel(1,branch[ix][1]*c[2][2],l, 3122,211);
    }

    // Sigma* -> Sigma + pion
    //....|S-*> -> 1/2|S0,pi-> + 1/2|S-,pi+>
    //....|S0*> -> 1/2|S+,pi->  + 1/2|S-,pi+>
    //....|S+*> -> 1/2|S+,pi0>  + 1/2|S0,pi+>
    if(branch[ix][2] > 0.0) {
      //brat(4)=rdec(3)*c(ich,4)
      //brat(5)=rdec(3)*c(ich,5)
        int l = angmom[ix][2]+3;
       	pa1->addChannel(1,branch[ix][2]*c[0][3],l, 3212,-211);
       	pa1->addChannel(1,branch[ix][2]*c[0][4],l, 3112,111);
       	pa2->addChannel(1,branch[ix][2]*c[1][3],l, 3222,-211);
       	pa2->addChannel(1,branch[ix][2]*c[1][4],l, 3112,211);
       	pa3->addChannel(1,branch[ix][2]*c[2][3],l, 3222,111);
       	pa3->addChannel(1,branch[ix][2]*c[2][4],l, 3212,211);
    }

    // Sigma* ->  Sigma + eta
    if(branch[ix][3] > 0.0) {
      //brat(6)=rdec(4)*c(ich,6)
        int l = angmom[ix][3]+3;
       	pa1->addChannel(1,branch[ix][3]*c[0][5],l, 3112,221);
       	pa2->addChannel(1,branch[ix][3]*c[1][5],l, 3212,221);
       	pa3->addChannel(1,branch[ix][3]*c[2][5],l, 3222,221);
    }
    // Sigma* ->  Xi +  Kaon
    //....|S-*> -> |Xi-,K0>
    //....|S0*> -> 1/2|Xi0,K0>  + 1/2|Xi-,K+>
    //....|S+*> -> |Xi0,K+>
    if(branch[ix][4] > 0.0) {
      //brat(7)=rdec(5)*c(ich,7)
      //brat(8)=rdec(5)*c(ich,8)
        int l = angmom[ix][4]+3;
       	pa1->addChannel(1,branch[ix][4]*c[0][6],l, 3312,311);
       	pa2->addChannel(1,branch[ix][4]*c[1][6],l, 3322,311);
       	pa2->addChannel(1,branch[ix][4]*c[1][7],l, 3312,321);
       	pa3->addChannel(1,branch[ix][4]*c[2][6],l, 3322,321);
    }

    // Sigma* ->  Sigma(1385) +  pion
    if(branch[ix][5] > 0.0) {
      //brat(9)=rdec(6)*c(ich,9)
      //brat(10)=rdec(6)*c(ich,10)
        int l = angmom[ix][5]+3;
       	pa1->addChannel(1,branch[ix][5]*c[0][8],l, 3214,-211);
       	pa1->addChannel(1,branch[ix][5]*c[0][9],l, 3114,111);
       	pa2->addChannel(1,branch[ix][5]*c[1][8],l, 3224,-211);
       	pa2->addChannel(1,branch[ix][5]*c[1][9],l, 3114,211);
       	pa3->addChannel(1,branch[ix][5]*c[2][8],l, 3224,111);
       	pa3->addChannel(1,branch[ix][5]*c[2][9],l, 3214,211);
    }
    // Sigma* ->  Lambda(1405) +  pion
    if(branch[ix][6] > 0.0) {
      //brat(11)=rdec(7)*c(ich,11)
        int l = angmom[ix][6]+3;
       	pa1->addChannel(1,branch[ix][6]*c[0][10],l, 13122,-211);
       	pa2->addChannel(1,branch[ix][6]*c[1][10],l, 13122,111);
       	pa3->addChannel(1,branch[ix][6]*c[2][10],l, 13122,211);
    }

    // Sigma* ->  Lambda(1520) +  pion
    if(branch[ix][7] > 0.0) {
      //brat(12)=rdec(8)*c(ich,12)
        int l = angmom[ix][7]+3;
       	pa1->addChannel(1,branch[ix][7]*c[0][11],l, 3124,-211);
       	pa2->addChannel(1,branch[ix][7]*c[1][11],l, 3124,111);
       	pa3->addChannel(1,branch[ix][7]*c[2][11],l, 3124,211);
    }

    // Sigma* ->  Delta + Kbar
    //....|S-*> -> 3/4|D-,Kb0>  + 1/4|D0,Kb->
    //....|S0*> -> 1/2|D0,Kb0>  + 1/2|D+,Kb->
    //....|S+*> -> 3/4|D++,Kb-> + 1/4|D+,Kb0>
    if(branch[ix][8] > 0.0) {
      //brat(13)=rdec(9)*c(ich,13)
      //brat(14)=rdec(9)*c(ich,14)
        int l = angmom[ix][8]+3;
       	pa1->addChannel(1,branch[ix][8]*c[0][12],l,1114,-311);
       	pa1->addChannel(1,branch[ix][8]*c[0][13],l,2114,-321);
       	pa2->addChannel(1,branch[ix][8]*c[1][12],l,2114,-311);
       	pa2->addChannel(1,branch[ix][8]*c[1][13],l,2214,-321);
       	pa3->addChannel(1,branch[ix][8]*c[2][12],l,2224,-321);
       	pa3->addChannel(1,branch[ix][8]*c[2][13],l,2214,-311);
    }

    // Sigma* ->  N + Kbar(892)
    if(branch[ix][9] > 0.0) {
      //brat(15)=rdec(10)*c(ich,15)
      //brat(16)=rdec(10)*c(ich,16)
        int l = angmom[ix][9]+3;
        pa1->addChannel(1,branch[ix][9]*c[0][14],l,2112,-323);
        pa2->addChannel(1,branch[ix][9]*c[1][14],l,2112,-313);
       	pa2->addChannel(1,branch[ix][9]*c[1][15],l,2212,-323);
        pa3->addChannel(1,branch[ix][9]*c[2][14],l,2212,-313);
    }
    // Sigma* ==> Xi(1530) + Kaon
    if(branch[ix][10] > 0.0) {
      //brat(17)=rdec(11)*c(ich,17)
      //brat(18)=rdec(11)*c(ich,18)
        int l = angmom[ix][10]+3;
       	pa1->addChannel(1,branch[ix][10]*c[0][16],l,3314,311);
       	pa2->addChannel(1,branch[ix][10]*c[1][16],l,3324,311);
       	pa2->addChannel(1,branch[ix][10]*c[1][17],l,3314,321);
       	pa3->addChannel(1,branch[ix][10]*c[2][16],l,3324,321);
    }

    // Sigma* ==> Lambda + rho
    if(branch[ix][11] > 0.0) {
      //brat(19)=rdec(12)*c(ich,19)
        int l = angmom[ix][11]+3;
       	pa1->addChannel(1,branch[ix][11]*c[0][18],l,3122,-213);
       	pa2->addChannel(1,branch[ix][11]*c[1][18],l,3122,113);
       	pa3->addChannel(1,branch[ix][11]*c[2][18],l,3122,213);
    }
};



void Baryons::addXiResonance(int ix,ParticleData* table,ParticleTable* xi,
	ParticleTable* xmstar, ParticleTable* x0star)
{
    using namespace xi_resonance;
    ParticleDataEntry *p = table->findParticle(pdg1[ix]);
    if(p) {
	xi->add(p);
	ParticleDataEntry *p2 = table->findParticle(pdg2[ix]);
        xi->add(p2);
	xmstar->add(p); x0star->add(p2);
	return;
    }

    //int id=ParticleTable::id_xis;

//....|X-*> -> 1/3|X-,pi0> + 2/3|X0,pi->
//....|X0*> -> 1/3|X0,pi0> + 2/3|X-,pi+>
//....|X-*> -> 1/3|S0,K->  + 2/3|S-,K0>
//....|X0*> -> 1/3|S0,K0>  + 2/3|S+,K->
    static const double c[7] = {0.333333,0.666666,  1.0,
              0.333333,0.666666,  0.333333,0.666666};

/*
      data ((kd(1,i,j),i=1,2),j=1,maxd)/
     $ 3312, 111,   ! X*- -> X-  + pi0
     $ 3322,-211,   ! X*- -> X0  + pi-
     $ 3122,-321,   ! X*- -> L   + K-
     $ 3212,-321,   ! X*- -> S0  + K-
     $ 3112,-311,   ! X*- -> S-  + K0
     $ 3314, 111,   ! X*- -> X*- + pi0
     $ 3324,-211/   ! X*- -> X*0 + pi-
      data ((kd(2,i,j),i=1,2),j=1,maxd)/
     $ 3322, 111,   ! X*0 -> X0  + pi0
     $ 3312, 211,   ! X*0 -> X-  + pi+
     $ 3122,-311,   ! X*0 -> L   + K0
     $ 3212,-311,   ! X*0 -> S0  + K0
     $ 3222,-321,   ! X*0 -> S+  + K-
     $ 3324, 111,   ! X*0 -> X*0 + pi0
     $ 3314, 211/   ! X*0 -> X*- + pi+
*/

    //double mMin=1.46, mMax=2.0;
    double mMin=mass[ix]-width[ix];
    double mMax=mass[ix]+width[ix];

    string nam = name[ix] + '0';
    string anti = name[ix] + "bar0"; 
    table->addParticle(pdg2[ix],nam,anti,spin[ix],0,0,
	    mass[ix],width[ix],mMin,mMax);

    string nam2  = name[ix] + '-';
    string anti2 = name[ix] + "bar+"; 
    table->addParticle(pdg1[ix],nam2,anti2,spin[ix],-3,0,
		mass[ix],width[ix],mMin,mMax);

    ParticleDataEntry* pa1 = table->particleDataEntryPtr(pdg2[ix]);
    ParticleDataEntry* pa2 = table->particleDataEntryPtr(pdg1[ix]);
    xi->add(pa1);
    xi->add(pa2);
    xmstar->add(pa2);
    x0star->add(pa1);

    // Xi* -> Xi + pi
    if(branch[ix][0] > 0.0) {
	double  brat1=branch[ix][0]*c[0]; // x+pi
        double  brat2=branch[ix][0]*c[1]; // x+pi
        int l = angmom[ix][0]+3;
       	pa1->addChannel(1,brat1,l,3322,111);
       	pa1->addChannel(1,brat2,l,3312,211);
       	pa2->addChannel(1,brat1,l,3312,111);
       	pa2->addChannel(1,brat2,l,3322,-211);
    }
    // Xi* -> Lambda + Kbar
    if(branch[ix][1] > 0.0) {
        int l = angmom[ix][1]+3;
       	pa1->addChannel(1,branch[ix][1],l,3122,-311);
       	pa2->addChannel(1,branch[ix][1],l,3122,-321);
    }

    // Xi* -> Sigma + Kbar
    if(branch[ix][2] > 0.0) {
	double brat1 = branch[ix][2]*c[3]; // S+akaon
	double brat2 = branch[ix][2]*c[4]; // S+akaon
        int l = angmom[ix][2]+3;
       	pa1->addChannel(1,brat1,l, 3212,-311);
       	pa1->addChannel(1,brat2,l, 3222,-321);
       	pa2->addChannel(1,brat1,l, 3212,-321);
       	pa2->addChannel(1,brat2,l, 3112,-311);
    }
    if(branch[ix][3] > 0.0) {
	double brat1 = branch[ix][3]*c[5]; // x(1530)+pion
	double brat2 = branch[ix][3]*c[6]; // x(1530)+pion
        int l = angmom[ix][3]+3;
       	pa1->addChannel(1,brat1,l, 3324 ,111);
       	pa1->addChannel(1,brat2,l, 3314 ,211);
       	pa2->addChannel(1,brat1,l, 3314 ,111);
       	pa2->addChannel(1,brat2,l, 3324 ,-211);
    }

};

void Baryons::setNuclResonance(ParticleData* table,ParticleTable* nucl,
	ParticleTable* nstar, ParticleTable* pstar)
{
    nucl->add(table->findParticle(2112));
    nucl->add(table->findParticle(2212));

    for(int i=0; i < nucleon_resonance::num; i++) {
	addNuclResonance(i,table,nucl,nstar,pstar);
    }
}

void Baryons::setDeltaResonance(ParticleData* table,ParticleTable* delta,
	ParticleTable* dmstar, ParticleTable* d0star, ParticleTable* dpstar,
	ParticleTable* dppstar)
{
    delta->add(table->findParticle(1114));
    delta->add(table->findParticle(2114));
    delta->add(table->findParticle(2214));
    delta->add(table->findParticle(2224));


    for(int i=0; i < delta_resonance::num; i++) {
	//if(i==0) continue;
	addDeltaResonance(i,table,delta,dmstar,d0star,dpstar,dppstar);
    }
}

void Baryons::setLambdaResonance(ParticleData* table,ParticleTable* lam)
{
    lam->add(table->findParticle(3122));

    for(int i=0; i < lambda_resonance::num; i++) {
	addLambdaResonance(i,table,lam);
    }

}

void Baryons::setSigmaResonance(ParticleData* table,ParticleTable* sigma,
	ParticleTable* smstar, ParticleTable* s0star, ParticleTable* spstar)
{
    sigma->add(table->findParticle(3112));
    sigma->add(table->findParticle(3212));
    sigma->add(table->findParticle(3222));

    for(int i=0; i < sigma_resonance::num; i++) {
	addSigmaResonance(i,table,sigma,smstar,s0star,spstar);
    }

}

void Baryons::setXiResonance(ParticleData* table,ParticleTable* xi,
	ParticleTable* xmstar, ParticleTable* x0star)
{
    xi->add(table->findParticle(3312));
    xi->add(table->findParticle(3322));

    for(int i=0; i < xi_resonance::num; i++) {
	addXiResonance(i,table,xi,xmstar,x0star);
    }
}

//void Baryons::addBaryons(Pythia8::ParticleData* table)
//{
    //Baryons::setNucleons(table);
    //Baryons::setDeltas(table);
    //Baryons::setLambdas(table);
    //Baryons::setSigmas(table);
    //Baryons::setXis(table);

    //p = new JParticle("Omega-","Omegabar+",3334,ParticleTable::id_omega,
//	    3,-3,3,-3,0, 0, 0, 0,
//           1.672    , 0.000    , 0.000    ,0.3363E-43,1, 0);
//    table->add(p);
//    p->addChannel( 1,0,1,  0.6760    , 3122, "K-");
//    p->addChannel( 1,0,1,  0.2340    , "Xi0", -211);
//    p->addChannel( 1,0,1,  0.0850    , "Xi-", 111);
//    p->addChannel( 1,0,0,  0.0050    , "nu_ebar", "e-", "Xi0");

//}

}
