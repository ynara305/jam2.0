#include <Pythia8/PythiaStdlib.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/hadrons/JamStdlib.h>
#include <jam2/hadrons/Baryons.h>
#include <jam2/hadrons/BaryonTable.h>
#include <jam2/hadrons/Mesons.h>

namespace jam2 {
using namespace std;

using Pythia8::ParticleData;

JamParticleData::JamParticleData(Pythia8::ParticleData* pd, Pythia8::Rndm* rnd)
{
  table=pd;
  rndm = rnd;

  Baryons* b=new Baryons();
  Mesons* m=new Mesons();

    //ElementaryParticle* e = new ElementaryParticle();
    //table = new ParticleData();

    nucleons = new ParticleTable();
    deltas = new ParticleTable();
    lambdas = new ParticleTable();
    sigmas = new ParticleTable();
    xis = new ParticleTable();
    mesons = new ParticleTable();

    nStar = new ParticleTable();
    pStar = new ParticleTable();
    dmStar = new ParticleTable();
    d0Star = new ParticleTable();
    dpStar = new ParticleTable();
    dppStar = new ParticleTable();
    smStar = new ParticleTable();
    s0Star = new ParticleTable();
    spStar = new ParticleTable();
    xmStar = new ParticleTable();
    x0Star = new ParticleTable();
    light0Meson = new ParticleTable();
    light1Meson0 = new ParticleTable();
    light1Mesonp = new ParticleTable();
    strMeson0 = new ParticleTable();
    strMesonp = new ParticleTable();

    // include Delta(1232).
    ParticleDataEntry* dm=table->findParticle(1114);
    ParticleDataEntry* d0=table->findParticle(2114);
    ParticleDataEntry* dp=table->findParticle(2214);
    ParticleDataEntry* dpp=table->findParticle(2224);
    /*
    dm->channel(0).meMode(3);
    d0->channel(0).meMode(3);
    d0->channel(1).meMode(3);
    dp->channel(0).meMode(3);
    dp->channel(1).meMode(3);
    dpp->channel(0).meMode(3);
    */

    dmStar->add(dm);
    d0Star->add(d0);
    dpStar->add(dp);
    dppStar->add(dpp);

    proton =  table->findParticle(2212);
    neutron = table->findParticle(2112);
    deltam  = table->findParticle(1114);
    delta0  = table->findParticle(2114);
    deltap  = table->findParticle(2214);
    deltapp = table->findParticle(2224);

    b->setNuclResonance(table,nucleons,nStar,pStar);
    b->setDeltaResonance(table,deltas,dmStar,d0Star,dpStar,dppStar);
    b->setLambdaResonance(table,lambdas);
    b->setSigmaResonance(table,sigmas,smStar,s0Star,spStar);
    b->setXiResonance(table,xis,xmStar,x0Star);

    //double eKinMin=0.004;
    double mpip=table->findParticle(211)->m0(); // 0.13957
    double mpi0=table->findParticle(111)->m0(); // 0.13498
    //double mk0=table->findParticle(311)->m0(); //  0.49761
    //double mkp=table->findParticle(321)->m0(); //  0.49368
    // change the minimum mass of the hadron resonances
    table->mMin(223,2*mpip+mpi0+eKinMin); // omega
    //table->mMin(313,2*mpip+mk0+eKinMin); // K*0
    //table->mMin(323,2*mpip+mk0+eKinMin); // K*0

    m->addMesons(table,mesons,light0Meson, light1Meson0, light1Mesonp,
	    strMeson0,strMesonp);

    // table for the hadrons which may produced from hydro.
    smHadrons = new ParticleTable();
    smHadrons->merge(mesons);

    /*
    for(auto it=smHadrons.begin(); it != smHadrons.end();) {
      if((*it)->id() == 130 || (*it)->id() == 310) {
        auto d=it;
        it=smHadrons.erase(it);
        continue;
      }
      it++;
    }
    */

    smHadrons->merge(nucleons);
    smHadrons->merge(deltas);
    smHadrons->merge(lambdas);
    smHadrons->merge(sigmas);
    smHadrons->merge(xis);
    smHadrons->add(table->findParticle(3334)); // add Omega-


    delete b;
    delete m;

    n1440   = table->findParticle(12112);
    p1440   = table->findParticle(12212);

    //nucleons->print(cout);
    //deltas->print(cout);
    //lambdas->print(cout);
    //sigmas->print(cout);
    //xis->print(cout);
    //mesons->print(cout);

    addDiffractiveParticle(table);

    setPidMap();

  // map of ID to n* d* number. 
  nstarID[2112] = 0;
  nstarID[2212] = 0;
  for(int i=0;i<nucleon_resonance::num;i++) {
    nstarID[nucleon_resonance::pdg1[i]] = i+1;
    nstarID[nucleon_resonance::pdg2[i]] = i+1;
  }
  dstarID[1114] = 0;
  dstarID[2114] = 0;
  dstarID[2214] = 0;
  dstarID[2224] = 0;
  for(int i=0;i<delta_resonance::num;i++) {
    dstarID[delta_resonance::pdg1[i]] = i+1;
    dstarID[delta_resonance::pdg2[i]] = i+1;
    dstarID[delta_resonance::pdg3[i]] = i+1;
    dstarID[delta_resonance::pdg4[i]] = i+1;
  }

  // make hadron decay table.
  decayWidth = new DecayWidth(table);

  int id=111;
  //cout << table->nextId(id) <<endl;
  while (table->nextId(id) != 0 ) {
    int id6 = (id / 1000000 ) % 10;
    int id9 = (id / 1000000000 ) % 10;
    if((id6 > 0 && id6 < 9) || id9 !=0) {id=table->nextId(id); continue;}
    ParticleDataEntry* p=table->findParticle(id);
    id=table->nextId(id);
    //cout << " id= "<< id << " name= "<< p->name() << endl;
    DecayWidthTable* dwt = new DecayWidthTable(p,decayWidth);
    decayWidthTable[p]=dwt;
  }

    //out.open("width.dat");
}

void JamParticleData::setPidMap()
{
    pID[11] = id_lept; // e+
    pID[12] = id_lept; // nu_e
    pID[13] = id_lept; // mu+
    pID[14] = id_lept; // nu_mu
    pID[15] = id_lept; // tau-
    pID[16] = id_lept; // nu_tau
    pID[17] = id_lept; // tau'-
    pID[18] = id_lept; // nu'_tau
    //pID[990] = id_pom;

    pID[1]=id_quark; pID[2]=id_quark; pID[3]=id_quark;
    pID[4]=id_quark; pID[5]=id_quark; pID[6]=id_quark;

    pID[1103]=id_diq;
    pID[2203]=id_diq;
    pID[2101]=id_diq; pID[2103]=id_diq;

    pID[3101]=id_diq; pID[3103]=id_diq;
    pID[3201]=id_diq; pID[3203]=id_diq;
    pID[3303]=id_diq;

    pID[4101]=id_diq; pID[4103]=id_diq;
    pID[4201]=id_diq; pID[4203]=id_diq;
    pID[4301]=id_diq; pID[4303]=id_diq;
    pID[4403]=id_diq;

    pID[5101]=id_diq; pID[5103]=id_diq;
    pID[5201]=id_diq; pID[5203]=id_diq;
    pID[5301]=id_diq; pID[5303]=id_diq;
    pID[5401]=id_diq; pID[5403]=id_diq;
    pID[5503]=id_diq;

    pID[22]=id_boson;
    pID[21]=id_boson;
    pID[23]=id_boson;
    pID[24]=id_boson;

    for(int i=0, N=nucleons->size(); i<N; i++) {
	int id=nucleons->getParticle(i)->id();
	pID[id]=id_nucls;
    }
    pID[2112]=id_nucl; pID[2212]=id_nucl;

    for(int i=0, N=deltas->size(); i<N; i++) {
	int id=deltas->getParticle(i)->id();
	pID[id]=id_delts;
    }
    pID[1114]=id_delt; pID[2114]=id_delt;
    pID[2214]=id_delt; pID[2224]=id_delt;

    for(int i=0, N=lambdas->size(); i<N; i++) {
	int id=lambdas->getParticle(i)->id();
	pID[id]=id_lambdas;
    }
    pID[3122]=id_lambda;

    for(int i=0, N=sigmas->size(); i<N; i++) {
	int id=sigmas->getParticle(i)->id();
	pID[id]=id_sigmas;
    }
    pID[3112]=id_sigma;
    pID[3212]=id_sigma;
    pID[3222]=id_sigma;

    for(int i=0, N=xis->size(); i<N; i++) {
	int id=xis->getParticle(i)->id();
	pID[id]=id_xis;
    }
    pID[3312]=id_xi;
    pID[3322]=id_xi;
    pID[3334]=id_omega;

    // Mesons.
    for(int i=0, N=mesons->size(); i<N; i++) {
	int id=mesons->getParticle(i)->id();
	int kflb=(id/100) % 10;
	int kflc=(id/10)  % 10;
	int id1= kflb*10 + kflc;
	if(id1==22 || id1==33) {
	    pID[id]=id_light0;
	} else if(id1==11 || id1==21) {
	    pID[id]=id_light1;
	} else if(id1==31 || id1==32)  {
	    pID[id]=id_str;
	} else {
	    cout << " meson id ? id= " << id << " id1= " << id1
		<< " i= " << i <<endl;
	    exit(1);
	}
    }
    pID[111]=id_pi; pID[211]=id_pi;
    pID[130]=id_str; pID[310]=id_str;


    // Charm and bottom hadrons have not set yet.
    const int charm_mes[30]={411,413,415,421,423,425,413,431,423,432,433,435,
	441,443,445,
	10411,10413,10421,10423,10431,10433,
	10441,10443,20443,20433,20423,20413,
        100441,100443, 30443};
    for(int i=0;i<30;i++) pID[charm_mes[i]]=id_charm;

    const int bottom_meson[24] = {511,513,515,521,523,525,
	531,533,535,541,543,545,20513,20523,20533,20543,
        10511,10513,10521,10523,10531,10533,10541,10543};
    for(int i=0;i<24;i++) pID[bottom_meson[i]]=id_bott;
    pID[551]=id_bb;
    pID[10551]=id_bb;
    pID[553]=id_bb;
    pID[10553]=id_bb;
    pID[20553]=id_bb;
    pID[100553]=id_bb;
    pID[200553]=id_bb;
    pID[555]=id_bb;

    pID[ 4122] = id_charmb; // Lambda_c
    pID[ 4124] = id_charmb; // Lambda_c(2625)
    pID[14122] = id_charmb; // Lambda_c(2593)+

    pID[4222]=id_charmb; // Sigma++_c
    pID[4212]=id_charmb; // Sigma+_c
    pID[4112]=id_charmb; // Sigma0_c

    pID[4224]=id_charmb; // Sigma*++_c
    pID[4214]=id_charmb; // Sigma*+_c
    pID[4114]=id_charmb; // Sigma*0_c

    pID[4232]=id_charmb; //Xi+_c
    pID[4132]=id_charmb; //Xi_c0
    pID[4314]=id_charmb; //Xi*_c0
    pID[4312]=id_charmb; //Xi'_c0
    pID[4322]=id_charmb; //Xi'_c+
    pID[4324]=id_charmb; //Xi*_c+
    pID[4332]=id_charmb; //Omega_c0
    pID[4334]=id_charmb; //Omega*_c0

    pID[4412] = id_charmb; // Xi_cc+
    pID[4414] = id_charmb; // Xi*_cc+
    pID[4422] = id_charmb; // Xi_cc++
    pID[4424] = id_charmb; // Xi*_cc++
    pID[4432] = id_charmb; // Omega_cc+
    pID[4434] = id_charmb; // Omega*_cc
    pID[4444] = id_charmb; // Omega*_ccc++

    const int bottom_baryon[]={5112,5114,5122,5132,5142,5212,5214,5222,5224,5232,
    5242,5312,5314,5322,5324,5332,5334,5342,
    5412,5414,5422,5424,5432,5434,5442,5444,
    5512,5514,5522,5524,5532,5534,5542,5544,5554};
    for(int i=0;i<35;i++) pID[bottom_baryon[i]]=id_bottb;
}

JamParticleData::~JamParticleData()
{
    delete nucleons;
    delete deltas;
    delete lambdas;
    delete sigmas;
    delete xis;
    delete mesons;
    delete nStar;
    delete pStar;
    delete d0Star;
    delete dmStar;
    delete dpStar;
    delete dppStar;
    delete spStar;
    delete s0Star;
    delete smStar;
    delete x0Star;
    delete xmStar;
    delete light0Meson;
    delete light1Meson0;
    delete light1Mesonp;
    delete strMeson0;
    delete strMesonp;
    delete smHadrons;

  for (auto pdt = decayWidthTable.begin(); pdt != decayWidthTable.end(); ++pdt) {
    delete pdt->second;
  }
  delete decayWidth;

  //out.close();

}

void JamParticleData::addDiffractiveParticle(ParticleData* pd)
{
  pd->addParticle(9901210,"N*_diff0","N*_diffbar0",0, 0);
  pd->addParticle(9901110,"Delta_diff-","Delta_diffbar+",0, -3);
  //pd.addParticle(9902110,"Delta_diff0","Delta_diffbar0",0,  0);
  //pd.addParticle(9902210,"Delta_diff+","Delta_diffbar-",0,  3);
  pd->addParticle(9902220,"Delta_diff++","Delta_diffbar--",0,6);

  pd->addParticle(9903120,"Lambda_diff0","Lambda_diffbar0",0,0);

  pd->addParticle(9903110,"Sigma_diff-","Sigma_diffbar-",0,-3);
  pd->addParticle(9903210,"Sigma_diff0","Sigma_diffbar0",0, 0);
  pd->addParticle(9903220,"Sigma_diff+","Sigma_diffbar+",0, 3);

  pd->addParticle(9903310,"Xi_diff-","Xi_diffbar-",0,-3);
  pd->addParticle(9903320,"Xi_diff0","Xi_diffbar0",0, 0);

  pd->addParticle(9903330,"Omega_diff-","Omega_diffbar-",0,-3);

  pd->addParticle(9900310,"K0_diff0","Kbar0_diff0",0,0);
  pd->addParticle(9900320,"K+_diff0","K-_diff0",0,3);
}

double JamParticleData::lifeTime(ParticleDataEntry* pd,double m, double e)
{
  if(m<1e-10) return 1.e+35;

  // constant width.
  //double wid0 = pd->mWidth()/HBARC;
  //double t0 = wid0 > 1e-8 ?  -log(max(rndm->flat(), 1.e-35) ) /wid0 : 1e+35;
  //return t0 * e/m;

  //if(m > 6.0) {
  if(m > 0.0) {
    //cout << "JamParticleData::lifeTime mass exceed table limit "<< m<<endl;
    double wid = decayWidth->getTotalWidth(pd,m)/HBARC;
    double t = wid > 1e-8 ?  -log(max(rndm->flat(), 1.e-35) ) /wid : 1e+35;

    //cout << pd->name() << " wid= " << wid << " wid0= "<< pd->mWidth()/HBARC <<" t= " << t <<  " gam= "<< e/m <<endl;
    //cin.get();

    /*
    if(pd->id()==1114 || pd->id()==2114 || pd->id()==2214 || pd->id()==2224) {
    out << setw(14) << m
        << setw(14) << t
        << setw(14) << t*e/m
        << setw(14) << 1.0/wid
	<<endl;
    }
    */

    return t * e/m;
  }
  double wid = decayWidthTable[pd] !=0 ?
       decayWidthTable[pd]->getWidth(m)/HBARC :
      decayWidth->getTotalWidth(pd,m)/HBARC;

  if(decayWidthTable[pd]==0)
      cout << "JamParticleData::lifetime id= "<< pd->id()<<endl;

  //double wid=decayWidthTable[pd]->getWidth(m)/HBARC;
  double t = wid > 1e-8 ?  -log(max(rndm->flat(), 1.e-35) ) /wid : 1e+35;
  return t * e/m;
}

// subroutine jambwmas(emmin,emmax,emr,wid,em,icon)
double JamParticleData::BWMass(double emmin,double emmax,double emr,double wid)
{
// generate mass according to the Breit Wigner distribution.
//==================================================================*
//  emmin : minimam mass           (input)
//  emmax : max. mass              (input)
//  emr   : resonance peak mass    (input)
//  wid   : resonance full width   (input)
//  em    : resonance mass         (output)
//==================================================================*

//...Check boundary.
  if(emmax <= emmin) return emmax;

  double const0=2.0*(emmin-emr)/wid;
  double const1=atan(const0);
  double const2=M_PI/2.0-const1;
  double xmax=(atan(2.0*(emmax-emr)/wid)-const1)/const2;
  xmax=min(1.0,xmax);
  double x=xmax*rndm->flat();
  double t=tan(x*const2);
  return emr+0.5*wid*(const0+t)/(1.0-const0*t);
}

}
