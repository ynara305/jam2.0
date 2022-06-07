#ifndef jam2_hadrons_ParticleTable_h
#define jam2_hadrons_ParticleTable_h

#include <iostream>
#include <vector>
#include <map>
#include <Pythia8/ParticleData.h>

namespace jam2 {

using Pythia8::ParticleDataEntry;

class ParticleTable
{
private:
    typedef std::map<int,ParticleDataEntry*>  MapII;
    typedef std::map<std::string,ParticleDataEntry*>  MapSI;

    std::vector<ParticleDataEntry*> allParticleTable;
    MapII kfTable;
    MapSI pTable,apTable;
    bool  antiTable;
public:
    ParticleTable(bool anti=false);
    ~ParticleTable();

/*
enum particleID {id_quark=1,id_lept=2,id_exc=3,
    id_boson=4,id_diq=5,
    id_tec=6,id_susy=7,id_special=6,
//...Light Mesons.
    id_pi=101,id_light1=100, id_light0=120,
    id_str=130,id_charm=140,id_bott=150,id_cc=160,id_bb=170,id_mdiff=199,
    id_nucl=11, id_nucls=12, id_delt=13,id_delts=14,
    id_lambda=21,id_lambdas=22, id_sigma=31,id_sigmas=32, id_xi=41,id_xis=42,
    id_omega=51, id_charmb=61,id_bottb=72, id_bdiff=99};
    */

    ParticleDataEntry* find(int kf);
    ParticleDataEntry* find(std::string name);
    ParticleDataEntry* findP(std::string name);
    ParticleDataEntry* findAntiP(std::string name);
    void add(ParticleDataEntry* p);
    void merge(ParticleTable* t);

    ParticleDataEntry* getParticle(const int i) {
	if(i >= 0 && i < (int)allParticleTable.size()) {
	    return allParticleTable[i];
	}else {
	    std::cerr << "(ParticleTable::) invalid i " << i << std::endl;
	    exit(1);
	}
    }
    int size()    {return (int)allParticleTable.size();}
    //void activateDecayTable(ParticleTable* table);
   // void activateDecayTable();
    void print(std::ostream& os);

};
}
#endif

