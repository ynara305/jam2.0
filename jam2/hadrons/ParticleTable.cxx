#include <iterator>
#include <algorithm>

#include <jam2/hadrons/ParticleTable.h>
#include <jam2/hadrons/Baryons.h>
#include <jam2/hadrons/Mesons.h>

using namespace std;

namespace jam2 {

ParticleTable::ParticleTable(bool anti)
{
    antiTable=anti;
}

ParticleTable::~ParticleTable()
{
    /*
    for(int i = 0; i < (int)allParticleTable.size(); i++) {
	if(allParticleTable[i] !=0)
	delete allParticleTable[i];
    }
    }*/
    allParticleTable.clear();

}

ParticleDataEntry* ParticleTable::find(const string name)
{
    MapSI::iterator it = pTable.find(name);
    MapSI::iterator at = apTable.find(name);
    if(it != pTable.end()) {
	return it->second;
    } else if(at !=apTable.end()) {
	return at->second;
    }else {
	cerr << "ParticleTable::find:Warning particle not found "
	     << name
	     << endl;
	return (ParticleDataEntry*)0;
    }

}

ParticleDataEntry* ParticleTable::findP(const string name)
{
    MapSI::iterator it = pTable.find(name);
    if(it != pTable.end()) {
	return it->second;
    }else {
	return (ParticleDataEntry*)0;
    }

}

ParticleDataEntry* ParticleTable::findAntiP(const string name)
{
    MapSI::iterator at = apTable.find(name);
    if(at !=apTable.end()) {
	return at->second;
    }else {
	return (ParticleDataEntry*)0;
    }

}

ParticleDataEntry* ParticleTable::find(const int kf)
{
    int kfa=kf;
    if(!antiTable) kfa=abs(kf);
    MapII::iterator it = kfTable.find(kfa);
    if(it != kfTable.end()) {
	return it->second;
    } else {
	return (ParticleDataEntry*)0;
    }

}


void ParticleTable::print(ostream& os)
{
    os << " Total Particle: " << allParticleTable.size() << endl;
    os << endl;
    for(int i=0; i< (int)allParticleTable.size(); i++) {
	//os << *allParticleTable[i] << endl;

       ParticleDataEntry p = *allParticleTable[i];

       os << " Name = " << p.name() << " id= " << p.id();
    if(p.hasAnti()) {
        os <<" antiparticle: " << p.hasAnti()
	   << " id= " << -p.id() << std::endl;
    }
        os << " Mass: " << p.m0() << " Width: " << p.mWidth()
           << std::endl;

        os << " color: " << p.colType();
        os << " Charge: " << p.chargeType()
	    << " Spin: " << p.spinType() << std::endl;

	os << endl;
    }
}

void ParticleTable::add(ParticleDataEntry* p)
{
    MapII::iterator it = kfTable.find(p->id());
    if(it == kfTable.end()) {
	//p->setKc(allParticleTable.size());
	allParticleTable.push_back(p);
	kfTable[p->id()]  = p;
	pTable[p->name()] = p;

	if(p->hasAnti()) apTable[p->name(-1)] = p;

	// create anti particle
	/*
	if(antiTable) {
	if(p->name(-1).length() > 0) {
	    ParticleDataEntry* ap = new ParticleDataEntry(*p,-1);
	    //ap->setKc(allParticleTable.size());
	    allParticleTable.push_back(ap);
	    kfTable[ap->getKF()]  = ap;
	    pTable[ap->getName()] = ap;
	}
	}
	*/

    } else {
	cerr << "Warning: ParticleTable::add:This particle already defined kf= "
	     << p->name()
	     << " id= " << p->id()
	     << " skip this particle "  << endl;
    }

}

void ParticleTable::merge(ParticleTable* t)
{
  for(int i=0;i<(int)t->size();i++) {
    ParticleDataEntry* p = t->getParticle(i);
    add(p);
  }
}
 
/*
// activate decay channel for all particle in table.
void ParticleTable::activateDecayTable(ParticleTable* table)
{
    //allParticleTable[0]->addDecayTable();
    for(int i = 0; i < (int)allParticleTable.size(); i++) {
	allParticleTable[i]->activateDecayParticle(table);
    }
}

// activate decay channel for all particle in table.
void ParticleTable::activateDecayTable()
{
    for(int i = 0; i < (int)allParticleTable.size(); i++) {
	allParticleTable[i]->activateDecayParticle(this);
    }
}
*/

}

