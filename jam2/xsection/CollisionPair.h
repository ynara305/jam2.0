#ifndef jam2_xsection_CollisionPair_h
#define jam2_xsection_CollisionPair_h

#include <vector>
#include <utility> // std::swap
#include <Pythia8/ParticleData.h>

//#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )
//#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
//#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

namespace jam2 {

class EventParticle;

class OutGoing
{
public:
    OutGoing(double sigin=0.0,int id01=0, int id02=0,int id03=0, int id04=0):
	sig(sigin), id1(id01), id2(id02), id3(id03), id4(id04)
    { m1=m2=m3=m4=0.0;}
    double sig;
    int id1, id2,id3,id4;
    double m1,m2,m3,m4;
    //ParticleDataEntry *pa1, *pa2, *pa3, *pa4;
    void swap() {
	if(id2 !=0) {
	  std::swap(id1,id2);
	  std::swap(m1,m2);
	}
    }
};

class CollisionPair
{
private:
    int mMode;
    bool anti;
    int iz[2],id[2],pid[2],ibar[2], str[2];
    double m[2], qFac[2];
    bool hasAnti[2];
    int collType;
    double sig,sigel,sigabs,sigt,sigelbw,sigch;
    double xSig;
    int iChannel;
    std::vector<OutGoing>  outgoing;
    double eCM,pCM;
    //Pythia8::ParticleDataEntry* pdata[2];

public:
    CollisionPair() { }
    CollisionPair(int icoltype, EventParticle* p1, EventParticle* p2,
	    double srt=0.0,double pr=0.0); 
    ~CollisionPair() {outgoing.clear(); }

    //Pythia8::ParticleDataEntry* getParticleDataEntry(int i) {return pdata[i];}

    void swap() {
	std::swap(iz[0],iz[1]);
	std::swap(id[0],id[1]);
	std::swap(pid[0],pid[1]);
	std::swap(ibar[0],ibar[1]);
	std::swap(str[0],str[1]);
	std::swap(m[0],m[1]);
	std::swap(hasAnti[0],hasAnti[1]);
	std::swap(qFac[0],qFac[1]);
    }
    void swapOutgoing() {
	for(int i=0;i<(int)outgoing.size();i++) outgoing[i].swap();
    }
    void setAnti() { anti=true;
	for(int jt=0; jt<2; jt++)
	if(hasAnti[jt]) {
	    id[jt] *= -1;
	    iz[jt] *= -1;
	    str[jt] *= -1;
	    ibar[jt] *= -1;
	}
    }
    bool isAnti() {return anti;}

    int mode() {return mMode;}
    void setMode(int i) {mMode=i;}
    void setSizeOutGoing(int i) {outgoing.resize(i);}
    void setOutGoing(int i,double sig0=0.0, int id1=0, int id2=0,
	    int id3=0,int id4=0) {
	outgoing[i].sig=sig0;
	outgoing[i].id1=id1;
	outgoing[i].id2=id2;
	outgoing[i].id3=id3;
	outgoing[i].id4=id4;
    }
    void setOutGoing(double sig0=0.0,int id1=0,int id2=0,int id3=0,int id4=0) {
	outgoing.push_back(OutGoing(sig,id1,id2,id3,id4));
    }
    void setOutGoingID(int id1=0,int id2=0,int id3=0,int id4=0) {
	outgoing[iChannel].id1=id1;
	outgoing[iChannel].id2=id2;
	outgoing[iChannel].id3=id3;
	outgoing[iChannel].id4=id4;
    }
    void setOutGoingMass(double m1,double m2=0.0,double m3=0.0,double m4=0.0) {
	outgoing[iChannel].m1=m1;
	outgoing[iChannel].m2=m2;
	outgoing[iChannel].m3=m3;
	outgoing[iChannel].m4=m4;
    }

    OutGoing&  getOutGoing() {return outgoing[iChannel];}
    OutGoing&  getOutGoing(int i) {return outgoing[i];}
    int outgoingSize() {return outgoing.size();}
    int baryon(int i)  {return ibar[i];}
    int getZ(int i) const {return iz[i];}
    int getCollType() const {return collType;}
    double getCMenergy() const {return eCM;}
    double getCMmomentum() const {return pCM;}
    void setECM(double e) {eCM=e;}
    void setPCM(double p) {pCM=p;}
    void setSigma(double a,double b) {sig=a;sigel=b;}
    void setXS(double a,double b,double c=0.0, double d=0.0) {
	sig=a; sigel=b;sigt=c; sigabs=d;}
    void setXSTotal(double a) {sig=a;}
    void setXSElastic(double a) {sigel=a;}
    void setXsig(double x) {xSig=x;}
    void setSigAbs(double x) {sigabs=x;}
    void setSigInT(double x) {sigt=x;}
    void setElasticBW(double a) {sigelbw=a;}
    void setChargeEx(double a) {sigch=a;}
    double getXsig()  {return xSig;}
    double getSigma() const {return sig;}
    double getSigmaElastic() const {return sigel;}
    double getSigAbs() const {return sigabs;}
    double getSigInT() const {return sigt;}
    double getElasticBW() const {return sigelbw;}
    double getChargeEx() const {return sigch;}
    double M(int i) {return m[i];}
    int getID(int i) {return id[i];}
    int getPID(int i) {return pid[i];}
    int getTotalStrangeness() {return getStrange(0)+getStrange(1);}
    int getStrange(int i) {return str[i];}
    void setID(int i,int k) {id[i]=k;}
    void setM(int i,double a) {m[i]=a;}
    double qFactor(int i) {return qFac[i];}
    void qFactor(double q1, double q2) {qFac[0]=q1; qFac[1]=q2;}

    inline int getPairID() {
      int idmin=std::min(pid[0],pid[1]);
      int idmax=std::max(pid[0],pid[1]);
      return (idmax*(idmax-1))/2+idmin;
    }
    static inline int Pair(int id1,int id2) {
      int idmin=std::min(id1,id2);
      int idmax=std::max(id1,id2);
      return (idmax*(idmax-1))/2+idmin;
    }
    int getChannel() {return iChannel;}
    void setChannel(int i) {iChannel=i;}
    int selectChannel();
    void checkChannel();

};
}
#endif

