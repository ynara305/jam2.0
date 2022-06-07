#ifndef jam2_collision_InterType_H
#define jam2_collision_InterType_H

//#include <jam2/collision/EventParticle.h>

namespace jam2 {

class InterType
{
private:
  int cType;
  //EventParticle* partner;
  double      cTime;
  double      oTime;
public:
  InterType() { };
  InterType(int ctype, double ctime, double tc) :
  cType(ctype), cTime(ctime), oTime(tc) { };
  // cType=0: no operation
  //      =1: two body collision
  //      =2: decay
  //      =3: wall collision

  int ctype() const {return cType;}
  //EventParticle* particle() const {return partner;}
  double collisionTime() const {return cTime;}
  double orderTime() const {return oTime;}
  void clear() {cType=0; cTime=1e+25; oTime=1e+25;}
  void reset(int ctype, double ctime, double tc) {
    cType=ctype; cTime=ctime; oTime=tc;
  }

};
}
#endif
