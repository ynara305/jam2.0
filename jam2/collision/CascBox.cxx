#include <jam2/collision/CascBox.h>
#include <jam2/collision/TimeOrder.h>
#include <jam2/collision/EventParticle.h>
#include <jam2/collision/InterList.h>
#include <Pythia8/Basics.h>
#include <jam2/hadrons/JamStdlib.h>

//! @def jam2_collision_CascBox_InterListContainerType
//! This macro selects the implementation of InterListContainer.
#define jam2_collision_CascBox_InterListContainerType 1 /* std::vector 46.55s (default) */
//#define jam2_collision_CascBox_InterListContainerType 2 /* std::multiset 63.93s */
//#define jam2_collision_CascBox_InterListContainerType 3 /* std::list 62.51s */

#include <cstddef>
#include <ostream>
#include <algorithm>
#if jam2_collision_CascBox_InterListContainerType == 2
# include <set>
#elif jam2_collision_CascBox_InterListContainerType == 3
# include <list>
#else /* jam2_collision_CascBox_InterListContainerType == 1 (default) */
# include <vector>
#endif

namespace jam2 {

template<typename InterListRange>
static bool checkNextCollisionTime(InterListRange& range, InterList const* inter, EventParticle const* p1, double dtau1, EventParticle const* p2, double dtau2) {
  if (p2) {
    for (InterList const* iw : range) {
      if(iw == 0 || iw == inter) continue;
      if(iw->getParticle(0) == p1 && iw->getCollisionTime(0) < dtau1) return true;
      if(iw->getParticle(1) == p1 && iw->getCollisionTime(1) < dtau1) return true;
      if(iw->getParticle(0) == p2 && iw->getCollisionTime(0) < dtau2) return true;
      if(iw->getParticle(1) == p2 && iw->getCollisionTime(1) < dtau2) return true;
    }
  } else {
    for (InterList const* iw : range) {
      if(iw == 0 || iw == inter) continue;
      if(iw->getParticle(0) == p1 && iw->getCollisionTime(0) < dtau1) return true;
      if(iw->getParticle(1) == p1 && iw->getCollisionTime(1) < dtau1) return true;
    }
  }
  return false;
}

#if jam2_collision_CascBox_InterListContainerType == 2

class InterListContainer {
  std::multiset<InterList*, TimeOrder> m_tree;

public:
  std::size_t size() const { return m_tree.size(); }

  void push(InterList* inter) { m_tree.insert(inter); }

  void print(std::ostream& ostr) const {
    for(auto& inter : m_tree)
      ostr << " i= "<< inter << " t= "<< inter->getCollisionOrderTime()<< std::endl;
  }

public:
  InterList* next() {
    if (m_tree.empty()) return nullptr;
    return *m_tree.begin();
  }

  void clean(EventParticle const* p1, EventParticle const* p2 = nullptr) {
    if (p2) {
      for (auto it = m_tree.begin(); it != m_tree.end(); ) {
	EventParticle* const ip1 = (*it)->getParticle(0);
	EventParticle* const ip2 = (*it)->getParticle(1);
	if (ip1 == p1 || ip2 == p1 || ip1 == p2 || ip2 == p2) {
	  delete *it;
	  it = m_tree.erase(it);
	} else {
	  ++it;
	}
      }
    } else {
      for (auto it = m_tree.begin(); it != m_tree.end(); ) {
	EventParticle* const ip1 = (*it)->getParticle(0);
	EventParticle* const ip2 = (*it)->getParticle(1);
	if (ip1 == p1 || ip2 == p1) {
	  delete *it;
	  it = m_tree.erase(it);
	} else {
	  ++it;
	}
      }
    }
  }

  bool remove(InterList* inter) {
    if (m_tree.empty()) return false;

    // "inter" is probably the first element because the interactions are
    // processed from eailier ones.
    if (inter == *m_tree.begin()) {
      m_tree.erase(m_tree.begin());
      return true;
    }

    // When "inter" is not found at the end, do a binary search.
#if __cplusplus >= 201402L
    auto it = m_tree.lower_bound(inter);
#else
    auto it = std::lower_bound(m_tree.begin(), m_tree.end(), inter, TimeOrderDescend());
#endif
    for (; it != m_tree.end(); ++it)
      {
	if (*it != inter) continue;
	delete inter;
	m_tree.erase(it);
	return true;
      }

    return false;
  }

  std::size_t clear() {
    std::size_t const original_size = m_tree.size();
    for(InterList const* inter : m_tree) delete inter;
    m_tree.clear();
    return original_size;
  }

  bool checkNextCollisionTime(InterList const* inter, EventParticle const* p1, double dtau1, EventParticle const* p2, double dtau2) const {
    return jam2::checkNextCollisionTime(m_tree, inter, p1, dtau1, p2, dtau2);
  }

  ~InterListContainer() { this->clear(); }
};

#elif jam2_collision_CascBox_InterListContainerType == 3

class InterListContainer {
  std::list<InterList*> m_list;
  bool m_unsorted = false;

public:
  std::size_t size() const { return m_list.size(); }

  void push(InterList* inter) { m_list.push_back(inter); m_unsorted = true; }

  void print(std::ostream& ostr) const {
    for(auto& inter : m_list)
      ostr << " i= "<< inter << " t= "<< inter->getCollisionOrderTime()<< std::endl;
  }

public:
  InterList* next() {
    if (m_list.empty()) return nullptr;
    if (m_unsorted) {
      m_list.sort(TimeOrder());
      m_unsorted = false;
    }
    return *m_list.begin();
  }

  void clean(EventParticle const* p1, EventParticle const* p2 = nullptr) {
    if (p2) {
      for (auto it = m_list.begin(); it != m_list.end(); ) {
	EventParticle* const ip1 = (*it)->getParticle(0);
	EventParticle* const ip2 = (*it)->getParticle(1);
	if (ip1 == p1 || ip2 == p1 || ip1 == p2 || ip2 == p2) {
	  delete *it;
	  it = m_list.erase(it);
	} else {
	  ++it;
	}
      }
    } else {
      for (auto it = m_list.begin(); it != m_list.end(); ) {
	EventParticle* const ip1 = (*it)->getParticle(0);
	EventParticle* const ip2 = (*it)->getParticle(1);
	if (ip1 == p1 || ip2 == p1) {
	  delete *it;
	  it = m_list.erase(it);
	} else {
	  ++it;
	}
      }
    }
  }

  bool remove(InterList* inter) {
    if (m_list.empty()) return false;

    // When "inter" is not found at the end, do a binary search.
    auto const it = std::find(m_list.begin(), m_list.end(), inter);
    if (it != m_list.end() && *it == inter)
      {
	delete inter;
	m_list.erase(it);
	return true;
      }

    return false;
  }

  std::size_t clear() {
    std::size_t const original_size = m_list.size();
    for(InterList const* inter : m_list) delete inter;
    m_list.clear();
    m_unsorted = false;
    return original_size;
  }

  bool checkNextCollisionTime(InterList const* inter, EventParticle const* p1, double dtau1, EventParticle const* p2, double dtau2) const {
    return jam2::checkNextCollisionTime(m_list, inter, p1, dtau1, p2, dtau2);
  }

  ~InterListContainer() { this->clear(); }
};

#else
// default (jam2_collision_CascBox_InterListContainerType == 1)

class InterListContainer {
  typedef std::vector<InterList*> buffer_type;
  buffer_type m_data;
  buffer_type m_buff;
  std::size_t m_sorted_size = 0;
  std::size_t stat_sort = 0, stat_ninter = 0, stat_rmmiss = 0;

public:
  std::size_t size() const {
    return m_data.size();
  }

  void push(InterList* inter) {
    stat_ninter++;
    m_data.push_back(inter);
  }

  void print(std::ostream& ostr) const {
    for(auto& inter : m_data)
      ostr << " i= "<< inter << " t= "<< inter->getCollisionOrderTime()<< std::endl;
  }

private:
  void sort() {
    if (m_sorted_size == m_data.size()) return;
    stat_sort++;
    auto const beg = m_data.begin();
    auto const mid = beg + m_sorted_size;
    auto const end = m_data.end();
    std::sort(mid, end, TimeOrderDescend());
    if (m_sorted_size) {
      m_buff.resize(m_data.size());
      std::merge(beg, mid, mid, end, m_buff.begin(), TimeOrderDescend());
      m_data.swap(m_buff);
    }
    m_sorted_size = m_data.size();
  }

public:
  InterList* next() {
    if (m_data.empty()) return nullptr;
    this->sort();
    return m_data.back();
  }

  void clean(EventParticle const* p1, EventParticle const* p2 = nullptr) {
    this->sort();
    if (p2) {
      m_data.erase(
	std::remove_if(m_data.begin(), m_data.end(), [p1, p2](InterList* const entry) {
	  EventParticle* const ip1 = entry->getParticle(0);
	  EventParticle* const ip2 = entry->getParticle(1);
	  if (ip1 == p1 || ip2 == p1 || ip1 == p2 || ip2 == p2) {
	    delete entry;
	    return true;
	  }
	  return false;
	}), m_data.end());
    } else {
      m_data.erase(
	std::remove_if(m_data.begin(), m_data.end(), [p1](InterList* const entry) {
	  EventParticle* const ip1 = entry->getParticle(0);
	  EventParticle* const ip2 = entry->getParticle(1);
	  if (ip1 == p1 || ip2 == p1) {
	    delete entry;
	    return true;
	  }
	  return false;
	}), m_data.end());
    }
    m_sorted_size = m_data.size();
  }

  bool remove(InterList* inter) {
    if (m_data.empty()) return false;

    // "inter" is probably the last element because the interactions are
    // processed from eailier ones.
    if (inter == m_data.back()) {
      delete inter;
      if (m_sorted_size == m_data.size()) m_sorted_size--;
      m_data.pop_back();
      return true;
    }

    // When "inter" is not found at the end, do a binary search.
    stat_rmmiss++;
    this->sort();
    for (auto it = std::lower_bound(m_data.begin(), m_data.end(), inter, TimeOrderDescend());
	 it != m_data.end(); ++it)
      {
	if (*it != inter) continue;
	delete inter;
	m_data.erase(it);
	m_sorted_size--;
	return true;
      }

    return false;
  }

  std::size_t clear() {
    std::size_t const original_size = m_data.size();
    for(InterList const* inter : m_data) delete inter;
    m_data.clear();
    m_sorted_size = 0;
    return original_size;
  }

  bool checkNextCollisionTime(InterList const* inter, EventParticle const* p1, double dtau1, EventParticle const* p2, double dtau2) const {
    return jam2::checkNextCollisionTime(m_data, inter, p1, dtau1, p2, dtau2);
  }

  ~InterListContainer() {
    this->clear();
    // if (stat_ninter)
    //   std::cerr << "ninter=" << stat_ninter << " sort=" << stat_sort << " nmiss=" << stat_rmmiss << std::endl;
  }
};

#endif

CascBox::CascBox(Pythia8::Settings* s,std::array<double,3>& xmin, std::array<double,3>& xmax, 
    std::array<int,3>& pos, bool ed): settings(s),xMin(xmin),xMax(xmax),position(pos),edgeBox(ed) {
  m_interList = std::make_unique<InterListContainer>();
  withBox=settings->mode("Cascade:box");
  xBox=settings->parm("Cascade:boxLx");
  yBox=settings->parm("Cascade:boxLy");
  zBox=settings->parm("Cascade:boxLz");
  gWidth = settings->parm("Cascade:GaussianWidth");
  if(settings->mode("Cascade:PauliBlocking")==1) {
    pauliR = 1.0/(2*gWidth);
    pauliP = 2.0*gWidth/(HBARC*HBARC);
  // Fusimi function.
  } else {
    pauliR = 1.0/(4*gWidth);
    pauliP = gWidth/(HBARC*HBARC);
  }
}

CascBox:: ~CascBox()
{
  for (auto& p : particles) p = 0;
  particles.clear();
  for (auto& p : neighbors1) p = 0;
  neighbors1.clear();
  for (auto& p : neighbors2) p = 0;
  neighbors2.clear();

  if (m_interList->size() > 0) {
    cout << "CascBox::~CascBox interlist size = "<< m_interList->size()
      << " xmin= "<< xMin[0] << " ymin= "<< xMin[1] << " zmin= "<< xMin[2]
      << " xmax= "<< xMax[0] << " ymax= "<< xMax[1] << " zmax= "<< xMax[2]
      <<endl;
  }
}

void CascBox::eraseParticle(EventParticle* ip) 
{
  delete ip;
  particles.erase(std::find(particles.begin(),particles.end(),ip));
}

void CascBox::clearParticle()
{
  for(auto& p : particles) {
    p->clearBox();
    //cout << " p= "<< p->box()<<endl;
  }
  particles.clear();
}

void CascBox::propagate(double ctime, int opt, int step)
{
  if(step==1) return;
  for(auto& i : particles) {
    if(ctime > i->getT()) i->updateR(ctime,opt);
  }
}

InterList* CascBox::sortInterList() { return m_interList->next(); }

int CascBox::interListSize() { return m_interList->size(); }

void CascBox::setInterList(InterList* inter) { m_interList->push(inter); }

void CascBox::printInterList() { m_interList->print(std::cout); }

void CascBox::cleanInterList(EventParticle* i1, EventParticle* i2) {
  m_interList->clean(i1, i2);
  //for(auto& b : neighbors2) b->cleanInterList(i1,i2); 
}

void CascBox::removeInterList(EventParticle* i1) { m_interList->clean(i1); }

void CascBox::cleanInterList(EventParticle* i1)
{
  for(auto& b : neighbors2)  b->removeInterList(i1);
}

int CascBox::cleanInterList() { return m_interList->clear(); }

void CascBox::removeInter(InterList* inter)
{
  if (!m_interList->remove(inter)) {
    std::cout << "CascBox not find inter "<< inter <<std::endl;
    std::exit(1);
  }
}

bool CascBox::checkNextCollisionTime(InterList const* inter, EventParticle const* p1, double dtau1, EventParticle const* p2, double dtau2) const {
  return m_interList->checkNextCollisionTime(inter, p1, dtau1, p2, dtau2);
}

bool CascBox::isNeighbor(const EventParticle& p)
{
  if(std::find(neighbors2.begin(),neighbors2.end(),p.box()) != neighbors2.end()) return true;
  return false;
}

bool CascBox::isNeighbor(const CascBox& box)
{
  if(std::find(neighbors2.begin(),neighbors2.end(),&box) != neighbors2.end()) return true;
  return false;
}

double CascBox::phaseSpace(EventParticle* i1, EventParticle* i2, EventParticle* ip,double ctime,int opt)
{
  Vec4 r = ip->getR();
  Vec4 p = ip->getP();
  int idp = ip->getID();

  // Loop over all particles. 
  double phase = 0.0;
  for(const auto& i3 : particles) {
      if(i1 == i3) continue; // exclude incoming particle 1.
      if(i2 == i3) continue; // exclude incoming particle 2.
      if(idp != i3->getID()) continue;  // not a nucleon.
      if(ctime < i3->getT()) continue; // not formed
      if(ctime < i3->getTf()) continue; // not formed
      Vec4 r3 = i3->propagate(ctime,opt);
      Vec4 p3 = i3->getP();
      Vec4 dr = r - r3;
      Vec4 dp = p - p3;
      dr[0]=0.0;

      if(withBox) {
        dr[1] = modulo(dr[1] + xBox/2, xBox) - xBox/2;
        dr[2] = modulo(dr[2] + yBox/2, yBox) - yBox/2;
        dr[3] = modulo(dr[3] + zBox/2, zBox) - zBox/2;
      }

      double s = m2(p,p3);
      Vec4 P = p+p3;
      double dr2 = dr.m2Calc() - pow2(dr*P)/s;
      double dp2 = dp.m2Calc() - pow2(dp*P)/s;
      phase += exp(pauliR*dr2 + pauliP*dp2);
  }

  return phase;

}

/*
void CascBox::searchCollisionPair()
{
//typedef std::list<EventParticle*>::iterator EventParIt;

  auto& i0 = ++particles.begin();
  for(auto& i1 = i0; i1 != particles.end(); ++i1)
  for(auto& i2 = particles.begin(); i2 != i1; ++i2) {
    if (TwoBodyInterList* it = hit(*i1,*i2))
      m_interList->push(it);
 }

  // search for collision in the neighbor grid.
  for(auto& neighbor : neighbors1) {
    list<EventParticle*> pl2=neighbor->getParticles();
    for(auto& p1:particles)
    for(auto& p2:pl2) {
      if (TwoBodyInterList* it = hit(p1,p2))
	m_interList->push(it);
    }
  }
}

// search collision between newly produced particle p1 and other old particles.
void CascBox::searchCollisionPair(EventParticle* p1)
{
  // find collision with particles including the neighbor grid.
  for(const auto& neighbor : box->neighbors2) {
    // loop over all particle in this box.
    for(const auto& p2:neighbor->getParticles()) {
      if(p1==p2) continue;
      if (TwoBodyInterList* it = hit(p1,p2))
	m_interList->push(it);
    }
  }

}
*/

} // end of namespace jam2
