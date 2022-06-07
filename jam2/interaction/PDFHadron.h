#ifndef jam2_interaction_PDFHadron_h
#define jam2_interaction_PDFHadron_h

#include <Pythia8/PartonDistributions.h>
//#include <Pythia8/Basics.h>
//#include <Pythia8/ParticleData.h>
//#include <string>

namespace jam2 {

using Pythia8::PDF;

class PDFHadron : public PDF
{
public:
    PDFHadron(int id=2212) : PDF(id) {}
    void setIdBeam(int id) {idBeam=id; idBeamAbs=abs(idBeam);}
    void setPDF(PDF* pdf1, PDF* pdf2) {pdfProton=pdf1; pdfPion=pdf2;}
private:
    PDF* pdfProton;
    PDF* pdfPion;
    double xf(int id, double x, double Q2);
    double xfVal(int id, double x, double Q2);
    double xfSea(int id, double x, double Q2);
    void xfUpdate(int id, double x, double Q2);
};

}
#endif


