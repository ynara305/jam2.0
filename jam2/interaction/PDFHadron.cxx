#include <jam2/interaction/PDFHadron.h>

namespace jam2 {

using namespace std;

double PDFHadron::xf(int id, double x, double Q2) {

  // Need to update if flavour, x or Q2 changed.
  // Use idSav = 9 to indicate that ALL flavours are up-to-date.
  // Assume that flavour and antiflavour always updated simultaneously.
  if ( (abs(idSav) != abs(id) && idSav != 9) || x != xSav || Q2 != Q2Sav) {
	idSav = id; xfUpdate(id, x, Q2); xSav = x; Q2Sav = Q2;}

//2017/2/15/ynara+
    int id1 = (idBeamAbs/1000)%10;
    int id2 = (idBeamAbs/100)%10;
    int id3 = (idBeamAbs/10)%10;
    int idprod = id1*id2*id3;
    //int idp=0;
    //if(idprod == 4) idp=1; // proton like
    //if(idprod == 2) idp=2; // neutron like
//2017/2/15/ynara-


  // Baryon and nondiagonal meson beams: only p, pbar, pi+, pi- for now.
  if (idBeamAbs == 2212 || idBeamAbs == 211) {
    int idNow = (idBeam > 0) ? id : -id;
    int idAbs = abs(id);
    if (idNow ==  0 || idAbs == 21) return max(0., xg);
    if (idNow ==  1) return max(0., xd);
    if (idNow == -1) return max(0., xdbar);
    if (idNow ==  2) return max(0., xu);
    if (idNow == -2) return max(0., xubar);
    if (idNow ==  3) return max(0., xs);
    if (idNow == -3) return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Baryon beams: n and nbar by isospin conjugation of p and pbar.
  } else if (idBeamAbs == 2112) {
    int idNow = (idBeam > 0) ? id : -id;
    int idAbs = abs(id);
    if (idNow ==  0 || idAbs == 21) return max(0., xg);
    if (idNow ==  1) return max(0., xu);
    if (idNow == -1) return max(0., xubar);
    if (idNow ==  2) return max(0., xd);
    if (idNow == -2) return max(0., xdbar);
    if (idNow ==  3) return max(0., xs);
    if (idNow == -3) return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

//2017/2/15/ynara begin
    // Simple recipes for baryons; average valence parton distribution.
  } else if (idprod != 0) {
      double xval = (xuVal + xdVal)/3.0;
      double xsea = (xuSea + xdSea)/2.0;
      enum Partons {G,D,U,S,C,B,T,Dbar,Ubar,Sbar,Cbar,Bbar,Tbar};
      double xp[13];
      for(int i=0;i<13;i++) xp[i]=0.0;
      xp[D] = xsea;
      xp[U] = xsea;
      xp[S] = xs;
      xp[C] = xc;
      xp[B] = xb;
      xp[Dbar] = xsea;
      xp[Ubar] = xsea;
      xp[Sbar] = xsbar;
      xp[id1] += xval;
      xp[id2] += xval;
      xp[id3] += xval;

    int idNow = (idBeam > 0) ? id : -id;
    int idAbs = abs(id);
    if (idNow ==  0 || idAbs == 21) return max(0., xg);
    if (idNow ==  1) return max(0., xp[D]);
    if (idNow == -1) return max(0., xp[Dbar]);
    if (idNow ==  2) return max(0., xp[U]);
    if (idNow == -2) return max(0., xp[Ubar]);
    if (idNow ==  3) return max(0., xp[S]);
    if (idNow == -3) return max(0., xp[Sbar]);
    if (idAbs ==  4) return max(0., xp[C]);
    if (idAbs ==  5) return max(0., xp[B]);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;
//2017/2/15/ynara end



  // Diagonal meson beams: only pi0, Pomeron for now.
  } else if (idBeam == 111 || idBeam == 990) {
    int idAbs = abs(id);
    if (id ==  0 || idAbs == 21) return max(0., xg);
    if (id == idVal1 || id == idVal2) return max(0., xu);
    if (idAbs <=  2) return max(0., xubar);
    if (idAbs ==  3) return max(0., xs);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Photon beam.
  } else if (idBeam == 22) {
    int idAbs = abs(id);
    if (id ==  0 || idAbs == 21) return max(0., xg);
    if (id ==  1)    return max(0., xd);
    if (id == -1)    return max(0., xdbar);
    if (id ==  2)    return max(0., xu);
    if (id == -2)    return max(0., xubar);
    if (id ==  3)    return max(0., xs);
    if (id == -3)    return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Photon beam inside lepton beam.
  } else if ( ( idBeamAbs == 11 || idBeamAbs == 13 || idBeamAbs == 15 )
    && hasGammaInLepton ) {
    int idAbs = abs(id);
    if (idAbs ==  0 || idAbs == 21) return max(0., xg);
    if (idAbs ==  1) return max(0., xd);
    if (idAbs ==  2) return max(0., xu);
    if (idAbs ==  3) return max(0., xs);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Lepton beam.
  } else {
    if (id == idBeam ) return max(0., xlepton);
    if (abs(id) == 22) return max(0., xgamma);

    //ynara
    cout << " PDF invalid particle " << idBeam
	<<endl;
    exit(1);
    return 0.;
  }

}

//--------------------------------------------------------------------------

// Only valence part of parton densities.

double PDFHadron::xfVal(int id, double x, double Q2) {

  // Need to update if flavour, x or Q2 changed.
  // Use idSav = 9 to indicate that ALL flavours are up-to-date.
  // Assume that flavour and antiflavour always updated simultaneously.
  if ( (abs(idSav) != abs(id) && idSav != 9) || x != xSav || Q2 != Q2Sav)
    {idSav = id; xfUpdate(id, x, Q2); xSav = x; Q2Sav = Q2;}

//2017/2/15/ynara+
    int id1 = (idBeamAbs/1000)%10;
    int id2 = (idBeamAbs/100)%10;
    int id3 = (idBeamAbs/10)%10;
    int idprod = id1*id2*id3;
    //int idp=0;
    //if(idprod == 4) idp=1; // proton like
    //if(idprod == 2) idp=2; // neutron like
//2017/2/15/ynara-


  // Baryon and nondiagonal meson beams: only p, pbar, n, nbar, pi+, pi-.
  if (idBeamAbs == 2212) {
    int idNow = (idBeam > 0) ? id : -id;
    if (idNow == 1) return max(0., xdVal);
    if (idNow == 2) return max(0., xuVal);
    return 0.;
  } else if (idBeamAbs == 2112) {
    int idNow = (idBeam > 0) ? id : -id;
    if (idNow == 1) return max(0., xuVal);
    if (idNow == 2) return max(0., xdVal);
    return 0.;

//2017/2/15/ynara+
  } else if (idprod != 0) {
    int idNow = (idBeam > 0) ? id : -id;
      double xval = (xuVal + xdVal)/3.0;
      //enum Partons {G,D,U,S,C,B};
      double xp[6];
      for(int i=0;i<6;i++) xp[i]=0.0;
      xp[id1] += xval;
      xp[id2] += xval;
      xp[id3] += xval;
      return xp[idNow];
//2017/2/15/ynara-


  } else if (idBeamAbs == 211) {
    int idNow = (idBeam > 0) ? id : -id;
    if (idNow == 2 || idNow == -1) return max(0., xuVal);
    return 0.;

  // Diagonal meson beams: only pi0, Pomeron for now.
  } else if (idBeam == 111 || idBeam == 990) {
    if (id == idVal1 || id == idVal2) return max(0., xuVal);
    return 0.;

  // Photon beam.
  } else if (idBeam == 22) {
    int idAbs = abs(id);
    if (id == idVal1 || id == idVal2) {
      if (idAbs == 1) return max(0., xdVal);
      if (idAbs == 2) return max(0., xuVal);
      if (idAbs == 3) return max(0., xsVal);
      if (idAbs == 4) return max(0., xcVal);
      if (idAbs == 5) return max(0., xbVal);
    }
    return 0.;

  // Lepton beam.
  } else {
    if (id == idBeam) return max(0., xlepton);
    return 0.;
  }

}

//--------------------------------------------------------------------------

// Only sea part of parton densities.

double PDFHadron::xfSea(int id, double x, double Q2) {

  // Need to update if flavour, x or Q2 changed.
  // Use idSav = 9 to indicate that ALL flavours are up-to-date.
  // Assume that flavour and antiflavour always updated simultaneously.
  if ( (abs(idSav) != abs(id) && idSav != 9) || x != xSav || Q2 != Q2Sav)
    {idSav = id; xfUpdate(id, x, Q2); xSav = x; Q2Sav = Q2;}

//2017/2/15/ynara+
    int id1 = (idBeamAbs/1000)%10;
    int id2 = (idBeamAbs/100)%10;
    int id3 = (idBeamAbs/10)%10;
    int idprod = id1*id2*id3;
    //int idp=0;
    //if(idprod == 4) idp=1; // proton like
    //if(idprod == 2) idp=2; // neutron like
//2017/2/15/ynara-

  // Hadron beams.
  if (idBeamAbs > 100) {
    int idNow = (idBeam > 0) ? id : -id;
    int idAbs = abs(id);
    if (idNow == 0 || idAbs == 21) return max(0., xg);
    if (idBeamAbs == 2212) {
      if (idNow ==  1) return max(0., xdSea);
      if (idNow == -1) return max(0., xdbar);
      if (idNow ==  2) return max(0., xuSea);
      if (idNow == -2) return max(0., xubar);
    } else if (idBeamAbs == 2112) {
      if (idNow ==  1) return max(0., xuSea);
      if (idNow == -1) return max(0., xubar);
      if (idNow ==  2) return max(0., xdSea);
      if (idNow == -2) return max(0., xdbar);
//2017/2/15/ynara+
    } else if(idprod != 0) {
      if (idAbs <=  2) return max(0., (xuSea+xdSea)/2.0);
//2017/2/15/ynara-

    } else {
      if (idAbs <=  2) return max(0., xuSea);
    }
    if (idNow ==  3) return max(0., xs);
    if (idNow == -3) return max(0., xsbar);
    if (idAbs ==  4) return max(0., xc);
    if (idAbs ==  5) return max(0., xb);
    if (idAbs == 22) return max(0., xgamma);
    return 0.;

  // Photon beam.
  } else if (idBeamAbs == 22) {
    int idAbs = abs(id);
    if ( id == 0 || idAbs == 21 ) return max(0., xg);
    if ( idAbs == 22 ) return max(0., xgamma);

    // If a valence parton return only the sea part.
    // Otherwise return the total PDF.
    if ( id == idVal1 || id == idVal2 ) {
      if (idAbs ==  1) return max(0., xdSea);
      if (idAbs ==  2) return max(0., xuSea);
      if (idAbs ==  3) return max(0., xsSea);
      if (idAbs ==  4) return max(0., xcSea);
      if (idAbs ==  5) return max(0., xbSea);
    } else {
      if (idAbs ==  1) return max(0., xd);
      if (idAbs ==  2) return max(0., xu);
      if (idAbs ==  3) return max(0., xs);
      if (idAbs ==  4) return max(0., xc);
      if (idAbs ==  5) return max(0., xb);
    }
    return 0.;

  // Lepton beam.
  } else {
    if (abs(id) == 22) return max(0., xgamma);
    return 0.;
  }

}

void PDFHadron::xfUpdate(int id, double x, double Q2)
{
    if ((idBeamAbs/1000)%10 == 0) {
	//pdfPion->xfUpdate(id,x,Q2);
	xg     = pdfPion->xfSea(21,x,Q2);
	xuVal  = pdfPion->xfSea(2,x,Q2);
	xuSea  = pdfPion->xfSea(1,x,Q2);

	xu    = pdfPion->xf(1,x,Q2);
	xubar = pdfPion->xf(-1,x,Q2);
	xs    = pdfPion->xf(3,x,Q2);
	xsbar = pdfPion->xf(-3,x,Q2);
	xc    = pdfPion->xf(4,x,Q2);
	xb    = pdfPion->xf(5,x,Q2);
        xgamma  = pdfPion->xf(22,x,Q2);
    } else {
	//pdfProton->xfUpdate(id,x,Q2);
	xdVal = pdfProton->xfVal(1,x,Q2);
	xuVal = pdfProton->xfVal(2,x,Q2);

	xg    = pdfProton->xfSea(21,x,Q2);
	xdSea = pdfProton->xfSea(1,x,Q2);
	xuSea = pdfProton->xfSea(2,x,Q2);

	xd    = pdfProton->xf(1,x,Q2);
	xdbar = pdfProton->xf(-1,x,Q2);
	xu    = pdfProton->xf(2,x,Q2);
	xubar = pdfProton->xf(-2,x,Q2);
	xs    = pdfProton->xf(3,x,Q2);
	xsbar = pdfProton->xf(-3,x,Q2);
	xc    = pdfProton->xf(4,x,Q2);
	xb    = pdfProton->xf(5,x,Q2);
        xgamma  = pdfProton->xf(22,x,Q2);
    }
}

} // end namespace jam2
