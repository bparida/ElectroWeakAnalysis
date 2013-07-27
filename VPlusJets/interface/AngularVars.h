// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

//
// find kinematic quantities for W+W -> u+v + j + j
//

// find the cross product of two vectors and sin(theta) between them
void dg_cross2(TVector3 &cross, double &sintheta, 
               TLorentzVector p1, TLorentzVector p2)
{
  TVector3 p1Vector = p1.Vect();
  TVector3 p2Vector = p2.Vect();

  cross = p1Vector.Cross(p2Vector);
  sintheta = cross.Mag() / ( p1Vector.Mag() * p2Vector.Mag() );
}



// inverse euler angle rotations on daughter b from direction of parent a
TLorentzVector dgieuler( TLorentzVector parent, TLorentzVector daughter) {


  TVector3 a = parent.Vect(); 
  TVector3 b = daughter.Vect();

  // Parent's: costheta, sintheta, cosphi, sinphi
  double ct = a.CosTheta();
  double st = sqrt(1. - ct*ct);
  double cp  = cos(a.Phi());
  double sp  = sin(a.Phi());
 
  //bx, by, pz
  double bx = b.Px();
  double by = b.Py();
  double bz = b.Pz();

  // Euler's angle
  double rx =  ct*cp*bx + ct*sp*by - st*bz;
  double ry =  -sp*bx + cp*by ;
  double rz =  st*cp*bx + st*sp*by + ct*bz;
  double rt =  daughter.E();

  TLorentzVector r ( rx, ry, rz, rt);
  return r;
}



////////////////////////////////

// does lorentz trans by beta, gamma (sense ikey)
// on p vector (px,py,pz,e)
TLorentzVector  dgloren( TLorentzVector p, double b, 
                         double g, double ikey) {

  double rx =  p.Px();
  double ry =  p.Py();
  double rz = g *( p.Pz() + ikey * b *p.E() );
  double rt = sqrt( p.M2() + rx*rx + ry*ry + rz*rz );

  TLorentzVector r ( rx, ry, rz, rt);
  return r;
}



/// routine to find cm decay angle of p1 in CM system p1+p2
// returns the cosine of the Jackson angle
double JacksonAngle( TLorentzVector p1, TLorentzVector p2) {

  TLorentzVector ppar = p1 + p2;

  // rotate so z axis is parent direction for p1
  TLorentzVector newp1 = dgieuler(ppar,p1);

  // boost to cm, along z axis
  TLorentzVector  pb = dgloren(newp1, ppar.Beta(), ppar.Gamma(), -1.0);

  // resulting angle is Jackson angle: return cosine theta
  return pb.Pz() / pb.P();
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


void dg_kin_Wuv_Wjj( TLorentzVector pu, TLorentzVector pv, 
                     TLorentzVector pj1, TLorentzVector pj2, 
		     float &cosphipl, float &ctuv, float &ctjj){

  // parent momentum is the sum of all four momenta
  TLorentzVector  ppar = pu + pv + pj1 + pj2;


  //rotate so z axis is parent direction for u,v,j1,j2
  TLorentzVector  pus  =   dgieuler( ppar, pu );
  TLorentzVector  pvs  =   dgieuler( ppar, pv );
  TLorentzVector  pj1s =   dgieuler( ppar, pj1 );
  TLorentzVector  pj2s =   dgieuler( ppar, pj2 );

  //

  //
  // boost to WW cm, along z axis = WW direction
  //
  double gamma = ppar.Gamma();
  double beta = ppar.Beta();

  double ikey = -1.;
  TLorentzVector  pust  = dgloren(pus,beta,gamma,ikey);
  TLorentzVector  pvst  = dgloren(pvs,beta,gamma,ikey);
  TLorentzVector  pj1st = dgloren(pj1s,beta,gamma,ikey);
  TLorentzVector  pj2st = dgloren(pj2s,beta,gamma,ikey);


  //
  // while in WW rest frame look at Wuv vs Wjj decay plane orientation
  // decay plane normal vectors
  //
  TVector3 cruv;
  TVector3 crjj;
  double sinuv, sinjj;
  dg_cross2( cruv, sinuv, pust, pvst);
  dg_cross2( crjj, sinjj, pj1st, pj2st);


  // cosine of angle between normals of the 2 decay planes  
  if( crjj.Mag()==0.0 || crjj.Mag()==0.0) cosphipl = -10.0;
  else cosphipl = cruv.Dot(crjj) / (cruv.Mag() * crjj.Mag() );


  // rotate so that u+v is z axis in WW rest frame
  ppar = pust + pvst;
  TLorentzVector  puro = dgieuler(ppar,pust);


  // boost u to u+v rest frame 
  gamma = ppar.Gamma();
  beta = ppar.Beta();
  TLorentzVector  pucm = dgloren(puro,beta,gamma,ikey);
  ctuv = pucm.CosTheta();


  // boost u to j1+j2 rest frame
  ppar = pj1st + pj2st;
  TLorentzVector  pj1ro = dgieuler(ppar, pj1st);
  gamma = ppar.Gamma();
  beta = ppar.Beta();
  TLorentzVector  pj1cm = dgloren(pj1ro,beta,gamma,ikey);
  ctjj = pj1cm.CosTheta();
}

