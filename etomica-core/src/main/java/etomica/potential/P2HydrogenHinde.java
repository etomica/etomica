/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
//This subroutine computes the H2-H2 potential energy surface as
//described in the Journal of Chemical Physics article "A six-
//dimensional H2-H2 potential energy surface for bound state
//spectroscopy" by Robert J. Hinde.
//
//Please direct comments to rhinde@utk.edu or hinde.robert@gmail.com
//
//This version of the subroutine is dated 1 June 2010.  An earlier
//version of the subroutine (dated 31 October 2007) contained an
//error in the computation of the the potential energy at fixed
//(r, t1, t2, phi) as a power series in (r1, r2).  Specifically,
//the term including the coefficient c20 contained a typographical
//error, and was given as c20*(r2-r0)**2.  This has been changed
//to the correct expression c20*(r1-r0)**2.  Piotr Jankowski first
//notified me of this problem and I appreciate his assistance in
//locating the typographical error.
//
//Before calling this subroutine, your program must call the vinit
//subroutine that performs some initialization work.  The variables
//used in this initialization work are stored in the common block
///vh2h2c/.
//
//The subroutine takes six input parameters:
//
//    r = H2-H2 distance in bohr
//    r1 = bond length of monomer 1 in bohr
//    r2 = bond length of monomer 2 in bohr
//    t1 = theta1 in radians
//    t2 = theta2 in radians
//    phi = phi in radians
//
//    The potential energy is returned to the calling program in the
//    variable potl, in units of wavenumbers.

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.space.Vector;
import etomica.atom.AtomHydrogen;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Degree;
import etomica.util.Constants;

public class P2HydrogenHinde implements IPotential {
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        P2HydrogenHinde pHin = new P2HydrogenHinde(space);
        P2HydrogenPatkowski pPat = new P2HydrogenPatkowski(space);        
        double r0 = BohrRadius.UNIT.toSim(1.448736);
        double t1 = Degree.UNIT.toSim(0);
        double t2 = Degree.UNIT.toSim(0);
        double p1 = Degree.UNIT.toSim(90);
//        System.out.println((Double.POSITIVE_INFINITY-Double.POSITIVE_INFINITY+5));
//        System.exit(1);
        for (int i=0; i<1000; i++) {
            double R = 0.60 + 9.75*i/1000;
            double eHin = pHin.vH2H2(R,r0,r0,t1,t1,p1);
            double ePat = pPat.vH2H2(R, t1, t2, p1);
            System.out.println(R+" "+ePat+" "+eHin);            
        }
    }    
    protected final double [][][][][] coef = new double[18][3][3][4][4];
    protected final double [][][] arep = new double [3][3][4];
    protected final double [][][] crep = new double [3][3][4];
    protected final double [][] c5 = new double [3][3];
    protected final double [][][] c6 = new double [3][3][4];
    protected final double [][][] c8 = new double [3][3][4];
    protected final double [][] cten = new double [3][3];
    protected Boundary boundary;
    protected final Vector dr,com0,com1,hh0,hh1,n0,n1;
    public P2HydrogenHinde(Space space) {
        dr = space.makeVector();
        com0 = space.makeVector();
        com1 = space.makeVector();  
        hh0 = space.makeVector();
        hh1 = space.makeVector();  
        n0 = space.makeVector();
        n1 = space.makeVector();        
        fillData();
    }
    public double vH2H2(double r,double r1, double r2,double th1,double th2,double phi) {
        double pr = BohrRadius.UNIT.fromSim(r);
        double pr1 = BohrRadius.UNIT.fromSim(r1);
        double pr2 = BohrRadius.UNIT.fromSim(r2);
        double[] v = new double[4];
        double[] vv = new double[9];
        double[][][] vij = new double[3][3][4];
        double r0 = 1.4;
        double dr1 = 0.3;
        double dr2 = dr1*dr1;
        double r5 = pr*pr*pr*pr*pr;
        double r6 = r5*pr;
        double r8 = r6*pr*pr;
        double r10 = r5*r5;
        boolean goto200 = false;


//  --- test for R outside spline endpoints.

        if (pr <= 4.210625) {            
//            System.out.println("Warning, R is less than recommended R");
           for (int i=0; i<3; i++) {
               for (int j=0; j<3; j++) {
                   for (int k=0; k<4; k++) {
                       vij[i][j][k] = arep[i][j][k]*Math.exp(pr*crep[i][j][k]);
                       if (Double.isInfinite(vij[i][j][k])) return Double.POSITIVE_INFINITY;
//                       System.out.println(i+" "+j+" "+k+" "+arep[i][j][k]+" "+crep[i][j][k]);
                   }
               }
           }           
           goto200 = true;
        }
        
        if (pr>=12.00) {
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    vij[i][j][0] = cten[i][j]/r10;
                    vij[i][j][3] = c5[i][j]/r5;
                    for (int k=1; k<3; k++) {
                        vij[i][j][k] = 0.0;
                    }
                }
            }

            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    for (int k=0; k<4; k++) {
                        vij[i][j][k] += c6[i][j][k]/r6 + c8[i][j][k]/r8;
                    }
                }
            }
            goto200 = true;
        }

        if (!goto200) {

            //  compute "old R" before shift.
            //  for R < 6.5, it's more convenient to use the old R.
            double rold = (pr+6.5*0.0175)/1.0175;
            //
            //  compute index into spline look-up tables.
            int idx = 0;
            if (rold < 6.5) {
                idx = (int)((rold-4.25)*4.0);
            }
            else if (pr >= 6.5 && pr < 10.0) {
                idx = (int)((rold -6.5)*2.0) + 9;
            }
            else if (pr >= 10.0 && pr < 12.0) {
                idx = (int)(rold -10.0) + 16;
            }
            else throw new RuntimeException("spline look up failed");        

            //  evaluate splines for all 9 (r1, r2) pairs.
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    for (int k=0; k<4; k++) {
                        vij[i][j][k] = 0.0;
                    }
                } 
            }
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    for (int k=0; k<4; k++) {
                        for (int l=0; l<4; l++) {
                            vij[i][j][k] = pr*vij[i][j][k] + coef[idx][i][j][k][l];
                        }
                    }
                } 
            }
        }        
        //  use power series in (r1, r2) to compute potential energy.
        for (int k=0; k<4; k++) {
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    vv[3*i+j] = vij[i][j][k];
//                    System.out.println((3*i+j)+" "+vv[3*i+j] +" "+vij[i][j][k]);
                }
            }

            //  perform symmetry checks.
            if (k == 0 || k == 3) {
                if (vv[1] != vv[3]) throw new RuntimeException("angular terms dont match");
                if (vv[2] != vv[6]) throw new RuntimeException("angular terms dont match");
                if (vv[5] != vv[7]) throw new RuntimeException("angular terms dont match");
            }
            //
            //  compute power series coefficients.
            double c00 = vv[4];
            double c01 = 0.5*(vv[5] - vv[3])/dr1;
            double c10 = 0.5*(vv[7] - vv[1])/dr1;
            double c02 = (0.5*(vv[5]+vv[3]) - vv[4])/dr2;
            double c20 = (0.5*(vv[7]+vv[1]) - vv[4])/dr2;
            double c11 = 0.25*(vv[8] - vv[2] - vv[6] + vv[0])/dr2;
            double aa = 0.5*(vv[2] + vv[0] - 2.0*vv[1])/dr2;
            double bb = 0.5*(vv[8] + vv[6] - 2.0*vv[7])/dr2;
            double c12 = 0.5*(bb-aa)/dr1;
            aa = 0.5*(vv[6] + vv[0] - 2.0*vv[3])/dr2;
            bb = 0.5*(vv[8] + vv[2] - 2.0*vv[5])/dr2;
            double c21 = 0.5*(bb-aa)/dr1;
            double cc = 0.5*(vv[7] + vv[1] - 2.0*vv[4])/dr2;
            double c22 = 0.5*(aa + bb - 2.0*cc)/dr2;

            //  evaluate power series.
            v[k] = c00 + c01*(pr2-r0) + c10*(pr1 - r0) + c02*(pr2 - r0)*(pr2 - r0) + c20*(pr1 - r0)*(pr1 - r0) + c11*(pr1 - r0)*(pr2 - r0) + c12*(pr1 - r0)*(pr2 - r0)*(pr2-r0) +c21*(pr1 - r0)*(pr1 - r0)*(pr2 - r0) + c22*(pr2 - r0)*(pr2-r0)*(pr1 - r0)*(pr1 - r0);
        }        
        //  compute the angular functions.
        double c1 = Math.cos(th1);
        double c2 = Math.cos(th2);
        double s1 = Math.sin(th1);
        double s2 = Math.sin(th2);
        double g000 = 1.0;
        double g202 = 2.5*(3*c1*c1 - 1.0);
        double g022 = 2.5*(3*c2*c2 - 1.0);
        double g224 = 45.0*(0.32*g022*g202 - 16.0*s1*c1*s2*c2*Math.cos(phi)+s1*s1*s2*s2*Math.cos(2*phi))/(4.0*Math.sqrt(70.0));
        //  compute the potential energy from the angular functions and the
        //  interpolated (r1, r2)-dependent coefficients.
        double vpot = v[0]*g000 + v[1]*g202 + v[2]*g022 + v[3]*g224;
//        System.out.println(vpot);
        return (vpot*1E-8*Constants.PLANCK_H*Constants.LIGHT_SPEED);
    }

    protected void fillData() {                
        FileReader fileReader = null;
        String fileName = "P2HydrogenHinde_allCoefs.dat";
        
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);            
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    for (int k=0; k<4; k++) {                        
                        for (int l=0; l<18; l++) {
                            String[] coefs = bufReader.readLine().trim().split(" +");
                            for (int m=0; m<4; m++) {                                
                                coef[l][i][j][k][m] = Double.valueOf(coefs[m]).doubleValue();
                            }
                        }
                    }
                }

            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }

        
        fileName = "P2HydrogenHinde_shortRange.dat";        
        try {
            fileReader = new FileReader(fileName);                              
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }
        try {            
            BufferedReader bufReader = new BufferedReader(fileReader);
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {                    
                    for (int k=0; k<4; k++) {
                        String[] str = bufReader.readLine().trim().split(" +");
                        arep[i][j][k] = Double.valueOf(str[0]).doubleValue();
                        crep[i][j][k] = Double.valueOf(str[1]).doubleValue();
                    }
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());
        }
        
        fileName = "P2HydrogenHinde_longRange.dat";
        
        try {
            fileReader = new FileReader(fileName);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
        }                    
        try {            
            BufferedReader bufReader = new BufferedReader(fileReader);
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    String[] str0 = bufReader.readLine().trim().split(" +");
                    c5[i][j] = Double.valueOf(str0[0]).doubleValue();
                    String[] str1 = bufReader.readLine().trim().split(" +");
                    c6[i][j][0] = Double.valueOf(str1[0]).doubleValue();
                    c8[i][j][0] = Double.valueOf(str1[1]).doubleValue();
                    cten[i][j] = Double.valueOf(str1[2]).doubleValue();                    
                    for (int k=1; k<4; k++) {
                        String[] str2 = bufReader.readLine().trim().split(" +");
                        c6[i][j][k] = Double.valueOf(str2[0]).doubleValue();
                        c8[i][j][k] = Double.valueOf(str2[1]).doubleValue();
                    }
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading from "+fileName+", caught IOException: " + e.getMessage());                    
        }
    }    
    public double getRange() {        
        return Double.POSITIVE_INFINITY;
    }
    
    public void setBox(Box box) {
        boundary = box.getBoundary();    
    }

    public int nBody() {
        return 2;
    }
    public static class P2HydrogenHindeMolecular extends P2HydrogenHinde implements IPotentialMolecular {
        public P2HydrogenHindeMolecular(Space space) {
            super(space);     
        }
        public double energy(IMoleculeList molecules) {
            IMolecule m0 = molecules.getMolecule(0);
            IMolecule m1 = molecules.getMolecule(1);
            IAtom a00 = m0.getChildList().getAtom(0);
            IAtom a01 = m0.getChildList().getAtom(1);
            IAtom a10 = m1.getChildList().getAtom(0);
            IAtom a11 = m1.getChildList().getAtom(1);
            Vector orient00 = a00.getPosition();
            Vector orient01 = a01.getPosition();
            Vector orient10 = a10.getPosition();
            Vector orient11 = a11.getPosition();

            com0.Ev1Pv2(orient00, orient01);
            com0.TE(0.5);        
            com1.Ev1Pv2(orient10, orient11);
            com1.TE(0.5);
            dr.Ev1Mv2(com1, com0);    
            boundary.nearestImage(dr);    
            double r01 = Math.sqrt(dr.squared());
            dr.normalize();
            hh0.Ev1Mv2(orient01, orient00);
            hh1.Ev1Mv2(orient11, orient10);
            double r0 = Math.sqrt(hh0.squared());            
            double s0 = Math.sqrt(hh1.squared());            
            hh0.normalize();            
            hh1.normalize();
            double th1 = Math.acos(dr.dot(hh0));
            double th2 = Math.acos(dr.dot(hh1));
            n0.E(hh0);
            n1.E(hh1);
            n0.XE(dr);
            n1.XE(dr);
            double phi = Math.acos(n0.dot(n1));
            if (th1 > (Math.PI/2.0)) {
                th1 = Math.PI - th1;
                phi = Math.PI - phi;
            } 
            if (th2 > (Math.PI/2.0)) {
                th2 = Math.PI - th2;
                phi = Math.PI - phi;
            }
            double E = vH2H2(r01,r0,s0,th1,th2,phi);
            return E;
        } 
    }
    public static class P2HydrogenHindeAtomic extends P2HydrogenHinde implements IPotentialAtomic {
        public P2HydrogenHindeAtomic(Space space) {
            super(space);     
        }

        public double energy(IAtomList atoms) {
            AtomHydrogen m0 = (AtomHydrogen) atoms.getAtom(0);
            AtomHydrogen m1 = (AtomHydrogen) atoms.getAtom(1);            
            Vector hh0 = m0.getOrientation().getDirection();
            Vector hh1 = m1.getOrientation().getDirection();
            Vector com0 = m0.getPosition();
            Vector com1 = m1.getPosition();

            dr.Ev1Mv2(com1, com0);    
            boundary.nearestImage(dr);    
            double r01 = Math.sqrt(dr.squared());        
            dr.normalize();        
            double th1 = Math.acos(dr.dot(hh0));
            double th2 = Math.acos(dr.dot(hh1));
            n0.E(hh0);
            n1.E(hh1);
            n0.XE(dr);
            n1.XE(dr);
            double phi = Math.acos(n0.dot(n1));
            if (th1 > (Math.PI/2.0)) {
                th1 = Math.PI - th1;
                phi = Math.PI - phi;
            } 
            if (th2 > (Math.PI/2.0)) {
                th2 = Math.PI - th2;
                phi = Math.PI - phi;
            }
            if (th2 == 0) phi = 0;
            double r0 = m0.getBondLength();
            double s0 = m1.getBondLength();
            double E = vH2H2(r01,r0,s0,th1,th2,phi);
            if (E < -20000) {
                vH2H2(r01,r0,s0,th1,th2,phi);
            }
            return E;
        } 
    }
}
