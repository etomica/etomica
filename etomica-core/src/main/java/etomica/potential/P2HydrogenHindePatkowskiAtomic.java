/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomHydrogen;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.BohrRadius;
import etomica.units.Degree;
import etomica.units.Kelvin;

import java.io.FileWriter;
import java.io.IOException;

public class P2HydrogenHindePatkowskiAtomic implements IPotentialAtomic {
    protected Boundary boundary;
    protected Vector dr,com0,com1,hh0,hh1,n0,n1;
    protected P2HydrogenHindeAtomic p2Hinde;
    protected P2HydrogenPatkowskiAtomic p2Patkowski;
    protected static final double r0 = BohrRadius.UNIT.toSim(1.448736);
    protected static final double rMin = 1.4;
    public static final double blMax = 1.2;
    protected boolean print = false;
    public FileWriter filePat = null;
    public P2HydrogenHindePatkowskiAtomic(Space space) {
        p2Hinde = new P2HydrogenHindeAtomic(space);
        p2Patkowski = new P2HydrogenPatkowskiAtomic(space);        
        dr = space.makeVector();
        com0 = space.makeVector();
        com1 = space.makeVector();  
        hh0 = space.makeVector();
        hh1 = space.makeVector();  
        n0 = space.makeVector();
        n1 = space.makeVector();
    }

    public double getRange() {    
        return Double.POSITIVE_INFINITY;
    }


    public void setBox(Box box) {
        boundary = box.getBoundary();
        p2Patkowski.setBox(box);
        p2Hinde.setBox(box);
    }


    public int nBody() {
        return 2;
    }
    

        public double energy(IAtomList atoms) {
            AtomHydrogen m0 = (AtomHydrogen) atoms.get(0);
            AtomHydrogen m1 = (AtomHydrogen) atoms.get(1);
            Vector hh0 = m0.getOrientation().getDirection();
            Vector hh1 = m1.getOrientation().getDirection();
            Vector com0 = m0.getPosition();
            Vector com1 = m1.getPosition();
            
            dr.Ev1Mv2(com1, com0);    
            boundary.nearestImage(dr);    
            double r01 = Math.sqrt(dr.squared());        
            if (r01 == 0) return Double.POSITIVE_INFINITY;
            dr.normalize();
            
            double cth1 = dr.dot(hh0);
            double cth2 = -dr.dot(hh1);
            if (cth1 > 1.0) cth1 = 1.0;
            if (cth1 < -1.0) cth1 = -1.0;
            double th1 = Math.acos(cth1);
            if (cth2 > 1.0) cth2 = 1.0;
            if (cth2 < -1.0) cth2 = -1.0;
            double th2 = Math.acos(cth2);
            
            n0.E(hh0);
            n0.PEa1Tv1(-cth1, dr);            
            n0.normalize();

            n1.E(hh1);
            n1.PEa1Tv1(cth2, dr);
            n1.normalize();
            
            if (n0.isNaN() || n1.isNaN()) throw new RuntimeException("oops");

            double cphi = n0.dot(n1);
            if (cphi > 1.0) cphi = 1.0;
            if (cphi < -1.0) cphi = -1.0;
            double phi = Math.acos(cphi);
            
            if (th1 > (Math.PI/2.0)) {
                th1 = Math.PI - th1;
                phi = Math.PI - phi;
            } 
            if (th2 > (Math.PI/2.0)) {
                th2 = Math.PI - th2;
                phi = Math.PI - phi;
            }
            if (th2 == 0) phi = 0;
            double r = m0.getBondLength();
            double s = m1.getBondLength();
            if (r01 < rMin) r01 = rMin; // rMin = 1.4 angstroms , not so quite hard core
            // same as Garberoglio
            // Max bond length criteria also same as garberoglio
            if (r > blMax) r = blMax;
            if (s > blMax) s = blMax;
            double ePat = p2Patkowski.vH2H2(r01,th1,th2,phi);
            double eHin1 = p2Hinde.vH2H2(r01, r, s, th1, th2, phi);
            double eHin2 = p2Hinde.vH2H2(r01, r0, r0, th1, th2, phi);
            double E = 0;
            if (Double.isInfinite(ePat) || Double.isInfinite(eHin1) || Double.isInfinite(eHin2)) return Double.POSITIVE_INFINITY;
            E = ePat + eHin1 - eHin2;
            if (Double.isInfinite(E) || Double.isNaN(E)) throw new RuntimeException("oops "+E);
            
            
            if (print && Math.random() < 0.0001) {
                try {
                    if (filePat == null) filePat = new FileWriter("patConfigurations.dat");                    
                    double rBohr = BohrRadius.UNIT.fromSim(r01);
                    double Ek = Kelvin.UNIT.fromSim(ePat);
                    double dth1 = Degree.UNIT.fromSim(th1);
                    double dth2 = Degree.UNIT.fromSim(th2);
                    double dphi = Degree.UNIT.fromSim(phi);
                    filePat.write(rBohr+" "+dth1+" "+dth2+" "+dphi+" "+Ek+"\n");
                    filePat.flush();                    
                }
                catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
            return E;
        }
}
