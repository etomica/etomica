/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.api.IPotentialAtomic;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;
import etomica.space.Space;
import etomica.units.BohrRadius;
import etomica.units.Electron;
import etomica.units.Hartree;
import etomica.units.Kelvin;

public class P3NitrogenHellmannNonAdditive implements IPotentialAtomic{
            
    protected Boundary boundary;
    protected Space space;
    protected double[] q,pos;    
    public boolean parametersB = true; // set this to false if using parameters for V_12^A potential
    public static final double blN2 = 1.1014;
    protected static final double C9iso = 619.93; // (Eh bohr^9)
    protected static final double C9isoSim = Hartree.UNIT.toSim(C9iso)*Math.pow(BohrRadius.UNIT.toSim(1), 9);
    protected static final double alphaIso = 11.74; // (e^2 bohr^2/Eh)
    protected static final double alphaIsoSim = alphaIso*BohrRadius.UNIT.toSim(1)*BohrRadius.UNIT.toSim(1)*Electron.UNIT.toSim(1)*Electron.UNIT.toSim(1)/Hartree.UNIT.toSim(1);
    protected static final int [] siteID = {0,1,2,1,0};    
    protected static final double[] sitePosA = {-0.682390307412, -0.434260425799, 0.00, 0.434260425799, 0.682390307412};
    protected static final double[] qA = {-794.51215612, 1621.86520764, -1654.70610304, 1621.86520764, -794.51215612};    
    protected static final double[] sitePosB = {-0.680065710389,-0.447763006688, 0.00, 0.447763006688, 0.680065710389};
    protected static final double[] qB = {-832.77884541,1601.24507755,-1536.93246428,1601.24507755, -832.77884541};    

    public P3NitrogenHellmannNonAdditive(Space space) {
        this.space = space;
        q = new double[3];
        pos = new double[7];
        fillData();
    }
    
    protected void fillData() {
        if (parametersB) {     
            for (int i = 0; i < q.length; i++) {
                q[i] = qB[i];
            }
            for (int i = 0; i< sitePosB.length; i++) {
                pos[i] = sitePosB[i];
            }
            
        }
        else {
            for (int i = 0; i < q.length; i++) {
                q[i] = qA[i];
            }
            for (int i = 0; i< sitePosA.length; i++) {
                pos[i] = sitePosA[i];
            }
        }
        pos[5] = -blN2/2.0;
        pos[6] = blN2/2.0;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }
    
    public int nBody() {
        return 3;
    }
    
    public double energy(IAtomList atoms) {
        final IAtomOriented atom0 = (IAtomOriented) atoms.getAtom(0);
        final IAtomOriented atom1 = (IAtomOriented) atoms.getAtom(1);
        final IAtomOriented atom2 = (IAtomOriented) atoms.getAtom(2);
        Vector cm0 = atom0.getPosition();
        Vector cm1 = atom1.getPosition();
        Vector cm2 = atom2.getPosition();
        Vector or0 = atom0.getOrientation().getDirection();
        Vector or1 = atom1.getOrientation().getDirection();
        Vector or2 = atom2.getOrientation().getDirection();

        Vector r0 = space.makeVector();
        Vector r1 = space.makeVector();
        Vector r2 = space.makeVector();
        
        double vDisp = 0;
        for (int i=5; i<7; i++){
            r0.E(cm0);
            r0.PEa1Tv1(pos[i], or0);
            for (int j=5; j<7; j++) {
                r1.E(cm1);
                r1.PEa1Tv1(pos[j], or1);                                
                for (int k=5; k<7; k++) {
                    r2.E(cm2);
                    r2.PEa1Tv1(pos[k], or2);
                    
                    Vector dr01 = space.makeVector();
                    Vector dr12 = space.makeVector();
                    Vector dr02 = space.makeVector();
                    
                    dr01.Ev1Mv2(r1, r0);
                    dr12.Ev1Mv2(r2, r1);
                    dr02.Ev1Mv2(r2, r0);
                    
                    double R01ij = Math.sqrt(dr01.squared());
                    if (R01ij == 0) return Double.POSITIVE_INFINITY;
                    double R12jk = Math.sqrt(dr12.squared());
                    if (R12jk == 0) return Double.POSITIVE_INFINITY;                    
                    double R20ik = Math.sqrt(dr02.squared());
                    if (R20ik == 0) return Double.POSITIVE_INFINITY;
                    double r9 = R01ij*R01ij*R01ij*R12jk*R12jk*R12jk*R20ik*R20ik*R20ik;                    
                    
                    if (dr02.isZero() || dr02.isNaN() || dr01.isNaN() || dr01.isZero() || dr12.isNaN() || dr12.isZero()) throw new RuntimeException("oops");
                    dr02.normalize();
                    dr12.normalize();
                    dr01.normalize();

                    double cth0 = dr01.dot(dr02);
                    dr01.TE(-1);
                    double cth1 = dr01.dot(dr12);
                    dr02.TE(-1);
                    dr12.TE(-1);
                    double cth2 = dr02.dot(dr12);
                    
                    vDisp += (1.0 + 3*cth0*cth1*cth2)/r9;
                }
            }
        }
        vDisp *= C9isoSim/8.00;
        
        double cmDist01 = cm0.Mv1Squared(cm1);
        if (cmDist01 > 900) return vDisp;
        double cmDist12 = cm1.Mv1Squared(cm2);
        if (cmDist12 > 900) return vDisp;
        double cmDist02 = cm0.Mv1Squared(cm2);
        if (cmDist02 > 900) return vDisp;

        double vInd0 = 0;
        for (int i=5; i<7; i++) {
            r0.E(cm0);
            r0.PEa1Tv1(pos[i], or0);
            Vector dr01 = space.makeVector();
            Vector dr02 = space.makeVector();
            dr01.E(0);
            dr02.E(0);
            for (int j=0; j<5; j++){
                int jj = siteID[j];
                r1.E(cm1);
                r1.PEa1Tv1(pos[j], or1);
                r1.ME(r0);
                r1.TE(-1.0);
                double r01ij = Math.sqrt(r1.squared());
                double r301ij = r01ij*r01ij*r01ij;
                dr01.PEa1Tv1(q[jj]/r301ij, r1);
                                
                r2.E(cm2);
                r2.PEa1Tv1(pos[j], or2);
                r2.ME(r0);
                r2.TE(-1.0);
                double r02ik = Math.sqrt(r2.squared());
                double r302ik = r02ik*r02ik*r02ik;
                dr02.PEa1Tv1(q[jj]/r302ik, r2);
            }
            vInd0 += Kelvin.UNIT.toSim(dr01.dot(dr02));
        }
        vInd0 *= -0.5*alphaIsoSim;
        
        double vInd1 = 0;
        for (int i=5; i<7; i++) {            
            r1.E(cm1);
            r1.PEa1Tv1(pos[i], or1);
            Vector dr11 = space.makeVector();
            Vector dr12 = space.makeVector();
            dr11.E(0);
            dr12.E(0);
            for (int j=0; j<5; j++){
                int jj = siteID[j];
                r0.E(cm0);
                r0.PEa1Tv1(pos[j], or0);
                r0.ME(r1);
                r0.TE(-1.0);
                double r01ij = Math.sqrt(r0.squared());
                double r301ij = r01ij*r01ij*r01ij;
                dr11.PEa1Tv1(q[jj]/r301ij, r0);
                                
                r2.E(cm2);
                r2.PEa1Tv1(pos[j], or2);
                r2.ME(r1);
                r2.TE(-1.0);
                double r12ik = Math.sqrt(r2.squared());
                double r312ik = r12ik*r12ik*r12ik;
                dr12.PEa1Tv1(q[jj]/r312ik, r2);
            }
            vInd1 += Kelvin.UNIT.toSim(dr11.dot(dr12));
        }
        vInd1 *= -0.5*alphaIsoSim;
        
        double vInd2 = 0;
        for (int i=5; i<7; i++) {
            r2.E(cm2);
            r2.PEa1Tv1(pos[i], or2);
            Vector dr21 = space.makeVector();
            Vector dr22 = space.makeVector();
            dr21.E(0);
            dr22.E(0);
            for (int j=0; j<5; j++){
                int jj = siteID[j];
                r0.E(cm0);
                r0.PEa1Tv1(pos[j], or0);
                r0.ME(r2);
                r0.TE(-1.0);
                double r02ij = Math.sqrt(r0.squared());
                double r302ij = r02ij*r02ij*r02ij;
                dr21.PEa1Tv1(q[jj]/r302ij, r0);
                                
                r1.E(cm1);
                r1.PEa1Tv1(pos[j], or1);
                r1.ME(r2);
                r1.TE(-1.0);
                double r12ik = Math.sqrt(r1.squared());
                double r312ik = r12ik*r12ik*r12ik;
                dr22.PEa1Tv1(q[jj]/r312ik, r1);
            }
            vInd2 += Kelvin.UNIT.toSim(dr21.dot(dr22));
        }
        vInd2 *= -0.5*alphaIsoSim;
        
        double v = vInd0 + vInd1 + vInd2 + vDisp;
        return v;
    }
}
