/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationSolidSuper implements PotentialCalculation {
        
    protected final CoordinateDefinition coordinateDefinition;
    protected final Vector drSite0, drSite1, drA, dr, drB, dr0;
    protected final Space space;
    protected double pSum, virialSum;
    protected final Vector pSumXYZ, pTmp;
    protected double energySum, dadbSum;
    protected double fac1, fac2;
    protected Box box;
    protected Boundary boundary;
    protected boolean doD2;
    protected double d2sum;
    protected double pHarmonic, temperature;
    
    public PotentialCalculationSolidSuper(Space space, CoordinateDefinition coordinateDefinition) {
        this.coordinateDefinition = coordinateDefinition;
        this.space = space;
        drSite0 = space.makeVector();
        drSite1 = space.makeVector();
        drA = space.makeVector();
        dr = space.makeVector();
        drB = space.makeVector();
        setBox(coordinateDefinition.getBox());
        pSumXYZ = space.makeVector();
        pTmp = space.makeVector();
        dr0 = space.makeVector();
    }
    
    public void setDoSecondDerivative(boolean doD2) {
        this.doD2 = doD2;
    }
    
    public void setPHarmonic(double pHarmonic, double temperature) {
        this.pHarmonic = pHarmonic;
        this.temperature = temperature;
    }

    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        IAtom atom0 = atoms.get(0);
        IAtom atom1 = atoms.get(1);
//        if (atom0.getLeafIndex() >0 || atom1.getLeafIndex() > 1) return;
        energySum += potential.energy(atoms);
//        int i0 = atoms.getAtom(0).getLeafIndex();
//        int i1 = atoms.getAtom(1).getLeafIndex();
//        if (i0*i1 > 0 || i0+i1 != 8) return;
//        dr.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
//        boundary.nearestImage(dr);
//        double r = Math.sqrt(dr.squared()); // distance between atoms
        Vector site0 = coordinateDefinition.getLatticePosition(atom0);
        Vector site1 = coordinateDefinition.getLatticePosition(atom1);
//        if (debug) {
//            System.out.println(site0+" "+atom0.getPosition());
//            System.out.println(site1+" "+atom1.getPosition());
//        }

        drB.Ev1Mv2(site1, site0);
        boundary.nearestImage(drB);
        double r2site = drB.squared();

        drSite0.Ev1Mv2(atom0.getPosition(), site0);
        drSite0.ME(dr0);
        boundary.nearestImage(drSite0);
        drSite1.Ev1Mv2(atom1.getPosition(), site1);
        drSite1.ME(dr0);
        boundary.nearestImage(drSite1);

        drA.Ev1Mv2(drSite1, drSite0);
        // dr = drB + drSite1 - drSite0
        dr.Ev1Pv2(drB, drA);
        boundary.nearestImage(dr);
        double r2 = dr.squared();

        Potential2SoftSpherical potentialSoft = (Potential2SoftSpherical)potential;
        
        double du = potentialSoft.du(r2);
        pSum -= du/r2*dr.dot(drB);
        pTmp.Ea1Tv1(-du/r2, dr);
        pTmp.TE(drB);
        pSumXYZ.PE(pTmp);
        virialSum += du;

//        dadbSum -= f[0].dot(drSite0);
//        dadbSum -= f[1].dot(drSite1);
//        dadbSum -= f[1].dot(drA);
//        System.out.println(f[1].getX(0)+" "+Math.sqrt(f[1].squared()/r2)*dr.getX(0)+" "+potentialSoft.du(r2)/r2*dr.getX(0));
        double dot = dr.dot(drA);
        dadbSum -= du/r2*dot;
        pSum -= fac2*du/r2*dot;
        pTmp.Ea1Tv1(-fac2*du/r2, dr);
        pTmp.TE(drA);
        pSumXYZ.PE(pTmp);
//        if (atoms.getAtom(0).getLeafIndex()==0 && atoms.getAtom(1).getLeafIndex()==8) {
//            System.out.println("new 0 8 "+(dr.dot(f[0])+drA.dot(f[0])*fac1));
//        }
        //dr.PEa1Tv1(fac2, drA);
        //if (debug) System.out.print(drA+"\n "+dr+"\n");
        
        //sum += dr.dot(f[0]);
        
        // dr is vector between i and j
        // drA is vector between sites i and j
        if (doD2) {
            double r2A = drA.squared();
            double d2u = potentialSoft.d2u(r2);
//            System.out.println("rij "+Math.sqrt(r2));
//            System.out.println("dudr "+du/Math.sqrt(r2));
//            System.out.println("dot "+dot);
//            System.out.println();
//            System.out.println("drdT "+dot/Math.sqrt(r2));
//            System.out.println("d2udrdT "+d2u*dot/(r2*Math.sqrt(r2)));
//            System.out.println("ddotdT "+(r2A+dot));
//            System.out.println();
            d2sum -= d2u*dot*dot/(r2*r2) - du*dot*dot/(r2*r2)
                   + du*(r2A + dot)/r2;
//            System.out.println("gradient");
//            System.out.println(du*dr.getX(0)/r2+" "+du*dr.getX(1)/r2+" "+du*dr.getX(2)/r2);
//            double dud2u = (d2u + -du)/(r2*r2);
//            for (int i=0; i<3; i++) {
//                d2sum += du*drA.getX(i)*drA.getX(i);
//                for (int j=0; j<3; j++) {
//                    d2sum += dud2u * dr.getX(i)*dr.getX(j)*drA.getX(i)*drA.getX(j);
//                }
//            }
        }
    }
    
    public double getPressureSum() {
        return pSum;
    }
    
    public Vector getPressureSumXYZ() {
        return pSumXYZ;
    }
    
    public double getVirialSum() {
        return virialSum;
    }

    public double getEnergySum() {
        return energySum;
    }
    
    /**
     * Returns sum of Fi dot ri
     */
    public double getDADBSum() {
      return dadbSum;
    }
    
    public double getD2Sum() {
        return d2sum;
    }
    
    public void zeroSum() {
        pSum = virialSum = energySum = dadbSum = d2sum = 0;
        pSumXYZ.E(0);
        boundary = box.getBoundary();
        int D = boundary.getBoxSize().getD();
        double vol = boundary.volume();
        fac1 = 1.0/(D*vol);
        int N = box.getMoleculeList().size();
        fac2 = (-1/vol + pHarmonic/temperature)/(D*N-D);

        IAtom atom0 = box.getLeafList().get(0);
        Vector site0 = coordinateDefinition.getLatticePosition(atom0);
        dr0.Ev1Mv2(atom0.getPosition(), site0);
    }
    
    public void setBox(Box box) {
        this.box = box;
        boundary = box.getBoundary();
    }
}
