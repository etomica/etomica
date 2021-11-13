/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationSolidSuperCut implements PotentialCallback {

    protected final CoordinateDefinition coordinateDefinition;
    protected final Vector drSite0, drSite1, drA, dr, drB, dr0;
    protected final Space space;
    protected double[] pzxySum;
    protected final Vector pTmp1;
    protected double[] sum1, virialSum;
    protected double[] energySum, dadbSum;
    protected double fac1;
    protected Box box;
    protected Boundary boundary;
    protected final double[] r2Cut;

    public PotentialCalculationSolidSuperCut(Space space, CoordinateDefinition coordinateDefinition, double[] cutoffs) {
        this.coordinateDefinition = coordinateDefinition;
        this.space = space;
        this.r2Cut = new double[cutoffs.length];
        for (int i=0; i<r2Cut.length; i++) {
            r2Cut[i] = cutoffs[i]*cutoffs[i];
        }
        drSite0 = space.makeVector();
        drSite1 = space.makeVector();
        drA = space.makeVector();
        dr = space.makeVector();
        drB = space.makeVector();
        dr0 = space.makeVector();
        setBox(coordinateDefinition.getBox());
        sum1 = new double[r2Cut.length];
        virialSum = new double[r2Cut.length];
        pzxySum = new double[r2Cut.length];
        energySum = new double[r2Cut.length];
        dadbSum = new double[r2Cut.length];
        pTmp1 = space.makeVector();
    }

    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {
        IAtomList atoms = coordinateDefinition.getBox().getLeafList();

        IAtom atom0 = atoms.get(i);
        IAtom atom1 = atoms.get(j);
        Vector site0 = coordinateDefinition.getLatticePosition(atom0);
        Vector site1 = coordinateDefinition.getLatticePosition(atom1);

        drSite0.Ev1Mv2(atom0.getPosition(), site0);
        drSite0.ME(dr0);
        boundary.nearestImage(drSite0);
        drSite1.Ev1Mv2(atom1.getPosition(), site1);
        drSite1.ME(dr0);
        boundary.nearestImage(drSite1);

        drA.Ev1Mv2(drSite1, drSite0);

        // dr = drB + drSite1 - drSite0
        drB.Ev1Mv2(dr, drA);
        double r2site = drB.squared();
        double r2 = dr.squared();

        double p1 = -u012[1]/r2*dr.dot(drB)*fac1;
        double dadb = -u012[1]/r2*dr.dot(drA);

        pTmp1.Ea1Tv1(-u012[1] / r2, dr);
        double pzxy = pTmp1.getX(2) * dr.getX(2) - (pTmp1.getX(0) * dr.getX(0) + pTmp1.getX(1) * dr.getX(1)) / 2;

        int n=r2Cut.length;
        for (int k=n-1; k>=0; k--) {
            if (r2Cut[k] < r2site) break;
            if (k == 0 && i*j == 0 && false) {
                System.out.println("PCSSCF "+(i+j)+" "+u012[0]);
            }
            energySum[k] += u012[0];
            sum1[k] += p1;
            virialSum[k] += u012[1];
            pzxySum[k] += pzxy;
            dadbSum[k] += dadb;
        }
    }
    
    public double[] getPressure1() {
        return sum1;
    }
    
    public double[] getVirialSum() {
        return virialSum;
    }

    public double[] getPzxy() {
        return pzxySum;
    }

    public double[] getEnergySum() {
      return energySum;
    }
    
    /**
     * Returns sum of Fi dot ri
     */
    public double[] getDADBSum() {
      return dadbSum;
    }
    
    public void zeroSum() {
        int n = r2Cut.length;
        for (int i=0; i<n; i++) {
            sum1[i] = virialSum[i] = energySum[i] = dadbSum[i] = pzxySum[i] = 0;
        }
        boundary = box.getBoundary();
        int D = boundary.getBoxSize().getD();
        double vol = boundary.volume();
        fac1 = 1.0/(D*vol);

        IAtom atom0 = box.getLeafList().get(0);
        Vector site0 = coordinateDefinition.getLatticePosition(atom0);
        dr0.Ev1Mv2(atom0.getPosition(), site0);
    }
    
    public void setBox(Box box) {
        this.box = box;
        boundary = box.getBoundary();
    }

}
