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
import java.lang.Math;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationSolidSuper implements PotentialCalculation {
        
    protected final CoordinateDefinition coordinateDefinition;
    protected final Vector atom0DistanceFromLatticeSite0, atom1DistanceFromLatticeSite1, distanceBetweenAtoms, dr, distanceBetweenLatticeSites, dr0;
    protected final Space space;
    protected double virialSum;
    protected double energySum, dadbSum;
    protected double fac1;
    protected Box box;
    protected Boundary boundary;
    protected boolean doD2;
    protected double d2sum, normalModeSum;

    public PotentialCalculationSolidSuper(Space space, CoordinateDefinition coordinateDefinition) {
        this.coordinateDefinition = coordinateDefinition;
        this.space = space;
        atom0DistanceFromLatticeSite0 = space.makeVector();
        atom1DistanceFromLatticeSite1 = space.makeVector();
        distanceBetweenAtoms = space.makeVector();
        dr = space.makeVector();
        distanceBetweenLatticeSites = space.makeVector();
        setBox(coordinateDefinition.getBox());
        dr0 = space.makeVector();
    }

    public void setDoSecondDerivative(boolean doD2) {
        this.doD2 = doD2;
    }
    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        IAtom atom0 = atoms.get(0);
        IAtom atom1 = atoms.get(1);
        energySum += potential.energy(atoms);
        Vector site0 = coordinateDefinition.getLatticePosition(atom0);
        Vector site1 = coordinateDefinition.getLatticePosition(atom1);

        // TODO: These names are not the clearest. What's 'drB'?
        distanceBetweenLatticeSites.Ev1Mv2(site1, site0);
        boundary.nearestImage(distanceBetweenLatticeSites);

        atom0DistanceFromLatticeSite0.Ev1Mv2(atom0.getPosition(), site0);
        atom0DistanceFromLatticeSite0.ME(dr0);
        boundary.nearestImage(atom0DistanceFromLatticeSite0);
        atom1DistanceFromLatticeSite1.Ev1Mv2(atom1.getPosition(), site1);
        atom1DistanceFromLatticeSite1.ME(dr0);        // Decrements the above vector by 'dr0'.
        boundary.nearestImage(atom1DistanceFromLatticeSite1);

        distanceBetweenAtoms.Ev1Mv2(atom1DistanceFromLatticeSite1, atom0DistanceFromLatticeSite0);
        // dr = drB + drSite1 - drSite0
        dr.Ev1Pv2(distanceBetweenLatticeSites, distanceBetweenAtoms);        // Adding together the components of the vector.
        boundary.nearestImage(dr);
        double r2 = dr.squared();

        Potential2SoftSpherical potentialSoft = (Potential2SoftSpherical)potential;
        
        double du = potentialSoft.du(r2);       // This returns the dudr * r2. If you divide this by r2, you will obtain the force.
        virialSum += du;

        double dot = dr.dot(distanceBetweenAtoms);       // dot product between dr and drA. dr is the vector describing the distance between the two atoms. drA is the distance between the sites.
        dadbSum -= du/r2*dot;       // NOTE: The du/r2 is equivalent to the force acting on the atom of interest.


        // dr is vector between i and j (This is the vector describing the distance between atoms.)
        // drA is vector between sites i and j (I believe that this is referring to the lattice sites and not the atoms).
        if (doD2) {
            // TODO: This code can be refactored into a function?
            double d2u = potentialSoft.d2u(r2);
            double dud2u = (du - d2u) / (r2 * r2);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double der2 = dr.getX(i) * dr.getX(j) * dud2u;
                    if (i == j) der2 -= du / r2;
                    d2sum += der2 * (2 * atom0DistanceFromLatticeSite0.getX(i) * atom1DistanceFromLatticeSite1.getX(j)
                            - atom0DistanceFromLatticeSite0.getX(i) * atom0DistanceFromLatticeSite0.getX(j) - atom1DistanceFromLatticeSite1.getX(i) * atom1DistanceFromLatticeSite1.getX(j));
                }
            }
        }
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

    public double getNormalModeSum() {
        return normalModeSum;
    }
    
    public double getD2Sum() {
        return d2sum;
    }

    public void zeroSum() {
        virialSum = energySum = dadbSum = d2sum = normalModeSum = 0;

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
