/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Force;

public class MeterWallForce extends DataSourceScalar {

    protected final PotentialMaster potentialMaster;
    protected final ISpecies topWall;
    protected final PotentialCalculationWallForce pc;
    protected final Space space;
    protected final Box box;
    protected final IteratorDirective id;
    
    public MeterWallForce(Space space, PotentialMaster potentialMaster, Box box, ISpecies topWall) {
        super("Force", Force.DIMENSION);
        this.space = space;
        this.box = box;
        this.potentialMaster = potentialMaster;
        this.topWall = topWall;
        pc = new PotentialCalculationWallForce();
        pc.setAtomType(topWall.getAtomType(0));
        id = new IteratorDirective(null);
    }

    public double getDataAsScalar() {
        IMoleculeList topWallMolecules = box.getMoleculeList(topWall);
        pc.reset();
        for (int i=0; i<topWallMolecules.getMoleculeCount(); i++) {
            IAtom wallAtom = topWallMolecules.getMolecule(i).getChildList().getAtom(0);
            id.setTargetAtom(wallAtom);
            potentialMaster.calculate(box, id, pc);
        }
        return pc.getSum();
    }

    /**
     * Sums the force on each iterated atom and adds it to the integrator agent
     * associated with the atom.
     */
    public static class PotentialCalculationWallForce implements PotentialCalculation {

        protected double sum;
        protected AtomType atomType;
    
        public double getSum() {
            return sum;
        }

        public void setAtomType(AtomType type) {
            atomType = type;
        }
        
        /**
         * Re-zeros the force vector.
         *
         */
        public void reset(){
            sum = 0;
        }
    
        /**
         * Adds forces due to given potential acting on the atoms produced by the iterator.
         * Implemented for only 1- and 2-body potentials.
         */
        public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
            if (atoms.getAtomCount()==1) return;
            PotentialSoft potentialSoft = (PotentialSoft)potential;
            Vector[] f = potentialSoft.gradient(atoms);
            if (atoms.getAtom(0).getType() == atomType) {
                sum -= f[0].getX(2);
            }
            else {
                sum -= f[1].getX(2);
            }
        }
    }


}
