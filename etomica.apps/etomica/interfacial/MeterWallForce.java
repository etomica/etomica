package etomica.interfacial;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMaster;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.units.Force;

public class MeterWallForce extends DataSourceScalar {

    protected final IPotentialMaster potentialMaster;
    protected final ISpecies topWall;
    protected final PotentialCalculationWallForce pc;
    protected final ISpace space;
    protected final IBox box;
    protected final IteratorDirective id;
    
    public MeterWallForce(ISpace space, IPotentialMaster potentialMaster, IBox box, ISpecies topWall) {
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
        protected IAtomType atomType;
    
        public double getSum() {
            return sum;
        }
        
        public void setAtomType(IAtomType type) {
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
            IVector[] f = potentialSoft.gradient(atoms);
            if (atoms.getAtom(0).getType() == atomType) {
                sum -= f[0].getX(2);
            }
            else {
                sum -= f[1].getX(2);
            }
        }
    }


}
