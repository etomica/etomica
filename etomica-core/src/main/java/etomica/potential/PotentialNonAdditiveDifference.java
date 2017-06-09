package etomica.potential;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePair;
import etomica.space.Space;

public class PotentialNonAdditiveDifference extends PotentialMolecular {

    protected final IPotentialMolecular p2, pFull;
    protected final MoleculePair pair;
    
    public PotentialNonAdditiveDifference(Space space, IPotentialMolecular p2, IPotentialMolecular pFull) {
        super(Integer.MAX_VALUE, space);
        this.p2 = p2;
        this.pFull = pFull;
        pair = new MoleculePair();
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(IMoleculeList molecules) {
        double u = pFull.energy(molecules);
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            pair.atom0 = molecules.getMolecule(i);
            for (int j=i+1; j<molecules.getMoleculeCount(); j++) {
                pair.atom1 = molecules.getMolecule(j);
                u -= p2.energy(pair);
            }
        }
        return u;
    }

    public void setBox(Box box) {
        pFull.setBox(box);
        p2.setBox(box);
    }

}
