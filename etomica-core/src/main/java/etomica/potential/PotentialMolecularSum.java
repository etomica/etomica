package etomica.potential;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;

/**
 * Molecular potential class that simply sums up contributions from multiple
 * (molecular) potentials.
 * 
 * @author Andrew Schultz
 */
public class PotentialMolecularSum implements IPotentialMolecular {

    protected final IPotentialMolecular[] p;
    
    public PotentialMolecularSum(IPotentialMolecular[] p) {
        this.p = p;
    }

    public double getRange() {
        double r = 0;
        for (int i=0; i<p.length; i++) {
            if (r < p[i].getRange()) r = p[i].getRange();
        }
        return r;
    }

    public void setBox(Box box) {
        for (int i=0; i<p.length; i++) {
            p[i].setBox(box);
        }
    }

    public int nBody() {
        return p[0].nBody();
    }

    public double energy(IMoleculeList molecules) {
        double sum = 0;
        for (int i=0; i<p.length; i++) {
            sum += p[i].energy(molecules);
        }
        return sum;
    }

}
