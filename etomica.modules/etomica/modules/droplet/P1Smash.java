package etomica.modules.droplet;

import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * Gravity-like potential that pushes the molecules toward the center.
 */
public class P1Smash implements PotentialSoft {

    public P1Smash(ISpace space) {
        gradient = new IVectorMutable[1];
        gradient[0] = space.makeVector();
        g = 1;
    }
    
    public void setBox(IBox newBox) {}
    
    public int nBody() {
        return 1;
    }
    
    public void setG(double newG) {
        g = newG;
    }
    
    public double getG() {
        return g;
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public IVector[] gradient(IAtomList atoms) {
        IAtomPositioned a = ((IAtomPositioned)atoms.getAtom(0));
        if (a.getPosition().getX(2) > 0) {
            gradient[0].setX(2, g);
        }
        else {
            gradient[0].setX(2,-g);
        }
        return gradient;
    }

    public double energy(IAtomList atoms) {
        IAtomPositioned a = ((IAtomPositioned)atoms.getAtom(0));
        return Math.abs(a.getPosition().getX(2))*g;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    protected final IVectorMutable[] gradient;
    protected double g;
}
