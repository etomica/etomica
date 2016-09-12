package etomica.interfacial;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;

public class P1Wall extends Potential1 implements PotentialSoft {

    protected final double spring, springPosition;
    protected final double gSat;
    protected final IVectorMutable[] grad;
    
    public P1Wall(ISpace space, double spring, double springPosition, double gSat) {
        super(space);
        this.spring = spring;
        this.springPosition = springPosition;
        this.gSat = gSat;
        grad = new IVectorMutable[]{space.makeVector()};
    }

    public double energy(IAtomList atoms) {
        double dz = atoms.getAtom(0).getPosition().getX(2)-springPosition;
        double uSpring = 0.5*spring*dz*dz;
        return uSpring + gSat*dz;
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    public IVector[] gradient(IAtomList atoms) {
        double dz = atoms.getAtom(0).getPosition().getX(2)-springPosition;
        grad[0].setX(2, gSat + spring*dz);
        return grad;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

}
