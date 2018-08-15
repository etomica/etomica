package etomica.mappedDensity;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class P1Sine implements IPotentialAtomic, PotentialSoft {

    private final int n;
    private final double temperature;
    private double L;
    private final Vector[] gradient;

    public P1Sine(Space space, int n, double temperature) {
        this.n = n;
        this.temperature = temperature;
        gradient = new Vector[]{space.makeVector()};
    }

    @Override
    public double virial(IAtomList atoms) {
        return 0;
    }

    @Override
    public Vector[] gradient(IAtomList atoms) {
        double z = atoms.get(0).getPosition().getX(2);
        // temperature/(2+sin())*cos()*(2 PI n / L)
        double arg = 2 * Math.PI * n / L;
        gradient[0].setX(2, -temperature / (2 + Math.sin(arg * z)) * Math.cos(arg * z) * arg);
//        System.out.println(z+" "+energy(atoms)+" "+gradient[0].getX(2));
        return gradient;
    }

    @Override
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        throw new RuntimeException("not implemented");
    }

    @Override
    public double energy(IAtomList atoms) {
        double z = atoms.get(0).getPosition().getX(2);
        return -temperature * Math.log(2 + Math.sin(2 * Math.PI * n * z / L));
    }

    @Override
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public void setBox(Box box) {
        L = box.getBoundary().getBoxSize().getX(2);
    }

    @Override
    public int nBody() {
        return 1;
    }
}
