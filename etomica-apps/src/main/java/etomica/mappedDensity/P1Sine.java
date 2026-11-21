package etomica.mappedDensity;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.potential.IPotential1;
import etomica.space.Vector;
 /**
 *
 * sinusoidal external field
 *
 */
public class P1Sine implements IPotential1 {
//defines the external potential that will be added
    private final int n;
    private final double temperature;
    private double L;

    public P1Sine(int n, double temperature, Box box) {
        this.n = n;
        this.temperature = temperature;
        L = box.getBoundary().getBoxSize().getX(2);
    }

    @Override
    public double u(IAtom atom) {
        double z = atom.getPosition().getX(2);
        return -temperature * Math.log(2 + Math.sin(2 * Math.PI * n * z / L));
    }

    @Override
    public double udu(IAtom atom, Vector f) {
        double z = atom.getPosition().getX(2);
        // temperature/(2+sin())*cos()*(2 PI n / L)
        double arg = 2 * Math.PI * n / L;
        f.setX(2, f.getX(2) - (temperature / (2 + Math.sin(arg * z))) * Math.cos(arg * z) * arg);
        return -temperature * Math.log(2 + Math.sin(2 * Math.PI * n * z / L));
    }
}
