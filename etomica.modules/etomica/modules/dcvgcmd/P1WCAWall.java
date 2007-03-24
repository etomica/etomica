package etomica.modules.dcvgcmd;

import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * 1-D potential that has a WCA form in the Z direction.
 */

public class P1WCAWall extends Potential1 implements PotentialSoft {

    private static final long serialVersionUID = 1L;
    private final IVector[] gradient;
    private double sigma;
    private double epsilon;
    private double cutoff;

    public P1WCAWall(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().atomSize, sim.getDefaults().potentialWell);
    }

    public P1WCAWall(Space space, double sigma, double epsilon) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        gradient = new IVector[1];
        gradient[0] = space.makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo(
                "WCA LJ Potential in the Z-Coordinate");
        return info;
    }
    
    public double getRange() {
        return cutoff;
    }

    public double energy(AtomSet atom) {
        AtomLeaf a = (AtomLeaf) atom;
        IVector dimensions = boundary.getDimensions();
        double rz = a.getCoord().getPosition().x(2);
        double dzHalf = 0.5 * dimensions.x(2);
        return energy(dzHalf + rz) + energy(dzHalf - rz);
    }//end of energy

    private double energy(double r) {
        if (r < cutoff) {
            double rr = sigma / r;
            double r2 = rr * rr;
            double r6 = r2 * r2 * r2;
            return 4 * epsilon * r6 * (r6 - 1.0) + epsilon;
        }
        return 0;
    }

    private double gradient(double r) {
        if (r < cutoff) {
            double rr = sigma / r;
            double r2 = rr * rr;
            double r6 = r2 * r2 * r2;
            return -48 * epsilon * r6 * (r6 - 0.5);
        }
        return 0;
    }

    public IVector[] gradient(AtomSet atom) {
        AtomLeaf a = (AtomLeaf) atom;
        IVector dimensions = boundary.getDimensions();
        double rz = a.getCoord().getPosition().x(2);
        double dzHalf = 0.5 * dimensions.x(2);
        double gradz = gradient(rz + dzHalf) - gradient(dzHalf - rz);
        gradient[0].setX(2, gradz);
        return gradient;
    }
    
    public IVector[] gradient(AtomSet atom, Tensor pressureTensor) {
        return gradient(atom);
    }
    
    public double virial(AtomSet atoms) {
        return 0.0;
    }

    /**
     * Returns the radius.
     * 
     * @return double
     */
    public double getSigma() {
        return sigma;
    }

    /**
     * Sets the radius.
     * 
     * @param radius
     *            The radius to set
     */
    public void setSigma(double radius) {
        this.sigma = radius;
        cutoff = radius * Math.pow(2, 1. / 6.);
    }

    /**
     * @return Returns the epsilon.
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * @param epsilon
     *            The epsilon to set.
     */
    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }
}