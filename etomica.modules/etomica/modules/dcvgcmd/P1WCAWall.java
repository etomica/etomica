package etomica.modules.dcvgcmd;

import etomica.api.IAtomSet;
import etomica.api.IAtomPositioned;
import etomica.api.IVector;

import etomica.EtomicaInfo;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * 1-D potential that has a WCA form in the Z direction.
 */

public class P1WCAWall extends Potential1 implements PotentialSoft {

    private static final long serialVersionUID = 1L;
    protected final IVector[] gradient;
    protected double sigma;
    protected double epsilon;
    protected double cutoff;
    protected int wallDim;

    public P1WCAWall(Space space, int wallDim) {
        this(space, wallDim, 1.0, 1.0);
    }

    public P1WCAWall(Space space, int wallDim, double sigma, double epsilon) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setWallDim(wallDim);
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

    public double energy(IAtomSet atom) {
        IVector dimensions = boundary.getDimensions();
        double rz = ((IAtomPositioned)atom.getAtom(0)).getPosition().x(wallDim);
        double dzHalf = 0.5 * dimensions.x(wallDim);
        return energy(sigma + dzHalf + rz) + energy(sigma + dzHalf - rz);
    }

    private double energy(double r) {
        if (r > cutoff) {
            return 0;
        }
        double rr = sigma / r;
        double r2 = rr * rr;
        double r6 = r2 * r2 * r2;
        return 4 * epsilon * r6 * (r6 - 1.0) + epsilon;
    }

    private double gradient(double r) {
        if (r > cutoff) {
            return 0;
        }
        double rr = sigma / r;
        double r2 = rr * rr;
        double r6 = r2 * r2 * r2;
        return -48 * epsilon * r6 * (r6 - 0.5);
    }

    public IVector[] gradient(IAtomSet atom) {
        IVector dimensions = boundary.getDimensions();
        double rz = ((IAtomPositioned)atom.getAtom(0)).getPosition().x(wallDim);
        double dzHalf = 0.5 * dimensions.x(wallDim);
        double gradz = gradient(sigma + rz + dzHalf) - gradient(sigma + dzHalf - rz);
        gradient[0].setX(wallDim, gradz);
        return gradient;
    }
    
    public IVector[] gradient(IAtomSet atom, Tensor pressureTensor) {
        return gradient(atom);
    }
    
    public double virial(IAtomSet atoms) {
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

    public int getWallDim() {
        return wallDim;
    }

    public void setWallDim(int wallDim) {
        this.wallDim = wallDim;
    }
}