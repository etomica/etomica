package etomica.modules.sam;

import etomica.EtomicaInfo;
import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IVectorMutable;
import etomica.api.IVector;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * 1-D potential that has a WCA form in the Z direction.
 */

public class P1WCAWall extends Potential1 implements PotentialSoft {

    private static final long serialVersionUID = 1L;
    protected final IVectorMutable[] gradient;
    protected double sigma, sigma2;
    protected double epsilon;
    protected double cutoff, cutoff2;
    protected int wallDim;
    protected double wallPosition;

    public P1WCAWall(ISpace space, int wallDim) {
        this(space, wallDim, 1.0, 1.0);
    }

    public P1WCAWall(ISpace space, int wallDim, double sigma, double epsilon) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
        setWallDim(wallDim);
        gradient = new IVectorMutable[1];
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

    public double energy(IAtomList atom) {
        double rz = ((IAtomPositioned)atom.getAtom(0)).getPosition().x(wallDim) - wallPosition;
        return energy(rz*rz);
    }

    private double energy(double r2) {
        if (r2 > cutoff2) {
            return 0;
        }
        double s2 = sigma2 / r2;
        double s6 = s2 * s2 * s2;
        return 4 * epsilon * s6 * (s6 - 1.0) + epsilon;
    }

    private double gradient(double r2) {
        if (r2 > cutoff2) {
            return 0;
        }
        double s2 = sigma2 / r2;
        double s6 = s2 * s2 * s2;
        return -48 * epsilon * s6 * (s6 - 0.5);
    }

    public IVector[] gradient(IAtomList atom) {
        double rz = ((IAtomPositioned)atom.getAtom(0)).getPosition().x(wallDim) - wallPosition;
        double gradz = gradient(rz*rz);
        gradient[0].setX(wallDim, rz > 0 ? gradz : -gradz);
        return gradient;
    }
    
    public IVector[] gradient(IAtomList atom, Tensor pressureTensor) {
        return gradient(atom);
    }
    
    public double virial(IAtomList atoms) {
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
        sigma2 = sigma*sigma;
        cutoff2 = sigma2 * Math.pow(2, 1. / 3.);
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

    public double getWallPosition() {
        return wallPosition;
    }

    public void setWallPosition(double newWallPosition) {
        wallPosition = newWallPosition;
    }
}