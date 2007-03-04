package etomica.space;

import etomica.util.IRandom;

/**
 * Vector interface for vectors having random number methods 
 */
public interface IVectorRandom extends IVector {

    /**
     * Assigns this vector to equal a point chosen randomly on the 
     * surface of a unit sphere.
     */
    public void setRandomSphere(IRandom random);

    /**
     * Assigns each component to (its own) random value between -0.5 and + 0.5.
     */
    public void setRandomCube(IRandom random);

    /**
     * Assigns this vector to equal a point chosen randomly in the volume
     * of a unit spheres.
     */
    public void setRandomInSphere(IRandom random);
}