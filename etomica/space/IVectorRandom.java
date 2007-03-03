package etomica.space;


/**
 * Vector interface for vectors having random number methods 
 */
public interface IVectorRandom extends IVector {

    /**
     * Assigns this vector to equal a point chosen randomly on the 
     * surface of a unit sphere.
     */
    public void setRandomSphere();

    /**
     * Assigns each component to (its own) random value between -0.5 and + 0.5.
     */
    public void setRandomCube();

    /**
     * Assigns this vector to equal a point chosen randomly in the volume
     * of a unit spheres.
     */
    public void setRandomInSphere();
}