package etomica.overlap;

/**
 * Class that can return values of alpha.
 * 
 * @author Andrew Schultz
 */
public interface AlphaSource {
    
    /**
     * Returns the number of possible alpha values.
     */
    public int getNumAlpha();

    /**
     * Returns the ith alpha value (i=0..(n-1))
     */
    public double getAlpha(int i);
}