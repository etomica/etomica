package etomica.util;

/**
 * Interface for random number generation.
 *
 * @author Andrew Schultz
 */
public interface IRandom {

    /**
     * Returns a pseudorandom double, uniformly distributed between 
     * 0.0 (inclusive) and 1.0 (exclusive).
     */
    public double nextDouble();

    /**
     * Returns a pseudorandom integer, uniformly distributed between 
     * 0 (inclusive) and maxInt (exclusive).
     */
    public int nextInt(int maxInt);

    /**
     * Returns a pseudorandom double, Gaussian ("normally") distributed
     * value with mean 0.0 and standard deviation 1.0.
     */
    public double nextGaussian();

}