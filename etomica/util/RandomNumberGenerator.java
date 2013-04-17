package etomica.util;

import java.util.Random;

import etomica.api.IRandom;

/**
 * Random number generator class.  This just wraps Java's own RNG.
 * 
 * @author Andrew Schultz
 */
public class RandomNumberGenerator implements IRandom, java.io.Serializable {

    public RandomNumberGenerator() {
        jRandom = new Random();
    }

    public RandomNumberGenerator(long seed) {
        jRandom = new Random(seed);
    }
    
    // XXX we're just wrapping Java's implementation, which will be fixed-point
    public double nextDouble() {
        return jRandom.nextDouble();
    }
    
    public double nextFixedDouble() {
        return jRandom.nextDouble();
    }

    public int nextInt(int maxInt) {
        return jRandom.nextInt(maxInt);
    }
    
    public double nextGaussian() {
        return jRandom.nextGaussian();
    }
    
    public Random getWrappedRandom() {
        return jRandom;
    }
    
    private static final long serialVersionUID = 1L;
    private final Random jRandom;
}
