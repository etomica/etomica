package etomica.util;

import java.util.Random;

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
    
    public double nextDouble() {
        return jRandom.nextDouble();
    }
    
    public int nextInt(int maxInt) {
        return jRandom.nextInt(maxInt);
    }
    
    public double nextGaussian() {
        return jRandom.nextGaussian();
    }
    
    private static final long serialVersionUID = 1L;
    private final Random jRandom;
}
