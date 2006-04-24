package etomica.species;

import java.lang.reflect.Constructor;

/**
 * Class containing minimal information on how to recreate a Species instance. 
 */
public class SpeciesSignature implements Comparable {

    /**
     * Creates an instance for the given name, constructor and constructor parameters.
     * The actual constructor should take a Simulation followed by the given constructor 
     * parameters.
     */
    public SpeciesSignature(String speciesName, Constructor speciesConstructor, Object[] constructorParameters) {
        name = speciesName;
        constructor = speciesConstructor;
        parameters = (Object[])constructorParameters.clone();
    }
    
    public int compareTo(Object signature) {
        if (constructor != ((SpeciesSignature)signature).constructor ||
                java.util.Arrays.equals(parameters,((SpeciesSignature)signature).parameters)) {
            return 1;
        }
        return 0;
    }
    
    public final String name;
    public final Constructor constructor;
    public final Object[] parameters;
}
