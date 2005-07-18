/*
 * Created on Jul 18, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica;

import java.lang.reflect.Constructor;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class SpeciesSignature implements Comparable {

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
