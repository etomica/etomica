package etomica.meta.wrappers;


import etomica.meta.properties.CollectionElementProperty;
import etomica.meta.SimulationModel;

import java.util.Collection;

/**
 * Wraps a collection of objects, such that the getChildProps gives property instances for each element of the collection.
 */
public class CollectionWrapper extends Wrapper<Collection> {

    public CollectionWrapper(Collection wrapped, SimulationModel simModel) {
        super(wrapped, simModel);

        for(Object el : wrapped) {
            properties.add(new CollectionElementProperty(el));
        }

        childProps.addAll(properties);
    }

}
