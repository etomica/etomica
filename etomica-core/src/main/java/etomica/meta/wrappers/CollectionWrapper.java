package etomica.meta.wrappers;


import etomica.meta.SimulationModel;
import etomica.meta.properties.CollectionElementProperty;

import java.util.Collection;

/**
 * Wraps a collection of objects, such that the getChildProps gives property instances for each element of the collection.
 */
public class CollectionWrapper extends Wrapper<Collection> {

    public CollectionWrapper(Collection wrapped, SimulationModel simModel, boolean doSerialize) {
        super(wrapped, simModel, doSerialize);

        for(Object el : wrapped) {
            properties.add(new CollectionElementProperty(el));
        }

        childProps.addAll(properties);
    }

}
