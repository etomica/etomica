package etomica.meta.wrappers;

import etomica.meta.SimulationModel;
import etomica.meta.properties.CollectionElementProperty;

/**
 * Wraps an array of objects, such that the getChildProperties gives property instances for each element of the collection.
 */
public class ArrayWrapper extends Wrapper<Object[]> {

    public ArrayWrapper(Object[] wrapped, SimulationModel simModel, boolean doSerialize) {
        super(wrapped, simModel, doSerialize);

        for(Object el : wrapped) {
            childProps.add(new CollectionElementProperty(el));
        }

    }
}
