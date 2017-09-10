package etomica.meta;

import etomica.meta.properties.Property;
import etomica.meta.wrappers.Wrapper;
import etomica.simulation.Simulation;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SimulationModel {

    private final Map<Class, List<Property>> classes;
    private final Map<Object, Wrapper> wrappers;
    private final Map<Long, Wrapper> wrappersById;
    private long objectCount = 0;
    private final Simulation simulation;

    public SimulationModel(Simulation sim) {
        this.classes = new HashMap<>();
        this.wrappers = new HashMap<>();
        wrappersById = new HashMap<>();
        makeWrapper(sim);
//        this.tree = new PropertyTree(base, base.getWrappedClass().getSimpleName(), classes);
        simulation = sim;
    }

    public Map<Class, List<Property>> getClasses() {
        return classes;
    }

    public long getId(Object o) {
        if(wrappers.containsKey(o)) {
            return wrappers.get(o).getWrappedId();
        }
        else return -1;
    }

    public long getNewId() {
        return objectCount++;
    }

    public Wrapper getWrapper(Object o) {
        return wrappers.get(o);
    }

    public Wrapper getWrapperById(long id) {
        return wrappersById.get(id);
    }

    /**
     * Ensures that a wrapper is defined for the given object and, recursively, all other objects referenced by it.
     * If wrapper for object exists, takes no action; otherwise, creates new wrapper for it and adds to map.
     * @param o the object to be wrapped
     */
    public void makeWrapper(Object o) {

        if(o == null || wrappers.containsKey(o)) return;

        Wrapper<?> newWrapper = WrapperIndex.getWrapper(o, this);
        if (newWrapper.doSerialize()) {
            wrappers.put(o, newWrapper);
            wrappersById.put(newWrapper.getWrappedId(), newWrapper);
        }

        List<Property> childProps = newWrapper.getChildProperties();
        for (Property prop : childProps) {
            if(prop.isIndexedProperty()) {
                if(!prop.canCount()) continue;
                int count = prop.invokeCount();
                for(int i=0; i<count; i++) {
                    makeWrapper(prop.invokeReader(i));
                }
            } else {
                makeWrapper(prop.invokeReader());
            }
        }

    }

    public Collection<Wrapper> allWrappers() {
        return wrappers.values();
    }

    public Simulation getSimulation() {
        return simulation;
    }
}
