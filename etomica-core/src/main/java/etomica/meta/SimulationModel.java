package etomica.meta;

import etomica.meta.properties.Property;
import etomica.meta.wrappers.Wrapper;
import etomica.simulation.Simulation;

import java.util.*;

public class SimulationModel {

    private final Map<Class, List<Property>> classes;
    private final Map<Object, Wrapper> wrappers;
    private long objectCount = 0;

    public SimulationModel(Simulation sim) {
        this.classes = new HashMap<>();
        this.wrappers = new HashMap<>();
        makeWrapper(sim);
//        this.tree = new PropertyTree(base, base.getWrappedClass().getSimpleName(), classes);
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

    /**
     * Ensures that a wrapper is defined for the given object and, recursively, all other objects referenced by it.
     * If wrapper for object exists, takes no action; otherwise, creates new wrapper for it and adds to map.
     * @param o the object to be wrapped
     */
    public void makeWrapper(Object o) {

        if(o == null || wrappers.containsKey(o)) return;

        Wrapper newWrapper = WrapperIndex.getWrapper(o, this);
        if (newWrapper.doSerialize()) wrappers.put(o, newWrapper);

        List<Property> childProps = newWrapper.getChildProps();
        for (Property prop : childProps) {
            if(prop.isIndexedProperty()) {
                if(!prop.canCount()) continue;
                ArrayList list = new ArrayList();
                int count = prop.invokeCount();
                for(int i=0; i<count; i++) {
                    list.add(prop.invokeReader(i));
                }
                makeWrapper(list);
            } else {
                makeWrapper(prop.invokeReader());
            }
        }

    }

    public Collection<Wrapper> allWrappers() {
        return wrappers.values();
    }

}
