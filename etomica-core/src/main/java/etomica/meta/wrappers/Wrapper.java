package etomica.meta.wrappers;

import etomica.meta.SimulationModel;
import etomica.meta.properties.Property;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Surrogate for an instance of an object that is part of a simulation, providing greater control over which
 * fields are made accessible in the simulation model.
 * @param <T> the class of the wrapped object
 */
public abstract class Wrapper<T> {
    protected final T wrapped;
    protected final Class wrappedClass;
    protected final long wrappedId;
    protected final SimulationModel simModel;

    protected final List<Property> properties = new ArrayList<>();
    protected final List<Property> childProps = new ArrayList<>();
    protected final List<Property> valueProps = new ArrayList<>();

    protected final boolean doSerialize;

    //SimulationModel is needed for information about wrapper id, both for this instance and regarding child instances
    public Wrapper(T wrapped, SimulationModel simModel, boolean doSerialize) {
        this.wrapped = wrapped;
        this.wrappedClass = wrapped.getClass();
        this.simModel = simModel;
        wrappedId = simModel.getNewId();
        this.doSerialize = doSerialize;
    }

    public boolean doSerialize() {
        return doSerialize;
    }

    /**
     *
     * @return the class of the object contained in this wrapper
     */
    public Class getWrappedClass() {
        return wrappedClass;
    }

    public long getWrappedId() {
        return wrappedId;
    }

    public List<Property> getProperties() {
        return properties;
    }

    /**
     *
     * @return a Map from the property name to the current value held in the wrapped instance for all properties.
     * Values of child properties are given by their wrapper id
     */
    public Map<String, Object> getValues() {
        Map<String, Object> values = new HashMap<>();
        for (Property prop : valueProps) {
            if (prop.isIndexedProperty()) {
                if(!prop.canCount()) {
                    continue;
                } else {
                    int count = prop.invokeCount();
                    Object[] propValues = new Object[count];
                    for (int i = 0; i < count; i++) {
                        propValues[i] = prop.invokeReader(i);
                    }

                    values.put(prop.getName(), propValues);
                }

            } else {
                values.put(prop.getName(), prop.invokeReader());
            }
        }

        for (Property prop : childProps) {
            if (prop.isIndexedProperty()) {
                if (!prop.canCount()) {
                    continue;
                }
                else {
                    int count = prop.invokeCount();
                    Object[] propValues = new Object[count];
                    for (int i = 0; i < count; i++) {
                        Wrapper wrapper = simModel.getWrapper(prop.invokeReader(i));
                        if (wrapped != null) propValues[i] = wrapper.getWrappedId();
                    }

                    values.put("#" + prop.getName(), propValues);
                }
            }
            else {
                Wrapper wrapper = simModel.getWrapper(prop.invokeReader());
                if (wrapper != null) {
                    values.put("#" + prop.getName(), wrapper.getWrappedId());
                }
            }
        }

        return values;
    }

    /**
     *
     * @return a list of all the properties that are themselves objects
     */
    public List<Property> getChildProps() {
        return childProps;
    }

/*
    @SuppressWarnings("unchecked")
    public Map<String, Wrapper<?>> getChildren() {
        Map<String, Wrapper<?>> children = new HashMap<>();
        for (InstanceProperty prop : childProps) {
            System.out.printf("%s -> %s (%s)\n", this.wrappedClass.getSimpleName(), prop.getName(), prop.getPropertyType().getSimpleName());

            if (prop.isIndexedProperty()) {
                if (!prop.canCount()) {
                    continue;
                }

                int count = prop.invokeCount();
                Object[] propValues = new Object[count];
                for (int i = 0; i < count; i++) {
                    propValues[i] = prop.invokeReader(i);
                }
                children.put(prop.getName(), new ArrayWrapper(propValues));
            } else if (List.class.isAssignableFrom(prop.getPropertyType())) {
                List childList = (List) prop.invokeReader();
                for (Object o : childList) {
                    children.put(prop.getName(), WrapperIndex.getWrapper(o));
                }
            } else {
                children.put(prop.getName(), WrapperIndex.getWrapper(prop.invokeReader()));
            }
        }

        return children;
    }
*/

}
