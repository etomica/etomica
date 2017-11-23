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

    public T getWrapped() {
        return wrapped;
    }

    public long getWrappedId() {
        return wrappedId;
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
                values.put((prop.canWrite() ? "$" : "") + prop.getName(), prop.invokeReader());
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
                        if (wrapper != null) propValues[i] = wrapper.getWrappedId();
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
    public List<Property> getChildProperties() {
        return childProps;
    }

    public List<Property> getValueProperties() {
        return valueProps;
    }

}
