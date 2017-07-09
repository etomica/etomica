package etomica.meta.wrappers;

import com.fasterxml.jackson.core.JsonProcessingException;
import etomica.meta.InstanceProperty;
import etomica.meta.WrapperIndex;
import etomica.meta.annotations.IgnoreProperty;
import etomica.simulation.prototypes.HSMD2D;

import java.beans.IndexedPropertyDescriptor;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.lang.reflect.Method;
import java.util.*;

public class Wrapper<T> {
    protected final T wrapped;
    protected final Class wrappedClass;

    protected final List<InstanceProperty> properties = new ArrayList<>();


    public Wrapper(T wrapped) {
        this.wrapped = wrapped;
        this.wrappedClass = wrapped.getClass();

        try {
            Arrays.stream(Introspector.getBeanInfo(wrappedClass).getPropertyDescriptors())
                    .filter(propertyDescriptor -> !propertyDescriptor.getName().equalsIgnoreCase("class"))
                    .filter(propertyDescriptor -> propertyDescriptorMethod(propertyDescriptor) != null)
                    .filter(propertyDescriptor -> !propertyDescriptorMethod(propertyDescriptor).isAnnotationPresent(IgnoreProperty.class))
                    .forEach(propertyDescriptor -> properties.add(new InstanceProperty(wrapped, propertyDescriptor)));
        } catch (IntrospectionException e) {
            e.printStackTrace();
        }
    }

    public List<InstanceProperty> getProperties() {
        return properties;
    }

    @SuppressWarnings("unchecked")
    public List<Wrapper> getChildren() {
        List<Wrapper> children = new ArrayList<>();
        for(InstanceProperty prop : properties) {
            if(prop.isIndexedProperty()) {
                for (int i = 0; i < prop.invokeCount(); i++) {
                    children.add(WrapperIndex.getWrapper(prop.invokeReader(i)));
                }
            } else {
                children.add(WrapperIndex.getWrapper(prop.invokeReader()));
            }
        }

        return children;
    }

    private static Method propertyDescriptorMethod(PropertyDescriptor pd) {
        if(pd instanceof IndexedPropertyDescriptor) {
            return ((IndexedPropertyDescriptor) pd).getIndexedReadMethod();
        } else {
            return pd.getReadMethod();
        }
    }

    public static void main(String[] args) throws JsonProcessingException {
        SimulationWrapper simWrapper = new SimulationWrapper(new HSMD2D());
        List<InstanceProperty> props = simWrapper.getProperties();
        List<Wrapper> tree = simWrapper.getChildren();
        System.out.println(props);
    }

}
