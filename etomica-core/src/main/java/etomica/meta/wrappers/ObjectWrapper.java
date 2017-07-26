package etomica.meta.wrappers;

import etomica.meta.SimulationModel;
import etomica.meta.annotations.IgnoreProperty;
import etomica.meta.properties.ArrayProperty;
import etomica.meta.properties.InstanceProperty;
import etomica.meta.properties.Property;
import etomica.meta.properties.VectorProperty;
import etomica.space.Vector;
import org.apache.commons.lang3.reflect.MethodUtils;

import java.beans.IndexedPropertyDescriptor;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.lang.annotation.Annotation;
import java.lang.reflect.Method;
import java.util.Arrays;

public class ObjectWrapper<T> extends Wrapper<T> {

    public ObjectWrapper(T wrapped, SimulationModel simModel, boolean doSerialize) {
        super(wrapped, simModel, doSerialize);

        try {
            Arrays.stream(Introspector.getBeanInfo(wrappedClass).getPropertyDescriptors())
                    .filter(propertyDescriptor -> !propertyDescriptor.getName().equalsIgnoreCase("class"))
                    .filter(propertyDescriptor -> propertyDescriptorMethod(propertyDescriptor) != null)
                    .filter(propertyDescriptor -> !hasAnnotation(propertyDescriptorMethod(propertyDescriptor), IgnoreProperty.class))
                    .filter(propertyDescriptor -> !propertyDescriptor.getName().toLowerCase().endsWith("dimension"))
                    .forEach(propertyDescriptor -> properties.add(makeProperty(wrapped, propertyDescriptor)));

            for (Property p : properties) {
                if (p.isValueProperty()) {
                    valueProps.add(p);
                } else {
                    childProps.add(p);
                }
            }
        } catch (IntrospectionException e) {
            e.printStackTrace();
        }
    }

    private static Property makeProperty(Object o, PropertyDescriptor propertyDescriptor) {
        Class propertyType = propertyDescriptor.getPropertyType();
        if (propertyType != null && Vector.class.isAssignableFrom(propertyType)) {
            return new VectorProperty(o, propertyDescriptor);
        }
        if (!(propertyDescriptor instanceof IndexedPropertyDescriptor) && propertyType.isArray() &&
                !(propertyType.getComponentType().isPrimitive() || propertyType.getComponentType().equals(String.class))) {
            return new ArrayProperty(o, propertyDescriptor);
        }
        return new InstanceProperty(o, propertyDescriptor);
    }

    private static Method propertyDescriptorMethod(PropertyDescriptor pd) {
        if (pd instanceof IndexedPropertyDescriptor) {
            return ((IndexedPropertyDescriptor) pd).getIndexedReadMethod();
        } else {
            return pd.getReadMethod();
        }
    }

    private static boolean hasAnnotation(Method method, Class<? extends Annotation> ann) {
        return MethodUtils.getAnnotation(method, ann, true, true) != null;
    }

}
