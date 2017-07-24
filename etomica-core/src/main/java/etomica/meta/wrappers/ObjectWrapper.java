package etomica.meta.wrappers;

import etomica.meta.properties.InstanceProperty;
import etomica.meta.properties.Property;
import etomica.meta.SimulationModel;
import etomica.meta.annotations.IgnoreProperty;
import org.apache.commons.lang3.reflect.MethodUtils;

import java.beans.IndexedPropertyDescriptor;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.lang.annotation.Annotation;
import java.lang.reflect.Method;
import java.util.Arrays;

public class ObjectWrapper<T> extends Wrapper<T> {

    public ObjectWrapper(T wrapped, SimulationModel simModel) {
        super(wrapped, simModel);

        try {
            Arrays.stream(Introspector.getBeanInfo(wrappedClass).getPropertyDescriptors())
                    .filter(propertyDescriptor -> !propertyDescriptor.getName().equalsIgnoreCase("class"))
                    .filter(propertyDescriptor -> propertyDescriptorMethod(propertyDescriptor) != null)
                    .filter(propertyDescriptor -> !hasAnnotation(propertyDescriptorMethod(propertyDescriptor), IgnoreProperty.class))
                    .filter(propertyDescriptor -> !propertyDescriptor.getName().toLowerCase().endsWith("dimension"))
                    .forEach(propertyDescriptor -> properties.add(new InstanceProperty(wrapped, propertyDescriptor)));

            for (Property p : properties) {
                Class<?> type = p.getPropertyType();
                if (type.isPrimitive() || type.equals(String.class) || type.equals(Class.class)
                        || (type.isArray() && (type.getComponentType().isPrimitive() || type.getComponentType().equals(String.class)))) {
                    valueProps.add(p);
                } else {
                    childProps.add(p);
                }
            }
        } catch (IntrospectionException e) {
            e.printStackTrace();
        }
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
