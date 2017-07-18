package etomica.meta.wrappers;

import com.fasterxml.jackson.core.JsonProcessingException;
import etomica.meta.InstanceProperty;
import etomica.meta.WrapperIndex;
import etomica.meta.annotations.IgnoreProperty;
import etomica.simulation.prototypes.HSMD2D;
import org.apache.commons.lang3.reflect.MethodUtils;

import java.beans.IndexedPropertyDescriptor;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.lang.annotation.Annotation;
import java.lang.reflect.Method;
import java.util.*;


public class Wrapper<T> {
    protected final T wrapped;
    protected final Class wrappedClass;

    protected final List<InstanceProperty> properties = new ArrayList<>();

    protected final List<InstanceProperty> childProps = new ArrayList<>();
    protected final List<InstanceProperty> valueProps = new ArrayList<>();


    public Wrapper(T wrapped) {
        this.wrapped = wrapped;
        this.wrappedClass = wrapped.getClass();

        try {
            Arrays.stream(Introspector.getBeanInfo(wrappedClass).getPropertyDescriptors())
                    .filter(propertyDescriptor -> !propertyDescriptor.getName().equalsIgnoreCase("class"))
                    .filter(propertyDescriptor -> propertyDescriptorMethod(propertyDescriptor) != null)
                    .filter(propertyDescriptor -> !hasAnnotation(propertyDescriptorMethod(propertyDescriptor), IgnoreProperty.class))
                    .filter(propertyDescriptor -> !propertyDescriptor.getName().toLowerCase().endsWith("dimension"))
                    .forEach(propertyDescriptor -> properties.add(new InstanceProperty(wrapped, propertyDescriptor)));

            for(InstanceProperty p : properties) {
                Class<?> type = p.getPropertyType();
                if(type.isPrimitive() || type.equals(String.class) || type.equals(Class.class)) {
                    valueProps.add(p);
                } else {
                    childProps.add(p);
                }
            }
        } catch (IntrospectionException e) {
            e.printStackTrace();
        }
    }

    public Class getWrappedClass() {
        return wrappedClass;
    }

    public List<InstanceProperty> getProperties() {
        return properties;
    }

    public Map<String, Object> getValues() {
        Map<String, Object> values = new HashMap<>();
        for(InstanceProperty prop : valueProps) {
            if(prop.isIndexedProperty()) {
                System.err.println("Indexed value prop");
                // TODO
            } else {
                values.put(prop.getName(), prop.invokeReader());
            }
        }

        return values;
    }

    @SuppressWarnings("unchecked")
    public Map<String, Wrapper<?>> getChildren() {
        Map<String, Wrapper<?>> children = new HashMap<>();
        for(InstanceProperty prop : childProps) {
            System.out.printf("%s -> %s (%s)\n", this.wrappedClass.getSimpleName(), prop.getName(), prop.getPropertyType().getSimpleName());

            if(prop.isIndexedProperty()) {
                if(!prop.canCount()) { continue; }

                for (int i = 0; i < prop.invokeCount(); i++) {
                    children.put(prop.getName(), WrapperIndex.getWrapper(prop.invokeReader(i)));
                }
            } else if(List.class.isAssignableFrom(prop.getPropertyType())) {
                List childList = (List) prop.invokeReader();
                for(Object o : childList) {
                    children.put(prop.getName(), WrapperIndex.getWrapper(o));
                }
//            } else if(prop.getPropertyType().isArray()) {
//                for(Object o : (Object[]) prop.invokeReader()) {
//                    children.add(WrapperIndex.getWrapper(o));
//                }
            } else {
                children.put(prop.getName(), WrapperIndex.getWrapper(prop.invokeReader()));
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

    private static boolean hasAnnotation(Method method, Class<? extends Annotation> ann) {
        return MethodUtils.getAnnotation(method, ann, true, true) != null;
    }

    public static void main(String[] args) throws JsonProcessingException {
        SimulationWrapper simWrapper = new SimulationWrapper(new HSMD2D());
        List<InstanceProperty> props = simWrapper.getProperties();
        System.out.println(props);
    }

}
