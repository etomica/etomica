package etomica.meta.properties;

import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Container class for mapping an instance to an array property.
 * \
 */
public class ArrayProperty implements Property {
    private final PropertyDescriptor descriptor;
    private final Class<?> propertyType;
    private final Object instance;
    private final Method reader;

    public ArrayProperty(Object instance, PropertyDescriptor descriptor) {
        this.descriptor = descriptor;
        this.instance = instance;
        this.propertyType = descriptor.getPropertyType().getComponentType();

        reader = descriptor.getReadMethod();
    }

    @Override
    public Object invokeReader() {
        throw new RuntimeException("not a simple property");
    }

    @Override
    public String getName() {
        return descriptor.getName();
    }

    @Override
    public Object invokeReader(int i) {
        try {
            return ((Object[]) reader.invoke(instance))[i];
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void invokeWriter(Object... params) {
        throw new RuntimeException("not a simple property");
    }

    @Override
    public void invokeWriter(int i, Object... params) {
        throw new RuntimeException("not a simple property");
    }

    @Override
    public void invokeAdder(Object o) {
        throw new RuntimeException("not a addable property");
    }

    @Override
    public void invokeRemover(Object o) {
        throw new RuntimeException("not a removable property");
    }

    @Override
    public int invokeCount() {
        try {
            return ((Object[]) reader.invoke(instance)).length;
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public final boolean isIndexedProperty() {
        return true;
    }

    @Override
    public boolean canRead() {
        return reader != null;
    }

    @Override
    public boolean canWrite() {
        return false;
    }

    @Override
    public boolean canAdd() {
        return false;
    }

    @Override
    public boolean canRemove() {
        return false;
    }

    @Override
    public boolean canCount() {
        return true;
    }

    @Override
    public Class<?> getPropertyType() {
        return propertyType;
    }
}
