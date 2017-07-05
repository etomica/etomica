package etomica.meta;

import java.beans.IndexedPropertyDescriptor;
import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;

public class InstanceProperty {
    private final PropertyDescriptor descriptor;
    private final Object instance;


    public InstanceProperty(Object instance, PropertyDescriptor descriptor) {
        this.descriptor = descriptor;
        this.instance = instance;
    }

    public Object invokeReader() throws InvocationTargetException, IllegalAccessException {
        return descriptor.getReadMethod().invoke(instance);
    }

    public Object invokeReader(int i) throws InvocationTargetException, IllegalAccessException {
        return ((IndexedPropertyDescriptor) descriptor).getIndexedReadMethod().invoke(instance, i);
    }

    public void invokeWriter(Object... params) throws InvocationTargetException, IllegalAccessException {
        descriptor.getWriteMethod().invoke(instance, params);
    }

    public void invokeWriter(int i, Object... params) throws InvocationTargetException, IllegalAccessException {
        ((IndexedPropertyDescriptor) descriptor).getIndexedWriteMethod().invoke(instance, i, params);
    }

    public boolean isIndexedProperty() {
        return descriptor instanceof IndexedPropertyDescriptor;
    }
}
