package etomica.meta;

import java.beans.IndexedPropertyDescriptor;
import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

public class InstanceProperty {
    private final PropertyDescriptor descriptor;
    private final Class<?> propertyType;
    private final Object instance;
    private final Method reader;
    private final Method writer;
    private final Method adder;
    private final Method remover;
    private final Method counter;


    public InstanceProperty(Object instance, PropertyDescriptor descriptor) {
        this.descriptor = descriptor;
        this.instance = instance;
        this.propertyType = isIndexedProperty() ? ((IndexedPropertyDescriptor) descriptor).getIndexedPropertyType() : descriptor.getPropertyType();

        reader = isIndexedProperty() ? ((IndexedPropertyDescriptor) descriptor).getIndexedReadMethod() : descriptor.getReadMethod();
        writer = isIndexedProperty() ? ((IndexedPropertyDescriptor) descriptor).getIndexedWriteMethod() : descriptor.getWriteMethod();

        String baseName = descriptor.getName().substring(0, 1).toUpperCase() + descriptor.getName().substring(1);

        adder = getMethod("add" + baseName, propertyType);
        remover = getMethod("remove" + baseName, propertyType);
        counter = getMethod("get" + baseName + "Count");
    }

    public Object invokeReader() {
        try {
            return reader.invoke(instance);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    public String getName() {
        return descriptor.getName();
    }

    public Object invokeReader(int i) {
        try {
            return reader.invoke(instance, i);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    public void invokeWriter(Object... params) {
        try {
            writer.invoke(instance, params);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    public void invokeWriter(int i, Object... params) {
        try {
            writer.invoke(instance, i, params);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    public void invokeAdder(Object o) {
        try {
            adder.invoke(instance, o);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    public void invokeRemover(Object o) {
        try {
            remover.invoke(instance, o);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    public int invokeCount() {
        try {
            return (int) counter.invoke(instance);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }

    }

    public boolean isIndexedProperty() {
        return descriptor instanceof IndexedPropertyDescriptor;
    }

    public boolean canRead() {
        return reader != null;
    }

    public boolean canWrite() {
        return reader != null;
    }

    public boolean canAdd() {
        return adder != null;
    }

    public boolean canRemove() {
        return remover != null;
    }

    public Class<?> getPropertyType() {
        return propertyType;
    }

    private Method getMethod(String name, Class<?>... paramTypes) {
        try {
            return instance.getClass().getMethod(name, paramTypes);
        } catch (NoSuchMethodException e) {
            return null;
        }
    }
}
