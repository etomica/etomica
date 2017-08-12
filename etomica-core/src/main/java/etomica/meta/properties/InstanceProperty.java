package etomica.meta.properties;

import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Dimensioned;

import java.beans.IndexedPropertyDescriptor;
import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Container class for mapping an instance to a property.
 * <p>
 * Also contains additional adder and remover methods for that property if they exist.
 */
public class InstanceProperty implements Property {
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

    @Override
    public Object invokeReader() {
        try {
            return reader.invoke(instance);
        } catch(IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public String getName() {
        return descriptor.getName();
    }

    @Override
    public Object invokeReader(int i) {
        try {
            return reader.invoke(instance, i);
        } catch(IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void invokeWriter(Object... params) {
        try {
            writer.invoke(instance, params);
        } catch(IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void invokeWriter(int i, Object... params) {
        try {
            writer.invoke(instance, i, params);
        } catch(IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void invokeAdder(Object o) {
        try {
            adder.invoke(instance, o);
        } catch(IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void invokeRemover(Object o) {
        try {
            remover.invoke(instance, o);
        } catch(IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public int invokeCount() {
        try {
            return (int) counter.invoke(instance);
        } catch(IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }

    }

    @Override
    public final boolean isIndexedProperty() {
        return descriptor instanceof IndexedPropertyDescriptor;
    }

    @Override
    public boolean canRead() {
        return reader != null;
    }

    @Override
    public boolean canWrite() {
        return writer != null;
    }

    @Override
    public boolean canAdd() {
        return adder != null;
    }

    @Override
    public boolean canRemove() {
        return remover != null;
    }

    @Override
    public boolean canCount() {
        return counter != null;
    }

    @Override
    public Class<?> getPropertyType() {
        return propertyType;
    }

    private Method getMethod(String name, Class<?>... paramTypes) {
        try {
            return instance.getClass().getMethod(name, paramTypes);
        } catch(NoSuchMethodException e) {
            return null;
        }
    }

    public boolean isValueProperty() {
        return (propertyType.isPrimitive() || propertyType.equals(String.class) || propertyType.equals(Class.class)
                || (propertyType.isArray() && (propertyType.getComponentType().isPrimitive() || propertyType.getComponentType().equals(String.class))));
    }

    private Class<? extends Dimension> getDimension() {
        Class<?> cls = instance.getClass();
        Dimensioned ann = null;

        while(cls != null) {
            ann = cls.getAnnotation(Dimensioned.class);
            if(ann != null) {
                break;
            }
            cls = cls.getSuperclass();
        }

        if(ann == null) {
            return null;
        }

        return ann.dimension();
    }
}
