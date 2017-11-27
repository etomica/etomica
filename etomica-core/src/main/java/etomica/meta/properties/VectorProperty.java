package etomica.meta.properties;

import etomica.space.Vector;

import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.List;

/**
 * Created by kofke on 7/24/17.
 */
public class VectorProperty implements Property {

    private final Object instance;
    private final String name;
    private final Method reader;
    private final Method writer;

    public VectorProperty(Object instance, PropertyDescriptor descriptor) {
        this.name = descriptor.getName();
        this.instance = instance;

        reader = descriptor.getReadMethod();
        writer = descriptor.getWriteMethod();
    }

    @Override
    public Object invokeReader() {
        throw new RuntimeException("not a simple property");
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public Object invokeReader(int i) {
        try {
            return ((Vector) reader.invoke(instance)).getX(i);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void invokeWriter(Object... params) {
        try {
            Vector v = (Vector) reader.invoke(instance);
            if(double[].class.isAssignableFrom(params[0].getClass())) {
                v.E((double[]) params[0]);
            } else if(List.class.isAssignableFrom(params[0].getClass())) {
                List<Double> l = (List<Double>) params[0];
                double[] d = new double[l.size()];
                for (int i = 0; i < l.size(); i++) {
                    d[i] = l.get(i);
                }
                v.E(d);

            }
        } catch (IllegalAccessException | InvocationTargetException e) {
            e.printStackTrace();
        }

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
            return ((Vector) reader.invoke(instance)).getD();
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public boolean isIndexedProperty() {
        return true;
    }

    @Override
    public boolean canRead() {
        return true;
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
        return Vector.class;
    }

    public boolean isValueProperty() {
        return true;
    }
}
