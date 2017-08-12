package etomica.meta.properties;

import etomica.atom.IAtomList;

import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Created by kofke on 7/24/17.
 */
public class AtomListProperty implements Property {

    private final Object instance;
    private final String name;
    private final Method reader;

    public AtomListProperty(Object instance, PropertyDescriptor descriptor) {
        this.name = descriptor.getName();
        this.instance = instance;

        reader = descriptor.getReadMethod();
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
            return ((IAtomList) reader.invoke(instance)).getAtom(i);
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
            return ((IAtomList) reader.invoke(instance)).getAtomCount();
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
        return IAtomList.class;
    }

    public boolean isValueProperty() {
        return false;
    }
}
