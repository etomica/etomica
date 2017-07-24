package etomica.meta.properties;

import etomica.meta.properties.Property;
import etomica.space.Vector;

/**
 * Created by kofke on 7/24/17.
 */
public class VectorProperty implements Property {

    private final Vector v;
    public VectorProperty(Vector v) {
        this.v = v;
    }

    @Override
    public Object invokeReader() {
        throw new RuntimeException();
    }

    @Override
    public String getName() {
        return "Vector name";  //don't know what to put here
    }

    @Override
    public Object invokeReader(int i) {
        return v.getX(i);
    }

    @Override
    public void invokeWriter(Object... params) {
        throw new RuntimeException();
    }

    @Override
    public void invokeWriter(int i, Object... params) {
        v.setX(i, ((Number)params[0]).doubleValue());//is this right?
    }

    @Override
    public void invokeAdder(Object o) {
        throw new RuntimeException();
    }

    @Override
    public void invokeRemover(Object o) {
        throw new RuntimeException();
    }

    @Override
    public int invokeCount() {
        return v.getD();
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
        return true;
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
}
