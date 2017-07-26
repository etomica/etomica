package etomica.meta.properties;

/**
 * Created by kofke on 7/23/17.
 */
public interface Property {
    Object invokeReader();

    String getName();

    Object invokeReader(int i);

    void invokeWriter(Object... params);

    void invokeWriter(int i, Object... params);

    void invokeAdder(Object o);

    void invokeRemover(Object o);

    int invokeCount();

    boolean isIndexedProperty();

    boolean canRead();

    boolean canWrite();

    boolean canAdd();

    boolean canRemove();

    boolean canCount();

    Class<?> getPropertyType();

    boolean isValueProperty();
}
