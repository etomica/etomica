package etomica.box.storage;

import etomica.space.Space;

public interface Token<T extends Storage> {
    T createStorage(Space space);

    default void init(int idx, T storage) {}
}