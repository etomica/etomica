package etomica.box.storage;

import etomica.space.Vector;
import etomica.space3d.Space3D;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Objects;

import static org.junit.jupiter.api.Assertions.*;

class VectorStorageTest {

    @BeforeEach
    void setUp() {

    }

    @Test
    void equalStandaloneVector() {
        VectorStorage storage = new VectorStorage(Space3D.getInstance(), 5);
        storage.create(2).E(3);
        Vector v = storage.get(2);
        assertEquals(Vector.of(3, 3, 3), v);
    }

    @Test
    void get() {
        VectorStorage storage = new VectorStorage(Space3D.getInstance(), 5);

        assertNull(storage.get(2));
        storage.create(2).E(new double[]{1, 2, 3});
        assertEquals(Vector.of(1, 2, 3), storage.get(2));
        storage.get(2).E(5);
        assertEquals(Vector.of(5, 5, 5), storage.get(2));
    }

    @Test
    void addNull() {
        VectorStorage storage = new VectorStorage(Space3D.getInstance(), 5);
        storage.addNull(2);
        assertNull(storage.get(5));
        assertNull(storage.get(6));
    }

    @Test
    void add() {
        VectorStorage storage = new VectorStorage(Space3D.getInstance(), 5);
        storage.add(2);
        assertEquals(Vector.of(0, 0, 0), storage.get(5));
        assertEquals(Vector.of(0, 0, 0), storage.get(6));
    }

    @Test
    void swapRemove() {
        VectorStorage storage = new VectorStorage(Space3D.getInstance(), 5);
        assertEquals(5, storage.size());
        storage.create(2).E(2);
        storage.create(3).E(3);
        storage.create(4).E(4);

        storage.swapRemove(2);
        assertEquals(Vector.of(4, 4, 4), storage.get(2));
    }
}