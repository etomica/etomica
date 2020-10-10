package etomica.space;

import etomica.box.storage.VectorStorage;
import etomica.space3d.Space3D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.StringJoiner;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

class VectorTest {

    Vector v1;
    Vector vv1;
    Vector base1;
    Vector base2;
    Vector vBase1;
    Vector vBase2;
    VectorStorage storage;

    private static BiConsumer<Vector, Vector> withString(BiConsumer<Vector, Vector> method, String name) {
        return new BiConsumer<Vector, Vector>() {
            @Override
            public void accept(Vector vector, Vector vector2) {
                method.accept(vector, vector2);
            }

            @Override
            public String toString() {
                return name;
            }
        };
    }

    private static Stream<BiConsumer<Vector, Vector>> makeMethods() {
        return Stream.of(
                withString(Vector::E, "E"),
                withString(Vector::PE, "PE"),
                withString((v1, v2) -> v1.PE(v2.x()), "PE double"),
                withString((v1, v2) -> v1.TE(v2.x()), "TE double"),
                withString(Vector::ME, "ME"),
                withString(Vector::TE, "TE"),
                withString(Vector::DE, "DE"),
                withString((v1, v2) -> v1.PEa1Tv1(3.1, v2), "PEa1Tv1"),
                withString((v1, v2) -> v1.Ev1Mv2(v2, v2), "Ev1Mv2"),
                withString((v1, v2) -> v1.Ev1Pv2(v2, v2), "Ev1Pv2"),
                withString(Vector::mod, "mod"),
                withString(Vector::XE, "XE"),
                withString(Vector::nearestImage, "nearestImage"),
                withString((v1, v2) -> v1.E(v1.dot(v2)), "dot"),
                withString((v1, v2) -> v1.E(v1.squared() + v2.squared()), "squared"),
                withString((v1, v2) -> {
                    v1.normalize();
                    v2.normalize();
                    v1.PE(v2);
                }, "normalize"),
                withString((v1, v2) -> v1.E(v1.Mv1Squared(v2)), "Mv1Squared")
        );
    }

    @BeforeEach
    void setUp() {
        v1 = Vector.of(1, 2, 3);
        storage = new VectorStorage(Space3D.getInstance(), 5);
        vv1 = storage.create(1);
        vv1.E(1, 2, 3);

        vBase1 = storage.create(3);
        vBase2 = storage.create(4);
        vBase1.E(10, 11, 12);
        vBase2.E(10, 11, 12);
        base1 = Vector.of(10, 11, 12);
        base2 = Vector.of(10, 11, 12);

    }

    @ParameterizedTest
    @MethodSource("makeMethods")
    void testCombos(BiConsumer<Vector, Vector> method) {
        method.accept(base1, v1);
        method.accept(base2, vv1);
        method.accept(vBase1, v1);
        method.accept(vBase2, vv1);

        assertAll(
                () -> assertEquals(base1, base2, "vec <op> viewVec"),
                () -> assertEquals(base1, vBase1, "viewVec <op> vec"),
                () -> assertEquals(base1, vBase2, "viewVec <op> viewVec")
        );

    }

    @Test
    void squared() {
    }

    @Test
    void dot() {
    }

    @Test
    void mv1Squared() {
    }

    @Test
    void PE() {
    }

    @Test
    void ME() {
    }

    @Test
    void TE() {
    }

    @Test
    void DE() {
    }

    @Test
    void ea1Tv1() {
    }

    @Test
    void PEa1Tv1() {
    }

    @Test
    void ev1Pv2() {
    }

    @Test
    void ev1Mv2() {
    }

    @Test
    void mod() {
    }

    @Test
    void normalize() {
    }

    @Test
    void nearestImage() {
    }
}