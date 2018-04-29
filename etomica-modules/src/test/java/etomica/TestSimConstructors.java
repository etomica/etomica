package etomica;

import etomica.normalmode.NormalModeAnalysisDisplay3D;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import io.github.lukehutch.fastclasspathscanner.FastClasspathScanner;
import io.github.lukehutch.fastclasspathscanner.scanner.ScanResult;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;

public class TestSimConstructors {
    public static final ScanResult SCAN = new FastClasspathScanner("etomica").scan();
    private static final List<Constructor<?>> constructors = new ArrayList<>();
    private static final List<Constructor<?>> spaceConstructors = new ArrayList<>();

    @BeforeAll
    static void setUp() {
        List<Class<?>> classes =  SCAN.classNamesToClassRefs(SCAN.getNamesOfSubclassesOf(Simulation.class));
        for (Class<?> cls : classes) {
            try {
                constructors.add(cls.getConstructor());
            } catch (NoSuchMethodException e) {
                continue;
            }
        }

        for (Class<?> cls : classes) {
            try {
                spaceConstructors.add(cls.getConstructor(Space.class));
            } catch (NoSuchMethodException e) {
                continue;
            }
        }

    }

    static List<Constructor<?>> getConstructors() {
        return constructors;
    }

    static List<Constructor<?>> getSpaceConstructors() {
        return spaceConstructors;
    }

    @ParameterizedTest
    @MethodSource("getConstructors")
    public void testNullaryConstructors(Constructor<?> constructor) {
        // Special case to skip this class. It requires files that don't exists. TODO delete if we can't find them
        Assumptions.assumeFalse(constructor.getDeclaringClass() == NormalModeAnalysisDisplay3D.class);

        assertDoesNotThrow(() -> {
            Simulation sim = (Simulation) constructor.newInstance();
        });
    }

    @ParameterizedTest
    @MethodSource("getSpaceConstructors")
    public void testSpaceConstructors(Constructor<?> constr) {
        assertDoesNotThrow(() -> {
            constr.newInstance(Space3D.getInstance());
            constr.newInstance(Space2D.getInstance());
        });
    }
}
