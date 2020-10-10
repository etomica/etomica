package etomica;

import etomica.integrator.Integrator;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.TestReporter;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

import static etomica.TestSimConstructors.SCAN;

public class TestSimShortRuns {
    private static final ExecutorService exec = Executors.newSingleThreadExecutor();

    /**
     * Sim classes that are excluded for various reasons causing them not to work with the reflection hacks here.
     */
    private static final Set<Class<?>> EXCLUDED = new HashSet<>(Arrays.asList(
    ));

    static List<Arguments> getIntegrators() throws IllegalAccessException, InstantiationException, InvocationTargetException {
        // This beautiful method will search through constructible simulations looking for ones
        // with one and only one field of type Integrator (getIntegrator isn't reliable, and with
        // more than one we don't know which one to use), instantiate the simulations, and return the integrators.
        List<Class<?>> classes =  SCAN.getSubclasses(Simulation.class.getName()).loadClasses();
        List<Arguments> args = new ArrayList<>();

        for (Class<?> cl : classes) {
            if (EXCLUDED.contains(cl)) {
                continue;
            }

            Field[] fields = cl.getDeclaredFields();
            List<Field> integratorFields = Arrays.stream(fields)
                    .filter(f -> Integrator.class.isAssignableFrom(f.getType()))
                    .collect(Collectors.toList());

            if (integratorFields.size() == 1) {
                Simulation sim = null;
                try {
                    sim = (Simulation) cl.getConstructor().newInstance();
                } catch (NoSuchMethodException e) {
                    try {
                        sim = (Simulation) cl.getConstructor(Space.class).newInstance(Space3D.getInstance());
                    } catch (NoSuchMethodException e1) {
                        continue;
                    }
                }

                integratorFields.get(0).setAccessible(true);
                args.add(Arguments.of(cl, integratorFields.get(0).get(sim)));
            }
        }

        return args;
    }

    @ParameterizedTest
    @MethodSource("getIntegrators")
    void testShortRuns(Class<?> simClass, Integrator simIntegrator, TestReporter reporter) {
        Assertions.assertDoesNotThrow(() -> {
            simIntegrator.reset();
            Future future = exec.submit(() -> {
                while (!Thread.interrupted()) {
                    simIntegrator.doStep();
                }
            });
            Thread.sleep(400);
            future.cancel(true);
        });
        reporter.publishEntry(simClass + " steps", Long.toString(simIntegrator.getStepCount()));
    }
}
