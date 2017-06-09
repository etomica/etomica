package etomica.simulation.prototypes;

import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMD;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.profile.GCProfiler;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.results.RunResult;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.Collection;

/**
 * Created by alex on 4/29/17.
 */
@State(Scope.Benchmark)
public class HSMD3DBenchmark {

    private Simulation sim;

    @Setup
    public void setUp() {
        sim = new HSMD3D(Space3D.getInstance());
    }

    @Benchmark
    public void hsmd3d() {
        sim.getIntegrator().doStep();
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(HSMD3DBenchmark.class.getSimpleName())
                .forks(0)
                .addProfiler(GCProfiler.class)
                .build();

        new Runner(opts).run();
    }
}
