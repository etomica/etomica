package etomica.simulation;

import etomica.modules.catalysis.CatalysisFasterer;
import etomica.space3d.Space3D;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimCatalysis {
    private CatalysisFasterer sim;

    @Param("20")
    private int nCellsZ;

    @Setup(Level.Iteration)
    public void setup() {
        sim = new CatalysisFasterer(Space3D.getInstance(), nCellsZ);
        sim.integrator.reset();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 3, iterations = 5)
    @Measurement(time = 10, timeUnit = TimeUnit.SECONDS, iterations = 3)
    public long integratorStep() {
        sim.integrator.doStep();
        return sim.integrator.getStepCount();
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(BenchSimCatalysis.class.getSimpleName())
                .addProfiler(StackProfiler.class)
                .build();

        new Runner(opts).run();
    }
}
