package etomica.simulation;

import etomica.experimental.LJMD3DVecSys;
import etomica.integrator.Integrator;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.prototypes.LJMD3DNbr;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.profile.LinuxPerfAsmProfiler;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJMD3DFast {


    private Integrator integrator;

    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include(BenchSimLJMD3DFast.class.getSimpleName())
//                .addProfiler(LinuxPerfAsmProfiler.class)
                .jvmArgs(
                        "-XX:+UnlockDiagnosticVMOptions"
//                        "-XX:+PrintAssembly",
//                        "-XX:CompileCommand=print etomica.experimental.VectorSystem3D::sub"
//                        ,"-XX:-UseSuperWord"
                )
                .build();

        new Runner(opts).run();
    }

    @Setup(Level.Iteration)
    public void setUp() {
        LJMD3DVecSys sim = new LJMD3DVecSys();
        integrator = sim.integrator;
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public long integratorStep() {
        integrator.doStep();
        return integrator.getStepCount();
    }
}
