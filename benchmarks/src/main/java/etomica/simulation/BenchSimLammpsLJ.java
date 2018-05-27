package etomica.simulation;

import etomica.action.ActionIntegrate;
import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.tests.TestLammpsLJ;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.profile.LinuxPerfAsmProfiler;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLammpsLJ {

    TestLammpsLJ sim;
    IntegratorVelocityVerlet integrator;
    ActionIntegrate ai;

    @Setup(Level.Iteration)
    public void setUp() {
        sim = new TestLammpsLJ();
        integrator = sim.integrator;
        integrator.reset();
        ai = new ActionIntegrate(sim.integrator);
        ai.setMaxSteps(100);
    }

    @Benchmark
    @Measurement(iterations = 5, time = 5)
    @Warmup(iterations = 3, time = 3)
    public long integratorStep() {
        integrator.doStep();
        return integrator.getStepCount();
    }

    @Benchmark
    @BenchmarkMode(Mode.SingleShotTime)
    @Measurement(iterations = 1)
    @Warmup(iterations = 1)
    public long fullRun() {
        ai.actionPerformed();
        return integrator.getStepCount();
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include("BenchSimLammpsLJ")
//                .jvmArgs(
//                        "-XX:+UnlockDiagnosticVMOptions",
//                        "-XX:+PrintAssembly",
//                        "-XX:PrintAssemblyOptions=intel",
//                        "-XX:CompileCommand=print,*BoxBench.bench*")
                .build();

        new Runner(opts).run();
    }
}
