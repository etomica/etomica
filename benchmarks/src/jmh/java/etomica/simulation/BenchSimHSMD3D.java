package etomica.simulation;

import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.data.meter.MeterPressureHard;
import etomica.space3d.Space3D;
import etomica.tests.TestHSMD3D;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimHSMD3D {

    @Param({"500", "4000", "32000"})
    private int numMolecules;

    // 20_000_000
    @Param({"20000000"})
    private int numSteps;

    private TestHSMD3D sim;
    private MeterPressureHard pMeter;

    @Setup(Level.Iteration)
    public void setUp() {

        Configuration config = new ConfigurationResourceFile(
                String.format("tests/HSMD3D%d.pos", numMolecules),
                TestHSMD3D.class
        );

        sim = new TestHSMD3D(Space3D.getInstance(), numMolecules, numSteps / numMolecules, config);

        pMeter = new MeterPressureHard(sim.space);
        pMeter.setIntegrator(sim.integrator);
        sim.integrator.reset();
    }

    @Benchmark
    @BenchmarkMode(Mode.SingleShotTime)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(iterations = 0)
    public double hsmd3d() {
        sim.getController().actionPerformed();
        return pMeter.getDataAsScalar() + sim.integrator.getTemperature();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 5, timeUnit = TimeUnit.SECONDS, iterations = 5)
    public double integratorStep() {
        sim.integrator.doStep();
        return pMeter.getDataAsScalar();
    }
}
