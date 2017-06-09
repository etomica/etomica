package etomica.simulation;

import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.data.meter.MeterPressure;
import etomica.tests.TestLJMC3D;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJMC3D {

    @Param({"500", "4000"})
    private int numMolecules;

    @Param({"200000"})
    private int numSteps;

    private TestLJMC3D sim;
    private MeterPressure pMeter;

    @Setup
    public void setUp() {

        Configuration config = new ConfigurationResourceFile(
                String.format("tests/LJMC3D%d.pos", numMolecules),
                TestLJMC3D.class
        );

        sim = new TestLJMC3D(numMolecules, numSteps / numMolecules, config);

        pMeter = new MeterPressure(sim.space);
        pMeter.setIntegrator(sim.integrator);
    }

    @Benchmark
    @BenchmarkMode(Mode.SingleShotTime)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(iterations = 0)
    public double ljmc3d() {
        sim.getController().actionPerformed();
        return pMeter.getDataAsScalar() + sim.integrator.getTemperature();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 1, timeUnit = TimeUnit.SECONDS)
    public long integratorStep() {
        sim.integrator.doStep();
        return sim.integrator.getStepCount();
    }
}

