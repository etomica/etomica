package etomica.simulation.prototypes.simulation;

import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.space3d.Space3D;
import etomica.tests.TestSWChain;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.SECONDS)
public class BenchSimSWChain {

    @Param({"500", "4000"})
    public int numMolecules;

    @Param({"100000"})
    public int numSteps;

    private TestSWChain sim;
    private MeterPressureHard pMeter;

    @Setup
    public void setUp() {
        double simTime = numSteps / numMolecules;
        Configuration config = new ConfigurationResourceFile(
                String.format("tests/SWChain%d.pos", numMolecules),
                TestSWChain.class
        );

        sim = new TestSWChain(Space3D.getInstance(), numMolecules, simTime, config);

        pMeter = new MeterPressureHard(sim.getSpace());
        pMeter.setIntegrator(sim.integrator);
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed();
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(energyPump);
    }

    @Benchmark
    public double swChain() {
        sim.getController().actionPerformed();
        return pMeter.getDataAsScalar() + sim.integrator.getTemperature();
    }
}
