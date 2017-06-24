package etomica.simulation;

import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterPressureTensor;
import etomica.space3d.Space3D;
import etomica.surfacetension.LJMC;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimSurfacetensionLJMC {

    @Param({"1000", "20000"})
    private int numAtoms;

    private LJMC sim;
    private MeterPressureTensor meter;

    @Setup(Level.Iteration)
    public void setUp() {
        sim = new LJMC(Space3D.getInstance(), numAtoms, 1.1, 6);

        meter = new MeterPressureTensor(sim.potentialMaster, sim.space);
        meter.setBox(sim.box);
        meter.setTemperature(1.1);
        DataPumpListener pumpListener = new DataPumpListener(meter, new AccumulatorAverageFixed(10), 2 * numAtoms);
        sim.integrator.getEventManager().addListener(pumpListener);
        sim.integrator.reset();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 5, timeUnit = TimeUnit.SECONDS, iterations = 5)
    public IData integratorStep() {
        sim.integrator.doStep();
        return meter.getData();
    }
}
