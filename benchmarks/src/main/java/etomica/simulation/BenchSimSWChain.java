/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

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
import org.openjdk.jmh.profile.GCProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimSWChain {

    @Param({"500", "4000"})
    public int numMolecules;

    @Param({"100000"})
    public int numSteps;

    private TestSWChain sim;
    private MeterPressureHard pMeter;

    @Setup(Level.Iteration)
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
        sim.integrator.reset();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 3, iterations = 3)
    @Measurement(iterations = 5, time = 5, timeUnit = TimeUnit.SECONDS)
    public long integratorStep() {
        sim.integrator.doStep();
        return sim.integrator.getStepCount();
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(BenchSimSWChain.class.getSimpleName())
                .build();

        new Runner(opts).run();
    }
}
