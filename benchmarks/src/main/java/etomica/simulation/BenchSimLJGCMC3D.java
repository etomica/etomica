/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPressure;
import etomica.tests.TestLJGCMC3D;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJGCMC3D {

    @Param({"500", "4000"})
    private int numMolecules;

    @Param({"200000"})
    private int numSteps;

    private TestLJGCMC3D sim;
    private MeterPressure pMeter;

    @Setup(Level.Iteration)
    public void setUp() {

        Configuration config = new ConfigurationResourceFile(
                // Uses the same position config as LJMC3D
                String.format("tests/LJMC3D%d.pos", numMolecules),
                TestLJGCMC3D.class
        );

        sim = new TestLJGCMC3D(numMolecules, numSteps, config);

        pMeter = new MeterPressure(sim.space);
        pMeter.setIntegrator(sim.integrator);
        DataPumpListener pumpListener = new DataPumpListener(pMeter, new AccumulatorAverageFixed(10), 2 * numMolecules);
        sim.integrator.getEventManager().addListener(pumpListener);
        sim.integrator.reset();
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

