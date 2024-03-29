/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPressure;
import etomica.tests.TestLJMC3D;
import etomica.tests.TestLJMC3DBrute;
import etomica.tests.TestLJMC3DNew;
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
    private TestLJMC3DNew simNew;
    private TestLJMC3DBrute simBrute;

    @Setup(Level.Iteration)
    public void setUp() {

        Configuration config = new ConfigurationResourceFile(
                String.format("LJMC3D%d.pos", numMolecules),
                TestLJMC3D.class
        );

//        {
//            simSlowBrute = new TestLJMC3DSlowBrute(numMolecules, numSteps, config);
//
//            MeterPressure pMeter = new MeterPressure(simSlowBrute.space);
//            pMeter.setIntegrator(simSlowBrute.integrator);
//            DataPumpListener pumpListener = new DataPumpListener(pMeter, new AccumulatorAverageFixed(10), 2 * numMolecules);
//            simSlowBrute.integrator.getEventManager().addListener(pumpListener);
//            simSlowBrute.integrator.reset();
//        }

//        {
//            simSlowerer = new TestLJMC3DSlowerer(numMolecules, numSteps, config);
//
//            MeterPressure pMeter = new MeterPressure(simSlowerer.space);
//            pMeter.setIntegrator(simSlowerer.integrator);
//            DataPumpListener pumpListener = new DataPumpListener(pMeter, new AccumulatorAverageFixed(10), 2 * numMolecules);
//            simSlowerer.integrator.getEventManager().addListener(pumpListener);
//            simSlowerer.integrator.reset();
//        }

        {
            simBrute = new TestLJMC3DBrute(numMolecules, config);

            MeterPressure pMeter = new MeterPressure(simBrute.box, simBrute.potentialMaster);
            pMeter.setTemperature(simBrute.integrator.getTemperature());
            DataPumpListener pumpListener = new DataPumpListener(pMeter, new AccumulatorAverageFixed(10), 2 * numMolecules);
            simBrute.integrator.getEventManager().addListener(pumpListener);
            simBrute.integrator.reset();
        }

        {
            sim = new TestLJMC3D(numMolecules, config);

            MeterPressure pMeter = new MeterPressure(sim.box, sim.potentialMaster);
            pMeter.setTemperature(sim.integrator.getTemperature());
            DataPumpListener pumpListener = new DataPumpListener(pMeter, new AccumulatorAverageFixed(10), 2 * numMolecules);
            sim.integrator.getEventManager().addListener(pumpListener);
            sim.integrator.reset();
        }
        {
            simNew = new TestLJMC3DNew(numMolecules, config);
            MeterPressure pMeter = new MeterPressure(simNew.box, simNew.integrator.getPotentialCompute());
            pMeter.setTemperature(simNew.integrator.getTemperature());
            DataPumpListener pumpListener = new DataPumpListener(pMeter, new AccumulatorAverageFixed(10), 2 * numMolecules);
            simNew.integrator.getEventManager().addListener(pumpListener);
            simNew.integrator.reset();
        }

    }

//    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 10, iterations = 3)
    public void integratorStepBrute() {
        simBrute.integrator.doStep();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 10, iterations = 3)
    public void integratorStep() {
        sim.integrator.doStep();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 10, iterations = 3)
    public void integratorStepNew() {
        simNew.integrator.doStep();
    }


}

