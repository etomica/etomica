/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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

    @Setup(Level.Iteration)
    public void setUp() {

        Configuration config = new ConfigurationResourceFile(
                String.format("HSMD3D%d.pos", numMolecules),
                TestHSMD3D.class
        );

        {
            sim = new TestHSMD3D(Space3D.getInstance(), numMolecules, config);

            MeterPressureHard pMeter = new MeterPressureHard(sim.integrator);
            sim.integrator.reset();
        }
    }

    //    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, timeUnit = TimeUnit.SECONDS, iterations = 5)
    public void integratorStep() {
        sim.integrator.doStep();
    }

}
