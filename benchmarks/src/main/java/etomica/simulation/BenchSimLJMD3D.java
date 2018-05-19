/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.simulation.prototypes.LJMD3D;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJMD3D {

    private LJMD3D sim;

    @Setup(Level.Iteration)
    public void setUp() {

        sim = new LJMD3D();
        sim.integrator.reset();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public long integratorStep() {
        sim.integrator.doStep();
        //sim.integrator.getPotentialMaster().calculate(sim.box, id, pc);
        return sim.integrator.getStepCount();
    }
}
