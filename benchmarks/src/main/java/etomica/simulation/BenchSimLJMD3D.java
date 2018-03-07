/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.nbr.list.PotentialMasterListFast;
import etomica.nbr.list.PotentialMasterListFaster;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.prototypes.LJMD3D;
import etomica.simulation.prototypes.LJMD3DFast;
import etomica.simulation.prototypes.LJMD3DFaster;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJMD3D {

    private LJMD3D sim;
    private LJMD3DFast simFast;
    private LJMD3DFaster simFaster;
    private PotentialCalculationForceSum pc, pcFast;
    private IteratorDirective id;

    @Setup(Level.Iteration)
    public void setUp() {

        sim = new LJMD3D();
        id = new IteratorDirective();
        id.setDirection(IteratorDirective.Direction.UP);

        simFast = new LJMD3DFast();
        simFaster = new LJMD3DFaster();
        pc = sim.integrator.getForceSum();
        pcFast = simFast.integrator.getForceSum();
        sim.integrator.reset();
        simFast.integrator.reset();
        simFaster.integrator.reset();
    }

//    @Setup(Level.Invocation)
//    public void setUpCall() {
//
//    }

//    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public long integratorStep() {
//        sim.integrator.doStep();
        sim.integrator.getPotentialMaster().calculate(sim.box, id, pc);
        return sim.integrator.getStepCount();
    }


    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public long integratorStepFast() {
//        simFast.integrator.doStep();
        ((PotentialMasterListFast) simFast.integrator.getPotentialMaster()).calculateFast(simFast.box, id, pcFast);
        return simFast.integrator.getStepCount();
    }

//    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public long integratorStepFaster() {
//        simFaster.integrator.doStep();
        simFaster.integrator.getPotentialMaster().calculate(simFast.box, id, pcFast);
        return simFaster.integrator.getStepCount();
    }
}

