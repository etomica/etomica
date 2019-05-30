/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.glass;

import etomica.data.*;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterPressureHardTensor;
import etomica.data.meter.MeterPressureTensorFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.util.ParseArgs;

public class GlassProd {
    public static void main(String[] args) {
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.doSwap = true;
            params.potential = SimGlass.PotentialChoice.LJ;
            params.nA = 256;
            params.nB = 0;
            params.density = 0.844;
            params.D = 3;
            params.temperature = 0.722;
            params.numStepsEq = 10000;
            params.numSteps =   100000;
        }

        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential);

        //Initialize
        sim.integrator.setIsothermal(true);
        sim.integrator.setIntegratorMC(sim.integratorMC, 10000);
        sim.activityIntegrate.setMaxSteps(params.numStepsEq);
        sim.getController().actionPerformed();
        sim.getController().reset();

        //Production
        sim.integrator.setIntegratorMC(null, 0);
        sim.integrator.setIsothermal(false);
        sim.activityIntegrate.setMaxSteps(params.numSteps);


        //G
        IDataSource pTensorMeter;
        if (sim.integrator instanceof IntegratorVelocityVerlet) {
            pTensorMeter = new MeterPressureTensorFromIntegrator(sim.getSpace());
            ((MeterPressureTensorFromIntegrator) pTensorMeter).setIntegrator((IntegratorVelocityVerlet) sim.integrator);
        } else {
            pTensorMeter = new MeterPressureHardTensor(sim.getSpace());
            ((MeterPressureHardTensor) pTensorMeter).setIntegrator((IntegratorHard) sim.integrator);
            new MeterPressureHard((IntegratorHard) sim.integrator);
        }

        //Viscosity
        DataFork pTensorFork = new DataFork();
        AccumulatorPTensor pTensorAccumVisc = new AccumulatorPTensor(sim.integrator, sim.integrator.getTimeStep());
        pTensorAccumVisc.setEnabled(true);
        pTensorFork.addDataSink(pTensorAccumVisc);
        DataPumpListener pTensorAccumViscPump = new DataPumpListener(pTensorMeter, pTensorFork);


        AccumulatorAverageFixed pTensorAccumulator = new AccumulatorAverageFixed(1);
        pTensorFork.addDataSink(pTensorAccumulator);

        sim.integrator.getEventManager().addListener(pTensorAccumViscPump);

        GlassGraphic.DataProcessorTensorTrace tracer = new GlassGraphic.DataProcessorTensorTrace();
        pTensorFork.addDataSink(tracer);
        AccumulatorAverageFixed pAccumulator = new AccumulatorAverageFixed(10);
        pTensorFork.addDataSink(pAccumulator);

        //Run
        sim.getController().actionPerformed();

        int[] pIndex;
        if(params.D==2){
            pIndex = new int[] {1};
        }else{
            pIndex = new int[] {1,2,5};
        }

        //0-Pressure
        DataGroup dataP = (DataGroup)pTensorAccumulator.getData();
        IData dataPAvg = dataP.getData(pTensorAccumulator.AVERAGE.index);
        IData dataPErr = dataP.getData(pTensorAccumulator.ERROR.index);

        double pAvg = dataPAvg.getValue(0);
        double pErr = dataPErr.getValue(0);


        //1-Pressure Tensor
        DataGroup dataPTensor = (DataGroup)pTensorAccumulator.getData();
        IData dataPTensorAvg = dataPTensor.getData(pTensorAccumulator.AVERAGE.index);
        IData dataPTensorSD = dataPTensor.getData(pTensorAccumulator.STANDARD_DEVIATION.index);


        double sd2PTensor = 0;
        for(int i=0; i<pIndex.length; i++){sd2PTensor += dataPTensorSD.getValue(pIndex[i])*dataPTensorSD.getValue(pIndex[i]);}
        sd2PTensor/=pIndex.length;



      //2-Viscosity
        DataDoubleArray x = pTensorAccumVisc.getIndependentData(0);
        for (int i=0; i<pTensorAccumVisc.getData().getLength(); i++){
            double yi = pTensorAccumVisc.getData().getValue(i);
            double xi = x.getValue(i);
            if(!Double.isNaN(yi)){
                System.out.println(xi +"  "+ yi);
            }
        }

        System.out.println(pAvg +" "+ pErr +" "+ sd2PTensor*sim.box.getBoundary().volume()/params.temperature);

    }

    public static class SimParams extends SimGlass.GlassParams {
        public int numStepsEq = 10000;
        public int numSteps =   1000000;

    }
}