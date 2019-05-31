/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.meter.MeterPressureHard;
import etomica.data.meter.MeterPressureHardTensor;
import etomica.data.meter.MeterPressureTensorFromIntegrator;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.util.ParseArgs;

import java.io.FileWriter;
import java.io.IOException;
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
            params.density = 0.95492965855;//0.8442;
            params.D = 3;
            params.temperature = 0.722;
            params.numStepsEq = 1000;
            params.numSteps =   10000;
        }

        SimGlass sim = new SimGlass(params.D, params.nA, params.nB, params.density, params.temperature, params.doSwap, params.potential);
        System.out.println(sim.potentialChoice);
        System.out.println( params.numSteps + " MD steps after " + params.numStepsEq + " equilibaration steps");

        //Initialize
        sim.integrator.setIsothermal(true);
        sim.integrator.setIntegratorMC(sim.integratorMC, 1000);
        sim.integrator.setTemperature(params.temperature);
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

        long blocksize = params.numSteps/100;
        GlassGraphic.DataProcessorTensorTrace tracer = new GlassGraphic.DataProcessorTensorTrace();
        pTensorFork.addDataSink(tracer);
        AccumulatorAverageFixed pAccumulator = new AccumulatorAverageFixed(blocksize);
        pTensorFork.addDataSink(pAccumulator);

        AccumulatorAverageFixed tAccumulator = null;
        if(sim.potentialChoice != SimGlass.PotentialChoice.HS){
            MeterTemperature tMeter = new MeterTemperature(sim, sim.box, params.D);
            tAccumulator = new AccumulatorAverageFixed(blocksize);
            DataPumpListener tPump = new DataPumpListener(tMeter, tAccumulator, 10);
            tAccumulator.addDataSink(pTensorAccumVisc.makeTemperatureSink(),  new AccumulatorAverage.StatType[]{tAccumulator.AVERAGE});
            sim.integrator.getEventManager().addListener(tPump);
        }


        //Run
        sim.getController().actionPerformed();

        int[] pIndex;
        if(params.D==2){
            pIndex = new int[] {1};
        }else{
            pIndex = new int[] {1,2,5};
        }

        //Pressure
        DataGroup dataP = (DataGroup)pAccumulator.getData();
        IData dataPAvg = dataP.getData(pAccumulator.AVERAGE.index);
        IData dataPErr = dataP.getData(pAccumulator.ERROR.index);
        IData dataPCorr = dataP.getData(pAccumulator.BLOCK_CORRELATION.index);
        double pAvg  = dataPAvg.getValue(0);
        double pErr  = dataPErr.getValue(0);
        double pCorr = dataPCorr.getValue(0);

        //Pressure Tensor (G_inf)
        DataGroup dataPTensor = (DataGroup)pTensorAccumulator.getData();
        IData dataPTensorAvg = dataPTensor.getData(pTensorAccumulator.AVERAGE.index);
        IData dataPTensorSD = dataPTensor.getData(pTensorAccumulator.STANDARD_DEVIATION.index);
        double sd2PTensor = 0;
        for(int i=0; i<pIndex.length; i++){sd2PTensor += dataPTensorSD.getValue(pIndex[i])*dataPTensorSD.getValue(pIndex[i]);}
        sd2PTensor/=pIndex.length;

        double volume = sim.box.getBoundary().volume();
        int numAtoms = params.nA + params.nB;
        double rho= numAtoms/volume;
        System.out.println("rho: " + rho);
        String fileName;

        if(sim.potentialChoice == SimGlass.PotentialChoice.HS){
            double phi;
            if(params.D == 2){
                phi = Math.PI/4*(params.nA+params.nB/(1.4*1.4))/volume;
            }else{
                phi = Math.PI/6*(params.nA+params.nB/(1.4*1.4*1.4))/volume;
            }
            System.out.println("phi: " + phi);
            fileName = String.format("visc%1dDRho%1.3f", params.D,rho);
        }else{
            DataGroup dataT = (DataGroup)tAccumulator.getData();
            IData dataTAvg = dataT.getData(tAccumulator.AVERAGE.index);
            IData dataTErr = dataT.getData(tAccumulator.ERROR.index);
            IData dataTCorr = dataT.getData(tAccumulator.BLOCK_CORRELATION.index);
            double tAvg  = dataTAvg.getValue(0);
            double tErr  = dataTErr.getValue(0);
            double tCorr = dataTCorr.getValue(0);
            System.out.println("T: " + tAvg +"  "+ tErr +"  cor: "+tCorr);
            fileName = String.format("visc%1dDPho%1.3fT%1.3f", params.D, rho, params.temperature);
            System.out.println("G: " + sim.box.getBoundary().volume()/tAvg*sd2PTensor+"\n");
        }
        System.out.println("P: " + pAvg +"  "+ pErr +"  cor: "+pCorr);

        //Viscosity
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter(fileName, false);
            DataDoubleArray x = pTensorAccumVisc.getIndependentData(0);
            for (int i=0; i<pTensorAccumVisc.getData().getLength(); i++){
                double yi = pTensorAccumVisc.getData().getValue(i);
                double xi = x.getValue(i);
                if(!Double.isNaN(yi)){
                    fileWriter.write(xi + " " + yi + "\n");
                }
            }
            fileWriter.close();
        } catch (IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
    }

    public static class SimParams extends SimGlass.GlassParams {
        public int numStepsEq = 10000;
        public int numSteps =   1000000;

    }
}