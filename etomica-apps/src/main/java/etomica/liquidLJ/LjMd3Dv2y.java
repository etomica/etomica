/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.action.WriteConfigurationBinary;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationFileBinary;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.space3d.Space3D;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.io.File;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LjMd3Dv2y {
    
    public static void main(String[] args) {

        // according to Mastny & de Pablo
        // http://link.aip.org/link/doi/10.1063/1.2753149
        // triple point
        // T = 0.694
        // liquid density = 0.845435
        
        // Agrawal and Kofke:
        //      small      large
        // T    0.698    0.687(4)
        // p    0.0013   0.0011
        // rho  0.854    0.850
        
        // Orkoulas, http://link.aip.org/link/doi/10.1063/1.4758698 says
        // T = 0.7085(5)
        // P = 0.002264(17)
        // rhoL = 0.8405(3)
        // rhoFCC = 0.9587(2)
        // rhoV = 0.002298(18)

        LjMd3DParams params = new LjMd3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = false;
            params.numAtoms = 500;
            params.steps = 200000;
            params.v2 = 0;
            params.rcShort = 2.5; //*Math.pow(params.v2, 1.0/6.0);
            params.y = 1.5;
            params.hybridInterval = 20;
        }

        final int numAtoms = params.numAtoms;
        final boolean ss = params.v2 == 0;
        boolean useRho = ss && params.density > 0;
        final double v2 = ss ? (useRho ? 1.0/(params.density*params.density) : 1) : params.v2;
        final double density = useRho ? params.density : 1.0/Math.sqrt(v2);
        final double y = useRho ? (4/(v2*v2)) : params.y;
        final double temperature = useRho ? 1 : 4/(y*v2*v2);
        long steps = params.steps;
        double tStep = 0.002*Math.sqrt(y)*v2 * params.tStepFac;
        double rcShort = params.rcShort;
        int nrcMax = params.nrcMax;
        boolean graphics = params.graphics;
        final int hybridInterval = params.hybridInterval;
        int nAccBlocks = params.nAccBlocks;

        int longInterval = (numAtoms / 200 / hybridInterval) * hybridInterval;
        if (longInterval == 0) longInterval = hybridInterval;
        int intervalLS = longInterval*5;

    	if (!graphics) {
            System.out.println("Running LJ MD with N="+numAtoms+(useRho?"":(" at y="+y+" v2="+params.v2)));
    	    System.out.println("  T="+temperature+" density="+density);
    	    System.out.println("time step: "+tStep);
    	    System.out.println(steps+" steps ("+(steps*tStep)+" time units)");
    	    System.out.println("short cutoff: "+rcShort);
    	    System.out.println("hybrid MC interval: "+hybridInterval);
            System.out.println("full energy interval: "+longInterval);
    	}

    	double L = Math.pow(numAtoms/density, 1.0/3.0);
        final LjMd3D sim = new LjMd3D(numAtoms, temperature, density, Double.NaN, tStep, rcShort, 0.494*L, hybridInterval, null, ss);

        String configFilename;
        if (useRho) {
	    	configFilename = String.format("configN%d_rho%4.2f", numAtoms, density);
	        File inConfigFile = new File(configFilename+".pos");
	        if (inConfigFile.exists()) {
	            ConfigurationFileBinary configFile = new ConfigurationFileBinary(configFilename);
	            configFile.initializeCoordinates(sim.box);
	            System.out.println("Continuing from previous configuration");
	        }
	        else {
	            // try 8x fewer atoms
	            String tmpConfigFilename = String.format("configN%d_rho%4.2f", numAtoms/8, density);
	            inConfigFile = new File(tmpConfigFilename+".pos");
	            if (inConfigFile.exists()) {
	                System.out.println("bringing configuration up from N="+numAtoms/8);
	                ConfigurationFileBinary configFile = new ConfigurationFileBinary(tmpConfigFilename);
	                ConfigurationFileBinary.replicate(configFile, sim.box, new int[]{2,2,2}, Space3D.getInstance());
	            }
	            else {
	                // try higher density, then lower density
	                for (int i=10; i>=-10; i--) {
	                    if (i==0) continue;
	                    double rho0 = i>0 ? (density+0.01*(11-i)) : (density+0.01*i);
	                    tmpConfigFilename = String.format("configN%d_rho%4.2f", numAtoms, rho0);
	                    inConfigFile = new File(tmpConfigFilename+".pos");
	                    if (inConfigFile.exists()) {
	                        System.out.printf("bringing configuration from rho=%4.2f\n",rho0);
	                        ConfigurationFileBinary configFile = new ConfigurationFileBinary(tmpConfigFilename);
	                        configFile.initializeCoordinates(sim.box);
	                        break;
	                    }
	                }
	            }
	        }
        }
        else {
	    	if (Double.parseDouble(String.format("%4.2f", params.v2)) != params.v2) {
	    	    throw new RuntimeException(String.format("you're just trying to cause trouble, use v2=%4.2f", params.v2));
	    	}
	    	configFilename = String.format("configN%d_V%4.2f_y%4.2f", numAtoms, params.v2, y);
	        File inConfigFile = new File(configFilename+".pos");
	        if (inConfigFile.exists()) {
	            ConfigurationFileBinary configFile = new ConfigurationFileBinary(configFilename);
	            configFile.initializeCoordinates(sim.box);
	            System.out.println("Continuing from previous configuration");
	        }
	        else {
	            // try 8x fewer atoms
	            String tmpConfigFilename = String.format("configN%d_V%4.2f_y%4.2f", numAtoms/8, params.v2, y);
	            inConfigFile = new File(tmpConfigFilename+".pos");
	            if (inConfigFile.exists()) {
	                System.out.println("bringing configuration up from N="+numAtoms/8);
	                ConfigurationFileBinary configFile = new ConfigurationFileBinary(tmpConfigFilename);
	                ConfigurationFileBinary.replicate(configFile, sim.box, new int[]{2,2,2}, Space3D.getInstance());
	            }
	            else {
	                // try lower temperature, then higher temperature
	                for (int i=10; i>=-10; i--) {
	                    if (i==0) continue;
	                    double y0 = i>0 ? (y+0.01*(11-i)) : (y+0.01*i);
	                    tmpConfigFilename = String.format("configN%d_V%4.2f_y%4.2f", numAtoms, params.v2, y0);
	                    inConfigFile = new File(tmpConfigFilename+".pos");
	                    if (inConfigFile.exists()) {
	                        System.out.printf("bringing configuration from y=%4.2f\n",y0);
	                        ConfigurationFileBinary configFile = new ConfigurationFileBinary(tmpConfigFilename);
	                        configFile.initializeCoordinates(sim.box);
	                        break;
	                    }
	                }
	            }
	        }
        }

        sim.potentialMasterList.init();

        if (!graphics) {
            long eqSteps = steps/10;
            if (eqSteps < 4000) {
                eqSteps = steps/4;
                if (eqSteps > 4000) eqSteps = 4000;
            }
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, eqSteps));

            System.out.println("equilibration finished ("+eqSteps+" steps)");
        }

        if (graphics) {
            final String APP_NAME = "LjMd3D";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3);
    
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            simGraphic.makeAndDisplayFrame(APP_NAME);

            return;
        }
        
        final MeterPotentialEnergyFromIntegrator meterEnergyFast = new MeterPotentialEnergyFromIntegrator(sim.integrator);

        double rcMax = 0.494*L;
        double fac = 1.2;
        int nCutoffs = 1 + (int)(Math.log(rcMax/rcShort)/Math.log(fac));
        if (nCutoffs-1>nrcMax) nCutoffs = 1+nrcMax;
        final double[] cutoffs = new double[nCutoffs];
        cutoffs[0] = rcShort;
        for (int i=1; i<cutoffs.length; i++) {
            cutoffs[i] = cutoffs[i-1]*1.2;
        }
        MeterPotentialEnergyCutoff meterEnergyCut = new MeterPotentialEnergyCutoff(sim.potentialMasterLongCut, sim.getSpace(), cutoffs);
        meterEnergyCut.setBox(sim.box);
        final double[] uFacCut = new double[cutoffs.length];
        IData uCut = meterEnergyCut.getData();
        double uFast0 = meterEnergyFast.getDataAsScalar();
        for (int i=0; i<cutoffs.length; i++) {
            uFacCut[i] = uCut.getValue(i) - uFast0;
        }

        final ValueCache energyFastCache = new ValueCache(meterEnergyFast, sim.integrator);

        final MeterPUCut meterPU;

        P2LennardJones p2LJ = null;

        if (params.v2 == 0 && !useRho) {
            PotentialMaster potentialMasterLJ = new PotentialMaster(sim.getSpeciesManager(), sim.box, BondingInfo.noBonding());
            p2LJ = new P2LennardJones();
            potentialMasterLJ.setPairPotential(sim.species().getLeafType(), sim.species.getLeafType(), p2LJ);
            meterPU = new MeterPUCut(sim.box, sim.potentialMasterLongCut, potentialMasterLJ, temperature, cutoffs);
        }
        else {
            meterPU = new MeterPUCut(sim.box, sim.potentialMasterLongCut, temperature, cutoffs);
        }

        long bs = steps/(longInterval*nAccBlocks);

        DataProcessorReweight puReweight = new DataProcessorReweight(temperature, energyFastCache, uFacCut, sim.box, nCutoffs);
        DataPumpListener pumpPU = new DataPumpListener(meterPU, puReweight, longInterval);
        sim.integrator.getEventManager().addListener(pumpPU);
        final AccumulatorAverageCovariance accPU = new AccumulatorAverageCovariance(bs == 0 ? 1 : bs);
        puReweight.setDataSink(accPU);

        DataProcessorReweightRatio puReweightRatio = new DataProcessorReweightRatio(nCutoffs);
        accPU.setBlockDataSink(puReweightRatio);
        AccumulatorAverageCovariance accPUBlocks = new AccumulatorAverageCovariance(1, true);
        puReweightRatio.setDataSink(accPUBlocks);


        rcMax *= 3;
        int nCutoffsLS = 1 + (int)(Math.log(rcMax/rcShort)/Math.log(fac)); 
        if (nCutoffs-1 >= nrcMax) nCutoffsLS=0;
        else if (nCutoffsLS-1>nrcMax) nCutoffsLS=1+nrcMax;
        final double[] cutoffsLS = new double[nCutoffsLS];
        PotentialMasterCell potentialMasterLS;
        PotentialMasterCell potentialMasterLJLS;
        AccumulatorAverageCovariance accPULS = null;
        AccumulatorAverageCovariance accPULSBlocks = new AccumulatorAverageCovariance(1, true);
        double[] uFacCutLS = null;
        if (nCutoffsLS>0) {
            cutoffsLS[0] = rcShort;
            for (int i=1; i<cutoffsLS.length; i++) {
                cutoffsLS[i] = cutoffsLS[i-1]*1.2;
            }

            double rcMaxLS = cutoffsLS[nCutoffsLS-1];
            int cellRange = (int)Math.ceil(2.01*rcMaxLS / sim.box.getBoundary().getBoxSize().getX(0));
            cellRange = Math.max(cellRange, 3);

            potentialMasterLS = new PotentialMasterCell(sim.getSpeciesManager(), sim.box, cellRange, BondingInfo.noBonding());
            P2SoftSphericalSumTruncated p2tLS = new P2SoftSphericalSumTruncated(rcMaxLS, sim.potential);
            potentialMasterLS.setPairPotential(sim.species.getLeafType(), sim.species.getLeafType(), p2tLS);

            final MeterPUCut meterPULS;
            potentialMasterLS.init();

            if (params.v2 == 0 && !useRho) {
                potentialMasterLJLS = new PotentialMasterCell(sim.getSpeciesManager(), sim.box, cellRange, BondingInfo.noBonding());
                P2SoftSphericalSumTruncated p2tLJLS = new P2SoftSphericalSumTruncated(rcMaxLS, p2LJ);
                potentialMasterLJLS.setPairPotential(sim.species.getLeafType(), sim.species.getLeafType(), p2tLJLS);
                potentialMasterLJLS.init();
                meterPULS = new MeterPUCut(sim.box, potentialMasterLS, potentialMasterLJLS, temperature, cutoffsLS);
            }
            else {
                meterPULS = new MeterPUCut(sim.box, potentialMasterLS, temperature, cutoffsLS);
            }
            uCut = meterPULS.getData();

            uFacCutLS = new double[cutoffsLS.length];
            for (int i=0; i<uFacCut.length; i++) {
                uFacCutLS[i] = uFacCut[i];
            }
            for (int i=uFacCut.length; i<uFacCutLS.length; i++) {
                uFacCutLS[i] = uCut.getValue(4*i)*numAtoms - uFast0;
            }
            DataProcessorReweight puLSReweight = new DataProcessorReweight(temperature, energyFastCache, uFacCutLS, sim.box, nCutoffsLS);
            DataPumpListener pumpPULS = new DataPumpListener(meterPULS, puLSReweight, intervalLS);
            sim.integrator.getEventManager().addListener(pumpPULS);
            bs = steps/(intervalLS*nAccBlocks);
            accPULS = new AccumulatorAverageCovariance(bs == 0 ? 1 : bs);
            puLSReweight.setDataSink(accPULS);

            DataProcessorReweightRatio puLSReweightRatio = new DataProcessorReweightRatio(nCutoffsLS, nCutoffs-1);
            accPULS.setBlockDataSink(puLSReweightRatio);
            puLSReweightRatio.setDataSink(accPULSBlocks);
        }


    	long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
        long t2 = System.currentTimeMillis();

        System.out.println();

        WriteConfigurationBinary writeConfig = new WriteConfigurationBinary(sim.getSpace());
        writeConfig.setFileName(configFilename+".pos");
        writeConfig.setBox(sim.box);
        writeConfig.actionPerformed();
        
        System.out.println("hybrid acceptance: "+sim.integrator.getHybridAcceptance());
        // FIXME
        long numNbrUpdates = 0; //sim.potentialMasterList.getNeighborManager(sim.box).getNumUpdates();
        System.out.printf("avg steps between nbr update: %3.1f\n",((double)steps)/numNbrUpdates);

        System.out.println();

        DataGroup dataPU = (DataGroup)accPU.getData();
        IData avgPU = dataPU.getData(accPU.AVERAGE.index);
        IData errPU = dataPU.getData(accPU.ERROR.index);
        IData covPU = dataPU.getData(accPU.BLOCK_COVARIANCE.index);
        IData corPU = dataPU.getData(accPU.BLOCK_CORRELATION.index);
        
        int j = 0;
        for (int i=0; i<cutoffs.length; i++) {

            double ulrc = uduCorrection(sim, cutoffs[i])[0];

            double avgW = avgPU.getValue(j+4);
            double errW = errPU.getValue(j+4);
            double corW = corPU.getValue(j+4);

            System.out.printf("rc: %d  A-Afast: % 22.15e  %10.4e  % 5.2f  %6.4f\n", i, (ulrc + uFacCut[i] - temperature*Math.log(avgW))/numAtoms, temperature*errW/avgW/numAtoms, corW, errW/avgW);

            j+=5;
        }
        System.out.println();
        
        if (nCutoffsLS > 0) {
            DataGroup dataPULS = (DataGroup)accPULS.getData();
            IData avgPULS = dataPULS.getData(accPULS.AVERAGE.index);
            IData errPULS = dataPULS.getData(accPULS.ERROR.index);
            IData corPULS = dataPULS.getData(accPULS.BLOCK_CORRELATION.index);

            j = 0;
            for (int i=0; i<cutoffsLS.length; i++) {

                double ulrc = uduCorrection(sim, cutoffsLS[i])[0];

                double avgW = avgPULS.getValue(j+4);
                double errW = errPULS.getValue(j+4);
                double corW = corPULS.getValue(j+4);

                System.out.printf("rcLS: %d  A-Afast: % 22.15e  %10.4e  % 5.2f  %6.4f\n", i, (ulrc + uFacCutLS[i] - temperature*Math.log(avgW))/numAtoms, temperature*errW/avgW/numAtoms, corW, errW/avgW);

                j+=5;
            }

            System.out.println();
        }
        
        DataGroup dataPU1 = (DataGroup)accPUBlocks.getData();
        IData avgPU1 = dataPU1.getData(accPUBlocks.AVERAGE.index);
        IData covPU1 = dataPU1.getData(accPUBlocks.BLOCK_COVARIANCE.index);
        IData corPU1 = dataPU1.getData(accPUBlocks.BLOCK_CORRELATION.index);
        
        int n = avgPU1.getLength();
        int nRaw = avgPU.getLength();

        j = 0;
        int jRaw = 0;
        for (int i=0; i<cutoffs.length; i++) {

            double avgW = avgPU.getValue(jRaw+4);
            double errW = errPU.getValue(jRaw+4);
            double ratioW2 = errW*errW/(avgW*avgW);

            double[] uduLRC = uduCorrection(sim, cutoffs[i]);
            double ulrc = uduLRC[0] / numAtoms;

            double avgU = avgPU.getValue(jRaw+0);
            double errU = errPU.getValue(jRaw+0);
            double corUW = covPU.getValue((jRaw+0)*nRaw+jRaw+4) / Math.sqrt(covPU.getValue((jRaw+0)*nRaw+jRaw+0) * covPU.getValue((jRaw+4)*nRaw+jRaw+4));
            double biasU = ratioW2 - errU*errW/(avgU*avgW)*corUW;
            errU = Math.abs(avgU/avgW)*Math.sqrt(errU*errU/(avgU*avgU) + ratioW2 - 2*errU*errW/(avgU*avgW)*corUW);
            avgU /= avgW;
            double corU = corPU1.getValue(j+0);
            if (!useRho) System.out.printf("rc: %d  U:       % 22.15e  %10.4e  % 10.4e  % 5.2f\n", i, ulrc + avgU, errU, avgU*biasU, corU);

            double avgP = avgPU.getValue(jRaw+1);
            double errP = errPU.getValue(jRaw+1);
            double corPW = covPU.getValue((jRaw+1)*nRaw+jRaw+4) / Math.sqrt(covPU.getValue((jRaw+1)*nRaw+jRaw+1) * covPU.getValue((jRaw+4)*nRaw+jRaw+4));
            double biasP = ratioW2 - errP*errW/(avgP*avgW)*corPW;
            errP = Math.abs(avgP/avgW)*Math.sqrt(errP*errP/(avgP*avgP) + ratioW2 - 2*errP*errW/(avgP*avgW)*corPW);
            avgP /= avgW;
            double corP = corPU1.getValue(j+1);
            double vol = sim.box.getBoundary().volume();
            double plrc = -uduLRC[1]/(3*vol);
            double puCor = covPU1.getValue((j+0)*n+j+1) / Math.sqrt(covPU1.getValue((j+0)*n+j+0) * covPU1.getValue((j+1)*n+j+1));
            if (useRho) System.out.printf("rc: %d  P:       % 22.15e  %10.4e  % 10.4e  % 5.2f\n", i, plrc + avgP, errP, avgP*biasP, corP);
            else System.out.printf("rc: %d  P:       % 22.15e  %10.4e  % 10.4e  % 5.2f  % 7.4f\n", i, plrc + avgP, errP, avgP*biasP, corP, puCor);

            if (!useRho) {
	            double avgDADy = avgPU.getValue(jRaw+2);
	            double errDADy = errPU.getValue(jRaw+2);
	            double corDADyW = covPU.getValue((jRaw+2)*nRaw+jRaw+4) / Math.sqrt(covPU.getValue((jRaw+2)*nRaw+jRaw+2) * covPU.getValue((jRaw+4)*nRaw+jRaw+4));
	            double biasDADy = ratioW2 - errDADy*errW/(avgDADy*avgW)*corDADyW;
	            errDADy = Math.abs(avgDADy/avgW)*Math.sqrt(errDADy*errDADy/(avgDADy*avgDADy) + ratioW2 - 2*errDADy*errW/(avgDADy*avgW)*corDADyW);
	            avgDADy /= avgW;
	            double corDADy = corPU1.getValue(j+2);
	            System.out.printf("rc: %d  DADy:    % 22.15e  %10.4e  % 10.4e  % 5.2f\n", i, ulrc*Math.pow(density,-4)/4 + avgDADy, errDADy, avgDADy*biasDADy, corDADy);
	
	            double avgDADv2 = avgPU.getValue(jRaw+3);
	            double errDADv2 = errPU.getValue(jRaw+3);
	            double corDADv2W = covPU.getValue((jRaw+3)*nRaw+jRaw+4) / Math.sqrt(covPU.getValue((jRaw+3)*nRaw+jRaw+3) * covPU.getValue((jRaw+4)*nRaw+jRaw+4));
	            double biasDADv2 = ratioW2 - errDADv2*errW/(avgDADv2*avgW)*corDADv2W;
	            errDADv2 = Math.abs(avgDADv2/avgW)*Math.sqrt(errDADv2*errDADv2/(avgDADv2*avgDADv2) + ratioW2 - 2*errDADv2*errW/(avgDADv2*avgW)*corDADv2W);
	            avgDADv2 /= avgW;
	            double corDADv2 = corPU1.getValue(j+3);
	
	            // -(P/(temperature*density) - 1 - 4 * U / (temperature))*density*density/2;
	
	            if (p2LJ!=null) {
                    double[] udulj = uduCorrection(sim, p2LJ, cutoffs[i]);
	                ulrc = udulj[0]/numAtoms;
	                plrc = -udulj[1]/(3*vol);
	            }
	
	            double DADv2LRC = (-plrc/(temperature*density) + 4*ulrc/temperature)*density*density/2;
	            double dadCor = covPU1.getValue((j+2)*n+j+3) / Math.sqrt(covPU1.getValue((j+2)*n+j+2) * covPU1.getValue((j+3)*n+j+3));
	            System.out.printf("rc: %d  DADv2:   % 22.15e  %10.4e  % 10.4e  % 5.2f  % 7.4f\n", i, DADv2LRC + avgDADv2, errDADv2, avgDADv2*biasDADv2, corDADv2, dadCor);
	            System.out.println();
            }

            j+=4;
            jRaw+=5;
        }

        if (nCutoffsLS > 0) {

            dataPU = (DataGroup)accPULS.getData();
            avgPU = dataPU.getData(accPULS.AVERAGE.index);

            dataPU1 = (DataGroup)accPULSBlocks.getData();
            avgPU1 = dataPU1.getData(accPULSBlocks.AVERAGE.index);
            IData errPU1 = dataPU1.getData(accPULSBlocks.ERROR.index);
            covPU1 = dataPU1.getData(accPULSBlocks.BLOCK_COVARIANCE.index);
            corPU1 = dataPU1.getData(accPULSBlocks.BLOCK_CORRELATION.index);

            n = avgPU.getLength();
            int n1 = avgPU1.getLength();

            j =  0;
            int j1 = 0;
            double ulrcRef = 0, ulrcRefLJ = 0;
            double plrcRef = 0, plrcRefLJ = 0;
            double avgUref = 0, avgPref = 0, avgDADyref = 0, avgDADv2ref = 0;
            for (int i=0; i<cutoffsLS.length; i++) {
                if (i<cutoffs.length-1) {
                    j1+=4;
                    j+=5;
                    continue;
                }

                double avgW = avgPU.getValue(j+4);

                PotentialMaster dummy = new PotentialMaster(sim.getSpeciesManager(), sim.box, BondingInfo.noBonding());
                dummy.doAllTruncationCorrection = true;
                P2SoftSphericalTruncated p2t = new P2SoftSphericalTruncated(sim.potential, cutoffsLS[i]);
                dummy.setPairPotential(sim.species.getLeafType(), sim.species.getLeafType(), p2t);
                double[] ulrca = {0}, dulrca = {0};
                dummy.computeAllTruncationCorrection(ulrca, dulrca);
                double ulrc = ulrca[0] / numAtoms;
                if (i==cutoffs.length-1) ulrcRef = ulrc;
                else ulrc -= ulrcRef;

                double avgU = avgPU.getValue(j+0)/avgW;
                double avgU1 = avgPU1.getValue(j1+0);
                double errU1 = errPU1.getValue(j1+0);
                double corU1 = corPU1.getValue(j1+0);

                double avgP = avgPU.getValue(j+1)/avgW;
                double avgP1 = avgPU1.getValue(j1+1);
                double errP1 = errPU1.getValue(j1+1);
                double corP1 = corPU1.getValue(j1+1);
                double vol = sim.box.getBoundary().volume();
                double plrc = -dulrca[0]/(3*vol);
                if (i==cutoffs.length-1) plrcRef = plrc;
                else plrc -= plrcRef;
                double puCor = covPU1.getValue((j1+0)*n1+j1+1) / Math.sqrt(covPU1.getValue((j1+0)*n1+j1+0) * covPU1.getValue((j1+1)*n1+j1+1));

                double avgDADy = avgPU.getValue(j+2)/avgW;
                double avgDADy1 = avgPU1.getValue(j1+2);
                double errDADy1 = errPU1.getValue(j1+2);
                double corDADy1 = corPU1.getValue(j1+2);
                if (i>cutoffs.length-1) {
                    avgU -= avgUref;
                    avgP -= avgPref;
                    avgDADy -= avgDADyref;
                    if (!useRho) System.out.printf("drcLS: %d  U:       % 22.15e  %10.4e  % 10.4e  % 5.2f\n", i, ulrc + avgU, errU1, (avgU-avgU1)/nAccBlocks, corU1);
                    if (useRho) System.out.printf("drcLS: %d  P:       % 22.15e  %10.4e  % 10.4e  % 5.2f\n", i, plrc + avgP, errP1, (avgP-avgP1)/nAccBlocks, corP1);
                    else System.out.printf("drcLS: %d  P:       % 22.15e  %10.4e  % 10.4e  % 5.2f  % 7.4f\n", i, plrc + avgP, errP1, (avgP-avgP1)/nAccBlocks, corP1, puCor);
                    if (!useRho) System.out.printf("drcLS: %d  DADy:    % 22.15e  %10.4e  % 10.4e  % 5.2f\n", i, ulrc*Math.pow(density,-4)/4 + avgDADy, errDADy1, (avgDADy-avgDADy1)/nAccBlocks, corDADy1);
                }

                double avgDADv2 = avgPU.getValue(j+3)/avgW;
                double avgDADv21 = avgPU1.getValue(j1+3);
                double errDADv21 = errPU1.getValue(j1+3);
                double corDADv21 = corPU1.getValue(j1+3);

                // -(P/(temperature*density) - 1 - 4 * U / (temperature))*density*density/2;
                if (p2LJ!=null) {
                    p2t = new P2SoftSphericalTruncated(p2LJ, cutoffsLS[i]);
                    dummy.setPairPotential(sim.species.getLeafType(), sim.species.getLeafType(), p2t);
                    ulrca[0] = 0;
                    dulrca[0] = 0;
                    dummy.computeAllTruncationCorrection(ulrca, dulrca);
                    ulrc = ulrca[0] / numAtoms;
                    plrc = -dulrca[0]/(3*vol);
                    if (i==cutoffs.length-1) ulrcRefLJ = ulrc;
                    else ulrc -= ulrcRefLJ;
                    if (i==cutoffs.length-1) plrcRefLJ = plrc;
                    else plrc -= plrcRefLJ;
                }

                if (i==cutoffs.length-1) {
                    avgUref = avgU;
                    avgPref = avgP;
                    avgDADyref = avgDADy;
                    avgDADv2ref = avgDADv2;
                    j1+=4;
                    j+=5;
                    continue;
                }
                avgDADv2 -= avgDADv2ref;

                if (!useRho) {
	                double DADv2LRC = (-plrc/(temperature*density) + 4*ulrc/temperature)*density*density/2;
	                double dadCor = covPU1.getValue((j1+2)*n1+j1+3) / Math.sqrt(covPU1.getValue((j1+2)*n1+j1+2) * covPU1.getValue((j1+3)*n1+j1+3));
	                System.out.printf("drcLS: %d  DADv2:   % 22.15e  %10.4e  % 10.4e  % 5.2f  % 7.4f\n", i, DADv2LRC + avgDADv2, errDADv21, (avgDADv21-avgDADv2)/nAccBlocks, corDADv21, dadCor);
	                System.out.println();
                }

                j1+=4;
                j+=5;
            }
        }

        System.out.println("time: "+(t2-t1)/1000.0+" seconds");
    }

    public static class LjMd3DParams extends ParameterBase {
        public int numAtoms = 2048;
        public double y = 1.3;
        public double v2 = 0.8;
        public double tStepFac = 1;
        public long steps = 40000;
        public int hybridInterval = 20;
        public double rcShort = 2.5;
        public int nrcMax = 100;
        public boolean graphics = false;
        public int nAccBlocks = 100;
        public double density = 0;
    }

    public static double[] uduCorrection(LjMd3D sim, double rc) {
        return uduCorrection(sim, sim.potential, rc);
    }

    public static double[] uduCorrection(LjMd3D sim, IPotential2 p, double rc) {
        P2SoftSphericalSumTruncated pt = new P2SoftSphericalSumTruncated(rc, p);
        int N = sim.box.getMoleculeList().size();
        double V = sim.box.getBoundary().volume();
        int numPairs = N * (N - 1) / 2;
        double pairDensity = numPairs / V;

        double[] u = new double[1];
        double[] du = new double[1];
        pt.u01TruncationCorrection(sim.getSpace(), u, du);
        return new double[]{u[0] * pairDensity, du[0]*pairDensity};
    }
}
