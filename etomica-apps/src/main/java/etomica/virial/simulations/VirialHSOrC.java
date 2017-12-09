/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import etomica.math.function.IFunction;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataDistributer;
import etomica.data.DataFork;
import etomica.data.DataSplitter;
import etomica.data.IData;
import etomica.data.IDataSink;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.CalcFFT;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterWheatleyHS;
import etomica.virial.MCMoveClusterAtomDiscrete;
import etomica.virial.MayerFunction;
import etomica.virial.MayerHardSphere;

/**
 * Overlap sampling simulation to compute c(r) for HS
 */
public class VirialHSOrC {


    public static void main(String[] args) {

        VirialLJrParam params = new VirialLJrParam();
        ParseArgs.doParseArgs(params, args);
        boolean isCommandline = args.length > 0;
        
        if (args.length == 0) {
            params.numSteps = 10000000;
            params.nPoints = 4;
            params.temperature = 1.0;
            params.dr = 0.01;
            params.rPow = 0;
        }

        final int nPoints = params.nPoints;  
        final double temperature = params.temperature;
        long steps = params.numSteps;
        final long maxBlockSize = params.maxBlockSize;
        final double dr = params.dr;
        final double refFrac = params.refFrac;
        final double rPow = params.rPow;

        System.out.println("Calculating HS "+nPoints+" distributions ");
        System.out.println("r01 weight power: "+rPow);
        
        Space space = Space3D.getInstance();
        
        MayerFunction fTarget = new MayerHardSphere();
        final MayerFunction fRef = fTarget;
        
        // we want to get to ~40sigma in steps of dr
        int N = (int)Math.round(40/dr);
        int power = (int)Math.ceil(Math.log(N)/Math.log(2));
        N = 1<<power;
        int powerOffset = 0;
        if (power < 20) {
            // we will compute any distribution functions we need with 2^20 points
            powerOffset = 20-power;
        }

        int[][] rBondList = new int[0][0];
        List<Integer> strandList = new ArrayList<Integer>();
        for (int i=2; i<nPoints; i++) {
            rBondList = (int[][])Arrays.addObject(rBondList, new int[]{0,i});
            rBondList = (int[][])Arrays.addObject(rBondList, new int[]{1,i});
            strandList.add(3);
        }
        CalcFFT calcFFTRef = new CalcFFT(new IFunction() {
            public double f(double x) {
                return fRef.f(null, x*x, 1.0/temperature);
            }
        }, dr/(1<<powerOffset), power+powerOffset);
        double[] refFFT = calcFFTRef.value(strandList, false)[0];
        final double[] rRef = new double[N];
        for (int i=0; i<N; i++) {
            rRef[i] = refFFT[i*(1<<powerOffset)];
        }
        calcFFTRef = null;
        refFFT = null;
        ClusterBonds refBonds = new ClusterBonds(nPoints, new int[][][]{rBondList});
        ClusterAbstract refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double[]{1}, new MayerFunction[]{fRef});
        
        ClusterAbstract targetCluster = new ClusterWheatleyHS(nPoints, fTarget);

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        System.out.println(steps+" steps");
		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("A")), temperature, refCluster, targetCluster);
        sim.integratorOS.setNumSubSteps(1000);
        
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        long blockSize = steps/100000L;
        if (blockSize == 0) blockSize = 1;
        else if (steps > maxBlockSize) blockSize = maxBlockSize;
        final long bs = blockSize;
        System.out.println("block size: "+blockSize);

        DataDistributer.Indexer indexer = null;
        DataSplitter.IDataSinkFactory accFac = new DataSplitter.IDataSinkFactory() {
            public IDataSink makeDataSink(int i) {
                AccumulatorAverageFixed a = new AccumulatorAverageFixed(bs);
                a.setDoStrictBlockData(true);
                return a;
            }
        };

        MCMoveClusterAtomDiscrete[] mcDiscrete = new MCMoveClusterAtomDiscrete[2];
        for (int i=0; i<2; i++) {
            mcDiscrete[i] = new MCMoveClusterAtomDiscrete(sim.getRandom(), space, dr);
            mcDiscrete[i].setRPow(rPow);
            sim.integrators[i].getMoveManager().addMCMove(mcDiscrete[i]);
            sim.integrators[i].getMoveManager().removeMCMove(sim.mcMoveTranslate[i]);
            sim.box[i].getLeafList().getAtom(1).getPosition().setX(0, dr*Math.round(0.5/dr));
            if (nPoints>2) sim.box[i].getLeafList().getAtom(2).getPosition().setX(1, dr*Math.round(0.5/dr));
            sim.box[i].trialNotify();
            sim.box[i].acceptNotify();

        }
        indexer = new DataDistributer.Indexer() {
            public int getIndex() {
                double x = sim.box[1].getLeafList().getAtom(1).getPosition().getX(0);
                return (int)Math.round(Math.abs(x)/dr);
            }
        };
        DataDistributer dataDistributer = new DataDistributer(indexer, accFac);
        DataFork targetFork = new DataFork();
        sim.accumulatorPumps[1].setDataSink(targetFork);
        targetFork.addDataSink(sim.dpVirialOverlap[1]);
        
        targetFork.addDataSink(dataDistributer);
        
        steps /= 1000;

        boolean doGraphics = false;
        if (doGraphics) {
            double size = 10;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox = simGraphic.getDisplayBox(sim.box[0]);
            displayBox.setPixelUnit(new Pixel(300.0/size));
            displayBox.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox.canvas).setBackgroundColor(Color.WHITE);
            simGraphic.makeAndDisplayFrame();
            return;
        }

        long t1 = System.currentTimeMillis();
        // if running interactively, don't use the file
        String refFileName = isCommandline ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/20);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/10);
        
        System.out.println("equilibration finished");
        
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        
        int nAcc = dataDistributer.getNumDataSinks();
        for (int i=0; i<nAcc; i++) {
            AccumulatorAverage acc = (AccumulatorAverage)dataDistributer.getDataSink(i);
            if (acc != null) {
                acc.reset();
            }
        }

        sim.integratorOS.setNumSubSteps((int)steps);
        sim.setAccumulatorBlockSize(steps);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+mcDiscrete[i].getStepSize());
        }
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

        int digits = (int)Math.ceil(-Math.log10(dr));
        String rfmt = "%"+(4+digits)+"."+digits+"f";
        nAcc = dataDistributer.getNumDataSinks();
        
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        System.out.println((long)(steps*1000*(1.0-sim.integratorOS.getRefStepFraction())/bs)+" target blocks");
        long targetSteps = (long)(steps*1000*(1.0-sim.integratorOS.getRefStepFraction()));
        double refIntegral = 0;
        for (int i=0; i<rRef.length; i++) {
            refIntegral -= rRef[i]*Math.pow(i*dr,rPow)*nPoints;
        }
        System.out.println("HSB: "+refIntegral);
        sim.printResults(refIntegral);
        System.out.println();
        DataGroup allYourBaseRef = (DataGroup)sim.accumulators[0].getData();
        IData refAverageData = allYourBaseRef.getData(sim.accumulators[0].AVERAGE.index);
        double ra = refAverageData.getValue(0);
        double roa = refAverageData.getValue(1);
        DataGroup allYourBaseTarget = (DataGroup)sim.accumulators[1].getData();
        IData targetAverageData = allYourBaseTarget.getData(sim.accumulators[1].AVERAGE.index);
        double toa = targetAverageData.getValue(1);
        for (int i=0; i<nAcc; i++) {
            AccumulatorAverageFixed accumulator = (AccumulatorAverageFixed)dataDistributer.getDataSink(i);
            if (accumulator == null || accumulator.getBlockCount() == 0) {
                if (dr*i > 1.1) break;
                System.out.print(String.format("r: "+rfmt, +dr*i));
                System.out.print("  0 blocks  0 steps  ");
                System.out.print(String.format("avg: % 20.15e  %9.4e   c(r): % 20.15e  %9.4e   cor: % 6.4f\n", 0.0, 0.0, 0.0, 0.0, 0.0));
                continue;
            }
            System.out.print(String.format("r: "+rfmt, +dr*i));
            long iSteps = accumulator.getSampleCount();
            System.out.print("  "+accumulator.getBlockCount()+" blocks  "+iSteps+" steps  ");
            DataGroup allYourBase = (DataGroup)accumulator.getData();
            IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
            IData errorData = allYourBase.getData(accumulator.ERROR.index);
            IData correlationData = allYourBase.getData(accumulator.BLOCK_CORRELATION.index);
            
            if (Double.isNaN(averageData.getValue(0))) {
                accumulator.getData();
                throw new RuntimeException("oops");
            }

            double avg = averageData.getValue(0);
            double err = errorData.getValue(0);
            double fac = iSteps*roa*refIntegral/(toa*ra*targetSteps);
            if (rPow>0||i>0) fac*=Math.pow(dr*i,-rPow);
            
            System.out.print(String.format("avg: % 20.15e  %9.4e   c(r): % 20.15e  %9.4e   cor: % 6.4f\n",
                    avg, err, avg*fac, err*Math.abs(fac), correlationData.getValue(0)));
        }

        System.out.println("time: "+(t2-t1)/1000.0);
    }

    public enum PotentialChoice { LJ, WCA, SS, HS };
    
    /**
     * Inner class for parameters
     */
    public static class VirialLJrParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 1;
        public long numSteps = 100000000L;
        public double dr = 0.01;
        public long maxBlockSize = 100;
        public double refFrac = -1;
        public double rPow = 2;
    }
}
