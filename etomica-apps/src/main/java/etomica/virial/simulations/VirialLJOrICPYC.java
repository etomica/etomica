/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.*;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.math.function.IFunction;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphere;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;
import etomica.util.Arrays;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.VirialDiagrams;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialLJOrICPYC {


    public static void main(String[] args) {

        VirialLJrParam params = new VirialLJrParam();
        ParseArgs.doParseArgs(params, args);
        boolean isCommandline = args.length > 0;
        
        if (args.length == 0) {
            params.numSteps = 1000000000;
            params.nPoints = 6;
        }

        final int nPoints = params.nPoints;  
        final double temperature = params.temperature;
        long steps = params.numSteps;
        final long maxBlockSize = params.maxBlockSize;
        final double dr = params.dr;
        final double refFrac = params.refFrac;
        final double rPow = params.rPow;

        System.out.println("Calculating LJ "+nPoints+" distributions ");
        System.out.println("r01 weight power: "+rPow);
        
        Space space = Space3D.getInstance();
        
        Potential2Spherical pTarget = new P2LennardJones(space,1.0,1.0);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        final MayerSphericalPlus fRef = new MayerSphericalPlus(pTarget, new P2SoftSphere(space, 1, 1, 20));
        
        // we want to get to ~40sigma in steps of dr
        int N = (int)Math.round(40/dr);
        int power = (int)Math.ceil(Math.log(N)/Math.log(2));
        N = 1<<power;
        int powerOffset = 0;
        if (power < 20) {
            // we will compute any distribution functions we need with 2^20 points
            powerOffset = 20-power;
        }

        final double[] rWeight = new double[N];
       
        // we just need a functioning VirialDiagrams object to make our target cluster
        VirialDiagrams diagrams = new VirialDiagrams(nPoints, false, false);
        // we need to poke diagrams enough so that it initializes itself properly
        diagrams.makeRhoDiagrams();
        for (int i=0; i<N; i++) {
            rWeight[i] = (fTarget.f(null, i*i*dr*dr, 1.0/temperature)+1);
        }
    
        MayerFunction f12 = new MayerFunction() {
            public void setBox(Box box) {}
            public IPotential getPotential() {return null;}
            
            public double f(IMoleculeList pair, double r2, double beta) {
                double r = Math.sqrt(r2);
                int ir = (int)Math.round(r/dr);
                return rWeight[ir];
            }
        };
        

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
            rRef[i] = refFFT[i*(1<<powerOffset)]*rWeight[i];
//            if (i*dr <= 10) System.out.println(i*dr+" "+rWeight[i]);
        }
//        System.exit(1);
        calcFFTRef = null;
        refFFT = null;
        ClusterBonds refBonds = new ClusterBonds(nPoints, new int[][][]{rBondList,new int[][]{{0,1}}});
        ClusterAbstract refCluster = new ClusterSum(new ClusterBonds[]{refBonds}, new double[]{1}, new MayerFunction[]{fRef,f12});
        
             
        ClusterAbstract targetCluster = diagrams.makeVirialCluster(PYGenerator.getICPYCorrection((byte)nPoints, true), fTarget, null);
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
            ((MCMoveClusterAtomMulti)sim.mcMoveTranslate[i]).setStartAtom(2);
            sim.box[i].getLeafList().getAtom(1).getPosition().setX(0, dr*Math.round(1.5/dr));
            sim.box[i].getLeafList().getAtom(2).getPosition().setX(0, dr*Math.round(1.5/dr));
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

        if (false) {
//            for (int i=0; i<2; i++) {
//                Box box = sim.box[i];
//                final IAtomList atoms = box.getLeafList();
//                IAtom atom1 = atoms.getAtom(1);
//                IVectorMutable p1 = atom1.getPosition();
//                p1.setX(0, r12);
//                ((MCMoveClusterAtomMulti)sim.mcMoveTranslate[i]).setStartAtom(2);
//                sim.box[i].trialNotify();
//                sim.box[i].acceptNotify();
//            }
        }
        
        if (false) {
            double size = 10;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
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
        //sim.setRefPref(1.0082398078547523);
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
        
//        IAction progressReport = new IAction() {
//            public void actionPerformed() {
//                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
//                double ratio = sim.dsvo.getDataAsScalar();
//                double error = sim.dsvo.getError();
//                System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
//            }
//        };
//        sim.integratorOS.addIntervalAction(progressReport);
//        sim.integratorOS.setActionInterval(progressReport, (int)(steps/10));

        sim.integratorOS.setNumSubSteps((int)steps);
        sim.setAccumulatorBlockSize(steps);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+mcDiscrete[i].getStepSize());
        }
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

        int digits = (int)Math.ceil(-Math.log10(dr));
        String rfmt = "%"+(4+digits)+"."+digits+"f";
        nAcc = dataDistributer.getNumDataSinks();
        
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        System.out.println((long)(steps*1000*(1.0-sim.integratorOS.getRefStepFraction())/bs)+" target blocks");
        double refIntegral = 0;
        for (int i=0; i<rRef.length; i++) {
            refIntegral += rRef[i]*Math.pow(i*dr,rPow);
        }
        System.out.println("HSB: "+refIntegral);
        sim.printResults(refIntegral);
        System.out.println();
        for (int i=0; i<nAcc; i++) {
            AccumulatorAverageFixed accumulator = (AccumulatorAverageFixed)dataDistributer.getDataSink(i);
            if (accumulator == null || accumulator.getBlockCount() == 0) {
                if (dr*i > 1.1) break;
                System.out.print(String.format("r: "+rfmt, +dr*i));
                System.out.print("  0 blocks  0 steps  ");
                System.out.print(String.format("average: % 20.15e  error: %9.4e  cor: % 6.4f\n", 0.0, 0.0, 0.0));
                continue;
            }
            System.out.print(String.format("r: "+rfmt, +dr*i));
            System.out.print("  "+accumulator.getBlockCount()+" blocks  "+accumulator.getSampleCount()+" steps  ");
            DataGroup allYourBase = (DataGroup)accumulator.getData();
            IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
            IData errorData = allYourBase.getData(accumulator.ERROR.index);
            IData correlationData = allYourBase.getData(accumulator.BLOCK_CORRELATION.index);
            
            if (Double.isNaN(averageData.getValue(0))) {
                accumulator.getData();
                throw new RuntimeException("oops");
            }
            
//            System.out.print(String.format("reference average: %20.15e  error: %9.4e  cor: %6.4f\n",
//                    averageData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));
//            if (Math.abs(averageData.getValue(0)-1) > 1e-10) throw new RuntimeException("oops");

            System.out.print(String.format("average: % 20.15e  error: %9.4e  cor: % 6.4f\n",
                    averageData.getValue(0), errorData.getValue(0), correlationData.getValue(0)));

//            int nData = averageData.getLength();
//            int nCovData = covarianceData.getLength();
//            if (nData*nData != nCovData) {
//                // we need to know these to grab the right elements of covarianceData
//                throw new RuntimeException("unexpected number of data values");
//            }
//            double correlationCoef = covarianceData.getValue(1)/Math.sqrt(covarianceData.getValue(0)*covarianceData.getValue(nData+1));
//            correlationCoef = (Double.isNaN(correlationCoef) || Double.isInfinite(correlationCoef)) ? 0 : correlationCoef;
//    
//            System.out.print(String.format("ratio average: %20.15e  error: %9.4e  cor: %6.4f\n", ratioData.getValue(1), ratioErrorData.getValue(1), correlationCoef));
//            System.out.print(String.format("abs average: %20.15e  error: %9.4e\n", ratioData.getValue(1)*refIntegral, ratioErrorData.getValue(1)*refIntegral));
        }

        System.out.println("time: "+(t2-t1)/1000.0);
    }

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
