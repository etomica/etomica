/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.vdependent;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.IntSet;
import etomica.virial.IntSet.PropertyBin;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerFunction;
import etomica.virial.MeterVirialEBinMultiThreaded;
import etomica.virial.MeterVirialEBinMultiThreaded.MyData;
import etomica.virial.MeterVirialEBinMultiThreaded.MyDataCov;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.MCMoveClusterAtomHSChain;
import etomica.virial.mcmove.MCMoveClusterAtomHSTree;
import etomica.virial.mcmove.RandomPositionLatticeSQW;
import etomica.virial.simulations.SimulationVirial;
import etomica.virial.simulations.hardsphere.VirialHSBinMultiThreaded.DooDad;
import etomica.virial.wheatley.ClusterWheatleyExtendSW;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Calculation for virial coefficients of hard spheres
 */
public class VirialLatticeSQWBin {


    public static void main(String[] args) {

        VirialHSBinParam params = new VirialHSBinParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.nPoints = 4;
            params.numSteps = 10000000L;
            params.chainFrac = 0.5;
            params.nThreads = 2;
            params.L = 5;
        }

        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        double targetTemp = params.targetTemp;
        final int L = params.L;
        final double Y = targetTemp == 0 ? 1 : (Math.exp(1 / targetTemp) - 1);
        final double chainFrac = params.chainFrac;
        final int nThreads = params.nThreads;
        final double w = params.w;
        final int[] allRandomSeeds = params.randomSeeds;
        final boolean doCov = params.doCov;
        boolean shareData = nThreads == 1 || params.shareData;

        final double vhs = L == 2 ? 3 : 5;

        System.out.println("SQW lattice sampling B"+nPoints);
        System.out.println("Box Length: " + L);

        final Space space = Space.getInstance(params.D);

        Boundary b = (L < Double.POSITIVE_INFINITY && L > 0) ? new BoundaryRectangularPeriodic(space, L) : new BoundaryRectangularNonperiodic(space);

        MayerFunction fTargete2 = new MayerEHardSphere(0.9);
        fTargete2 = new MayerVDependent(fTargete2, b);
        MayerFunction fTargetf1 = new MayerFunction() {
            final double sigma2 = 0.9;
            final double well2 = 1.2*1.2;
            
            public void setBox(Box box) {}

            public double f(IMoleculeList pair, double r2, double beta) {
                if (r2 < sigma2 || r2 > well2) return 0;
                return Y;
            }
        };
        fTargetf1 = new MayerVDependent(fTargetf1, b);
        MayerFunction fRefPos = new MayerFunction() {
            
            public void setBox(Box box) {}

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < 0.9 ? 1 : (r2 < 1.2 ? Y : 0);
            }
        };
        fRefPos = new MayerVDependent(fRefPos, b);


        Map<IntSet, MyData> allMyData = new HashMap<IntSet, MyData>();
        boolean doReweight = true;
        if (w < 0) {
            // meter0 exists only to allow us to read in the weights file (if it exists) and to populate allMyData
            // this needs to be done (finished) before any thread gets started.
        	MeterVirialEBinMultiThreaded meter0 = new MeterVirialEBinMultiThreaded(null, null, null, new long[1], allMyData, 0, true, nPoints);
        	meter0.setDoCov(doCov);
            meter0.readWeights(params.runName+nPoints+"_weights.dat", nPoints);
            doReweight = allMyData.size() == 0;
        }
        if (!shareData) {
            if (doReweight) {
                // reweighting is needed, we have to share data
                shareData = true;
                System.out.println("Sharing data between threads");
            }
            else {
                // reweighting not needed, we don't want to share data
                allMyData = null;
            }
        }
        if (!shareData) {
            System.out.println("Not sharing data between threads");
        }
        System.out.println("using a ring/chain/tree reference");
        System.out.println("  chains: "+chainFrac);
        System.out.println("  trees: "+(1-chainFrac));

        double refIntegral = 1;

        // times fit to data from B4 up to B10
        // time to generate 10^ configurations of each type
        double tsChain = 0.110096 * Math.pow(nPoints, 1.76928);
        double tsTree = 0.117435 * Math.pow(nPoints, 1.77962);
        double ts = (1-chainFrac)*tsTree + chainFrac*tsChain;
        // tc, time to compute SQW fB for 10^6 configurations
        double tc = 0.480327 + 0.000238842*Math.pow(nPoints,3.15833)*Math.pow(3.56909,nPoints);
        double tRatio = tc/ts;
        System.out.println("tRatio: "+tRatio+"   tc: "+tc+"  ts: "+ts);

//        MeterVirialBDBinMultiThreadedOld.setTRatio(tRatio);


        long t1 = System.currentTimeMillis();
        SimulationWorker[] sw = new SimulationWorker[nThreads];
        long[] totalCount = new long[nThreads];
        for (int it=0; it<nThreads; it++) {
            int[] mySeeds = allRandomSeeds.length > 0 ? new int[allRandomSeeds.length/nThreads] : null;
            if (mySeeds != null) {
                System.arraycopy(allRandomSeeds, allRandomSeeds.length/nThreads*it, mySeeds, 0, mySeeds.length);
            }
            sw[it] = new SimulationWorker(it, nPoints, fTargetf1, fTargete2, fRefPos, vhs, chainFrac, steps, space, params.runName, tRatio, allMyData, w, totalCount, doReweight, mySeeds, doCov, targetTemp, L);
        }
        for (int it=0; it<nThreads; it++) {
            sw[it].start();
        }
        for (int it=0; it<nThreads; it++) {
            try {
                sw[it].join();
            }
            catch (InterruptedException ex) {
                System.out.println("interrupted "+it+"?  really?");
            }
        }

        long t2 = System.currentTimeMillis();
        System.out.println("total steps: "+(nThreads*steps));

        File out = null;
        String filename = null;
        int rn = 1;
        while (true) {
            filename = params.runName+nPoints+"_run"+rn+"_raw.dat";
            out = new File(filename);
            try {
                if (out.createNewFile()) break;
            }
            catch (IOException e) {
                throw new RuntimeException("couldn't create raw output");
            }
            rn++;
        }
        if (allMyData == null) {
            // data is not shared.  merge it now.
            for (int i=1; i<nThreads; i++) {
                sw[0].meter.mergeData(sw[i].meter.getAllMyData());
            }
            allMyData = sw[0].meter.getAllMyData();
        }
        System.out.println(allMyData.size()+" sets");
        MeterVirialEBinMultiThreaded.writeData(filename, allMyData, nThreads*steps, nPoints);
        
        if (doReweight) {
            System.out.println();
            
            MeterVirialEBinMultiThreaded.setQuiet(!doReweight);
            MeterVirialEBinMultiThreaded.recomputeWeights(allMyData, nThreads*steps, nPoints);
        }

        List<IntSet> pvs = new ArrayList<IntSet>();
        pvs.addAll(allMyData.keySet());
        Collections.sort(pvs);
        
        long totalSampleCount = 0;
        long totalNotScreenedCount = 0;
        System.out.println();
        int nn = 1+nPoints*(nPoints-1)/2;
        double[][] cov = new double[nn][nn];
        double[] isum = new double[nn];

        for (int i=0; i<nn; i++) {
            double sum = 0;
            double E0a2 = 0;
            double sumErrStdev = 0;
            for (IntSet pv : pvs) {
                MyData amd = allMyData.get(pv);
                long c = amd.unscreenedCount;

                if (i==0) totalNotScreenedCount += c;
                long sc = amd.sampleCount;
                if (sc == 0) {
                    continue;
                }
    
                if (i==0) totalSampleCount += sc;
    
                double avg = amd.getAvg(i);
                if (i==1 && c>0 && avg!=0) throw new RuntimeException("oops "+i+" "+c+" "+avg);
                double var = amd.getVar(i);
                sum += c*avg;
                E0a2 += c*avg*avg;
                sumErrStdev += var/sc*c*c;

                if (doCov) {
                    for (int j=0; j<nn; j++) {
                        cov[i][j] += c*((MyDataCov)amd).getCov(i,j)
                                   + c*avg*amd.getAvg(j);
                    }
                }
            }
            if (sum == 0) continue;
            if (i==1) throw new RuntimeException("oops "+i+" "+sum);
            sum /= nThreads*steps;
            isum[i] = sum;

            if (doReweight) {
	            double sumErrNum = E0a2 - sum*sum*nThreads*steps;
	            double finalErr = Math.sqrt(sumErrStdev + sumErrNum)*Math.abs(refIntegral)/(nThreads*steps);
                sum *= refIntegral;
                sum /= Math.pow(Y, i);
                finalErr /= Math.pow(Y, i);

                System.out.print(String.format("%2d average: %21.14e   error: %11.5e   # var frac: %5.3f\n", i, sum, finalErr, sumErrNum/(sumErrStdev + sumErrNum)));
            }
        }
//	            System.out.println("Difficulty: "+(finalErr*Math.sqrt(t2-t1)));
        if (doReweight && doCov) {
            System.out.println("\nCorrelations:");
            for (int j=0; j<nn; j++) {
                for (int k=0; k<nn; k++) {
                    cov[j][k] -= isum[j]*isum[k]*(nThreads*steps);
                }
            }
            for (int j=0; j<nn; j++) {
                System.out.print(String.format("%2d ", j));
                for (int k=0; k<nn; k++) {
                    double cor = cov[j][k];
                    double d = Math.sqrt(cov[j][j]*cov[k][k]);
                    if (cor!=0) { // && Math.abs(cor) < 10000*d) {
                        cor /= d;
                    }
                    else {
                        cor = 0;
                    }
                    System.out.print(String.format(" % 6.4f", cor));
                }
                System.out.print("\n");
            }
            System.out.println();
        }

        System.out.println("number time fraction: "+(nThreads*steps)/(nThreads*steps + totalSampleCount*tRatio));
        System.out.println("fraction not screened: "+((double)totalNotScreenedCount)/(nThreads*steps));
        System.out.println("fraction measured: "+((double)totalSampleCount)/totalNotScreenedCount);
        
        System.out.println(String.format("expected time: %d\n",(int)((steps*ts+totalSampleCount/nThreads*tc)/1e6)));
        System.out.println("time: "+(t2-t1)/1000.0);
//        }
    }

    public static class SimulationWorker extends Thread {
        
        protected final int nPoints;
        protected final MayerFunction fTargetf1, fTargete2;
        protected final MayerFunction fRefPos;
        protected final double vhs;
        protected final double chainFrac;
        protected final long steps;
        protected final Space space;
        protected final String runName;
        protected final double tRatio;
        protected final Map<IntSet, MyData> allMyData;
        protected final int iThread;
        protected final double w;
        protected final long[] totalCount;
        protected final boolean doReweight;
        protected final int[] mySeeds;
        protected final boolean doCov;
        protected final double targetTemp;
        public MeterVirialEBinMultiThreaded meter;
        protected final int L;
        
        public SimulationWorker(int iThread, int nPoints, MayerFunction fTargetf1, MayerFunction fTargete2,
                                MayerFunction fRefPos, double vhs, double chainFrac,
                                long steps, Space space, String runName, double tRatio,
                                Map<IntSet, MyData> allMyData, double w, long[] totalCount,
                                boolean doReweight, int[] mySeeds, boolean doCov, double targetTemp,
                                int L) {
            this.iThread = iThread;
            this.nPoints = nPoints;
            this.fTargetf1 = fTargetf1;
            this.fTargete2 = fTargete2;
            this.fRefPos = fRefPos;
            this.vhs = vhs;
            this.chainFrac = chainFrac;
            this.steps = steps;
            this.space = space;
            this.runName = runName;
            this.tRatio = tRatio;
            this.allMyData = allMyData;
            this.w = w;
            this.totalCount = totalCount;
            this.doReweight = doReweight;
            this.mySeeds = mySeeds;
            this.doCov = doCov;
            this.targetTemp = targetTemp;
            this.L = L;
        }
        
        public void run() {
            long t1 = System.currentTimeMillis();
            double Y = targetTemp == 0 ? 1 : (Math.exp(1 / targetTemp) - 1);
            final ClusterWheatleyExtendSW targetCluster = new ClusterWheatleyExtendSW(nPoints, fTargetf1, fTargete2);
            targetCluster.setTemperature(1.0);
            
            double vCore = space.sphereVolume(1);
            double vIntegral = vCore + (vhs - vCore) * Y;
            double chainIntegral = (SpecialFunctions.factorial(nPoints) / 2) * Math.pow(vIntegral, nPoints - 1);
            ClusterChainHS crc = new ClusterChainHS(nPoints, fRefPos, chainFrac/chainIntegral, 0);
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            ClusterAbstract refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{crc, ct});
            long numTreeDiagrams = 1;
            for (int i=0; i<nPoints-2; i++) {
                numTreeDiagrams *= nPoints;
            }

            double treeIntegral = numTreeDiagrams * Math.pow(vIntegral, nPoints - 1);

            // weighting for chain and ring are handled internally
            ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{1,(1-chainFrac)/treeIntegral});
            refCluster.setTemperature(1.0);
            
            ClusterAbstract[] targetDiagrams = new ClusterAbstract[]{targetCluster};


            System.out.println("thread "+iThread+", "+steps+" steps");

            ISpecies species = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("A")));
            final SimulationVirial sim = new SimulationVirial(space, new ISpecies[]{species}, new int[]{nPoints}, 1.0, ClusterWeightAbs.makeWeightCluster(refCluster),refCluster, targetDiagrams);
            sim.setDoWiggle(false);
            sim.setSeeds(mySeeds);
            sim.init();
            int[] randSeed = ((RandomMersenneTwister)sim.getRandom()).getSeedArray();
            if (randSeed == null) {
                System.out.println(iThread+" Random seed: "+((RandomMersenneTwister)sim.getRandom()).getSeed());
            }
            else {
                System.out.println(iThread+" Random seeds: "+Arrays.toString(randSeed));
            }
            sim.setMeter(null);
            PropertyBin pod0 = new PropertyBin() {
                final IntSet pv = new IntSet(new int[1]);
                public IntSet value() {
                    return pv;
                }
            };
            PropertyBin pod = new PropertyBin() {
                final IntSet pv = new IntSet(new int[2]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getCoreEdgeCount();
                    pv.v[1] = targetCluster.getWellEdgeCount();
                    return pv;
                }
            };
            
            PropertyBin podOD4 = new PropertyBin() {
                final IntSet pv = new IntSet(new int[4]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getCoreEdgeCount();
                    pv.v[1] = targetCluster.getWellEdgeCount();
                    pv.v[2] = pv.v[3] = 0;
                    int[] odc = targetCluster.getOutDegreeCore();
                    int[] odw = targetCluster.getOutDegreeWell();
                    for (int i=0; i<odc.length; i++) {
                        pv.v[2] += odc[i]*odc[i];
                        pv.v[3] += odw[i]*odw[i];
                    }
                    
                    return pv;
                }
            };
            PropertyBin podOD5 = new PropertyBin() {
                final IntSet pv = new IntSet(new int[5]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getCoreEdgeCount();
                    pv.v[1] = targetCluster.getWellEdgeCount();
                    pv.v[2] = pv.v[3] = pv.v[4] = 0;
                    int[] odc = targetCluster.getOutDegreeCore();
                    int[] odw = targetCluster.getOutDegreeWell();
                    for (int i=0; i<odc.length; i++) {
                        pv.v[2] += odc[i]*odc[i];
                        pv.v[3] += odw[i]*odw[i];
                        pv.v[4] += odc[i]*odw[i];
                    }
                    
                    return pv;
                }
            };
            PropertyBin podOD8 = new PropertyBin() {
                final IntSet pv = new IntSet(new int[8]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getCoreEdgeCount();
                    pv.v[1] = targetCluster.getWellEdgeCount();
                    pv.v[2] = pv.v[3] = pv.v[4] = pv.v[5] = pv.v[6] = pv.v[7] = 0;
                    int[] odc = targetCluster.getOutDegreeCore();
                    int[] odw = targetCluster.getOutDegreeWell();
                    for (int i=0; i<nPoints; i++) {
                        pv.v[2] += odc[i]*odc[i];
                        pv.v[3] += odw[i]*odw[i];
                        pv.v[4] += odc[i]*odw[i];
                        int not = (nPoints-1 - odc[i] - odw[i]);
                        pv.v[5] += not*not;
                        pv.v[6] += odc[i]*not;
                        pv.v[7] += odw[i]*not;
                    }
                    
                    return pv;
                }
            };
            PropertyBin podODCliq = new PropertyBin() {
                final IntSet pv = new IntSet(new int[9]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getCoreEdgeCount();
                    pv.v[1] = targetCluster.getWellEdgeCount();
                    pv.v[2] = pv.v[3] = pv.v[4] = 0;
                    int[] odc = targetCluster.getOutDegreeCore();
                    int[] odw = targetCluster.getOutDegreeWell();
                    for (int i=0; i<nPoints; i++) {
                        pv.v[2] += odc[i]*odc[i];
                        pv.v[3] += odw[i]*odw[i];
                        pv.v[4] += odc[i]*odw[i];
                    }
                    pv.v[5] = targetCluster.getF1CliqueCount();
                    pv.v[6] = targetCluster.getE2CliqueCount();
                    pv.v[7] = targetCluster.getEFCliqueCount();
//                    pv.v[8] = targetCluster.getNoneCliqueCount();
                    
                    return pv;
                }
            };
            PropertyBin podODCliq2 = new PropertyBin() {
                final IntSet pv = new IntSet(new int[9]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getCoreEdgeCount();
                    pv.v[1] = targetCluster.getWellEdgeCount();
                    pv.v[2] = pv.v[3] = pv.v[4] = 0;
                    int[] odc = targetCluster.getOutDegreeCore();
                    int[] odw = targetCluster.getOutDegreeWell();
                    for (int i=0; i<nPoints; i++) {
                        pv.v[2] += odc[i]*odc[i];
                        pv.v[3] += odw[i]*odw[i];
                        pv.v[4] += odc[i]*odw[i];
                    }
                    pv.v[5] = targetCluster.getF1CliqueCount();
                    pv.v[6] = targetCluster.getE2CliqueCount();
                    pv.v[7] = targetCluster.getEFCliqueCount();
                    pv.v[8] = targetCluster.getNoneCliqueCount();
                    
                    return pv;
                }
            };
            final DooDad dooDad = new DooDad(nPoints);
            PropertyBin podODCliqDoodad = new PropertyBin() {
                final IntSet pv = new IntSet(new int[13]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getCoreEdgeCount();
                    pv.v[1] = targetCluster.getWellEdgeCount();
                    pv.v[2] = pv.v[3] = pv.v[4] = 0;
                    int[] odc = targetCluster.getOutDegreeCore();
                    int[] odw = targetCluster.getOutDegreeWell();
                    for (int i=0; i<nPoints; i++) {
                        pv.v[2] += odc[i]*odc[i];
                        pv.v[3] += odw[i]*odw[i];
                        pv.v[4] += odc[i]*odw[i];
                    }
                    pv.v[5] = targetCluster.getF1CliqueCount();
                    pv.v[6] = targetCluster.getE2CliqueCount();
                    pv.v[7] = targetCluster.getEFCliqueCount();
                    pv.v[8] = targetCluster.getNoneCliqueCount();
                    pv.v[9] = dooDad.value(pv.v[5], targetCluster.getF1Cliques());
                    pv.v[10] = dooDad.value(pv.v[6], targetCluster.getE2Cliques());
                    pv.v[11] = dooDad.value(pv.v[7], targetCluster.getEFCliques());
                    pv.v[12] = dooDad.value(pv.v[8], targetCluster.getNoneCliques());
                    
                    return pv;
                }
            };
            PropertyBin[] myPODs= new PropertyBin[11];
            myPODs[2] = pod;
            myPODs[3] = pod;
            myPODs[4] = pod; //OD5;
            myPODs[5] = podOD5; //ODCliq;
            myPODs[6] = podODCliqDoodad;
            myPODs[7] = podODCliqDoodad;
            // podODCliqDoodad yields too many bins for n>=8.  even if our memory could hold it, it would take
            // a long time to get any benefit from binning
            myPODs[8] = podODCliq2;
            myPODs[9] = podODCliqDoodad;
            myPODs[10] = podOD5;
            meter = new MeterVirialEBinMultiThreaded(targetCluster, sim.getRandom(), myPODs[nPoints], totalCount, allMyData, iThread, doReweight, nPoints);
            meter.setDoCov(doCov);
            meter.setBox(sim.box);
            if (w>=0) {
                meter.setWeight(w);
            }
            else if (allMyData == null) {
                meter.readWeights(runName+nPoints+"_weights.dat", nPoints);
            }
//            double tRatio = nPoints == 4 ? 30.0 : 0.154*Math.exp(1.0683*nPoints);

            // based on fit of data using weight=1, weight=0, each for 10^9 steps.
//            double tRatio = 0.44869 * Math.exp(0.64714 * nPoints);
            //double tRatio = 0.292077*nPoints*nPoints + 0.00378375*Math.pow(3, nPoints);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meter));
            
            sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);

            RandomPositionLatticeSQW positionSource = new RandomPositionLatticeSQW(space, sim.getRandom(), new double[]{1,Y}, L);

            MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), sim.box, 1);
            mcMoveHSC.setPositionSource(positionSource);
            sim.integrator.getMoveManager().addMCMove(mcMoveHSC);
            sim.integrator.getMoveManager().setFrequency(mcMoveHSC, chainFrac);
            MCMoveClusterAtomHSTree mcMoveHST = new MCMoveClusterAtomHSTree(sim.getRandom(), sim.box, 1);
            mcMoveHST.setPositionSource(positionSource);
            sim.integrator.getMoveManager().addMCMove(mcMoveHST);
            sim.integrator.getMoveManager().setFrequency(mcMoveHST, 1-chainFrac);
            MeterVirialEBinMultiThreaded.setTRatio(tRatio);

            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
            long t2 = System.currentTimeMillis();
            System.out.println("thread "+iThread+" time: "+(t2-t1)*0.001);
        }
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHSBinParam extends ParameterBase {
        public int nPoints = 6;
        public long numSteps = 100000000;
        public double chainFrac = 0.5;
        public String runName = "sqwlat";
        public int nThreads = 1;
        public double w = -1;
        public int[] randomSeeds = new int[0];
        public boolean shareData = true;
        public boolean doCov = false;
        public double targetTemp = 0;
        public int L = 0;
        public int D = 2;
    }
    
}
