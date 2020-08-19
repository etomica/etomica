/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.*;
import etomica.virial.IntSet.PropertyBin;
import etomica.virial.MeterVirialBDBinMultiThreaded.MyData;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Calculation for virial coefficients of hard spheres
 */
public class VirialHSBinMultiThreaded {


    public static void main(String[] args) {

        VirialHSBinParam params = new VirialHSBinParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.nPoints = 3;
            params.numSteps = 100000L;
            params.ref = VirialHSBinParam.TREE_CHAIN_RING;
            params.chainFrac = 0;
            params.ringFrac = 1;
            params.nPtsTabulated = 0;
            params.nThreads = 1;
//            params.oldTRatio = true;
//            params.doWheatley = false;
//            params.w = 1;
            /*params.randomSeeds = new int[]{735039559, -1020604306, 1247967731, 1429728360,
                                           50165459, -1444698221, -2093020123, -1788518319,
                                           -1969314700, 130773570, 790544604, -1245785173,
                                           765517301, -139633052, -218420304, 1328715880};*/
//                                           1637746417, 657165813, 1400225750, 58695843,
//                                           -845441424, -835797795, 785941010, -1362915253,
//                                           1209810111, -1814202639, -431560078, 2084780809,
//                                           46116748, 1060935712, -478171462, 559810669};
        }
        
//        double[] cf = new double[]{0,0,0,0,0,0.2,0.2,0.3,0.5,0.6,0.6,0.7};
//        double[] rf = new double[]{0,0,0,0,0,0.8,0.7,0.6,0.4,0.3,0.3,0.2};
//        for (int nn=5; nn<=11; nn++) {
//        params.nPoints = nn;
//        params.chainFrac = cf[nn];
//        params.ringFrac = rf[nn];

        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final int ref = params.ref;
        final double chainFrac = params.chainFrac;
        final double ringFrac = params.ringFrac;
        final double sigmaHS = 1.0;
        final int nPtsTabulated = params.nPtsTabulated;
        final int nThreads = params.nThreads;
        final boolean doWheatley = params.doWheatley;
        final double w = doWheatley ? params.w : 1;
        final int[] allRandomSeeds = params.randomSeeds;

        double litHSB = Double.NaN;
        try {
            litHSB = Standard.BHS(nPoints, sigmaHS);
        }
        catch (RuntimeException e) {}
        final double vhs = (4.0/3.0)*Math.PI*sigmaHS*sigmaHS*sigmaHS;

        System.out.println("HS singly-connected sampling B"+nPoints);
		
        final Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHS);
        MayerFunction fRefPos = new MayerFunction() {
            
            public void setBox(Box box) {
            }
            
            public IPotential getPotential() {
                return null;
            }
            
            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHS*sigmaHS ? 1 : 0;
            }
        };

        Map<IntSet,MeterVirialBDBinMultiThreaded.MyData> allMyData = new HashMap<IntSet,MeterVirialBDBinMultiThreaded.MyData>();
        boolean doReweight = doWheatley;
        if (w < 0) {
            // meter0 exists only to allow us to read in the weights file (if it exists) and to populate allMyData
            // this needs to be done (finished) before any thread gets started.
            MeterVirialBDBinMultiThreaded meter0 = new MeterVirialBDBinMultiThreaded(null, null, null, new long[1], allMyData, 0, true);
            meter0.readWeights(params.runName+nPoints+"_weights.dat");
            doReweight = allMyData.size() == 0;
        }
        if (!doReweight && !params.shareData) {
            // reweighting not needed, we don't want to share data
            allMyData = null;
        }
        else if (doReweight && !params.shareData) {
            // reweighting is needed, we have to share data
            params.shareData = true;
        }
        if (ref == VirialHSBinParam.TREE) {
            System.out.println("using a tree reference");
        }
        else if (ref == VirialHSBinParam.TREE_CHAIN) {
            System.out.println("using a chain/tree reference ("+chainFrac+" chains)");
        }
        else if (ref == VirialHSBinParam.TREE_CHAIN_RING) {
            System.out.println("using a ring/chain/tree reference ("+ringFrac+" rings, "+chainFrac+" chains)");
        }
        // (nPoints-1)! is simply not included by ClusterWheatley, so do that here.
        double refIntegral = 1;

        double tRatio = 0;
        double t0Tree = -3.18286 + nPoints*(1.66621 + nPoints*(-0.199946 + nPoints*0.0113384));
        double t0Chain = -1.10143 + nPoints*(0.707278 + nPoints*(-0.0830664 + nPoints*0.00582492));
        double t0Ring = -3.10337 + nPoints*(1.81457 + nPoints*(-0.243239 + nPoints*0.0191835));
        t0Tree = Math.exp(-8.727299006700907e-01 + 2.437041849455276e-01*nPoints);
        t0Chain = Math.exp(-1.076280422298466e+00 + 2.319522918620784e-01*nPoints);
        t0Ring = Math.exp(-6.048590247634202e-01 + 2.898441951777069e-01*nPoints);
        double tc = Math.exp(8.644763200709733e-01 + 3.227819802929632e-03*nPoints +  5.355934080103668e-02*nPoints*nPoints);
        double fUnscreenedTree = 0.00446301 + 224.919*Math.pow(nPoints, -5.29);
        double fUnscreenedChain = 0.00521601 + 90.2381*Math.pow(nPoints, -4.73159);
        double fUnscreenedRing = 0.0409468 + 1530.97*Math.pow(nPoints, -5.5381);
        double tRatioTree = 5.858 + 0.0190644*Math.pow(2.65435, nPoints);
        double tRatioChain = 7.22624 + 0.0431276*Math.pow(2.53604, nPoints);
        double tRatioRing = 3.95055 + 0.0249246*Math.pow(2.39758, nPoints);
        double t0 = 0;
        if (nPtsTabulated > 0) {
            t0Chain = 0.625204 + nPoints*(-0.0439293 + nPoints*(0.0194865 + nPoints*0.00137849));
            t0Ring = 3.53511 + nPoints*(-1.22024 + nPoints*(0.191016 + nPoints*(-0.000160057)));
            t0Tree = -0.92244 + nPoints*(0.611055 + nPoints*(-0.0532705 + nPoints*0.00494932));
            tRatioChain = 7.50451 + 0.00972982*Math.pow(2.61644, nPoints);
            tRatioRing = 3.84716 + 0.00467202*Math.pow(2.53048, nPoints);
            tRatioTree = 6.21487 + 0.0125424*Math.pow(2.47438, nPoints);
        }
        if  (ref == VirialHSBinParam.TREE) {
            tRatio = tRatioTree;
        }
        else if (ref == VirialHSBinParam.TREE_CHAIN) {
            double t1Tree = t0Tree*(1-chainFrac);
            double t2Tree = tRatioTree*fUnscreenedTree*t1Tree;
            double t1Chain = t0Chain*chainFrac;
            double t2Chain = tRatioChain*fUnscreenedChain*t1Chain;
            t0 = t1Tree + t1Chain;
            double fUnscreened = fUnscreenedTree*(1-chainFrac) + fUnscreenedChain*chainFrac;
            tRatio = (t2Tree + t2Chain)/(fUnscreened*(t1Tree + t1Chain));
        }
        else if (ref == VirialHSBinParam.TREE_CHAIN_RING) {
            double t1Tree = t0Tree*(1-chainFrac-ringFrac);  // total time sampling
            double t2Tree = tRatioTree*fUnscreenedTree*t1Tree; // total time computing
            double t1Chain = t0Chain*chainFrac;
            double t2Chain = tRatioChain*fUnscreenedChain*t1Chain;
            double t1Ring = t0Ring*ringFrac;
            double t2Ring = tRatioRing*fUnscreenedRing*t1Ring;
            t0 = t1Tree + t1Chain + t1Ring;
            double fUnscreened = fUnscreenedTree*(1-chainFrac-ringFrac) + fUnscreenedChain*chainFrac + fUnscreenedRing*ringFrac; // overall fraction unscreened
//            System.out.println(t0Tree*tRatioTree+" "+t0Chain*tRatioChain+" "+t0Ring*tRatioRing);
//            System.out.println((1-chainFrac-ringFrac)*t0Tree*tRatioTree + chainFrac*t0Chain*tRatioChain + ringFrac*t0Ring*tRatioRing);
            tRatio = (t2Tree + t2Chain + t2Ring)/(fUnscreened*(t1Tree + t1Chain + t1Ring));
            tc = tRatio*t0;
//            System.out.println("old tRatio: "+tRatio);
//            System.out.println((t2Tree+t2Chain+t2Ring)/((1-chainFrac-ringFrac)*t0Tree*tRatioTree + chainFrac*t0Chain*tRatioChain + ringFrac*t0Ring*tRatioRing)+" "+fUnscreened);
//            System.out.println(((1-chainFrac-ringFrac)*t0Tree*tRatioTree + chainFrac*t0Chain*tRatioChain + ringFrac*t0Ring*tRatioRing)/(t1Tree + t1Chain + t1Ring));
//            System.out.println("new tRatio: "+tc/t0);
            tRatio = tc/t0;
        }
        System.out.println("tRatio: "+tRatio+"   tc: "+tc+"  ts: "+t0);

//        MeterVirialBDBinMultiThreadedOld.setTRatio(tRatio);


        long t1 = System.currentTimeMillis();
        SimulationWorker[] sw = new SimulationWorker[nThreads];
        long[] totalCount = new long[nThreads];
        for (int it=0; it<nThreads; it++) {
            int[] mySeeds = allRandomSeeds.length > 0 ? new int[allRandomSeeds.length/nThreads] : null;
            if (mySeeds != null) {
                System.arraycopy(allRandomSeeds, allRandomSeeds.length/nThreads*it, mySeeds, 0, mySeeds.length);
            }
            sw[it] = new SimulationWorker(it, nPtsTabulated, nPoints, fRef, fRefPos, ref, vhs, chainFrac, ringFrac, steps, space, params.runName, tRatio, allMyData, w, totalCount, doReweight, mySeeds, doWheatley);
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
        MeterVirialBDBinMultiThreaded.writeData(filename, allMyData, nThreads*steps);
        
        if (doReweight || !doWheatley) {
            if (!Double.isNaN(litHSB)) System.out.println("lit value "+litHSB);
            System.out.println();
            
            MeterVirialBDBinMultiThreaded.recomputeWeights(allMyData, nThreads*steps);

            List<IntSet> pvs = new ArrayList<IntSet>();
            pvs.addAll(allMyData.keySet());
            Collections.sort(pvs);
            double sum = 0, dsum = 0;
            double sumErrNum = 0;
            double sumErrStdev = 0, dsumErrStdev = 0;
            long totalSampleCount = 0;
            long totalNotScreenedCount = 0;
            int nSets = 0;
            /*FileWriter fw = null;
            try {
                fw = new FileWriter(params.runName+nPoints+"_run"+rn+".dat");
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }*/
            for (IntSet pv : pvs) {
                MyData amd = allMyData.get(pv);
                long c = amd.unscreenedCount;
    
                nSets++;
                totalNotScreenedCount += c;
                long sc = amd.sampleCount;
                if (sc == 0) {
                    /*try {
                        fw.write(Arrays.toString(pv.v)+" 0/"+c+"  average: 0  error: 0\n");
                    }
                    catch (IOException e) {
                        throw new RuntimeException(e);
                    }*/
                    continue;
                }
    
                totalSampleCount += sc;
    
                double avg = amd.getAvgDouble();
                double davg = amd.getDoubleAvg();
                double var = amd.getVarDouble();
                double dvar = amd.getDoubleVar();
                /*double err = Math.sqrt(var/sc);
                double derr = Math.sqrt(dvar/sc);
                try {
                    fw.write(Arrays.toString(pv.v)+" "+sc+"/"+c+"  average: "+avg+"    error: "+err+"\n");
                }
                catch (IOException e) {
                    throw new RuntimeException(e);
                }*/
                sum += c*avg;
                dsum += c*davg;
                sumErrStdev += var/sc*c*c;
                dsumErrStdev += dvar/sc*c*c;
                sumErrNum += c*((double)(nThreads*steps - c))/(nThreads*steps)*avg*avg;
            }
            /*try {
                fw.close();
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }*/
            System.out.println(nSets+" sets");
            sum *= refIntegral/(nThreads*steps);
            dsum *= refIntegral/(nThreads*steps);
            double finalErr = Math.sqrt(sumErrStdev + sumErrNum)*Math.abs(refIntegral)/(nThreads*steps);
            double dfinalErr = Math.sqrt(dsumErrStdev + sumErrNum)*Math.abs(refIntegral)/(nThreads*steps);
    
            System.out.println();
    
            System.out.println("abs average: "+sum+"  error: "+finalErr);
            System.out.println("dabs average: "+dsum+"  error: "+dfinalErr);
            if (sumErrNum > 0) System.out.println("number variance fraction: "+sumErrNum/(sumErrStdev + sumErrNum));
            System.out.println("number time fraction: "+(nThreads*steps)/(nThreads*steps + totalSampleCount*tRatio));
            System.out.println("fraction not screened: "+((double)totalNotScreenedCount)/(nThreads*steps));
            System.out.println("fraction measured: "+((double)totalSampleCount)/totalNotScreenedCount);
    
            System.out.println("Difficulty: "+(finalErr*Math.sqrt(t2-t1)));
        }
        System.out.println("time: "+(t2-t1)/1000.0);
//        }
    }

    public static class SimulationWorker extends Thread {
        
        protected final int nPtsTabulated;
        protected final int nPoints;
        protected final MayerFunction fRef;
        protected final int ref;
        protected final MayerFunction fRefPos;
        protected final double vhs;
        protected final double chainFrac, ringFrac;
        protected final long steps;
        protected final Space space;
        protected final String runName;
        protected final double tRatio;
        protected final Map<IntSet,MeterVirialBDBinMultiThreaded.MyData> allMyData;
        protected final int iThread;
        protected final double w;
        protected final long[] totalCount;
        protected final boolean doReweight;
        protected final int[] mySeeds;
        public MeterVirialBDBinMultiThreaded meter;
        protected final boolean doWheatley;
        
        public SimulationWorker(int iThread, int nPtsTabulated, int nPoints, MayerFunction fRef,
                                MayerFunction fRefPos, int ref, double vhs, double chainFrac, double ringFrac,
                                long steps, Space space, String runName, double tRatio,
                                Map<IntSet,MeterVirialBDBinMultiThreaded.MyData> allMyData, double w, long[] totalCount,
                                boolean doReweight, int[] mySeeds, boolean doWheatley) {
            this.iThread = iThread;
            this.nPtsTabulated = nPtsTabulated;
            this.nPoints = nPoints;
            this.fRef = fRef;
            this.fRefPos = fRefPos;
            this.ref = ref;
            this.vhs = vhs;
            this.chainFrac = chainFrac;
            this.ringFrac = ringFrac;
            this.steps = steps;
            this.space = space;
            this.runName = runName;
            this.tRatio = tRatio;
            this.allMyData = allMyData;
            this.w = w;
            this.totalCount = totalCount;
            this.doReweight = doReweight;
            this.mySeeds = mySeeds;
            this.doWheatley = doWheatley;
        }
        
        public void run() {
            long t1 = System.currentTimeMillis();
            final ClusterWheatley targetCluster;
            if (doWheatley) {
                targetCluster = nPtsTabulated > 0 ? new ClusterWheatleyPartitionScreening(nPoints, fRef, nPtsTabulated) : new ClusterWheatleyHS(nPoints, fRef);
            }
            else {
                VirialDiagrams v = new VirialDiagrams(nPoints, false, false);
                v.setDoReeHoover(true);
                v.setDoShortcut(true);
                targetCluster = v.makeVirialClusterHS(v.getMSMCGraphs(true,false), fRef);
            }
            targetCluster.setTemperature(1.0);
            
            ClusterAbstract refCluster = null;
            if (ref == VirialHSBinParam.TREE) {
                refCluster = new ClusterSinglyConnected(nPoints, fRefPos);
            }
            else if (ref == VirialHSBinParam.TREE_CHAIN) {
                ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);
                ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
                refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc, ct});
                long numTreeDiagrams = 1;
                for (int i=0; i<nPoints-2; i++) {
                    numTreeDiagrams *= nPoints;
                }
                ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{chainFrac/(SpecialFunctions.factorial(nPoints)/2),(1-chainFrac)/numTreeDiagrams});
            }
            else if (ref == VirialHSBinParam.TREE_CHAIN_RING) {
                ClusterChainHS cr = new ClusterChainHS(nPoints, fRefPos, true);
                long numRingDiagrams = cr.numDiagrams();
                double ringIntegral = numRingDiagrams*Standard.ringHS(nPoints);
                double chainIntegral = (SpecialFunctions.factorial(nPoints)/2)*Math.pow(vhs, nPoints-1);
                ClusterChainHS crc = new ClusterChainHS(nPoints, fRefPos, chainFrac/chainIntegral, ringFrac/ringIntegral);
                ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
                refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{crc, ct});
                long numTreeDiagrams = 1;
                for (int i=0; i<nPoints-2; i++) {
                    numTreeDiagrams *= nPoints;
                }

                double treeIntegral = numTreeDiagrams*Math.pow(vhs, nPoints-1);

                // weighting for chain and ring are handled internally
                ((ClusterWeightUmbrella)refCluster).setWeightCoefficients(new double[]{1,(1-ringFrac-chainFrac)/treeIntegral});
            }
            refCluster.setTemperature(1.0);
            
            ClusterAbstract[] targetDiagrams = new ClusterAbstract[]{targetCluster};


            System.out.println("thread "+iThread+", "+steps+" steps");
            
            final SimulationVirial sim = new SimulationVirial(space,new SpeciesSpheresMono(space, new ElementSimple("A")), 1.0,ClusterWeightAbs.makeWeightCluster(refCluster),refCluster, targetDiagrams, false, mySeeds);
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
                final IntSet pv = new IntSet(new int[1]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getEdgeCount();
                    return pv;
                }
            };
            PropertyBin pclique = new PropertyBin() {
                final IntSet pv = new IntSet(new int[2]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getEdgeCount();
                    pv.v[1] = targetCluster.getCliqueCount();
                    return pv;
                }
            };
            PropertyBin pefclique = new PropertyBin() {
                final IntSet pv = new IntSet(new int[3]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getEdgeCount();
                    pv.v[1] = targetCluster.getCliqueCount();
                    pv.v[2] = targetCluster.getECliqueCount();
                    return pv;
                }
            };
            final DooDad dooDad = new DooDad(nPoints);
            PropertyBin pefbclique = new PropertyBin() {
                final IntSet pv = new IntSet(new int[4]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getEdgeCount();
                    pv.v[1] = targetCluster.getCliqueCount();
                    pv.v[2] = targetCluster.getECliqueCount();
                    int nb = 0;
                    byte[] od = targetCluster.getOutDegree();
                    
                    for (int i=0; i<nPoints; i++) {
                        nb += od[i]*od[i];
                    }
                    pv.v[3] = nb;
                    return pv;
                }
            };

            PropertyBin pefclique2 = new PropertyBin() {

                final IntSet pv = new IntSet(new int[4]);
                public IntSet value() {
                    pv.v[0] = targetCluster.getEdgeCount();
                    pv.v[1] = targetCluster.getCliqueCount();
                    pv.v[2] = targetCluster.getECliqueCount();
                    pv.v[3] = pv.v[0] < nPoints*(nPoints-1)/2 ? dooDad.value(pv.v[1], targetCluster.getCliques()) : 0;
                    
                    return pv;
                }
            };
            PropertyBin pefcliqueEF = new PropertyBin() {
                final IntSet pv = new IntSet(new int[5]);
                public synchronized IntSet value() {
                    pv.v[0] = targetCluster.getEdgeCount();
                    pv.v[1] = targetCluster.getCliqueCount();
                    pv.v[2] = targetCluster.getECliqueCount();
                    pv.v[3] = pv.v[0] < nPoints*(nPoints-1)/2 ? dooDad.value(pv.v[1], targetCluster.getCliques()) : 0;
                    pv.v[4] = dooDad.value(pv.v[2], targetCluster.getECliques());
                    
                    return pv;
                }
            };
            meter = new MeterVirialBDBinMultiThreaded(targetCluster, sim.getRandom(), doWheatley ? (nPoints<6 ? pod : (nPoints<12 ? pefcliqueEF : pefclique2)) : pod0, totalCount, allMyData, iThread, doReweight);
            meter.setBox(sim.box);
            if (w>=0) {
                meter.setWeight(w);
            }
            else if (allMyData == null) {
                meter.readWeights(runName+nPoints+"_weights.dat");
            }
//            double tRatio = nPoints == 4 ? 30.0 : 0.154*Math.exp(1.0683*nPoints);

            // based on fit of data using weight=1, weight=0, each for 10^9 steps.
//            double tRatio = 0.44869 * Math.exp(0.64714 * nPoints);
            //double tRatio = 0.292077*nPoints*nPoints + 0.00378375*Math.pow(3, nPoints);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meter));
            
            sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);

            if  (ref == VirialHSBinParam.TREE) {
                MCMoveClusterAtomHSTree mcMoveHS = new MCMoveClusterAtomHSTree(sim.getRandom(), space, 1);
                sim.integrator.getMoveManager().addMCMove(mcMoveHS);
            }
            else if (ref == VirialHSBinParam.TREE_CHAIN) {
                MCMoveClusterAtomHSTree mcMoveHST = new MCMoveClusterAtomHSTree(sim.getRandom(), space, 1);
                sim.integrator.getMoveManager().addMCMove(mcMoveHST);
                sim.integrator.getMoveManager().setFrequency(mcMoveHST, 1-chainFrac);
                MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), space, 1);
                sim.integrator.getMoveManager().addMCMove(mcMoveHSC);
                sim.integrator.getMoveManager().setFrequency(mcMoveHSC, chainFrac);
            }
            else if (ref == VirialHSBinParam.TREE_CHAIN_RING) {
                MCMoveClusterAtomHSRing mcMoveHSR = new MCMoveClusterAtomHSRing(sim.getRandom(), space, 1);
                sim.integrator.getMoveManager().addMCMove(mcMoveHSR);
                sim.integrator.getMoveManager().setFrequency(mcMoveHSR, ringFrac);
                MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), space, 1);
                sim.integrator.getMoveManager().addMCMove(mcMoveHSC);
                sim.integrator.getMoveManager().setFrequency(mcMoveHSC, chainFrac);
                MCMoveClusterAtomHSTree mcMoveHST = new MCMoveClusterAtomHSTree(sim.getRandom(), space, 1);
                sim.integrator.getMoveManager().addMCMove(mcMoveHST);
                sim.integrator.getMoveManager().setFrequency(mcMoveHST, 1-ringFrac-chainFrac);
            }
            MeterVirialBDBinMultiThreaded.setTRatio(tRatio);

            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
            long t2 = System.currentTimeMillis();
            System.out.println("thread "+iThread+" time: "+(t2-t1)*0.001);
        }
    }
    
    public static class DooDad {
        boolean[] excludedCliques;
        public DooDad(int nPoints) {
            excludedCliques = new boolean[1<<nPoints];
        }
        public int value(int cliqueCount, int[] cliqueList) {
            int ncbonds = 0;
            for (int i=0; i<cliqueCount; i++) {
                excludedCliques[i] = false;
            }
            
            for (int i=cliqueCount-1; i>=2; i--) {
                if (excludedCliques[i]) continue;
                int ci = cliqueList[i];
                for (int j=i-1; j>=0; j--) {
                    if (excludedCliques[j]) continue;
                    if ((ci|cliqueList[j]) == ci) {
                        // j is a subset of i; ignore j
                        excludedCliques[j] = true;
                        continue;
                    }
                }
            }
            for (int i=cliqueCount-1; i>=1; i--) {
                if (excludedCliques[i]) continue;
                int ci = cliqueList[i];
                for (int j=i-1; j>=0; j--) {
                    if (excludedCliques[j]) continue;
                    int cj = cliqueList[j];
                    if ((ci|cj) != ci+cj) {
                        // i and j share some points in common
                        ncbonds++;
                    }
                }
            }
            return ncbonds;
        }
    }
    
    /**
     * Inner class for parameters
     */
    public static class VirialHSBinParam extends ParameterBase {
        public int nPoints = 6;
        public long numSteps = 100000000;
        public static final int TREE = 0, TREE_CHAIN = 1, TREE_CHAIN_RING = 2;
        public int ref = TREE;
        public double chainFrac = 0.4;
        public double ringFrac = 0.4;
        public String runName = "hs";
        public int nPtsTabulated = 0;
        public int nThreads = 1;
        public double w = -1;
        public int[] randomSeeds = new int[0];
//        public boolean oldTRatio = false;
        public boolean shareData = true;
        public boolean doWheatley = true;
    }
    
}
