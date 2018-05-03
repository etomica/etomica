/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.MoleculeActionTranslateTo;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.IData;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.histogram.HistogramSimple;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.models.water.PNWaterGCPM;
import etomica.models.water.PNWaterGCPM.Component;
import etomica.models.water.PNWaterGCPM.PNWaterGCPMCached;
import etomica.models.water.SpeciesWater4PCOM;
import etomica.potential.PotentialNonAdditiveDifference;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.CompoundUnit;
import etomica.units.Kelvin;
import etomica.units.Unit;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;
import java.util.Arrays;

/**
 * Computes virial coefficients and its temperature derivatives for GCPM Water
 * 
 * @author Andrew Schultz
 */
public class VirialH2OGCPMD {


    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize parameters here
            params.nPoints = 5;
            params.nDer = 20;
            params.temperature = 800;
            params.numSteps = 1000000;
            params.sigmaHSRef = 5;
            params.nonAdditive = Nonadditive.TOTAL;
            params.seed = null;
            params.doHist = false;
            isCommandline = true;
            params.dorefpref = false;
            params.tol = 1e-12;
        }

        final int nPoints = params.nPoints;
        final int nDer = params.nDer;
        final double temperatureK = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        int[] seed = params.seed;
        boolean dorefpref = params.dorefpref;
        final double tol = params.tol;

        final double refFrac = params.refFrac;
        final Nonadditive nonAdditive = nPoints < 3 ? Nonadditive.NONE : params.nonAdditive;

        final double HSB = Standard.BHS(nPoints, sigmaHSRef);

        System.out.println("Overlap sampling for H2O GCPM at " + temperatureK + " K " + "for B"+nPoints+" and "+nDer+" derivatives");
        if (nonAdditive != Nonadditive.NONE) {            
            System.out.println("Including induction");
        }
        
        double temperature = Kelvin.UNIT.toSim(temperatureK);

        final int precision = -3*(int)Math.log10(tol);

        final double BDAccFrac = 0.001;
        
        System.out.println("Reference diagram: B"+nPoints+" for hard spheres with diameter " + sigmaHSRef + " Angstroms");
        
        System.out.println("  B"+nPoints+"HS: "+HSB);
		
        final Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        SpeciesWater4PCOM speciesWater = new SpeciesWater4PCOM(space);
        
        final PNWaterGCPM pTarget = new PNWaterGCPM(space);

        MayerGeneral fTarget = new MayerGeneral(pTarget);

        ClusterAbstractMultivalue targetCluster = new ClusterWheatleySoftDerivatives(nPoints, fTarget, tol,nDer);
        ClusterAbstractMultivalue targetClusterBD = new ClusterWheatleySoftDerivativesBD(nPoints, fTarget, precision,nDer);

        if (nPoints==2) {
            // pure B2 for water.  we need flipping.
            // additive B3 for water should be fine and biconnectivity will help with mixture coefficients.
            ((ClusterWheatleySoftDerivatives)targetCluster).setTolerance(0);
            ((ClusterWheatleySoftDerivatives)targetCluster).setDoCaching(false);
            ((ClusterWheatleySoftDerivativesBD)targetClusterBD).setDoCaching(false);
            targetCluster = new ClusterCoupledFlippedMultivalue(targetCluster, targetClusterBD, space, 20, nDer, tol);
        }

        if (nonAdditive == Nonadditive.FULL || nonAdditive == Nonadditive.TOTAL) {
            PNWaterGCPMCached p2 = pTarget.makeCachedPairPolarization();
            PNWaterGCPM pFull = new PNWaterGCPM(space);
            pFull.setComponent(Component.INDUCTION);
            PotentialNonAdditiveDifference pnad = new PotentialNonAdditiveDifference(space, p2, pFull);
            MayerFunctionNonAdditiveFull fnad = new MayerFunctionNonAdditiveFull(pnad);            
            targetCluster = new ClusterWheatleyMultibodyDerivatives(nPoints, fTarget,fnad, 0, nDer, nonAdditive == Nonadditive.TOTAL);
            targetClusterBD = new ClusterWheatleyMultibodyDerivativesBD(nPoints, fTarget,fnad,new MayerFunctionNonAdditive[0], precision, nDer, nonAdditive == Nonadditive.TOTAL);
            ((ClusterWheatleyMultibodyDerivatives)targetCluster).setRCut(100);
            ((ClusterWheatleyMultibodyDerivativesBD)targetClusterBD).setRCut(100);
            // water induction requires flipping
            ((ClusterWheatleyMultibodyDerivatives)targetCluster).setDoCaching(false);
            ((ClusterWheatleyMultibodyDerivativesBD)targetClusterBD).setDoCaching(false);
            targetCluster = new ClusterCoupledFlippedMultivalue(targetCluster, targetClusterBD, space, 20, nDer, tol);
        }

        ClusterMultiToSingle[] primes = new ClusterMultiToSingle[nDer];
        for(int m=0;m<primes.length;m++){
            primes[m]= new ClusterMultiToSingle(targetCluster, m+1);
        }

        targetCluster.setTemperature(temperature);

        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);
                        
        System.out.println(steps+" steps (1000 IntegratorOverlap steps of "+(steps/1000)+")");
 		
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, speciesWater, nPoints, temperature, refCluster, targetCluster);
        if(seed!=null)sim.setRandom(new RandomMersenneTwister(seed));
        if(targetCluster instanceof ClusterCoupledFlippedMultivalue) {
            ((ClusterCoupledFlippedMultivalue) targetCluster).setBDAccFrac(BDAccFrac,sim.getRandom());
        }
        sim.setExtraTargetClusters(primes);

        //No weighting for BD flipping
        if(nPoints > 4){
            double r=1000;
            double w=1;
            if(nPoints >5){
                r = 1000;
                w = 1;
            }
            ClusterWeight[] sampleclusters = sim.getSampleClusters();
            sampleclusters[1] = new Clusterfoo(targetCluster, r ,w);
            sim.setSampleClusters(sampleclusters);
        }
        
        sim.init();

        System.out.println("random seeds: "+Arrays.toString(seed==null?sim.getRandomSeeds():seed));
        System.out.println("Big Decimal Tolerance: " + tol);
        System.out.println("Big Decimal Acceptance Fraction: " + BDAccFrac);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        
        if (nonAdditive != Nonadditive.NONE) {
            MoleculeActionTranslateTo act = new MoleculeActionTranslateTo(space);
            Vector pos = space.makeVector();
            double r = 4;
            for (int i=1; i<nPoints; i++) {
                double theta = 2*i*Math.PI/nPoints;
                pos.setX(0, r*(1-Math.cos(theta)));
                pos.setX(1, r*Math.sin(theta));
                act.setDestination(pos);
                act.actionPerformed(sim.box[1].getMoleculeList().getMolecule(i));
            }
            sim.box[1].trialNotify();
            sim.box[1].acceptNotify();
        }
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{40,40,40}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{40,40,40}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            ((ColorSchemeByType)displayBox1.getColorScheme()).setColor(sim.species[0].getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)displayBox1.getColorScheme()).setColor(sim.species[0].getAtomType(1), Color.RED);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
//            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);

//            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
//            displayBox0.setColorScheme(colorScheme);
//            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
//            displayBox1.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 10);
                    sim.equilibrate(null, 20);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            return;
        }


        
        long t1 = System.currentTimeMillis();
        
        sim.integratorOS.setNumSubSteps(1000);
        
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        steps /= 1000;
        sim.setAccumulatorBlockSize(steps);
        
        System.out.println();
        String refFileName = null;
        if (isCommandline) {
            // if running interactively, don't use the file
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints+"_"+(nonAdditive==Nonadditive.NONE?"none_":(nonAdditive==Nonadditive.FULL?"full_":""))+tempString;
            refFileName += "K";
        }



        final HistogramNotSoSimple targHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple targPiHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        
        final HistogramSimple targHistr = new HistogramSimple(70, new DoubleRange(-1, 8));
        final HistogramSimple targHistBD = new HistogramSimple(70, new DoubleRange(-1, 8));
        
        int nBins = 100;
        double dx = sigmaHSRef/nBins;
        final HistogramNotSoSimple hist = new HistogramNotSoSimple(nBins, new DoubleRange(dx*0.5, sigmaHSRef+dx*0.5));
        final HistogramNotSoSimple piHist = new HistogramNotSoSimple(nBins, new DoubleRange(dx*0.5, sigmaHSRef+dx*0.5));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();
        IntegratorListener histListenerRef = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                CoordinatePairSet cPairs = sim.box[0].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }
                double v = finalTargetCluster.value(sim.box[0]);
                hist.addValue(Math.sqrt(r2Max), v);
                piHist.addValue(Math.sqrt(r2Max), Math.abs(v));
            }
            
            public void integratorInitialized(IntegratorEvent e) {
            }
        };
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}
            
            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                double r2Min = Double.POSITIVE_INFINITY;
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i=0; i<nPoints; i++) {
                    for (int j=i+1; j<nPoints; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij < r2Min) r2Min = r2ij;
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }

                double v = finalTargetCluster.value(sim.box[1]);
                double r = Math.sqrt(r2Max);
                if (r > 1) {
                    r = Math.log(r);
                }
                else {
                    r -= 1;
                }
                targHist.addValue(r, v);
                targPiHist.addValue(r, Math.abs(v));
                
                targHistr.addValue(r);
                if( Math.abs(v)<tol){
                targHistBD.addValue(r);
                }
            }

            public void integratorInitialized(IntegratorEvent e) {}
        };

        if (params.doHist) {

            final ClusterAbstractMultivalue tempcluster = targetCluster;
            long t11 = System.currentTimeMillis();

            IntegratorListener histReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*100) % sim.ai.getMaxSteps() != 0) return;
                    System.out.println("**** reference ****");
                    double[] xValues = hist.xValues();
                    double[] h = hist.getHistogram();
                    double[] piH = piHist.getHistogram();
                    for (int i=0; i<xValues.length; i++) {
                        if (!Double.isNaN(h[i])) {
                            System.out.println(xValues[i]+" "+h[i]+" "+piH[i]);
                        }
                    }
                    System.out.println("**** target ****");
                    xValues = targHist.xValues();
                    h = targHist.getHistogram();
                    piH = targPiHist.getHistogram();
                    double[] hr = targHistr.getHistogram();
                    double [] hBD = targHistBD.getHistogram();
                    for (int i=0; i<xValues.length; i++) {
                        if (!Double.isNaN(h[i])) {
                            double r = xValues[i];
                            if (r < 0) r += 1;
                            else r = Math.exp(r);
                            System.out.println(r+" "+h[i]+" "+piH[i]+" "+ hr[i]+" "+hBD[i]);
                        }
                    }

                    if(nPoints!=2&&nonAdditive==Nonadditive.NONE ){
                        System.out.println("SoftBDcount: " + ((ClusterWheatleySoftDerivatives)tempcluster).getSoftBDcount() + " SoftBDfrac: " + ((ClusterWheatleySoftDerivatives)tempcluster).getSoftBDfrac() + " Softcount: " + ((ClusterWheatleySoftDerivatives)tempcluster).getSoftcount());
                    }
                    else{
                        ClusterCoupledFlippedMultivalue foo = (ClusterCoupledFlippedMultivalue)tempcluster;
                        System.out.println("BDcount: " + foo.getBDcount() + " BDfrac: " + foo.getBDfrac() + " totBDcount: " + foo.getBDtotcount());
                        System.out.println("FlipCount: " + foo.getflipcount() + " Flipfrac: " + foo.getflipfrac() + " FlipTotcount: " + foo.gettotcount());
                        xValues= foo.histe.xValues();
                        h=foo.histe.getHistogram();
                        for (int i=0; i<xValues.length; i++) {
                            if (h[i]!=0) {
                                System.out.println(Math.exp(xValues[i]) + " " + h[i]);
                            }
                        }
                    }
                    System.out.println("time: "+(System.currentTimeMillis()-t11)/1000.0);
                }
            };
            sim.integratorOS.getEventManager().addListener(histReport);

            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[0].getEventManager().addListener(histListenerRef);
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }

        sim.initRefPref(refFileName, steps/20);
        sim.equilibrate(refFileName, steps/10);

        System.out.println("equilibration finished");

        if(dorefpref){
            long t2 = System.currentTimeMillis();
            System.out.println("time: "+(t2-t1)/1000.0);
            return;
        }

        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(1000);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }

        sim.getController().actionPerformed();
        
        if (params.doHist) {
            double[] xValues = hist.xValues();
            double[] h = hist.getHistogram();
            
            System.out.println("final ref histogram");
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
//                    System.out.println(xValues[i]+" "+(-2*h[i]+1)+" "+Math.exp(-u/temperature));
                    System.out.println(xValues[i]+" "+(-2*h[i]+1));
                }
            }
        }


        System.out.println();
        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        
        System.out.println();
        
        sim.printResults(HSB);

        //printing conversion factor for simulation temperature to K
        Unit conv = new CompoundUnit(new Unit[]{Kelvin.UNIT}, new double[]{1});
        System.out.println("K/SimT "+ conv.fromSim(1));
        
        //For printing derivatives w.r.t. beta in K
        boolean derprint = false;
        if(derprint){
            DataGroup allData = (DataGroup)sim.accumulators[1].getData();
            IData dataAvg = allData.getData(AccumulatorAverage.AVERAGE.index);
            IData dataErr = allData.getData(AccumulatorAverage.ERROR.index);
            IData dataCov = allData.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
            // we'll ignore block correlation -- whatever effects are here should be in the full target results
            
            int nTotal = (primes.length+2);
            double oVar = dataCov.getValue(nTotal*nTotal-1);
            for (int i=0; i<primes.length; i++) {
                System.out.print("derivative "+(i+1));
                // average is vi/|v| average, error is the uncertainty on that average
                // ocor is the correlation coefficient for the average and overlap values (vi/|v| and o/|v|)
                double ivar = dataCov.getValue((i+1)*nTotal+(i+1));
                double ocor = ivar*oVar == 0 ? 0 : dataCov.getValue(nTotal*(i+1)+nTotal-1)/Math.sqrt(ivar*oVar);
                
                Unit u = new CompoundUnit(new Unit[]{Kelvin.UNIT}, new double[]{i+1});
                double v = u.fromSim(dataAvg.getValue(i+1));
                double e = u.fromSim(dataErr.getValue(i+1));
                System.out.print(String.format(" average: %20.15e  error: %10.15e  ocor: %7.5f", v, e, ocor));
                if (primes.length > 0) {
                    System.out.print("  dcor:");
                    for (int j=0; j<primes.length; j++) {
                        if (i==j) continue;
                        double jvar = dataCov.getValue((j+1)*nTotal+(j+1));
                        double dcor = ivar*jvar == 0 ? 0 : dataCov.getValue((i+1)*nTotal+(j+1))/Math.sqrt(ivar*jvar);
                        System.out.print(String.format(" %8.6f", dcor));
                    }
                }
                System.out.println();
            }
        }
        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)/1000.0);

        if(nPoints!=2&&nonAdditive==Nonadditive.NONE ){
            System.out.println("SoftBDcount: " + ((ClusterWheatleySoftDerivatives)targetCluster).getSoftBDcount() + " SoftBDfrac: " + ((ClusterWheatleySoftDerivatives)targetCluster).getSoftBDfrac() + " Softcount: " + ((ClusterWheatleySoftDerivatives)targetCluster).getSoftcount());
        }
        else{
            ClusterCoupledFlippedMultivalue foo = (ClusterCoupledFlippedMultivalue)targetCluster;
            System.out.println("BDcount: " + foo.getBDcount() + " BDfrac: " + foo.getBDfrac() + " totBDcount: " + foo.getBDtotcount());
            System.out.println("FlipCount: " + foo.getflipcount() + " Flipfrac: " + foo.getflipfrac() + " FlipTotcount: " + foo.gettotcount());
        }
                
    }
    
    enum Nonadditive {
        NONE,  FULL, TOTAL
    }
    
    public static class Clusterfoo implements ClusterWeight{
        
        public final ClusterAbstract cluster;
        public final double r;
        public final double weight;
        
        public Clusterfoo(ClusterAbstract cluster, double r, double weight){
            this.r = r;
            this.cluster = cluster;
            this.weight = weight;
        }
        
        @Override
        public ClusterAbstract makeCopy() {
            return new Clusterfoo(cluster.makeCopy(),r,weight);
        }

        @Override
        public int pointCount() {
            return cluster.pointCount();
        }

        @Override
        public double value(BoxCluster box) {
            double val = Math.abs(cluster.value(box));
            int npoints = cluster.pointCount();
            double r2 = r*r;
            for(int i=0;i<npoints-1;i++){
                for(int j=i+1;j<npoints;j++){
                    double r2c = box.getCPairSet().getr2(i, j);
                    if(r2c>r2){
                        return val/weight;                   
                    }
                }
            }
            return val;            
        }

        @Override
        public void setTemperature(double temperature) {            
            cluster.setTemperature(temperature);
        }
        
        
    }
    
    
    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        // don't change these
        public int nPoints = 2;
        public int nDer = 2;
        public double temperature = 100;
        public long numSteps = 10000000;
        public double refFrac = -1;
        public double sigmaHSRef = 5;
        public boolean doHist = false;
        public Nonadditive nonAdditive = Nonadditive.NONE;
        public int[] seed = null;
        public boolean dorefpref = false;
        public double tol = 1e-12;

    }
}
