package etomica.virial.simulations;


import java.awt.Color;

import etomica.api.IAtomList;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.potential.P2HePCKLJS;
import etomica.potential.P3CPSNonAdditiveHe;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSumNonAdditiveTrimerEnergy;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;

// Adapted by Kate from VirialGCPM

public class VirialCPSHeliumNonAdditive {

    public static void main(String[] args) {

        VirialParam params = new VirialParam();
        
        double temperature; final int nPoints; double sigmaHSRef;
        long steps;
        if (args.length == 0) {
        	
        	nPoints = params.nPoints;
            temperature = params.temperature;
            steps = params.numSteps;
            sigmaHSRef = params.sigmaHSRef;
            
            // number of overlap sampling steps
            // for each overlap sampling step, the simulation boxes are allotted
            // 1000 attempts for MC moves, total
            
        } else if (args.length == 4) {
            //ReadParameters paramReader = new ReadParameters(args[0], params);
            //paramReader.readParameters();
        	nPoints = Integer.parseInt(args[0]);
        	temperature = Double.parseDouble(args[1]);
            steps = Integer.parseInt(args[2]);
            sigmaHSRef = Double.parseDouble(args[3]);
            params.writeRefPref = true;
        	
        } else {
        	throw new IllegalArgumentException("Incorrect number of arguments passed.");
        }
        

        int numSubSteps = 1000;

        final double[] HSB = new double[7];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
        System.out.println("Helium overlap sampling B"+nPoints+"NonAdd at T="+temperature+ " K");
        
        //Next line not needed because energy in Kelvin
        //temperature = Kelvin.UNIT.toSim(temperature);

        System.out.println(steps+" steps ("+steps/1000+" blocks of 1000)");
        steps /= 1000;

        Space space = Space3D.getInstance();

        
        
        P2HePCKLJS p2 = new P2HePCKLJS(space);
        P3CPSNonAdditiveHe p3NonAdd = new P3CPSNonAdditiveHe(space);
    	MayerGeneralSpherical fTarget = new MayerGeneralSpherical(p2);
    	ClusterSumNonAdditiveTrimerEnergy targetCluster = Standard.virialNonAdditiveTrimerEnergy(nPoints, fTarget, nPoints>3, false);
    	targetCluster.setNo72B2B3NonAdd(true);
    	
    	MayerHardSphere fRef1 = new MayerHardSphere(sigmaHSRef);
    	MayerHardSphere fRef2 = new MayerHardSphere(sigmaHSRef*3);
    	//MayerEHardSphere eRef1 = new MayerEHardSphere(sigmaHSRef);
    	//MayerEHardSphere eRef2 = new MayerEHardSphere(sigmaHSRef*3);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, (MayerFunction)fRef1, nPoints>3, null, false);
        //ClusterAbstract refCluster = Standard.virialCluster2(nPoints, (MayerFunction)fRef1, (MayerFunction)fRef2, nPoints>3, null, null, false);
        
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);


        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactorySpheres(), 
                temperature, refCluster,targetCluster, false);
        
        /////////////////////////////
        // Initialize non-overlapped configuration
        IAtomList atoms = sim.box[1].getLeafList();
        if (nPoints == 3) {
	        for (int i=1;i<atoms.getAtomCount();i++) {
	        	atoms.getAtom(i).getPosition().setX(0, i*sigmaHSRef);
	        }
        } else if (nPoints == 4) {
	        
	        atoms.getAtom(1).getPosition().setX(0, sigmaHSRef);
	        
	        atoms.getAtom(2).getPosition().setX(0, sigmaHSRef);
	        atoms.getAtom(2).getPosition().setX(1, sigmaHSRef);
	        
	        atoms.getAtom(3).getPosition().setX(1, sigmaHSRef);
	        
        } else if (nPoints == 5) {
        	
        	atoms.getAtom(1).getPosition().setX(0, sigmaHSRef);
        	atoms.getAtom(1).getPosition().setX(1, sigmaHSRef);
        	
	        atoms.getAtom(2).getPosition().setX(0, sigmaHSRef);
	        atoms.getAtom(2).getPosition().setX(1, -sigmaHSRef);
	        
	        atoms.getAtom(3).getPosition().setX(0, -sigmaHSRef);
	        atoms.getAtom(3).getPosition().setX(1, sigmaHSRef);
	        
	        atoms.getAtom(4).getPosition().setX(0, -sigmaHSRef);
	        atoms.getAtom(4).getPosition().setX(1, -sigmaHSRef);
	        
        } else {
        	throw new RuntimeException("Wrong number of points");
        }
        
        /*
        IAtomList atoms0 = sim.box[0].getLeafList();
        for (int i=1;i<atoms0.getAtomCount();i++) {
        	atoms0.getAtom(i).getPosition().setX(0, i*sigmaHSRef*1.3);
        } */ 
        
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            SpeciesSpheresMono species = (SpeciesSpheresMono)sim.getSpecies(0);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            ((ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
            simGraphic.makeAndDisplayFrame();
    
            sim.integratorOS.setNumSubSteps(numSubSteps);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
//            sim.getController().addAction(new IAction() {
//                public void actionPerformed() {
//                    sim.initRefPref(null, 0);
//                    sim.equilibrate(null,0);
//                    sim.ai.setMaxSteps(Long.MAX_VALUE);
//                }
//            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            return;
        }

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+params.temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }
        
        sim.setAccumulatorBlockSize((int)steps);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());
        
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

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());
        
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio average: "+ratio+", error: "+error);
        System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
        IData ratioData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index);
        IData ratioErrorData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index);
        IData averageData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index);
        IData stdevData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.STANDARD_DEVIATION.index);
        IData errorData = ((DataGroup)sim.accumulators[0].getData()).getData(AccumulatorRatioAverage.StatType.ERROR.index);
        System.out.println("reference ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("reference   average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("reference overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
        
        ratioData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.RATIO.index);
        ratioErrorData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index);
        averageData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.AVERAGE.index);
        stdevData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.STANDARD_DEVIATION.index);
        errorData = ((DataGroup)sim.accumulators[1].getData()).getData(AccumulatorRatioAverage.StatType.ERROR.index);
        System.out.println("target ratio average: "+ratioData.getValue(1)+" error: "+ratioErrorData.getValue(1));
        System.out.println("target average: "+averageData.getValue(0)
                          +" stdev: "+stdevData.getValue(0)
                          +" error: "+errorData.getValue(0));
        System.out.println("target overlap average: "+averageData.getValue(1)
                          +" stdev: "+stdevData.getValue(1)
                          +" error: "+errorData.getValue(1));
        
        System.out.println();
        System.out.println("cm"+((nPoints-1)*3)+"/mol"+(nPoints-1)+": ");
        System.out.println("abs average: "+ratio*HSB[nPoints]*Math.pow(0.60221415,nPoints-1)+", error: "+error*HSB[nPoints]*Math.pow(0.60221415,nPoints-1));
	}



    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 25.0;   // Kelvin
        public long numSteps = 100000;
        public double sigmaHSRef = 5;
        public boolean writeRefPref;
    }
}
