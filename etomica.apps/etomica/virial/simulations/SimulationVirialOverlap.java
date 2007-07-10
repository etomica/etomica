package etomica.virial.simulations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeLeaf;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataPump;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ConfigurationCluster;
import etomica.virial.MCMoveClusterAtomMulti;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterRotateMoleculeMulti;
import etomica.virial.MCMoveClusterWiggleMulti;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.MeterVirial;
import etomica.virial.P0Cluster;
import etomica.virial.BoxCluster;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

/**
 * @author kofke
 *
 * Simulation implementing the overlap-sampling approach to evaluating a cluster
 * diagram.
 */
public class SimulationVirialOverlap extends Simulation {


    public SimulationVirialOverlap(Space aSpace, SpeciesFactory speciesFactory, 
			double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster) {
		this(aSpace,speciesFactory,temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster.makeCopy()),ClusterWeightAbs.makeWeightCluster(targetCluster.makeCopy())});
	}
	
	public SimulationVirialOverlap(Space aSpace, SpeciesFactory speciesFactory, 
            double temperature, final ClusterAbstract[] aValueClusters, final ClusterWeight[] aSampleClusters) {
		super(aSpace,false);
		PotentialMaster potentialMaster = new PotentialMaster(space);
        sampleClusters = aSampleClusters;
        int nMolecules = sampleClusters[0].pointCount();
        species = speciesFactory.makeSpecies(this);
        getSpeciesManager().addSpecies(species);
        accumulators = new AccumulatorVirialOverlapSingleAverage[sampleClusters.length];
        accumulatorPumps = new DataPump[sampleClusters.length];
        box = new BoxCluster[sampleClusters.length];
        integrators = new IntegratorMC[sampleClusters.length];
        meters = new MeterVirial[sampleClusters.length];
        mcMoveTranslate = new MCMoveBoxStep[sampleClusters.length];
        if (species.getMoleculeType() instanceof AtomTypeGroup) {
            mcMoveRotate = new MCMoveBoxStep[sampleClusters.length];
        }
        if (species instanceof SpeciesSpheres) {
            mcMoveWiggle = new MCMoveClusterWiggleMulti[sampleClusters.length];
        }
        
        P0Cluster p0 = new P0Cluster(space);
        potentialMaster.addPotential(p0,new Species[]{});
        
        for (int iBox=0; iBox<sampleClusters.length; iBox++) {
            // integrator for iBox samples based on iBox cluster
            box[iBox] = new BoxCluster(this,sampleClusters[iBox]);
            box[iBox].setNMolecules(species, nMolecules);
            
            integrators[iBox] = new IntegratorMC(this, potentialMaster);
            integrators[iBox].setTemperature(temperature);
            integrators[iBox].setBox(box[iBox]);
            integrators[iBox].getMoveManager().setEquilibrating(false);
            
            MCMoveManager moveManager = integrators[iBox].getMoveManager();
            
            if (species.getMoleculeType() instanceof AtomTypeLeaf) {
                mcMoveTranslate[iBox] = new MCMoveClusterAtomMulti(this, potentialMaster, nMolecules-1);
                mcMoveTranslate[iBox].setStepSize(0.41);
                moveManager.addMCMove(mcMoveTranslate[iBox]);
            }
            else {
                mcMoveRotate[iBox] = new MCMoveClusterRotateMoleculeMulti(potentialMaster,getRandom(),nMolecules-1);
                mcMoveRotate[iBox].setStepSize(Math.PI);
                moveManager.addMCMove(mcMoveRotate[iBox]);
                mcMoveTranslate[iBox] = new MCMoveClusterMoleculeMulti(potentialMaster, getRandom(), 0.41, nMolecules-1);
                moveManager.addMCMove(mcMoveTranslate[iBox]);
                if (species instanceof SpeciesSpheres) {
                    if (species.getMoleculeFactory().getNumChildAtoms() > 2) {
                        mcMoveWiggle[iBox] = new MCMoveClusterWiggleMulti(this, potentialMaster, nMolecules);
                        moveManager.addMCMove(mcMoveWiggle[iBox]);
                    }
                }
            }
            
            ConfigurationCluster configuration = new ConfigurationCluster(getRandom());
            configuration.initializeCoordinates(box[iBox]);
            MeterVirial meter = new MeterVirial(new ClusterAbstract[]{aValueClusters[iBox],aSampleClusters[1-iBox]});
            setMeter(meter,iBox);
            AccumulatorVirialOverlapSingleAverage acc = new AccumulatorVirialOverlapSingleAverage(11, iBox==0);
            setAccumulator(acc,iBox);
              
        }
        
        setRefPref(1,5);
        integratorOS = new IntegratorOverlap(getRandom(), integrators);
        integratorOS.setNumSubSteps(1000);
        integratorOS.setEventInterval(1);
        ai = new ActivityIntegrate(integratorOS);
        getController().addAction(ai);
		
		dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
        integratorOS.setDSVO(dsvo);
	}

    public void setRefPref(double refPrefCenter, double span) {
        refPref = refPrefCenter;
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(refPrefCenter,span);
    }
    
    public void setMeter(MeterVirial newMeter, int iBox) {
        if (accumulators[iBox] != null) {
            // we need a new accumulator so nuke the old one now.
            if (accumulatorPumps[iBox] != null) {
                integrators[iBox].removeIntervalAction(accumulatorPumps[iBox]);
                accumulatorPumps[iBox] = null;
            }
            accumulators[iBox] = null;
        }
        meters[iBox] = newMeter;
        newMeter.setBox(box[iBox]);
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iBox) {
        accumulators[iBox] = newAccumulator;
        if (accumulatorPumps[iBox] == null) {
            accumulatorPumps[iBox] = new DataPump(meters[iBox],newAccumulator);
            integrators[iBox].addIntervalAction(accumulatorPumps[iBox]);
        }
        else {
            accumulatorPumps[iBox].setDataSink(newAccumulator);
        }
        if (integratorOS != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
            integratorOS.setDSVO(dsvo);
        }
    }

    public void setRefPref(double newRefPref) {
        System.out.println("setting ref pref (explicitly) to "+newRefPref);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
        setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
        setRefPref(newRefPref,1);
    }
    
    public void initRefPref(String fileName, long initSteps) {
        // refPref = -1 indicates we are searching for an appropriate value
        refPref = -1.0;
        if (fileName != null) {
            try { 
                FileReader fileReader = new FileReader(fileName);
                BufferedReader bufReader = new BufferedReader(fileReader);
                String refPrefString = bufReader.readLine();
                refPref = Double.parseDouble(refPrefString);
                bufReader.close();
                fileReader.close();
                System.out.println("setting ref pref (from file) to "+refPref);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
                setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
                setRefPref(refPref,1);
            }
            catch (IOException e) {
                // file not there, which is ok.
            }
        }
        
        if (refPref == -1) {
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(21,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(21,false),1);
            setRefPref(10000,15);
            ai.setMaxSteps(initSteps);
            getController().actionPerformed();

            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(11,false),1);
            setRefPref(refPref,0.2);
            for (int i=0; i<2; i++) {
                try {
                    integrators[i].reset();
                }
                catch (ConfigurationOverlapException e) {}
            }
            // set refPref back to -1 so that later on we know that we've been looking for
            // the appropriate value
            refPref = -1;
            getController().reset();
        }

    }
    
    public void equilibrate(String fileName, long initSteps) {
        // run a short simulation to get reasonable MC Move step sizes and
        // (if needed) narrow in on a reference preference
        ai.setMaxSteps(initSteps);

        for (int i=0; i<2; i++) {
            integrators[i].getMoveManager().setEquilibrating(true);
        }
        getController().actionPerformed();
        getController().reset();

        if (refPref == -1) {
            int newMinDiffLoc = dsvo.minDiffLocation();
            refPref = accumulators[0].getBennetAverage(newMinDiffLoc)
                /accumulators[1].getBennetAverage(newMinDiffLoc);
            System.out.println("setting ref pref to "+refPref+" ("+newMinDiffLoc+")");
//            int n = sim.accumulators[0].getNBennetPoints();
//            for (int i=0; i<n; i++) {
//                System.out.println(i+" "+sim.accumulators[0].getBennetBias(i)+" "+sim.accumulators[0].getBennetAverage(i)/sim.accumulators[1].getBennetAverage(i));
//            }
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,true),0);
            setAccumulator(new AccumulatorVirialOverlapSingleAverage(1,false),1);
            setRefPref(refPref,1);
            if (fileName != null) {
                try {
                    FileWriter fileWriter = new FileWriter(fileName);
                    BufferedWriter bufWriter = new BufferedWriter(fileWriter);
                    bufWriter.write(String.valueOf(refPref)+"\n");
                    bufWriter.close();
                    fileWriter.close();
                }
                catch (IOException e) {
                    throw new RuntimeException("couldn't write to refpref file");
                }
            }
        }
        else {
            dsvo.reset();
        }
    }
    
    private static final long serialVersionUID = 1L;
	protected DisplayPlot plot;
	public DataSourceVirialOverlap dsvo;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    protected DataPump[] accumulatorPumps;
	protected final ClusterWeight[] sampleClusters;
    public BoxCluster[] box;
    public Species species;
    public IntegratorMC[] integrators;
    public MCMoveBoxStep[] mcMoveRotate;
    public MCMoveBoxStep[] mcMoveTranslate;
    public MCMoveClusterWiggleMulti[] mcMoveWiggle;
    public MeterVirial[] meters;
    public ActivityIntegrate ai;
    public IntegratorOverlap integratorOS;
    public double refPref;

	public static void main(String[] args) {

		int nPoints = 5;
		double temperature = 1.3; //temperature governing sampling of configurations
//		double sigmaHSRef = 1.2*sigmaLJ1B(1.0/temperature);  //diameter of reference HS system
		double sigmaHSRef = 1.6;
        double sigmaLJ = 1.0;
		double b0 = 2*Math.PI/3. * Math.pow(sigmaHSRef,3);
		System.out.println("sigmaHSRef: "+sigmaHSRef);
		System.out.println("B2HS: "+b0);
		System.out.println("B3HS: "+(5./8.*b0*b0)+" = "+(5.0/8.0)+" B2HS^2");
		System.out.println("B4HS: "+(b0*b0*b0*(219.0*Math.sqrt(2.0)/2240.0/Math.PI-89.0/280.0+4131.0/2240.0/Math.PI*Math.atan(Math.sqrt(2.0))))+" = "
				+(219.0*Math.sqrt(2.0)/2240.0/Math.PI-89.0/280.0+4131.0/2240.0/Math.PI*Math.atan(Math.sqrt(2.0)))+" B2HS^3");
		
        Space3D space = Space3D.getInstance();
		MayerHardSphere fRef = new MayerHardSphere(space,sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);
        P2LennardJones p2LJ = new P2LennardJones(space,sigmaLJ,1.0);
		MayerGeneralSpherical fTarget = new MayerGeneralSpherical(space,p2LJ);
        MayerESpherical eTarget = new MayerESpherical(space,p2LJ);
		
        ClusterAbstract refCluster = Standard.virialCluster(nPoints,fRef,true,eRef,true,true);
        refCluster.setTemperature(temperature);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints,fTarget,true,eTarget,true,true);
        targetCluster.setTemperature(temperature);

		int maxSteps = 100000;
		
//		while (true) {
			SimulationVirialOverlap sim = new SimulationVirialOverlap(space, new SpeciesFactorySpheres(), temperature, refCluster, targetCluster);
			sim.ai.setMaxSteps(maxSteps);
//            sim.integratorOS.setEquilibrating(true);
//            sim.integratorOS.setAdjustStepFreq(true);
            sim.integratorOS.setStepFreq0(0.0001);
			sim.ai.actionPerformed();
			System.out.println("average: "+sim.dsvo.getDataAsScalar()+", error: "+sim.dsvo.getError());
            DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
            System.out.println("hard sphere ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
            System.out.println("hard sphere   average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
            System.out.println("hard sphere overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
            allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
            System.out.println("lennard jones ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
            System.out.println("lennard jones average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
            System.out.println("lennard jones overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
//		}
		
	}//end of main
}
