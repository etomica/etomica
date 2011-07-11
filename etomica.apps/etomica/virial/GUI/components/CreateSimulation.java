package etomica.virial.GUI.components;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.atom.DiameterHashByType;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P22CLJQ;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.CompoundDimension;
import etomica.units.CompoundUnit;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Kelvin;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.Volume;
import etomica.util.Constants.CompassDirection;
import etomica.virial.ClusterAbstract;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactoryTangentSpheres;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlapRejected;

public class CreateSimulation {
	
	
	private String CustomName = null;
	
	private int[] MoleculeCount = new int[2];
	private SimulationGraphic simGraphic;
	private ParameterMapping[] potential = new ParameterMapping[2];
	
	private int nPoints;
	private boolean InterNonBondedPotentialFlag;
	private boolean Species1IntraNonBondedPotentialFlag;
	private boolean Species2IntraNonBondedPotentialFlag;
	private boolean Species1IntraBondedPotentialFlag;
	private boolean Species2IntraBondedPotentialFlag;
	
	private double temperature;
	private double SigmaHSRef;
	private int steps;
	
	private JFrame frame;
	
	public CreateSimulation(ParameterMapping[] Potential, int Molecule1Count, int Molecule2Count){
		this.potential[0] = Potential[0];
		if(Potential.length>1){
			this.potential[1] = Potential[1];
		}
		this.MoleculeCount[0]=Molecule1Count;
		this.MoleculeCount[1]=Molecule2Count;
		nPoints = Molecule1Count+Molecule2Count;
	}

	public void runSimulation(){
		
		System.out.println(potential[0].getClass().getName());
		
		if(potential[1] != null){
			InterNonBondedPotentialFlag = true;
			
			if(potential[1].getPotentialSites().length > 1){
				Species2IntraNonBondedPotentialFlag = true;
				Species2IntraBondedPotentialFlag = true;
			}
		}
		else{
			InterNonBondedPotentialFlag = false;
			Species2IntraNonBondedPotentialFlag = false;
			Species2IntraBondedPotentialFlag = false;
		}
		
		
		if(potential[0].getPotentialSites().length > 1){
			Species1IntraNonBondedPotentialFlag = true;
			Species1IntraBondedPotentialFlag = true;
		}
		else{
			Species1IntraNonBondedPotentialFlag = false;
			Species1IntraBondedPotentialFlag = false;
		}
		
		//All Environment variables set first
		temperature = potential[0].getDoubleDefaultParameters("TEMPERATURE");
		SigmaHSRef = potential[0].getDoubleDefaultParameters("SIGMAHSREF");
		String NoOfSteps = Double.toString(potential[0].getDoubleDefaultParameters("STEPS"));
		String[] IntSteps= NoOfSteps.split("\\.");
		steps = Integer.parseInt(IntSteps[0]);
		
		//InterNonBondedPotentialFlag
		
		/*
		//Get the simulation environment variables
		double sigmaHSRef = potential[0].getDoubleDefaultParameters("SIGMAHSREF");
		double temperature = potential[0].getDoubleDefaultParameters("TEMPERATURE");;
		String NoOfSteps = Double.toString(potential[0].getDoubleDefaultParameters("STEPS"));
		String[] IntSteps= NoOfSteps.split("\\.");
		int steps = Integer.parseInt(IntSteps[0]);
		
		
		
		
		double epsilon = Kelvin.UNIT.toSim(125.317);// for CO2
        double sigma = 3.0354; // for CO2
        double moment = 3.0255*epsilon*Math.pow(sigma,5);    // moment=Q^2/(epsilon*sigma^5), 3.0255 for CO2
        double bondL = 0.699*sigma;   // 0.699sigma for CO2;
        sigmaHSRef = sigmaHSRef*sigma;
        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        P22CLJQ pTarget = new P22CLJQ(space);
        
        pTarget.setEpsilon(epsilon);
        pTarget.setSigma(sigma);
        pTarget.setQuadrupolarMomentSquare(moment);
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        MayerEGeneral eTarget = new MayerEGeneral(pTarget);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
		
        ConformationLinear conformation = new ConformationLinear(space, bondL);
        final SimulationVirialOverlapRejected sim = new SimulationVirialOverlapRejected(space,new SpeciesFactoryTangentSpheres(2, conformation), temperature,refCluster,targetCluster);
        sim.integratorOS.setNumSubSteps(1000);
        
        if (true) {
            double size = (sigma*(1+bondL))*2;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            if(simGraphic == null){
            simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());}
            else{
            	if(frame != null){
            		frame.setVisible(false);
            	}
            	simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            }
            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0/size));
            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
            
            DiameterHashByType diameterManager = (DiameterHashByType)simGraphic.getDisplayBox(sim.box[0]).getDiameterHash();
            diameterManager.setDiameter(sim.species.getAtomType(0), sigma);
            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameterManager);
            frame = simGraphic.makeAndDisplayFrame();

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
            if (Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0) {
                throw new RuntimeException("Oops");
            }
            
            return;
        }
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }

        System.out.println("equilibration finished");
        sim.setAccumulatorBlockSize((int)steps);

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }
        
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());

        double ratio = sim.dsvo.getDataAsScalar();
        double error = sim.dsvo.getError();
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
        */
	}
	

}
