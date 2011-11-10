package etomica.virial.GUI.components;

import java.awt.Color;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiIntergroup;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Spherical;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.units.CompoundDimension;
import etomica.units.CompoundUnit;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
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
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap;
import etomica.virial.simulations.SimulationVirialOverlap2;

public class CreateSimulation {
	
	
	
	private SimulationEnvironment SimEnv;
	private Space space;
	
	private int[] MoleculeCount = new int[2];
	
	private SimulationGraphic simGraphic;
	private ParameterMapping[] potential = new ParameterMapping[2];
	
	private int nPoints;

	
	private double temperature;
	private double SigmaHSRef;
	private int steps;
	

	
	private MayerGeneral fTarget;
	private MayerEGeneral eTarget;
	private MayerGeneralSpherical fTarget1; 
    private MayerESpherical eTarget1;
	
	
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

	
	@SuppressWarnings("unchecked")
	public void runSimulation(SimulationEnvironment simenv, PotentialObject PObject) throws NoSuchMethodException{
		
		SimEnv = simenv;
		
		//All Environment variables set first
		temperature = SimEnv.getTemperature();
		steps = SimEnv.getNoOfSteps();
		

		//Are we having a mixture? 
		if(potential[1] != null){
			
		}
		else{
			
		}

		
		//SigmaHSRef will vary according to mixing rules
		SigmaHSRef = SimEnv.getSigmaHSRef();
		/*if(!InterNonBondedPotentialFlag){
			
		}
		//For Alkane mixtures
		
		final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(SigmaHSRef);
        HSB[3] = Standard.B3HS(SigmaHSRef);
        HSB[4] = Standard.B4HS(SigmaHSRef);
        HSB[5] = Standard.B5HS(SigmaHSRef);
        HSB[6] = Standard.B6HS(SigmaHSRef);
        HSB[7] = Standard.B7HS(SigmaHSRef);
        HSB[8] = Standard.B8HS(SigmaHSRef);
		
		 MayerHardSphere fRef = new MayerHardSphere(SigmaHSRef);
	     MayerEHardSphere eRef = new MayerEHardSphere(SigmaHSRef);
	     
	     
	     
	     space = potential[0].getSpace();
	     
	     ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
	     targetCluster.setTemperature(temperature);
	     ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
	     refCluster.setTemperature(temperature);
	     
	     System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
	     
	     final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,potential[0].createSpeciesFactory(), temperature,refCluster,targetCluster);
	     
	     sim.integratorOS.setNumSubSteps(1000);
	     int blocksize = 100;
	     sim.setAccumulatorBlockSize(blocksize);
	        
	        IAction progressReport = new IAction() {
	            public void actionPerformed() {
	                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
	                double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
	                System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
	            }
	        };
	        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
	        progressReportListener.setInterval((int)(steps/10));
	        sim.integratorOS.getEventManager().addListener(progressReportListener);

	        sim.ai.setMaxSteps(steps);
	        for (int i=0; i<2; i++) {
	            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
	        }
	        sim.getController().actionPerformed();

	        
		
	}
	public int Factorial(int n)
	{
		if (n == 0)
			return 1;
		else
			return n * Factorial(n-1);
	}
*/
	
	}

}
