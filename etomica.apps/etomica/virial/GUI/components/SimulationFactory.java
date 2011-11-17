package etomica.virial.GUI.components;

import java.awt.Color;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.api.IAtomType;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.api.ISpecies;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
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
import etomica.space3d.Space3D;
import etomica.species.Species;
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
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.SpeciesFactoryTangentSpheres;
import etomica.virial.GUI.models.CreateSpeciesDM_IFactory;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialMultiOverlap;
import etomica.virial.simulations.SimulationVirialOverlap;
import etomica.virial.simulations.SimulationVirialOverlap2;

public class SimulationFactory {
	
	
	
	private SimulationEnvironmentObject SimEnv;
	private Space space;
	
	
	private SimulationGraphic simGraphic;
	private CreateSpeciesDM_IFactory potential1;
	private CreateSpeciesDM_IFactory potential2;

	private MayerGeneralSpherical fTarget1; 
    private MayerESpherical eTarget1;
	
	
	private JFrame frame;
	
	public SimulationFactory(CreateSpeciesDM_IFactory Potential1, CreateSpeciesDM_IFactory Potential2){
		this.potential1 = Potential1;
		if(Potential2 != null){
			this.potential2 = Potential2;
		}
	}

	
	@SuppressWarnings("unchecked")
	public void runSimulation(SimulationEnvironmentObject simEnv, PotentialCollections PObject) throws NoSuchMethodException{
		
		SimEnv = simEnv;
		
		final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(SimEnv.getSigmaHSRef());
        HSB[3] = Standard.B3HS(SimEnv.getSigmaHSRef());
        HSB[4] = Standard.B4HS(SimEnv.getSigmaHSRef());
        HSB[5] = Standard.B5HS(SimEnv.getSigmaHSRef());
        HSB[6] = Standard.B6HS(SimEnv.getSigmaHSRef());
        HSB[7] = Standard.B7HS(SimEnv.getSigmaHSRef());
        
        System.out.println("sigmaHSRef: "+SimEnv.getSigmaHSRef());
        System.out.println("B"+SimEnv.getnPoints()+"HS: "+HSB[SimEnv.getnPoints()]);
        
        int[] nTypes = new int[]{SimEnv.getSpeciesA(),SimEnv.getSpeciesB()}; 
        double temperature = SimEnv.getTemperature();
        Space space = Space3D.getInstance();
        
		MayerHardSphere fRef = new MayerHardSphere(SimEnv.getSigmaHSRef());
	    MayerEHardSphere eRef = new MayerEHardSphere(SimEnv.getSigmaHSRef());
	    
	    
	    ClusterAbstract refCluster = Standard.virialCluster(SimEnv.getnPoints(), fRef, SimEnv.getnPoints()>3, eRef, true);
        refCluster.setTemperature(temperature);

	    
	    if(PObject instanceof AtomicPotentialsCollection){
	    	//Set up the Potential groups!!
	    	PObject.setInterPotentialGroupII(1, new PotentialGroup(2));
	    	PObject.setInterPotentialGroupII(2, new PotentialGroup(2));
	    	PObject.setInterPotentialGroupIJ(new PotentialGroup(2));
	    	
	    
	    	
	    	MayerGeneral fTargetII = new MayerGeneral(PObject.getInterPotentialGroupII(1));  
	        MayerGeneral fTargetJJ = new MayerGeneral(PObject.getInterPotentialGroupII(2));
	        MayerGeneral fTargetIJ = new MayerGeneral(PObject.getInterPotentialGroupIJ());
	        
	        MayerEGeneral eTargetII = new MayerEGeneral(PObject.getInterPotentialGroupII(1));
	        MayerEGeneral eTargetJJ = new MayerEGeneral(PObject.getInterPotentialGroupII(2));
	        MayerEGeneral eTargetIJ = new MayerEGeneral(PObject.getInterPotentialGroupIJ());
	        

	        ClusterAbstract targetCluster = Standard.virialClusterMixture(SimEnv.getnPoints(), new MayerFunction[][]{{fTargetII,fTargetIJ},{fTargetIJ,fTargetJJ}},
	                new MayerFunction[][]{{eTargetII,eTargetIJ},{eTargetIJ,eTargetJJ}}, nTypes);
	        targetCluster.setTemperature(temperature);
	        
	       
	        
	        SpeciesFactory[] speciesFactory = new SpeciesFactory[2];
	        
	        speciesFactory[0] = potential1.createSpeciesFactory();
	        speciesFactory[1] = potential2.createSpeciesFactory();
	        final SimulationVirialMultiOverlap sim = new SimulationVirialMultiOverlap(space, speciesFactory,temperature,refCluster,targetCluster,nTypes);
	        
	      /*  final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new Species[]{(Species)potential1.createSpeciesFactory(),(Species)potential2.createSpeciesFactory()}, nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},
	                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},SimEnv.isDoWiggle());
	        */
	        sim.integratorOS.setNumSubSteps(1000);
	        if(PObject.getBondedPotentialSets(1) != null){
	        PObject.setIntraPotentialGroup(1, sim.integrators[1].getPotentialMaster().makePotentialGroup(1));
	        
	        }
	        
	        if(PObject.getBondedPotentialSets(2) != null){
		        PObject.setIntraPotentialGroup(2, sim.integrators[1].getPotentialMaster().makePotentialGroup(1));
		        
	        
	        }
	        
	        Map AtomSetMapPureII = PObject.getAtomSetsPure(1);
	        Map AtomSetMapPureJJ = PObject.getAtomSetsPure(2);
	        Map AtomSetMapMix = PObject.getAtomSetsMix();
	        
	        AddPotentials(AtomSetMapPureII,PObject.getInterPotentialGroupII(1),PObject.getPotentialSets());
	        AddPotentials(AtomSetMapPureJJ,PObject.getInterPotentialGroupII(1),PObject.getPotentialSets());
	        AddPotentials(AtomSetMapMix,PObject.getInterPotentialGroupII(1),PObject.getPotentialSets());
	        
	        
	        
	        
	        
	        
	       
	        
	        
	        
	        
	    	
	    }else if(PObject instanceof AtomicPotentialCollection){

	    	PObject.setInterPotentialGroupII(new PotentialGroup(2));
	    	
	    	MayerGeneral fTargetII = new MayerGeneral(PObject.getInterPotentialGroupII());
	        MayerEGeneral eTargetII = new MayerEGeneral(PObject.getInterPotentialGroupII());
	        
	        ClusterAbstract targetCluster = Standard.virialCluster(SimEnv.getnPoints(), fTargetII, SimEnv.getnPoints()>3, eTargetII, SimEnv.isDoWiggle());
	        targetCluster.setTemperature(temperature);
	        
	        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,(Species)potential1.createSpecies(),
	        		temperature,refCluster,targetCluster, SimEnv.isDoWiggle());
	        
	        sim.integratorOS.setNumSubSteps(1000);
	        
	        if(PObject.getBondedPotentialSets() != null){
		        PObject.setIntraPotentialGroup(sim.integrators[1].getPotentialMaster().makePotentialGroup(1));}
	        

	    	
	    }else if(PObject instanceof MixedPotentialsCollection){
	    	
	    	MayerGeneral fTargetII = new MayerGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure());
	        MayerEGeneral eTargetII = new MayerEGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure());
	    	
	    	PObject.setInterPotentialGroupII(new PotentialGroup(2));
	        MayerGeneral fTargetJJ= new MayerGeneral(PObject.getInterPotentialGroupII());
	        MayerEGeneral eTargetJJ = new MayerEGeneral(PObject.getInterPotentialGroupII());
	        
	        
	    	PObject.setInterPotentialGroupIJ(new PotentialGroup(2));
	        MayerGeneral fTargetIJ= new MayerGeneral(PObject.getInterPotentialGroupIJ());
	        MayerEGeneral eTargetIJ = new MayerEGeneral(PObject.getInterPotentialGroupIJ());
	        
	    
	        
	        ClusterAbstract targetCluster = Standard.virialClusterMixture(SimEnv.getnPoints(), new MayerFunction[][]{{fTargetII,fTargetIJ},{fTargetIJ,fTargetJJ}},
	                new MayerFunction[][]{{eTargetII,eTargetIJ},{eTargetIJ,eTargetJJ}}, nTypes);
	        targetCluster.setTemperature(temperature);
	        
	        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,
	        		(etomica.api.IPotentialMolecular.class.isAssignableFrom(potential1.getPotential())) ? 
	        				new ISpecies[]{(Species)potential1.createSpecies(),(Species)potential2.createSpecies()} : 
	        				new ISpecies[]{(Species)potential2.createSpecies(),(Species)potential1.createSpecies()}, 
	        					nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},
	        					new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},SimEnv.isDoWiggle());
	        				
	        sim.integratorOS.setNumSubSteps(1000);
	        
	        if(PObject.getBondedPotentialSets() != null){
		        PObject.setIntraPotentialGroup(sim.integrators[1].getPotentialMaster().makePotentialGroup(1));}
	        
	        
	    	
	    	
	    }else if(PObject instanceof MolecularPotentialsCollection){
	    	
	    	MayerGeneral fTargetII = new MayerGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure(1));
	        MayerEGeneral eTargetII = new MayerEGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure(1));
	        
	        MayerGeneral fTargetJJ = new MayerGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure(2));
	        MayerEGeneral eTargetJJ = new MayerEGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure(2));
	        
	        MayerGeneral fTargetIJ = new MayerGeneral((IPotentialMolecular) PObject.getMolecularPotentialCross());
	        MayerEGeneral eTargetIJ = new MayerEGeneral((IPotentialMolecular) PObject.getMolecularPotentialCross());
	        


	        
	        ClusterAbstract targetCluster = Standard.virialClusterMixture(SimEnv.getnPoints(), new MayerFunction[][]{{fTargetII,fTargetIJ},{fTargetIJ,fTargetJJ}},
	                new MayerFunction[][]{{eTargetII,eTargetIJ},{eTargetIJ,eTargetJJ}}, nTypes);
	        targetCluster.setTemperature(temperature);
	        
	        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,(Species)potential1.createSpecies(), 
	        					temperature,new ClusterAbstract[]{refCluster,targetCluster},
	        					new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},SimEnv.isDoWiggle());
	        sim.integratorOS.setNumSubSteps(1000);
	        
	    	
	    }else if(PObject instanceof MolecularPotentials2Collection){
	    	
	    	MayerGeneral fTargetII = new MayerGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure(1));
	        MayerEGeneral eTargetII = new MayerEGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure(1));
	        
	        MayerGeneral fTargetJJ = new MayerGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure(2));
	        MayerEGeneral eTargetJJ = new MayerEGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure(2));
	        
	        PObject.setInterPotentialGroupIJ(new PotentialGroup(2));
	        MayerGeneral fTargetIJ= new MayerGeneral(PObject.getInterPotentialGroupIJ());
	        MayerEGeneral eTargetIJ = new MayerEGeneral(PObject.getInterPotentialGroupIJ());
	        
	        ClusterAbstract targetCluster = Standard.virialClusterMixture(SimEnv.getnPoints(), new MayerFunction[][]{{fTargetII,fTargetIJ},{fTargetIJ,fTargetJJ}},
	                new MayerFunction[][]{{eTargetII,eTargetIJ},{eTargetIJ,eTargetJJ}}, nTypes);
	        targetCluster.setTemperature(temperature);
	        
	        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,
	        				new ISpecies[]{(Species)potential1.createSpecies(),(Species)potential2.createSpecies()}, 
	        					nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},
	        					new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},SimEnv.isDoWiggle());
	        sim.integratorOS.setNumSubSteps(1000);
	        
	        

	    	
	    }else if(PObject instanceof MolecularPotentialCollection){
	    	

	    	MayerGeneral fTargetII = new MayerGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure());
	        MayerEGeneral eTargetII = new MayerEGeneral((IPotentialMolecular) PObject.getMolecularPotentialPure());
	        
	        ClusterAbstract targetCluster = Standard.virialCluster(SimEnv.getnPoints(), fTargetII, SimEnv.getnPoints()>3, eTargetII, true);
	        targetCluster.setTemperature(temperature);
	        
	        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,(Species)potential1.createSpecies(),temperature,refCluster,targetCluster,SimEnv.isDoWiggle());
	        sim.integratorOS.setNumSubSteps(1000);
	        
	        
	    	
	    }
	     
	     /*
	     
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

	        
		
	}*/

	
	}
	
	private void AddPotentials(Map atomSet,
			PotentialGroup interPotentialGroup, HashMap<String[], IPotential> PotentialSet) {
		// TODO Auto-generated method stub
		
		 	Map AtomSet = atomSet;
			Set AtomEntries = AtomSet.entrySet();
			Iterator AtomItr = AtomEntries.iterator();
			
			
			Map PotentialSetsMap = PotentialSet;
			Set PotentialEntries = PotentialSetsMap.entrySet();
			Iterator PotentialItr = PotentialEntries.iterator();
			
			while(PotentialItr.hasNext()){
				Map.Entry PotentialEntry = (Map.Entry) PotentialItr.next();
				
				String[] PotentialKey= (String[]) PotentialEntry.getKey();
				while (AtomItr.hasNext()) {
					Map.Entry AtomEntry = (Map.Entry) AtomItr.next();
					String[] AtomKey= (String[]) AtomEntry.getKey();
					if(AtomKey[0] == PotentialKey[0] && AtomKey[1] == PotentialKey[1] || AtomKey[0] == PotentialKey[1] && AtomKey[1] == PotentialKey[0]){
						Object[] AtomObjects = (Object[]) AtomEntry.getValue();
						
						if(PotentialEntry.getValue() instanceof P2LennardJones){
							interPotentialGroup.addPotential((P2LennardJones)PotentialEntry.getValue(),
									ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{(IAtomType)AtomObjects[0],(IAtomType)AtomObjects[1]}));
						}
						
					}
				}
				if(!AtomItr.hasNext()){
					AtomItr = AtomEntries.iterator();
				}
         }
		
	}


	public void runSimulation(SimulationEnvironmentObject simenv){
		
	}

}
