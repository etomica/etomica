/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;


import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.atom.IAtomType;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;

import etomica.api.IPotentialMolecular;
import etomica.api.ISpecies;

import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIntergroup;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
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
import etomica.potential.P4BondTorsion;
import etomica.potential.Potential2Spherical;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.units.CompoundDimension;
import etomica.units.CompoundUnit;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.Volume;
import etomica.util.Constants.CompassDirection;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.MCMoveClusterTorsionMulti;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.GUI.models.ModelSimulationEnvironment;
import etomica.virial.GUI.models.ModelSelectedSpecies;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap2;

public class SimulationRunner {
	
	
	private ModelSelectedSpecies speciesDataModel;
	private ModelSimulationEnvironment SimEnv;
	private Space space;
	
	public SimulationRunner(ModelSelectedSpecies speciesDataModel){
		this.speciesDataModel = speciesDataModel;
	}

	
	@SuppressWarnings("unchecked")
	public void runSimulation(ModelSimulationEnvironment simEnv, ArrayList<ICollectionPotential> pArrayListPotentialCollection) throws NoSuchMethodException{
		
		SimEnv = simEnv;
		
		final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(SimEnv.getSystemSigmaHSRef());
        HSB[3] = Standard.B3HS(SimEnv.getSystemSigmaHSRef());
        HSB[4] = Standard.B4HS(SimEnv.getSystemSigmaHSRef());
        HSB[5] = Standard.B5HS(SimEnv.getSystemSigmaHSRef());
        HSB[6] = Standard.B6HS(SimEnv.getSystemSigmaHSRef());
        HSB[7] = Standard.B7HS(SimEnv.getSystemSigmaHSRef());
        
        System.out.println("sigmaHSRef: "+SimEnv.getSystemSigmaHSRef());
        System.out.println("B"+SimEnv.getnPoints()+"HS: "+HSB[SimEnv.getnPoints()]);
        
        final int nPoints = SimEnv.getnPoints();
        int[] tempNTypes = speciesDataModel.getSpeciesMoleculeCount();
        int[] nTypes = new int[speciesDataModel.getNthSpeciesAdded()];
        for(int i=0;i<speciesDataModel.getNthSpeciesAdded();i++){
        	nTypes[i] = tempNTypes[i];
        }
        
        double temperature = SimEnv.getTemperature();
        space = Space3D.getInstance();
        
		MayerHardSphere fRef = new MayerHardSphere(SimEnv.getSystemSigmaHSRef());
	    MayerEHardSphere eRef = new MayerEHardSphere(SimEnv.getSystemSigmaHSRef());
	    
	    ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);
        //hashmap collection
        
	 	
        if(pArrayListPotentialCollection.size() == 1){
        		ClusterAbstract targetCluster = null;
        		
        		
        		if(pArrayListPotentialCollection.get(0) instanceof CollectionPotentialAtomicLike){
        			
        			CollectionPotentialAtomicLike tempCollectionPotential = (CollectionPotentialAtomicLike) pArrayListPotentialCollection.get(0);
        			if(tempCollectionPotential.getHashMapPotentialNonBonded().getHashMapPotentialNonBonded().size() > 1){
        				tempCollectionPotential.setPotentialGroupInterNonBondedLike(new PotentialGroup(2));
        				
        				MayerGeneral fTarget = new MayerGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedLike());
				        MayerEGeneral eTarget = new MayerEGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedLike());
				        
				        targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
				        targetCluster.setTemperature(temperature);
				        
				        
				        
				        AddPotentials(tempCollectionPotential.getHashMapAtomTypeLikePairs(),tempCollectionPotential.getPotentialGroupInterNonBondedLike(),
				        				tempCollectionPotential.getHashMapPotentialNonBonded().getHashMapPotentialNonBonded());
				        
				        ISpecies species =  speciesDataModel.getSpecies(0);
				        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species,temperature,refCluster,targetCluster);
				        
		    			
		    			if(tempCollectionPotential.getHashMapPotentialIntraBonded() != null){
		    				tempCollectionPotential.setPotentialGroupIntraBonded(sim.integrators[1].getPotentialMaster().makePotentialGroup(1));
		    				AddBondedPotentials(tempCollectionPotential.getPotentialGroupIntraBonded(),tempCollectionPotential.getHashMapPotentialIntraBonded(),
		    										tempCollectionPotential.getHashmapAtomsetIterator(),tempCollectionPotential.getHashMapPotentialNonBonded().getHashMapPotentialNonBonded());
		    	        	if(tempCollectionPotential.getMCMoveClusterTorsionMulti() != null){
		    	        		tempCollectionPotential.getMCMoveClusterTorsionMulti()[0] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, (P4BondTorsion) tempCollectionPotential.getPotentialIntraBondedTorsion(), 40);
		    	        		tempCollectionPotential.getMCMoveClusterTorsionMulti()[0].setTemperature(temperature);
		    	        		sim.integrators[0].getMoveManager().addMCMove(tempCollectionPotential.getMCMoveClusterTorsionMulti()[0]);
		    	        	}
		    			}
		    			
		    			
		    			new IIntegratorListener() {
		    		    	public void integratorInitialized(IIntegratorEvent e) {}
		    		    	public void integratorStepStarted(IIntegratorEvent e) {}
		    		    	public void integratorStepFinished(IIntegratorEvent e) {
		    		    		if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
		    		    		
		    		    		System.out.print(sim.integratorOS.getStepCount()+" steps: ");
		    		            double[] ratioAndError = sim.dvo.getAverageAndError();
		    		            double ratio = ratioAndError[0];
		    		            double error = ratioAndError[1];
		    		            System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
		    		        }
		    		    	
		    		    };
		    		    
		    		    sim.getController().addAction(new IAction() {
		    	            public void actionPerformed() {
		    	                sim.initRefPref(null, 10);
		    	                sim.equilibrate(null, 20);
		    	                sim.ai.setMaxSteps(Long.MAX_VALUE);
		    	            }
		    	        });
		    		    
		    		    final DisplayTextBox averageBox = new DisplayTextBox();
        		        averageBox.setLabel("Average");
        		        final DisplayTextBox errorBox = new DisplayTextBox();
        		        errorBox.setLabel("Error");
        			    IAction pushAnswer = new IAction() {
        		            public void actionPerformed() {
        		                double[] ratioAndError = sim.dvo.getAverageAndError();
        		                double ratio = ratioAndError[0];
        		                double error = ratioAndError[1];
        		                data.x = ratio;
        		                averageBox.putData(data);
        		                data.x = error;
        		                errorBox.putData(data);
        		            }
        		            
        		            DataDouble data = new DataDouble();
        		        };
		    		    
        		        runGraphic(sim,nPoints,averageBox,errorBox,pushAnswer);
		    		    //showReport(sim,progressReport);
		    		    sim.printResults(HSB[nPoints]);
        			}
        			else{
        				//for a single site and atomic potential, forexample LJ
        				Map<String[],IPotential> potentialMap = ((CollectionPotentialAtomicLike) pArrayListPotentialCollection.get(0)).
								getHashMapPotentialNonBonded().getHashMapPotentialNonBonded();
        				Set<Entry<String[],IPotential>> potentialMapEntries = potentialMap.entrySet();
        				Iterator<Entry<String[],IPotential>> potentialMapItr = potentialMapEntries.iterator();
        				
        				
        				Map.Entry potentialEntry = potentialMapItr.next();
            			potentialEntry.getKey();
            			IPotentialAtomic potential = (IPotentialAtomic) potentialEntry.getValue();
            			
            			
        				if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(potential.getClass()) ){
    						MayerGeneralSpherical fTarget = new MayerGeneralSpherical((Potential2Spherical) potential);
    						MayerESpherical eTarget = new MayerESpherical((Potential2Spherical) potential);
    						targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
    						targetCluster.setTemperature(temperature);
    						
    					}else if(etomica.potential.PotentialSoft.class.isAssignableFrom(potential.getClass())){
    						PotentialGroup pTargetGroup = new PotentialGroup(2);
    						MayerGeneral fTarget = new MayerGeneral(pTargetGroup);
    				        MayerEGeneral eTarget = new MayerEGeneral(pTargetGroup);
    				        targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
    				        targetCluster.setTemperature(temperature);
    				        pTargetGroup.addPotential(potential, new ApiIntergroup());
    						
    					}
        				ISpecies species =  speciesDataModel.getSpecies(0);
        				final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species,temperature,refCluster,targetCluster);
        				new IIntegratorListener() {
        			    	public void integratorInitialized(IIntegratorEvent e) {}
        			    	public void integratorStepStarted(IIntegratorEvent e) {}
        			    	public void integratorStepFinished(IIntegratorEvent e) {
        			    		if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
        			    		
        			    		System.out.print(sim.integratorOS.getStepCount()+" steps: ");
        			            double[] ratioAndError = sim.dvo.getAverageAndError();
        			            double ratio = ratioAndError[0];
        			            double error = ratioAndError[1];
        			            System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
        			        }
        			    };
        			    
        			    sim.getController().addAction(new IAction() {
        		            public void actionPerformed() {
        		                sim.initRefPref(null, 10);
        		                sim.equilibrate(null, 20);
        		                sim.ai.setMaxSteps(Long.MAX_VALUE);
        		            }
        		        });
        			    
        			    final DisplayTextBox averageBox = new DisplayTextBox();
        		        averageBox.setLabel("Average");
        		        final DisplayTextBox errorBox = new DisplayTextBox();
        		        errorBox.setLabel("Error");
        			    IAction pushAnswer = new IAction() {
        		            public void actionPerformed() {
        		                double[] ratioAndError = sim.dvo.getAverageAndError();
        		                double ratio = ratioAndError[0];
        		                double error = ratioAndError[1];
        		                data.x = ratio;
        		                averageBox.putData(data);
        		                data.x = error;
        		                errorBox.putData(data);
        		            }
        		            
        		            DataDouble data = new DataDouble();
        		        };
        		        runGraphic(sim,nPoints,averageBox,errorBox,pushAnswer);
        			    //showReport(sim,progressReport);
        			    
        				sim.printResults(HSB[nPoints]);
        			}
        		
        		}
        		else if(pArrayListPotentialCollection.get(0) instanceof CollectionPotentialMolecularLike){
        			IPotentialMolecular targetPotential = ((CollectionPotentialMolecularLike) pArrayListPotentialCollection.get(0)).getPotentialMolecularNonBondedLike();
					MayerGeneral fTarget = new MayerGeneral(targetPotential);
					MayerEGeneral eTarget = new MayerEGeneral(targetPotential);
					targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
					targetCluster.setTemperature(temperature);
					
					ISpecies species =  speciesDataModel.getSpecies(0);
					final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species,temperature,refCluster,targetCluster);
					
					new IIntegratorListener() {
				    	public void integratorInitialized(IIntegratorEvent e) {}
				    	public void integratorStepStarted(IIntegratorEvent e) {}
				    	public void integratorStepFinished(IIntegratorEvent e) {
				    		if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
				    		
				    		System.out.print(sim.integratorOS.getStepCount()+" steps: ");
				            double[] ratioAndError = sim.dvo.getAverageAndError();
				            double ratio = ratioAndError[0];
				            double error = ratioAndError[1];
				            System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
				        }
				    };
				    
				    sim.getController().addAction(new IAction() {
			            public void actionPerformed() {
			                sim.initRefPref(null, 10);
			                sim.equilibrate(null, 20);
			                sim.ai.setMaxSteps(Long.MAX_VALUE);
			            }
			        });
				    
				    final DisplayTextBox averageBox = new DisplayTextBox();
    		        averageBox.setLabel("Average");
    		        final DisplayTextBox errorBox = new DisplayTextBox();
    		        errorBox.setLabel("Error");
    			    IAction pushAnswer = new IAction() {
    		            public void actionPerformed() {
    		                double[] ratioAndError = sim.dvo.getAverageAndError();
    		                double ratio = ratioAndError[0];
    		                double error = ratioAndError[1];
    		                data.x = ratio;
    		                averageBox.putData(data);
    		                data.x = error;
    		                errorBox.putData(data);
    		            }
    		            
    		            DataDouble data = new DataDouble();
    		        };
    		        
				    //showReport(sim,progressReport);
				    runGraphic(sim,nPoints,averageBox,errorBox,pushAnswer);
				    sim.printResults(HSB[nPoints]);
        		}
	
       }else{
    	  ClusterAbstract targetCluster = null;
    	  int noOfSpecies = speciesDataModel.getNthSpeciesAdded();
    	  boolean[] bondedInteractionIndex = new boolean[noOfSpecies];
    	  for(int i=0;i<bondedInteractionIndex.length;i++){
    		  bondedInteractionIndex[i] = false;
    	  }
    	  MayerFunction[][] f = new MayerFunction[noOfSpecies][noOfSpecies];
    	  MayerFunction[][] e = new MayerFunction[noOfSpecies][noOfSpecies];
    	   
    	  speciesDataModel.getSpeciesInteractionIndex();
    	   
    	 
    		  
    		  for(int j=0;j<pArrayListPotentialCollection.size();j++){
    			  if(pArrayListPotentialCollection.get(j) instanceof CollectionPotentialAtomicLike){
    				  CollectionPotentialAtomicLike tempCollectionPotential = (CollectionPotentialAtomicLike) pArrayListPotentialCollection.get(j);
    				  int speciesIndex = tempCollectionPotential.getSpeciesIndex();
    				  if(speciesDataModel.getSpeciesDataModel(speciesIndex).getPotentialSites().length == 1 && 
    						  speciesDataModel.getSpeciesDataModel(speciesIndex).getPotentialSiteAtIndex(0).contains("LJ")){
    					  
    					  
    					  if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom
       						   (speciesDataModel.getSpeciesDataModel(speciesIndex).getPotential()) ){
    						  IPotential potential = null;
    						  Map<String[],IPotential> potentialMap = ((CollectionPotentialAtomicLike) pArrayListPotentialCollection.get(j)).getHashMapPotentialNonBonded().getHashMapPotentialNonBonded();
    	        			  Set<Entry<String[],IPotential>> potentialMapEntries = potentialMap.entrySet();
    	        			  Iterator<Entry<String[],IPotential>> potentialMapItr = potentialMapEntries.iterator();
    	        			  if(potentialMap.size() > 1){
    	        				  while(potentialMapItr.hasNext()){
    	        					  Map.Entry potentialEntry = potentialMapItr.next();
    	        					  String[] potentialMapKey= (String[]) potentialEntry.getKey();
    	        					  potentialEntry.getValue();
    	        				  
    	        					  if(potentialMapKey[0] == speciesDataModel.getSpeciesDataModel(speciesIndex).getPotentialSiteAtIndex(0) &&
    	        						  potentialMapKey[1] == speciesDataModel.getSpeciesDataModel(speciesIndex).getPotentialSiteAtIndex(0)){
    	        						  potential = (IPotential) potentialEntry.getValue();
    	        					  }
    	        				  }
    	        			  }else{
    	        				  Map.Entry potentialEntry = potentialMapItr.next();
	        					  potentialEntry.getKey();
	        					  potential= (IPotential) potentialEntry.getValue();
    	        			  }

    						  MayerGeneralSpherical fTarget = new MayerGeneralSpherical((Potential2Spherical) potential);
    						  MayerESpherical eTarget = new MayerESpherical((Potential2Spherical) potential);

      				          f[speciesIndex][speciesIndex] = fTarget;
      				          e[speciesIndex][speciesIndex] = eTarget;
    						  
    					  }else if( etomica.potential.PotentialSoft.class.isAssignableFrom
          						   (speciesDataModel.getSpeciesDataModel(speciesIndex).getPotential()) ){
    						  tempCollectionPotential.setPotentialGroupInterNonBondedLike(new PotentialGroup(2));
        					  MayerGeneral fTarget = new MayerGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedLike());
      				          MayerEGeneral eTarget = new MayerEGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedLike());
      				          
      				          AddPotentials(tempCollectionPotential.getHashMapAtomTypeLikePairs(),tempCollectionPotential.getPotentialGroupInterNonBondedLike(),
      				        		  			tempCollectionPotential.getHashMapPotentialNonBonded().getHashMapPotentialNonBonded());
      				          
      				          f[speciesIndex][speciesIndex] = fTarget;
      				          e[speciesIndex][speciesIndex] = eTarget;
    					  }
    				  }else{
    					  tempCollectionPotential.setPotentialGroupInterNonBondedLike(new PotentialGroup(2));
    					  MayerGeneral fTarget = new MayerGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedLike());
  				          MayerEGeneral eTarget = new MayerEGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedLike());
  				          
  				          AddPotentials(tempCollectionPotential.getHashMapAtomTypeLikePairs(),tempCollectionPotential.getPotentialGroupInterNonBondedLike(),
  				        		  			tempCollectionPotential.getHashMapPotentialNonBonded().getHashMapPotentialNonBonded());
  				          
  				          f[speciesIndex][speciesIndex] = fTarget;
  				          e[speciesIndex][speciesIndex] = eTarget;
  				          
  				        if(tempCollectionPotential.getHashMapPotentialIntraBonded() != null){
  				        	bondedInteractionIndex[j]= true;
  				        }
    					  
    				  }
    				  
    			  }else if(pArrayListPotentialCollection.get(j) instanceof CollectionPotentialMolecularLike){
    				  CollectionPotentialMolecularLike tempCollectionPotential = (CollectionPotentialMolecularLike) pArrayListPotentialCollection.get(j);
    				  int speciesIndex = tempCollectionPotential.getSpeciesIndex();
					  MayerGeneral fTarget = new MayerGeneral(tempCollectionPotential.getPotentialMolecularNonBondedLike());
				      MayerEGeneral eTarget = new MayerEGeneral(tempCollectionPotential.getPotentialMolecularNonBondedLike());
				          
				      f[speciesIndex][speciesIndex] = fTarget;
				      e[speciesIndex][speciesIndex] = eTarget;
    			  
    			  }else if(pArrayListPotentialCollection.get(j) instanceof CollectionPotentialAtomicUnlike){
    				  CollectionPotentialAtomicUnlike tempCollectionPotential = (CollectionPotentialAtomicUnlike) pArrayListPotentialCollection.get(j);
    				  int[] speciesIndex = tempCollectionPotential.getSpeciesIndex();
    				  
    				  if(speciesDataModel.getSpeciesDataModel(speciesIndex[0]).getPotentialSites().length == 1 && 
    						  speciesDataModel.getSpeciesDataModel(speciesIndex[1]).getPotentialSites().length == 1 &&
    						  speciesDataModel.getSpeciesDataModel(speciesIndex[0]).getPotentialSiteAtIndex(0).contains("LJ") &&
    						  speciesDataModel.getSpeciesDataModel(speciesIndex[1]).getPotentialSiteAtIndex(0).contains("LJ")){
    					  if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom
          						   (speciesDataModel.getSpeciesDataModel(speciesIndex[0]).getPotential()) ){
    						  IPotential potential = null;
    						  Map<String[],IPotential> potentialMap = ((CollectionPotentialAtomicUnlike) pArrayListPotentialCollection.get(j)).getHashMapPotentialNonBonded().getHashMapPotentialNonBonded();
    	        			  Set<Entry<String[],IPotential>> potentialMapEntries = potentialMap.entrySet();
    	        			  Iterator<Entry<String[],IPotential>> potentialMapItr = potentialMapEntries.iterator();
    	        			  
    	        			  if(potentialMap.size() > 1){
    	        				  while(potentialMapItr.hasNext()){
    	        					  Map.Entry potentialEntry = potentialMapItr.next();
    	        					  String[] potentialMapKey= (String[]) potentialEntry.getKey();
    	        					  potentialEntry.getValue();
    	        				  
    	        					  if(potentialMapKey[0] == speciesDataModel.getSpeciesDataModel(speciesIndex[0]).getPotentialSiteAtIndex(0) &&
    	        						  potentialMapKey[1] == speciesDataModel.getSpeciesDataModel(speciesIndex[1]).getPotentialSiteAtIndex(0)){
    	        						  potential = (IPotential) potentialEntry.getValue();
    	        					  }
    	        				  }
    	        			  }else{
    	        				  Map.Entry potentialEntry = potentialMapItr.next();
	        					  potentialEntry.getKey();
	        					  potential= (IPotential) potentialEntry.getValue();
    	        			  }
    	        			  
    	        			  MayerGeneralSpherical fTarget = new MayerGeneralSpherical((Potential2Spherical) potential);
    						  MayerESpherical eTarget = new MayerESpherical((Potential2Spherical) potential);
    						  
    						  f[speciesIndex[0]][speciesIndex[1]] = fTarget;
    					      e[speciesIndex[0]][speciesIndex[1]] = eTarget;
    					      
    					      f[speciesIndex[1]][speciesIndex[0]] = fTarget;
    					      e[speciesIndex[1]][speciesIndex[0]] = eTarget;
    	    			  
    						  
    					  }else if(etomica.potential.PotentialSoft.class.isAssignableFrom
         						   (speciesDataModel.getSpeciesDataModel(speciesIndex[0]).getPotential())){
    						  tempCollectionPotential.setPotentialGroupInterNonBondedUnlike(new PotentialGroup(2));
    	    				  
    	    				  
        					  MayerGeneral fTarget = new MayerGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedUnlike());
        					  MayerEGeneral eTarget = new MayerEGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedUnlike());
        				  
        					  AddPotentials(tempCollectionPotential.getHashMapAtomTypesUnlikePairs(),tempCollectionPotential.getPotentialGroupInterNonBondedUnlike(),
    	        		  			tempCollectionPotential.getHashMapPotentialNonBonded().getHashMapPotentialNonBonded());
    				      
    				
        					  f[speciesIndex[0]][speciesIndex[1]] = fTarget;
        					  e[speciesIndex[0]][speciesIndex[1]] = eTarget;
    				      
        					  f[speciesIndex[1]][speciesIndex[0]] = fTarget;
        					  e[speciesIndex[1]][speciesIndex[0]] = eTarget;
    					  }
    				  }else {
    					  
    				  
    					  tempCollectionPotential.setPotentialGroupInterNonBondedUnlike(new PotentialGroup(2));
    				  
    				  
    					  MayerGeneral fTarget = new MayerGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedUnlike());
    					  MayerEGeneral eTarget = new MayerEGeneral(tempCollectionPotential.getPotentialGroupInterNonBondedUnlike());
    				  
    					  AddPotentials(tempCollectionPotential.getHashMapAtomTypesUnlikePairs(),tempCollectionPotential.getPotentialGroupInterNonBondedUnlike(),
	        		  			tempCollectionPotential.getHashMapPotentialNonBonded().getHashMapPotentialNonBonded());
				      
				
    					  f[speciesIndex[0]][speciesIndex[1]] = fTarget;
    					  e[speciesIndex[0]][speciesIndex[1]] = eTarget;
				      
    					  f[speciesIndex[1]][speciesIndex[0]] = fTarget;
    					  e[speciesIndex[1]][speciesIndex[0]] = eTarget;
    				  }
    			  
    			  }else if(pArrayListPotentialCollection.get(j) instanceof CollectionPotentialMolecularUnlike){
    				  CollectionPotentialMolecularUnlike tempCollectionPotential = (CollectionPotentialMolecularUnlike) pArrayListPotentialCollection.get(j);
    				  int[] speciesIndex = tempCollectionPotential.getSpeciesIndex();
    				  
    				  MayerGeneral fTarget = new MayerGeneral(tempCollectionPotential.getPotentialMolecularUnlike());
				      MayerEGeneral eTarget = new MayerEGeneral(tempCollectionPotential.getPotentialMolecularUnlike());
				          
				      f[speciesIndex[0]][speciesIndex[1]] = fTarget;
				      e[speciesIndex[0]][speciesIndex[1]] = eTarget;
				      
				      f[speciesIndex[1]][speciesIndex[0]] = fTarget;
				      e[speciesIndex[1]][speciesIndex[0]] = eTarget;
    			  }
    		  }
    		  targetCluster = Standard.virialClusterMixture(nPoints,f,e,nTypes);
    	      targetCluster.setTemperature(temperature);  
    	      
    	      ISpecies[] species =  new Species[noOfSpecies];
    	      for(int i=0;i<noOfSpecies;i++){
    	    	  species[i] = speciesDataModel.getSpecies(i);
    	      }
    	      
			  final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,species,nTypes,temperature,new ClusterAbstract[]{refCluster,targetCluster},
		                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)},true);
				
			  for(int k=0;k<bondedInteractionIndex.length;k++){
    			  if(bondedInteractionIndex[k] == true){
    				  CollectionPotentialAtomicLike tempCollectionPotential = (CollectionPotentialAtomicLike) pArrayListPotentialCollection.get(k);
    				  if(tempCollectionPotential.getHashMapPotentialIntraBonded() != null){
    					  tempCollectionPotential.setPotentialGroupIntraBonded(sim.integrators[1].getPotentialMaster().makePotentialGroup(1));
    					  AddBondedPotentials(tempCollectionPotential.getPotentialGroupIntraBonded(),tempCollectionPotential.getHashMapPotentialIntraBonded(),
		    										tempCollectionPotential.getHashmapAtomsetIterator(),tempCollectionPotential.getHashMapPotentialNonBonded().getHashMapPotentialNonBonded());
    					  if(tempCollectionPotential.getMCMoveClusterTorsionMulti() != null){
    						  tempCollectionPotential.getMCMoveClusterTorsionMulti()[0] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, (P4BondTorsion) tempCollectionPotential.getPotentialIntraBondedTorsion(), 40);
    						  tempCollectionPotential.getMCMoveClusterTorsionMulti()[0].setTemperature(temperature);
    						  sim.integrators[0].getMoveManager().addMCMove(tempCollectionPotential.getMCMoveClusterTorsionMulti()[0]);
    					  }
    				  }
    			  }
			  }
			  
			  IIntegratorListener progressReport = new IIntegratorListener() {
  		    	public void integratorInitialized(IIntegratorEvent e) {}
  		    	public void integratorStepStarted(IIntegratorEvent e) {}
  		    	public void integratorStepFinished(IIntegratorEvent e) {
  		    		if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
  		    		
  		    		System.out.print(sim.integratorOS.getStepCount()+" steps: ");
  		            double[] ratioAndError = sim.dvo.getAverageAndError();
  		            double ratio = ratioAndError[0];
  		            double error = ratioAndError[1];
  		            System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
  		        }
  		    	
  		    };
  		    
  		    sim.getController().addAction(new IAction() {
	            public void actionPerformed() {
	                sim.initRefPref(null, 10);
	                sim.equilibrate(null, 20);
	                sim.ai.setMaxSteps(Long.MAX_VALUE);
	            }
	        });
  		    
  		  final DisplayTextBox averageBox = new DisplayTextBox();
	        averageBox.setLabel("Average");
	        final DisplayTextBox errorBox = new DisplayTextBox();
	        errorBox.setLabel("Error");
		    new IAction() {
	            public void actionPerformed() {
	                double[] ratioAndError = sim.dvo.getAverageAndError();
	                double ratio = ratioAndError[0];
	                double error = ratioAndError[1];
	                data.x = ratio;
	                averageBox.putData(data);
	                data.x = error;
	                errorBox.putData(data);
	            }
	            
	            DataDouble data = new DataDouble();
	        };
	        
		    //showReport(sim,progressReport);
		    //runGraphic(sim,nPoints,averageBox,errorBox,pushAnswer);
  		    
  		    showReport(sim,progressReport);
			  
    		  
       }
	}
	
	private void AddBondedPotentials(PotentialGroup intraPotentialGroupII,
			HashMap<String[], IPotential> bondedPotentialSets,
			HashMap<String[], AtomsetIteratorBasisDependent> iteratorSets, HashMap<String[], IPotential> potentialSets) {
		// TODO Auto-generated method stub
		Map<String[], IPotential> bondedPotentialsMap = bondedPotentialSets;
		Set<Entry<String[], IPotential>> bondedPotentialEntries = bondedPotentialsMap.entrySet();
		Iterator<Entry<String[], IPotential>> bondedPotentialsItr = bondedPotentialEntries.iterator();
		
		
		Map<String[], AtomsetIteratorBasisDependent> atomGroupsIteratorsMap = iteratorSets;
		Set<Entry<String[], AtomsetIteratorBasisDependent>> atomGroupsIteratorEntries = atomGroupsIteratorsMap.entrySet();
		Iterator<Entry<String[], AtomsetIteratorBasisDependent>> atomGroupsItr = atomGroupsIteratorEntries.iterator();
		
		Map<String[], IPotential> potentialsMap = potentialSets;
		Set<Entry<String[], IPotential>> potentialsEntries = potentialsMap.entrySet();
		Iterator<Entry<String[], IPotential>> potentialsItr = potentialsEntries.iterator();
		
		while(atomGroupsItr.hasNext()){
			Map.Entry atomGroupsEntry = atomGroupsItr.next();
			String[] atomGroupsMapKey= (String[]) atomGroupsEntry.getKey();
			
			while (bondedPotentialsItr.hasNext()) {
				Map.Entry bondedPotentialsEntry = bondedPotentialsItr.next();
				String[] bondedPotentialsMapKey= (String[]) bondedPotentialsEntry.getKey();
				if(bondedPotentialsMapKey[0] == atomGroupsMapKey[0] && bondedPotentialsMapKey[1] == atomGroupsMapKey[1] 
						|| bondedPotentialsMapKey[0] == atomGroupsMapKey[1] && bondedPotentialsMapKey[1] == atomGroupsMapKey[0]){
					
					intraPotentialGroupII.addPotential((IPotentialAtomic)bondedPotentialsEntry.getValue(), (AtomsetIteratorBasisDependent)atomGroupsEntry.getValue());
					
				}
			}
			if(!bondedPotentialsItr.hasNext()){
				bondedPotentialsItr = bondedPotentialEntries.iterator();
			}
			
		}
		
		if(!atomGroupsItr.hasNext()){
			atomGroupsItr = atomGroupsIteratorEntries.iterator();
		}
		
		while(atomGroupsItr.hasNext()){
			Map.Entry atomGroupsEntry = atomGroupsItr.next();
			String[] atomGroupsMapKey= (String[]) atomGroupsEntry.getKey();
			
			while (potentialsItr.hasNext()) {
				Map.Entry potentialsEntry = potentialsItr.next();
				String[] potentialsMapKey= (String[]) potentialsEntry.getKey();
				if(potentialsMapKey[0] == atomGroupsMapKey[0] && potentialsMapKey[1] == atomGroupsMapKey[1] 
						|| potentialsMapKey[0] == atomGroupsMapKey[1] && potentialsMapKey[1] == atomGroupsMapKey[0]){
					
					intraPotentialGroupII.addPotential((IPotentialAtomic)potentialsEntry.getValue(), (AtomsetIteratorBasisDependent)atomGroupsEntry.getValue());
				}
			}
			if(!bondedPotentialsItr.hasNext()){
				bondedPotentialsItr = bondedPotentialEntries.iterator();
			}
			
		}
		
		
		
	}


	private void AddPotentials(HashMap<String[], IAtomType[]> atomSet,
			PotentialGroup interPotentialGroup, HashMap<String[], IPotential> PotentialSet) {
		// TODO Auto-generated method stub
		
		 	HashMap<String[], IAtomType[]> atomSetMap = atomSet;
			Set<Entry<String[], IAtomType[]>> AtomEntries = atomSetMap.entrySet();
			Iterator<Entry<String[], IAtomType[]>> AtomItr = AtomEntries.iterator();
			
			
			Map<String[], IPotential> PotentialSetsMap = PotentialSet;
			Set<Entry<String[], IPotential>> PotentialEntries = PotentialSetsMap.entrySet();
			Iterator<Entry<String[], IPotential>> PotentialItr = PotentialEntries.iterator();
			
			while(PotentialItr.hasNext()){
				Map.Entry PotentialEntry = PotentialItr.next();
				
				String[] PotentialKey= (String[]) PotentialEntry.getKey();
				while (AtomItr.hasNext()) {
					Map.Entry AtomEntry = AtomItr.next();
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
	
	public void showReport(SimulationVirialOverlap2 sim, IIntegratorListener progressReport){
		sim.integratorOS.setNumSubSteps(1000);
		sim.initRefPref(null, SimEnv.getNoOfSteps()/100);
	        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
	        // if it does continue looking for a pref, it will write the value to the file
	    sim.equilibrate(null, SimEnv.getNoOfSteps()/40);
	    if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
	       throw new RuntimeException("oops");
	    }
	    System.out.println("equilibration finished");
	    sim.setAccumulatorBlockSize(SimEnv.getNoOfSteps());

	    sim.integratorOS.getEventManager().addListener(progressReport);
	    sim.integratorOS.getMoveManager().setEquilibrating(false);
	    sim.ai.setMaxSteps(SimEnv.getNoOfSteps());
	    for (int i=0; i<2; i++) {
	    	System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
	    }
	    

	    sim.getController().actionPerformed();

	    System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
	    System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
	    	
		
	}
	
	public void runGraphic(SimulationVirialOverlap2 sim, int nPoints,DisplayTextBox averageBox, DisplayTextBox errorBox, IAction pushAnswer){
		sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
        sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
        DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
        DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
//        displayBox0.setPixelUnit(new Pixel(300.0/size));
//        displayBox1.setPixelUnit(new Pixel(300.0/size));
        displayBox0.setShowBoundary(false);
        displayBox1.setShowBoundary(false);
        ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
        ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);
        
        
        ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
        displayBox0.setColorScheme(colorScheme);
        colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
        displayBox1.setColorScheme(colorScheme);
        simGraphic.makeAndDisplayFrame();

        sim.integratorOS.setNumSubSteps(1000);
        sim.setAccumulatorBlockSize(1000);
            
        // if running interactively, set filename to null so that it doens't read
        // (or write) to a refpref file
        sim.getController().removeAction(sim.ai);
        
        sim.getController().addAction(sim.ai);
        if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
            throw new RuntimeException("Oops");
        }
       
        
        JLabel jLabelPanelParentGroup = new JLabel("B"+nPoints+" (L/mol)^"+(nPoints-1));
        final JPanel panelParentGroup = new JPanel(new java.awt.BorderLayout());
        panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
        panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.WEST);
        panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
        simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());
        
        
        IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
        Unit unit = new CompoundUnit(new Unit[]{new UnitRatio(Liter.UNIT, Mole.UNIT)}, new double[]{nPoints-1});
        averageBox.putDataInfo(dataInfo);
        averageBox.setLabel("average");
        averageBox.setUnit(unit);
        errorBox.putDataInfo(dataInfo);
        errorBox.setLabel("error");
        errorBox.setPrecision(2);
        errorBox.setUnit(unit);
        sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));
        
        return;
	    	
		
	}
	
}
