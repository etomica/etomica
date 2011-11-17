package etomica.virial.GUI.models;

import etomica.virial.GUI.components.PotentialCollectionsBuilder;
import etomica.virial.GUI.components.SimulationFactory;
import etomica.virial.GUI.components.PotentialCollections;
import etomica.virial.GUI.components.SimulationEnvironment;
import etomica.virial.GUI.components.SimulationEnvironmentObject;

public class SimulationConstructorDM {

	
	private SimulationEnvironmentObject simEnvironmentObject;
	private SimulationEnvironment simEnvironment;
	private PotentialCollections potentialCollectionFactory;
	private SimulationFactory factory;
	private PotentialCollectionsBuilder  potentialsCollectionBuilder;
	
	
	SimulationConstructorDM(){
		reset();
	}


	private void reset() {
		// TODO Auto-generated method stub
		
		
	}


	public SimulationEnvironmentObject getSimEnvironmentObject() {
		return simEnvironmentObject;
	}


	public void setSimEnvironmentObject(
			SimulationEnvironmentObject simEnvironmentObject) {
		this.simEnvironmentObject = simEnvironmentObject;
	}


	public SimulationEnvironment getSimEnvironment() {
		return simEnvironment;
	}


	public void setSimEnvironment(SimulationEnvironment simEnvironment) {
		this.simEnvironment = simEnvironment;
	}


	public PotentialCollections getPotentialCollectionFactory() {
		return potentialCollectionFactory;
	}


	public void setPotentialCollectionFactory(
			PotentialCollections potentialCollectionFactory) {
		this.potentialCollectionFactory = potentialCollectionFactory;
	}


	public SimulationFactory getFactory() {
		return factory;
	}


	public void setFactory(SimulationFactory factory) {
		this.factory = factory;
	}


	public PotentialCollectionsBuilder getPotentialsCollectionBuilder() {
		return potentialsCollectionBuilder;
	}


	public void setPotentialsCollectionBuilder(
			PotentialCollectionsBuilder potentialsCollectionBuilder) {
		this.potentialsCollectionBuilder = potentialsCollectionBuilder;
	}

	
	
}
