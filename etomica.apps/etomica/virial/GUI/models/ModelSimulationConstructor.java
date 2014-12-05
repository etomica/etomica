/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import java.util.ArrayList;
import java.util.HashMap;

import etomica.virial.GUI.components.ICollectionPotential;
import etomica.virial.GUI.components.BuilderCollectionPotential;
import etomica.virial.GUI.components.SimulationRunner;


public class ModelSimulationConstructor {

	
	private ICollectionPotential potentialCollectionFactory;
	private SimulationRunner factory;
	private BuilderCollectionPotential  potentialsCollectionBuilder;
	private ArrayList<ICollectionPotential> arrayListPotentialCollection;
	
	public ModelSimulationConstructor(){
		
		reset();
	}


	private void reset() {
		// TODO Auto-generated method stub
		
		
	}


	public ICollectionPotential getPotentialCollectionFactory() {
		return potentialCollectionFactory;
	}


	public void setPotentialCollectionFactory(
			ICollectionPotential potentialCollectionFactory) {
		this.potentialCollectionFactory = potentialCollectionFactory;
	}


	public SimulationRunner getFactory() {
		return factory;
	}


	public void setFactory(SimulationRunner factory) {
		this.factory = factory;
	}


	public BuilderCollectionPotential getPotentialsCollectionBuilder() {
		return potentialsCollectionBuilder;
	}


	public void setPotentialsCollectionBuilder(
			BuilderCollectionPotential potentialsCollectionBuilder) {
		this.potentialsCollectionBuilder = potentialsCollectionBuilder;
	}


	public ArrayList<ICollectionPotential> getArrayListPotentialCollection() {
		return arrayListPotentialCollection;
	}


	
	public void setArrayListPotentialCollection(
			ArrayList<ICollectionPotential> hashMapPotentialCollection) {
		this.arrayListPotentialCollection = hashMapPotentialCollection;
	}

	
	
}
