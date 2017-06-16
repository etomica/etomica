/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import etomica.species.ISpecies;

public class ModelSelectedSpecies {
	
	private IMolecularModel_SpeciesFactory[] speciesDataModel;
	//private IMolecularModel_SpeciesFactory species2DataModel;
	//private IMolecularModel_SpeciesFactory species3DataModel;
	
	private ISpecies[] species;
	//private ISpecies species2;
	//private ISpecies species3;
	
	private int[] speciesMoleculeCount;
	//private int species2MoleculeCount;
	//private int species3MoleculeCount;

	private int nthSpeciesAdded;
	
	private int[][] speciesInteractionIndex;
	
	public ModelSelectedSpecies(){
		reset();
	}
	
	public IMolecularModel_SpeciesFactory getSpeciesDataModel(int index) {
		return speciesDataModel[index];
	}
	public IMolecularModel_SpeciesFactory[] getSpeciesDataModel() {
		return speciesDataModel;
	}

	public void setSpeciesDataModel(IMolecularModel_SpeciesFactory species1DataModel,int index) {
		this.speciesDataModel[index] = species1DataModel;
	}
/*
	public IMolecularModel_SpeciesFactory getSpecies2DataModel() {
		return species2DataModel;
	}

	public void setSpecies2DataModel(IMolecularModel_SpeciesFactory species2DataModel) {
		this.species2DataModel = species2DataModel;
	}
	
	public IMolecularModel_SpeciesFactory getSpecies3DataModel() {
		return species3DataModel;
	}

	public void setSpecies3DataModel(IMolecularModel_SpeciesFactory species3DataModel) {
		this.species3DataModel = species3DataModel;
	}
*/
	
	
	public int getSpeciesMoleculeCount(int index) {
		return speciesMoleculeCount[index];
	}
	
	public int[] getSpeciesMoleculeCount() {
		return speciesMoleculeCount;
	}

	public void setSpeciesMoleculeCount(int species1MoleculeCount,int index) {
		this.speciesMoleculeCount[index] = species1MoleculeCount;
	}

	/*
	public int getSpecies2MoleculeCount() {
		return species2MoleculeCount;
	}

	public void setSpecies2MoleculeCount(int species2MoleculeCount) {
		this.species2MoleculeCount = species2MoleculeCount;
	}

	public int getSpecies3MoleculeCount() {
		return species3MoleculeCount;
	}

	public void setSpecies3MoleculeCount(int species3MoleculeCount) {
		this.species3MoleculeCount = species3MoleculeCount;
	}
	
	*/
	//Methods to increment Species 1 and 2 Count
	public void addSpeciesMolecule(int index){
		this.speciesMoleculeCount[index] = speciesMoleculeCount[index] + 1;
	}
/*
	public void addSpecies2Molecule(){
		this.species2MoleculeCount = this.species2MoleculeCount + 1;
	}
	
	public void addSpecies3Molecule(){
		this.species3MoleculeCount = this.species3MoleculeCount + 1;
	}
*/	
	public void removeSpeciesMolecule(int index){
		this.speciesMoleculeCount[index] = speciesMoleculeCount[index] - 1;
	}
/*
	public void removeSpecies2Molecule(){
		this.species2MoleculeCount = this.species2MoleculeCount - 1;
	}
	
	public void removeSpecies3Molecule(){
		this.species3MoleculeCount = this.species3MoleculeCount - 1;
	}
*/	
	public ISpecies getSpecies(int index) {
		return species[index];
	}

	public void setSpecies(ISpecies species1,int index) {
		this.species[index] = species1;
	}
/*
	public ISpecies getSpecies2() {
		return species2;
	}

	public void setSpecies2(ISpecies species2) {
		this.species2 = species2;
	}

	public ISpecies getSpecies3() {
		return species3;
	}

	public void setSpecies3(ISpecies species3) {
		this.species3 = species3;
	}
	
	*/
	
	
	/** Reset to initial value. */
    public void reset() {
    	speciesDataModel = new IMolecularModel_SpeciesFactory[10];
    	speciesMoleculeCount = new int[10];
    	species = new ISpecies[10];
    	nthSpeciesAdded = 0;
    	//species2DataModel = null;
    	//species3DataModel = null;
    	//speciesMoleculeCount = 0;
    	//species2MoleculeCount = 0;
    	//species3MoleculeCount = 0;
    	
    }
    
    
    
	public int[][] getSpeciesInteractionIndex() {
		return speciesInteractionIndex;
	}

	public void setSpeciesInteractionIndex(int[][] speciesInteractionIndex) {
		this.speciesInteractionIndex = speciesInteractionIndex;
	}

	public int getNthSpeciesAdded() {
		return nthSpeciesAdded;
	}

	public void setNthSpeciesAdded(int nthSpeciesAdded) {
		this.nthSpeciesAdded = nthSpeciesAdded;
	}
}
