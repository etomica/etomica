/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerHardSphere;

public class ModelSimulationEnvironment {
	

	private double temperature;
	private int noOfSteps;
	
	
	private double[] sigmaHSRefSpecies;
	private int[] lengthAlkaneChain;
	private int[] countSpecies;

	
	
	private int nPoints;
	
	private double systemSigmaHSRef;

	public boolean doWiggle = false;
	public boolean stopSimulation = false;
	
	

	public ModelSimulationEnvironment(){
		reset();
	}
	
	public void reset() {
		// TODO Auto-generated method stub
		
		sigmaHSRefSpecies = new double[10];
		lengthAlkaneChain = new int[10];
		doWiggle = false;
		stopSimulation = false;
	}

	public ModelSimulationEnvironment(double Temperature,int NoOfSteps,ModelSelectedSpecies modelSelectedSpecies ){
		temperature = Temperature;
		noOfSteps = NoOfSteps;
		sigmaHSRefSpecies = new double[10];
		lengthAlkaneChain = new int[10];
		
		doWiggle = false;
		stopSimulation = false;
		
		for (int j=0;j<10;j++){
			if(modelSelectedSpecies.getSpeciesDataModel(j) != null){
				setSigmaHSRef(getSigmaHSRefSpecies(modelSelectedSpecies.getSpeciesDataModel(j)),j);
				setLengthAlkaneChain(getAlkaneSpheresSpecies(modelSelectedSpecies.getSpeciesDataModel(j)),j);
			}
		}

	}

	public ModelSimulationEnvironment(double Temperature,int NoOfSteps){
		temperature = Temperature;
		noOfSteps = NoOfSteps;
		sigmaHSRefSpecies = new double[10];
		lengthAlkaneChain = new int[10];
		countSpecies = new int[10];
		doWiggle = false;
		stopSimulation = false;
	}
	
	
	public double getTemperature() {
		return temperature;
	}

	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	public int getNoOfSteps() {
		return noOfSteps;
	}

	public void setNoOfSteps(int noOfSteps) {
		this.noOfSteps = noOfSteps;
	}

	public double getSigmaHSRef(int index) {
		return sigmaHSRefSpecies[index];
	}

	public void setSigmaHSRef(double sigmaHSRef,int index) {
		this.sigmaHSRefSpecies[index] = sigmaHSRef;
	}


	public int getLengthAlkaneChain(int index) {
		return lengthAlkaneChain[index];
	}

	public void setLengthAlkaneChain(int lengthAlkaneChain1,int index) {
		this.lengthAlkaneChain[index] = lengthAlkaneChain1;
	}



	public double getSigmaHSRefSpecies(IMolecularModel_SpeciesFactory Potential){
		
		double SigmaHSRefSpecies = 0.0;
		if(Potential.getClass().getName().contains("Alkane")){
			if(Potential.getDoubleDefaultParameters("NUMBER") >= 2){
				SigmaHSRefSpecies = Potential.getDoubleDefaultParameters("SIGMACH3") + (0.5 * Potential.getDoubleDefaultParameters("NUMBER"));
			}
			else{
				SigmaHSRefSpecies = Potential.getDoubleDefaultParameters("SIGMACH3") + (0.5 * Potential.getDoubleDefaultParameters("NUMBER"));
			}
		}else{
			SigmaHSRefSpecies = Potential.getDoubleDefaultParameters("SIGMAHSREF");
		}

		return SigmaHSRefSpecies;
		
	}
	

	public int getCountSpecies(int index) {
		return countSpecies[index];
	}


	public void setCountSpecies(int species1,int index) {
		countSpecies[index] = species1;
	}


	
	public int getnPoints() {
		return nPoints;
	}


	public void setnPoints(int[] countOfSpecies) {
		for(int i=0;i<10;i++){
			nPoints = nPoints + countOfSpecies[i];
		}
		//this.nPoints = this.countSpecies + this.countSpecies2 + this.countSpecies3;
	}

	
	public int getAlkaneSpheresSpecies(IMolecularModel_SpeciesFactory potential1) {
		// TODO Auto-generated method stub
		int AlkaneSphere = 0;
		if(potential1.getClass().getName().contains("Alkane")){
			String tempAlkaneSpheres = Double.toString(potential1.getDoubleDefaultParameters("NUMBER"));
			String[] AlkaneSpheres = tempAlkaneSpheres.split("\\.");
			AlkaneSphere = Integer.parseInt(AlkaneSpheres[0]);
		
		}
		return AlkaneSphere;
	}
	
	public void calculateSystemSigmaHSRef(ModelSelectedSpecies modelSelectedSpecies){
		/*IMolecularModel_SpeciesFactory potential1 = modelSelectedSpecies.getSpeciesDataModel();
		IMolecularModel_SpeciesFactory potential2 = null;
		IMolecularModel_SpeciesFactory potential3 = null;
		if(modelSelectedSpecies.getSpecies2DataModel() != null){
			potential2 = modelSelectedSpecies.getSpecies2DataModel();
			if(modelSelectedSpecies.getSpecies3DataModel() != null){
				potential3 = modelSelectedSpecies.getSpecies3DataModel();
				this.SigmaHSRef23 = calculateSigmaHSRef(potential2,potential3);
				this.SigmaHSRef13 = calculateSigmaHSRef(potential1, potential3);
			}else{
				this.SigmaHSRef23 = 0.0;
				this.SigmaHSRef13 = 0.0;
			}
			this.SigmaHSRef12 = calculateSigmaHSRef(potential1,potential2);
			
		}else{
			if(potential1.getClass().getName().contains("Alkane") && potential1.getDoubleDefaultParameters("NUMBER") >= 2){
				this.SigmaHSRef12 = potential1.getDoubleDefaultParameters("SIGMACH3") + (0.5 * potential1.getDoubleDefaultParameters("NUMBER"));
				this.doWiggle = true;
			}else if(potential1.getClass().getName().contains("Alkane") && potential1.getDoubleDefaultParameters("NUMBER") == 1){
				this.SigmaHSRef12 = potential1.getDoubleDefaultParameters("SIGMACH3") + (0.5 * potential1.getDoubleDefaultParameters("NUMBER"));
			}
			else{
				this.SigmaHSRef12 = this.getSigmaHSRef();
			}
			this.SigmaHSRef23 = 0.0;
			this.SigmaHSRef13 = 0.0;
		}*/
		
		
	}

	

	public double calculateSystemSigmaHSRef(IMolecularModel_SpeciesFactory[] speciesDataModel,int[] speciesMoleculeCount){
		
		double SigmaHSRef = 0.0;
		
		if(speciesDataModel.length == 2){
			if(speciesDataModel[0].getClass().getName().contains("Alkane") 
				&& speciesDataModel[1].getClass().getName().contains("Alkane") 
				&& speciesDataModel[0].getDoubleDefaultParameters("NUMBER") >= 2
				&& speciesDataModel[1].getDoubleDefaultParameters("NUMBER") >= 2){
			
				SigmaHSRef = (0.5 * (  speciesDataModel[0].getDoubleDefaultParameters("SIGMACH3") + (0.5 * speciesDataModel[0].getDoubleDefaultParameters("NUMBER")) +
								speciesDataModel[1].getDoubleDefaultParameters("SIGMACH3") + (0.5 * speciesDataModel[1].getDoubleDefaultParameters("NUMBER")) ));
				
			}
		
			if(speciesDataModel[0].getClass().getName().contains("Alkane") 
				&& speciesDataModel[1].getClass().getName().contains("Alkane") && speciesDataModel[0].getDoubleDefaultParameters("NUMBER") == 1){
			
				SigmaHSRef = (0.7 * ( speciesDataModel[0].getDoubleDefaultParameters("SIGMACH3")  + speciesDataModel[1].getDoubleDefaultParameters("SIGMACH3")));
				
			}
		
			if(speciesDataModel[0].getClass().getName().contains("Alkane") 
				&& speciesDataModel[1].getClass().getName().contains("Alkane") && speciesDataModel[1].getDoubleDefaultParameters("NUMBER") == 1){
			
				SigmaHSRef = (0.7 * ( speciesDataModel[0].getDoubleDefaultParameters("SIGMACH3")  + speciesDataModel[1].getDoubleDefaultParameters("SIGMACH3")));
				
			}
		
			if(speciesDataModel[0].getClass().getName().contains("Alkane") && speciesDataModel[1].getClass().getName().contains("CO2") || 
					speciesDataModel[1] .getClass().getName().contains("Alkane") && speciesDataModel[0].getClass().getName().contains("CO2")){
				SigmaHSRef = 7;
				/*if(lengthAlkaneChain1 > 2 || lengthAlkaneChain2 > 2){
					this.doWiggle = true;
				}*/
				
			}
			if(speciesDataModel[0].getClass().getName().contains("CO2") && speciesDataModel[1].getClass().getName().contains("Naphthalene") ||
					speciesDataModel[1].getClass().getName().contains("CO2") && speciesDataModel[0].getClass().getName().contains("Naphthalene")){
			
				SigmaHSRef = 7;
				
			}
			if(speciesDataModel[0].getClass().getName().contains("CO2") && speciesDataModel[1].getClass().getName().contains("Ethane") ||
					speciesDataModel[1].getClass().getName().contains("CO2") && speciesDataModel[0].getClass().getName().contains("Ethane")){
			
				SigmaHSRef = 7;
				
			}
			if(speciesDataModel[0].getClass().getName().contains("CO2") && speciesDataModel[1].getClass().getName().contains("H2O") ||
					speciesDataModel[1].getClass().getName().contains("H2O") && speciesDataModel[0].getClass().getName().contains("CO2")){
				
				//SigmaHSRef = (this.countSpecies * (potential1.getDoubleDefaultParameters("SIGMAHSREF"))) + (this.countSpecies2 * (potential2.getDoubleDefaultParameters("SIGMAHSREF")))/ (countSpecies + countSpecies2);
				SigmaHSRef = 3.2;
			}
			return SigmaHSRef;
		}else if(speciesDataModel.length == 1){
			SigmaHSRef = getSigmaHSRefSpecies(speciesDataModel[0]);
		}else {
			double sumSigmaHSRef = 0;
			int totalMoleculeCount = 0;
			for(int i=0;i<10;i++){
				if(speciesDataModel[i] != null){
				sumSigmaHSRef = sumSigmaHSRef + (getSigmaHSRefSpecies(speciesDataModel[i]) * speciesMoleculeCount[i]);
				totalMoleculeCount = totalMoleculeCount + speciesMoleculeCount[i];
				}
			}
			SigmaHSRef = sumSigmaHSRef/totalMoleculeCount;
		}
		return SigmaHSRef;
	}
	
	public boolean isDoWiggle() {
		return doWiggle;
	}


	public double getSystemSigmaHSRef() {
		return systemSigmaHSRef;
	}

	public void setSystemSigmaHSRef(double sigmaHSRef) {
		this.systemSigmaHSRef = sigmaHSRef;
	}
	
	public boolean isStopSimulation() {
		return stopSimulation;
	}

	

	public void setStopSimulation(boolean b) {
		// TODO Auto-generated method stub
		this.stopSimulation = b;
	}

}
