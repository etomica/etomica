package etomica.virial.GUI.components;

import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerHardSphere;

public class SimulationEnvironmentObject {
	

	private double temperature;
	private int noOfSteps;
	private double sigmaHSRefA;
	private double sigmaHSRefB;
	private int Alkane1Spheres;
	private int Alkane2Spheres;
	
	private int SpeciesA;
	private int SpeciesB;
	private int nPoints;

	public double SigmaHSRef;
	public boolean doWiggle = false;
	
	
	public SimulationEnvironmentObject(double Temperature,int NoOfSteps,MixtureBuilderSpeciesFactory potential1, MixtureBuilderSpeciesFactory potential2){
		temperature = Temperature;
		noOfSteps = NoOfSteps;
		sigmaHSRefA = getSigmaHSRefSpecies(potential1);
		if(potential2 != null){
			sigmaHSRefB = getSigmaHSRefSpecies(potential2);}
		else{
			sigmaHSRefB = 0;
		}
		Alkane1Spheres = getAlkaneSpheresSpecies(potential1);
		if(potential2 != null){
			Alkane2Spheres = getAlkaneSpheresSpecies(potential2);
		}
		else{
			Alkane2Spheres = 0;
		}
		
		
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

	public double getSigmaHSRefA() {
		return sigmaHSRefA;
	}

	public void setSigmaHSRefA(double sigmaHSRef) {
		this.sigmaHSRefA = sigmaHSRef;
	}
	
	public double getSigmaHSRefB() {
		return sigmaHSRefB;
	}

	public void setSigmaHSRefB(double sigmaHSRef) {
		this.sigmaHSRefB = sigmaHSRef;
	}
	
	public int getAlkane1Spheres() {
		return Alkane1Spheres;
	}

	public void setAlkane1Spheres(int alkane1Spheres) {
		Alkane1Spheres = alkane1Spheres;
	}

	public int getAlkane2Spheres() {
		return Alkane2Spheres;
	}

	public void setAlkane2Spheres(int alkane2Spheres) {
		Alkane2Spheres = alkane2Spheres;
	}

	public double getSigmaHSRefSpecies(MixtureBuilderSpeciesFactory Potential){
		
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
	

	public int getSpeciesA() {
		return SpeciesA;
	}


	public void setSpeciesA(int speciesA) {
		SpeciesA = speciesA;
	}


	public int getSpeciesB() {
		return SpeciesB;
	}


	public void setSpeciesB(int speciesB) {
		SpeciesB = speciesB;
	}
	

	public int getnPoints() {
		return nPoints;
	}


	public void setnPoints() {
		this.nPoints = this.SpeciesA + this.SpeciesB;
	}

	
	private int getAlkaneSpheresSpecies(MixtureBuilderSpeciesFactory potential1) {
		// TODO Auto-generated method stub
		int AlkaneSphere = 0;
		if(potential1.getClass().getName().contains("Alkane")){
			String tempAlkaneSpheres = Double.toString(potential1.getDoubleDefaultParameters("NUMBER"));
			String[] AlkaneSpheres = tempAlkaneSpheres.split("\\.");
			AlkaneSphere = Integer.parseInt(AlkaneSpheres[0]);
		
		}
		return AlkaneSphere;
	}
	
	public void calculateSystemSigmaHSRef(MixtureBuilderSpeciesFactory potential1, MixtureBuilderSpeciesFactory potential2){
		
		if(potential2 != null){
			if(potential1.getClass().getName().contains("Alkane") 
				&& potential2.getClass().getName().contains("Alkane") 
				&& potential1.getDoubleDefaultParameters("NUMBER") >= 2
				&& potential2.getDoubleDefaultParameters("NUMBER") >= 2){
			
				this.SigmaHSRef=  0.5 * (  potential1.getDoubleDefaultParameters("SIGMACH3") + (0.5 * potential1.getDoubleDefaultParameters("NUMBER")) +
										potential2.getDoubleDefaultParameters("SIGMACH3") + (0.5 * potential2.getDoubleDefaultParameters("NUMBER")) );
				this.doWiggle = true;
			}
		
			if(potential1.getClass().getName().contains("Alkane") 
				&& potential2.getClass().getName().contains("Alkane") && potential1.getDoubleDefaultParameters("NUMBER") == 1){
			
				this.SigmaHSRef= 0.7 * ( potential1.getDoubleDefaultParameters("SIGMACH3")  + potential2.getDoubleDefaultParameters("SIGMACH3"));
				this.doWiggle = true;
			}
		
			if(potential1.getClass().getName().contains("Alkane") 
				&& potential2.getClass().getName().contains("Alkane") && potential2.getDoubleDefaultParameters("NUMBER") == 1){
			
				this.SigmaHSRef= 0.7 * ( potential1.getDoubleDefaultParameters("SIGMACH3")  + potential2.getDoubleDefaultParameters("SIGMACH3"));
				this.doWiggle = true;
			}
		
			if(potential1.getClass().getName().contains("Alkane") && potential2.getClass().getName().contains("CO2") || 
					potential2 .getClass().getName().contains("Alkane") && potential1.getClass().getName().contains("CO2")){
				this.SigmaHSRef = 7;
				if(Alkane1Spheres > 2 || Alkane2Spheres > 2){
					this.doWiggle = true;
				}
				
			}
			if(potential1.getClass().getName().contains("CO2") && potential2.getClass().getName().contains("Naphthalene") ||
				potential2.getClass().getName().contains("CO2") && potential2.getClass().getName().contains("Naphthalene")){
			
				this.SigmaHSRef = 7;
				
			}
			if(potential1.getClass().getName().contains("CO2") && potential2.getClass().getName().contains("Ethane") ||
				potential2.getClass().getName().contains("CO2") && potential2.getClass().getName().contains("Ethane")){
			
				this.SigmaHSRef = 7;
				this.doWiggle = true;
			}
			if(potential1.getClass().getName().contains("CO2") && potential2.getClass().getName().contains("Water") ||
					potential2.getClass().getName().contains("Water") && potential2.getClass().getName().contains("CO2")){
				
				this.SigmaHSRef = (this.SpeciesA * (potential1.getDoubleDefaultParameters("SIGMAHSREF"))) + (this.SpeciesB * (potential2.getDoubleDefaultParameters("SIGMAHSREF")))/ (SpeciesA + SpeciesB);
				
			}
		}
		else{
			if(potential1.getClass().getName().contains("Alkane") && potential1.getDoubleDefaultParameters("NUMBER") >= 2){
				this.SigmaHSRef = potential1.getDoubleDefaultParameters("SIGMACH3") + (0.5 * potential1.getDoubleDefaultParameters("NUMBER"));
				this.doWiggle = true;
			}else if(potential1.getClass().getName().contains("Alkane") && potential1.getDoubleDefaultParameters("NUMBER") == 1){
				this.SigmaHSRef = potential1.getDoubleDefaultParameters("SIGMACH3") + (0.5 * potential1.getDoubleDefaultParameters("NUMBER"));
			}
			else{
				this.SigmaHSRef = this.getSigmaHSRefA();
			}

		}
	}
	

	public boolean isDoWiggle() {
		return doWiggle;
	}


	public double getSigmaHSRef() {
		return SigmaHSRef;
	}



}
