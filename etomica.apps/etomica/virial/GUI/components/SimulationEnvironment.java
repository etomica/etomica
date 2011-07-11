package etomica.virial.GUI.components;

public class SimulationEnvironment {
	
	private double temperature;
	private int noOfSteps;
	private double sigmaHSRef;
	
	public SimulationEnvironment(double Temperature,int NoOfSteps,double SigmaHSRef){
		temperature = Temperature;
		noOfSteps = NoOfSteps;
		sigmaHSRef = SigmaHSRef;
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

	public double getSigmaHSRef() {
		return sigmaHSRef;
	}

	public void setSigmaHSRef(double sigmaHSRef) {
		this.sigmaHSRef = sigmaHSRef;
	}
	
	

}
