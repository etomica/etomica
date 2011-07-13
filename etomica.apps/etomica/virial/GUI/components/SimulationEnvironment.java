package etomica.virial.GUI.components;

public class SimulationEnvironment {
	
	private double temperature;
	private int noOfSteps;
	private double sigmaHSRef;
	private int Alkane1Spheres;
	private int Alkane2Spheres;
	
	
	public SimulationEnvironment(double Temperature,int NoOfSteps,double SigmaHSRef){
		temperature = Temperature;
		noOfSteps = NoOfSteps;
		sigmaHSRef = SigmaHSRef;
		Alkane1Spheres = 0;
		Alkane2Spheres = 0;
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


}
