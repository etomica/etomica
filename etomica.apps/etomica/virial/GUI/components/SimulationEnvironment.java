package etomica.virial.GUI.components;

public class SimulationEnvironment {
	
	public int noOfSteps;
	public double temperature;
	
	

	private SimulationEnvironment() {
		this.noOfSteps = 100000;
		this.temperature = 450.0;
		
	}
	 
    /**
    * SingletonHolder is loaded on the first execution of Singleton.getInstance() 
    * or the first access to SingletonHolder.INSTANCE, not before.
    */
    private static class SimulationEnvironmentHolder { 
            public static final SimulationEnvironment instance = new SimulationEnvironment();
    }

    public static SimulationEnvironment getInstance() {
            return SimulationEnvironmentHolder.instance;
    }

	public int getNoOfSteps() {
		return noOfSteps;
	}

	public void setNoOfSteps(int noOfSteps) {
		this.noOfSteps = noOfSteps;
	}

	

	public double getTemperature() {
		return temperature;
	}

	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}
    

}
