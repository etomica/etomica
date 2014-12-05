/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

public class ModelTemperatureAndSteps {
	
	public int noOfSteps;
	public double temperature;
	
	

	private ModelTemperatureAndSteps() {
		this.noOfSteps = 100000;
		this.temperature = 450.0;
		
	}
	 
    /**
    * SingletonHolder is loaded on the first execution of Singleton.getInstance() 
    * or the first access to SingletonHolder.INSTANCE, not before.
    */
    private static class SimulationEnvironmentHolder { 
            public static final ModelTemperatureAndSteps instance = new ModelTemperatureAndSteps();
    }

    public static ModelTemperatureAndSteps getInstance() {
            return SimulationEnvironmentHolder.instance;
    }

	public int getNoOfSteps() {
		return noOfSteps;
	}

	public void setNoOfSteps(int noOfSteps) {
		this.noOfSteps = noOfSteps;
	}

	public void reset(){
		this.noOfSteps = 100000;
		this.temperature = 450.0;
	}
	

	public double getTemperature() {
		return temperature;
	}

	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}
    

}
