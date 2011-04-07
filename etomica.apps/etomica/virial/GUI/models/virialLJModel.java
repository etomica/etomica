package etomica.virial.GUI.models;


public class virialLJModel {
	
	/*
	Initial Values to be set when User runs the Interface for the 
	first time/reset or default values
	*/
	private static final int INITIAL_VALUE_nPoints= 3;
	private static final double INITIAL_VALUE_temperature= 2.0;
	private static final long INITIAL_VALUE_steps= 100000;
	private static final double INITIAL_VALUE_sigmaHSRef = 1.5;
	
	private int nPoints;
	private double temperature;  
	private long steps;
	private double sigmaHSRef;
	
	//Constructor
	public virialLJModel() {
        reset();
    }

	/** Reset to initial value. */
    public void reset() {
    	nPoints = INITIAL_VALUE_nPoints;
    	temperature = INITIAL_VALUE_temperature;
    	steps = INITIAL_VALUE_steps;
    	sigmaHSRef = INITIAL_VALUE_sigmaHSRef;
    }
    
    //Methods to set the value of parameter variables defined as part of model
    
    public void setValueNPoints(String value) {
    	nPoints = Integer.parseInt(value);
    }
    
    public void setValueTemperature(String value) {
    	temperature = Double.parseDouble(value);
    }

    public void setValueSteps(String value) {
    	steps = Long.parseLong(value);
    }

    public void setValueSigmaHSRef(String value) {
    	sigmaHSRef = Double.parseDouble(value);
    	
    }
    
  //Methods to get the value of parameter variables defined as part of model
    
    public String getValueNPoints() {
    	return Integer.toString(nPoints);
    }
    
    public String getValueTemperature() {
    	return Double.toString(temperature);
    }

    public String getValueSteps() {
    	return Long.toString(steps);
    }

    public String getValueSigmaHSRef() {
    	return Double.toString(sigmaHSRef);
    	
    }
}
