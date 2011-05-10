package etomica.virial.GUI.models;

import etomica.units.Kelvin;

public enum CopyOfParametersDouble {
	
	
	//Array Order so far
	/*
	 * LJ
	 * CO2
	 * C in EPM2
	 * O in EPM2
	 * CO in EPM2
	 * C in Trappe
	 * O in Trappe
	 * CO in Trappe
	 * 
	 * 
	 * 
	 * 
	 * 
	 */
	SIGMA(new Double[]{1.0,3.0354}, "Potential Well-Depth"),
	EPSILON(new Double[]{1.0,(Kelvin.UNIT.toSim(125.317))},"Distance at which Interatomic potential between particles is zero"),
	MOMENT(new Double[]{1.0,(3.0255*Kelvin.UNIT.toSim(125.317)*Math.pow(3.0354,5))},"Moment"),
	BONDL(new Double[]{1.0,(0.699*3.0354)},"Bond Length"),
	MOMENTSQR(new Double[]{1.0,0.0},"Moment square");

	private Double[] DefaultValues;
	private String Description;
	
	CopyOfParametersDouble(Double[] Value, String description){
		this.DefaultValues = Value;
		this.Description = description;
	}
	
	public Double DefaultValue(int i){
		return DefaultValues[i];
	}
	
	public String Description(){
		return Description;
	}
}
