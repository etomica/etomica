package etomica.virial.GUI.models;

import etomica.units.Kelvin;

public enum PotentialParamDM_Description {
	
	
	//Array Order so far
	/*
	 * 
	 * 
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
	 */
	SIGMA(0.0,"Potential Well-Depth"),
	
	EPSILON(0.0, "Distance at which Interatomic potential between particles is zero"),
	
	MOMENT(0.0, "Moment"),
	
	BONDL(0.0, "Bond Length"),
	
	MOMENTSQR(0.0,"Moment square"),
	
	CHARGE(0.0,"Quadrapole Charge"),
	
	NominalbondL(0.0,"Fixed bond length between neighboring pseudo-atoms"),
	
	theta(0.0,"Equilibrium bond angle"),
	
	forceconstant(0.0,"Force Constant(k0/kB)"),
	
	TEMPERATURE(0.0,"Temperature of Simulation"),
	
	STEPS(0.0,"No of Steps"),
	
	SIGMAHSREF(0.0,"Hard Sphere Reference"),
	
	NUMBER(0.0,"No of Spheres");
	

	private Double DefaultValues;
	private String Description;
	
	PotentialParamDM_Description(Double Value, String description){
		this.DefaultValues = Value;
		this.Description = description;
	}
	
	public Double DefaultValue(){
		return DefaultValues;
	}
	
	public String Description(){
		return Description;
	}
}
