package etomica.virial.GUI.components;

public enum ParametersDouble {
	
	SIGMA(1.0, "Potential Well-Depth"),
	EPSILON(1.0,"Distance at which Interatomic potential between particles is zero"),
	MOMENT(1.0,"Moment"),
	BONDLENGTH(1.0,"Bond Length"),
	MOMENTSQUARED(1.0,"Moment square");

	private Double DefaultValues;
	private String Description;
	
	ParametersDouble(Double Value, String description){
		this.DefaultValues = Value;
		this.Description = description;
	}
	
	public double DefaultValue(){
		return DefaultValues;
	}
	
	public String Description(){
		return Description;
	}
}
