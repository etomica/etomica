package etomica.virial.GUI.components;

import etomica.virial.SpeciesFactory;

public interface ParameterMapping {
	
	
	
	public int getParameterCount();
	
	public void setParameter(String Parameter,String ParameterValue);
	
	public String getDescription(String Parameter);
	
	public Double getDoubleDefaultParameters(String Parameter);

	// Added on April 27, 2011
	public String[] getParametersArray();

	public SpeciesFactory createSpeciesFactory();

	public Object clone();
	
}
