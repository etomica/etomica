package etomica.virial.GUI.components;

public interface ParameterMapping {
	
	
	
	public int getParameterCount();
	
	public void setParameter(String Parameter,String ParameterValue);
	
	public String getDescription(String Parameter);
	
	public Double getDoubleDefaultParameters(String Parameter);
}
