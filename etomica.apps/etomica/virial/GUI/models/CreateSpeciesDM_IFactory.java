package etomica.virial.GUI.models;

import etomica.api.ISpecies;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.virial.SpeciesFactory;

public interface CreateSpeciesDM_IFactory {
	
	
	public Space getSpace();
	
	public int getParameterCount();
	
	public void setParameter(String Parameter,String ParameterValue);
	
	public String getDescription(String Parameter);
	
	public Double getDoubleDefaultParameters(String Parameter);
	
	public String getMoleculeDisplayName();

	// Added on April 27, 2011
	public String[] getParametersArray();
	
	
	//Added June 20, 2011
	public String[][] getParamAndValues();

	public SpeciesFactory createSpeciesFactory();
	
	public ISpecies createSpecies();

	public Object clone();
	
	public String getCustomName();
	
	public String[] getPotentialSites();
	
	public String getPotentialSiteAtIndex(int index);
	
	@SuppressWarnings("rawtypes")
	public Class getPotential();
	
	public boolean hasElectrostaticInteraction();
	
	public String getNonBondedInteractionModel();
}
