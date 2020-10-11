/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.virial.SpeciesFactory;

public interface IMolecularModel_SpeciesFactory {
	
	
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

	public ISpecies createSpecies();

	public Object clone();
	
	public String getMolecularModelDisplayName();
	
	public String[] getPotentialSites();
	
	public String getPotentialSiteAtIndex(int index);
	
	@SuppressWarnings("rawtypes")
	public Class getPotential();
	
	public boolean hasElectrostaticInteraction();
	
	public String getNonBondedInteractionModel();
}
