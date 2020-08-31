/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import etomica.potential.P2LJQ;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.Species;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;
import etomica.virial.GUI.components.SimpleElementForSimilarSpecies;
import etomica.virial.SpeciesFactory;

public class MolecularModelLJQ_SpeciesLJ implements IMolecularModel_SpeciesFactory,Cloneable{
	
	private static String MoleculeDisplayName = "LJ with Quad";
	private Space space;
	private double sigma[];
	private double epsilon[];
	private double momentSquare[];
	private SimpleElementForSimilarSpecies simpleElement;
	
	private double sigmaHSRef;
	
	private int id;
	private static int numberOfInstances = 0;
	
	
	
	private String CustomClassName = "Spherical-2-Body-With-Quad";
	
	private String[] ComponentParameters  = {"SIGMA","EPSILON","MOMENTSQR"};
	
	private String[] SharedComponentParameters =null;
	
	private String[] PotentialSites = {"LJ"};
	
	private String[][] ComponentValues = {{"1.0","1.0","1.0"}};
	
	private String[] SharedComponentValues = null;
	
	private String[][] ParamAndValues; 
	
	private String[] SimEnvParameters = {"SIGMAHSREF"};
	
	private String[] SimEnvValues = {"1.5"};
	
	
	
	//Constructors for different Instantiations
	
	public MolecularModelLJQ_SpeciesLJ(){
		space = Space3D.getInstance();
		sigma = new double[PotentialSites.length];
		epsilon = new double[PotentialSites.length];
		momentSquare = new double[PotentialSites.length];
		ParamAndValues = setParameterValues();
		//p2LJQ = new P2LJQ[PotentialSites.length];
		id = ++numberOfInstances;
	}
	
private String[][] setParameterValues() {
		
		int NoOfParam = ComponentParameters.length;
		
		int NoOfSites = PotentialSites.length;
		int totalNoOfParam = NoOfParam*NoOfSites;
		String[][] ReturnArray = new String[totalNoOfParam][2];
		int index = 0;
		for(int i=0;i<NoOfSites;i++){
			for(int j=0;j<NoOfParam;j++){
				if(ComponentParameters[j]=="SIGMA"){
					setSigma(Double.parseDouble(ComponentValues[i][j]),i);
				}
				if(ComponentParameters[j]=="EPSILON"){
					setEpsilon(Double.parseDouble(ComponentValues[i][j]),i);
				}
				if(ComponentParameters[j]=="MOMENTSQR"){
					setMomentSquare(Double.parseDouble(ComponentValues[i][j]),i);
				}
				ReturnArray[index][0] = ComponentParameters[j]+PotentialSites[i];
				ReturnArray[index][1] = ComponentValues[i][j];
				index++;
				
			}
		}
		
		int NoOfSimEnvParam = 1;
		for(int l = 0;l<NoOfSimEnvParam;l++){
			
			
			if(SimEnvParameters[l]=="SIGMAHSREF"){
				setSigmaHSRef(Double.parseDouble(SimEnvValues[l]));
			}
		}
		return ReturnArray;
	}
	
	

	public int getId() {
		return id;
	}
	

	public double getSigmaHSRef() {
		return sigmaHSRef;
	}

	public void setSigmaHSRef(double sigmaHSRef) {
		this.sigmaHSRef = sigmaHSRef;
	}

		public String[][] getComponentValues() {
		return ComponentValues;
	}

	public void setComponentValues(String[][] componentValues) {
		ComponentValues = componentValues;
	}

	public String[][] getParamAndValues() {
		return ParamAndValues;
	}
	
	
	 public Object clone(){
		 try{
			 MolecularModelLJQ_SpeciesLJ cloned = (MolecularModelLJQ_SpeciesLJ)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
	 }
	
	//Creates the LJAtom Species
	public ISpecies createSpecies(){
		simpleElement = SimpleElementForSimilarSpecies.getInstance();
		SpeciesGeneral speciesLJQ = SpeciesSpheresRotating.create(space, simpleElement.getOelement());
        return speciesLJQ;
	}
	
	//Testing Class
	public static void main(String[] args){
		MolecularModelLJQ_SpeciesLJ lj = new MolecularModelLJQ_SpeciesLJ();
		
		System.out.println(lj.getDescription("EPSILONLJ"));
		System.out.println(lj.getDoubleDefaultParameters("EPSILONLJ"));
		for(int j=0;j<lj.ComponentParameters.length*lj.PotentialSites.length;j++){
			
			System.out.println(lj.ParamAndValues[j][0]+"\n");
			System.out.println(lj.ParamAndValues[j][1]+"\n");
	}
		//lj.setParameter("epsilon", "1.5");
		//System.out.println(lj.getEpsilon(i));
	}

	public double getSigma(int index) {
		return sigma[index];
	}

	public void setSigma(double sigma,int index) {
		this.sigma[index] = sigma;
	}
	public double getEpsilon(int index) {
		return epsilon[index];
	}

	public void setEpsilon(double epsilon,int index) {
		this.epsilon[index] = epsilon;
	}

	public double getMomentSquare(int index) {
		return momentSquare[index];
	}

	public void setMomentSquare(double momentSquare,int index) {
		this.momentSquare[index] = momentSquare;
	}
	
	

	public int getParameterCount() {
		return 3;
	}

	
	public void setParameter(String Parameter, String ParameterValue) {
		// TODO Auto-generated method stub
		
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMA.toString()+PotentialSites[i])){
				setSigma(Double.parseDouble(ParameterValue),i); 
				ParamAndValues[0][1] = Double.toString(getSigma(i));
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.EPSILON.toString()+PotentialSites[i])){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
				ParamAndValues[1][1] = Double.toString(getEpsilon(i));
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.MOMENTSQR.toString()+PotentialSites[i])){
				setMomentSquare(Double.parseDouble(ParameterValue),i); 
				ParamAndValues[2][1] = Double.toString(getMomentSquare(i));
			}
			
		}

		
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMAHSREF.toString())){
			setSigmaHSRef(Double.parseDouble(ParameterValue)); 
		}
		
	}


	public String getDescription(String Parameter) {
		String Description = null;
		for(int i = 0;i <PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMA.toString()+PotentialSites[i])){
				Description = EnumPotentialParamDescription.SIGMA.description();
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.EPSILON.toString()+PotentialSites[i])){
				Description = EnumPotentialParamDescription.EPSILON.description();
			}
		
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.MOMENTSQR.toString()+PotentialSites[i])){
				Description = EnumPotentialParamDescription.MOMENTSQR.description();
			}
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.TEMPERATURE.toString())){
			Description = EnumPotentialParamDescription.TEMPERATURE.description();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.STEPS.toString())){
			Description = EnumPotentialParamDescription.STEPS.description();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMAHSREF.toString())){
			Description = EnumPotentialParamDescription.SIGMAHSREF.description();
		}
		return Description;
	}


	public Double getDoubleDefaultParameters(String Parameter) {
		// TODO Auto-generated method stub
		
		Double parameterValue = null;
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMA.toString()+PotentialSites[i])){
				parameterValue = getSigma(i);
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.EPSILON.toString()+PotentialSites[i])){
				parameterValue = getEpsilon(i);
			}
		
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.MOMENTSQR.toString()+PotentialSites[i])){
				parameterValue = getMomentSquare(i);
			}
		}
		
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMAHSREF.toString())){
			parameterValue = getSigmaHSRef();
		}
		
		return parameterValue;
	}
	
	public String[] getParametersArray() {
		return ComponentParameters;
	}

	@Override
	public String getMolecularModelDisplayName() {
		// TODO Auto-generated method stub
		return "Spherical-2-Body-With-Q";
	}

	

	@Override
	public String[] getPotentialSites() {
		// TODO Auto-generated method stub
		return PotentialSites;
	}

	@Override
	public String getPotentialSiteAtIndex(int index) {
		
		return PotentialSites[index];
	
	}

	@Override
	public String getMoleculeDisplayName() {
		// TODO Auto-generated method stub
		return MoleculeDisplayName;
	}

	@SuppressWarnings("rawtypes")

	public Class getPotential() {
		// TODO Auto-generated method stub
		return P2LJQ.class;
	}

	@Override
	public Space getSpace() {
		// TODO Auto-generated method stub
		return this.space;
	}

	@Override
	public boolean hasElectrostaticInteraction() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String getNonBondedInteractionModel() {
		// TODO Auto-generated method stub
		return "LennardJonesWithQuadrapole";
	}

	
}
