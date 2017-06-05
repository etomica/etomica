/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;


import etomica.atom.IAtomList;
import etomica.api.IElement;
import etomica.api.ISpecies;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.config.IConformation;
import etomica.potential.P2CO2EMP2;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.virial.SpeciesFactory;

public class MolecularModelEMP2_SpeciesCO2 implements IMolecularModel_SpeciesFactory,Cloneable{
	private static String MoleculeDisplayName = "CO2 - EPM";
	private Space space;
	private double[] sigma;
	private double[] epsilon;
	private double[] charge;
	private double bondL;
	//private IConformation conformation;


	private double sigmaHSRef;
	
	private int id;
	private static int numberOfInstances = 0;
	
	
	private String[] ComponentParameters  =  {"SIGMA","EPSILON","CHARGE"};
			
	private String[] SharedComponentParameters ={"BONDL"};
	
	private String[][] ParamAndValues; 
	
	private String[] PotentialSites = {"C","O"};
	
	private String[] SimEnvParameters = {"SIGMAHSREF"};
	
	private String[] SimEnvValues = {Double.toString(5.0*1.149)};
	
	private String[][] ComponentValues = {
			{"2.7570",Double.toString(Kelvin.UNIT.toSim(28.129)),Double.toString(Electron.UNIT.toSim(0.6512))},
			{"3.0330",Double.toString(Kelvin.UNIT.toSim(80.507)),Double.toString((-0.5)*Electron.UNIT.toSim(0.6512))},
	};

	private String[] SharedComponentValues = {"1.149"};
	
	public String getPotentialSiteAtIndex(int index) {
		
		return PotentialSites[index];
	}
	
	public boolean IsPotentialSiteMoreThanOne(){
		return true;
	}
	

	public MolecularModelEMP2_SpeciesCO2(){
		space = Space3D.getInstance();
		sigma = new double[PotentialSites.length];
		epsilon = new double[PotentialSites.length];
		charge = new double[PotentialSites.length];
		ParamAndValues=setParameterValues();
		//setConformation();
		id=++numberOfInstances;
	}
	
private String[][] setParameterValues() {
		
		int NoOfParam = ComponentParameters.length;
		int NoOfCommonParam = SharedComponentParameters.length;
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
				if(ComponentParameters[j]=="CHARGE"){
					setCharge(Double.parseDouble(ComponentValues[i][j]),i);
				}
				
				ReturnArray[index][0] = ComponentParameters[j]+PotentialSites[i];
				ReturnArray[index][1] = ComponentValues[i][j];
				index++;
			}
		}
		for(int k = 0;k<NoOfCommonParam;k++){
			if(SharedComponentParameters[k]=="BONDL"){
				setBondL(Double.parseDouble(SharedComponentValues[k]));
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
	
	
	
public double getSigmaHSRef() {
	return sigmaHSRef;
}

public void setSigmaHSRef(double sigmaHSRef) {
	this.sigmaHSRef = sigmaHSRef;
}

	public double getCharge(int index) {
		return charge[index];
	}

	public void setCharge(double charge, int index) {
		this.charge[index] = charge;
	}

	public double getBondL() {
		return bondL;
	}

	public void setBondL(double bondL) {
		this.bondL = bondL;
	}
	
	public int getId() {
		return id;
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

	public void setEpsilon(double epsilon, int index) {
		this.epsilon[index] = epsilon;
	}

	/*
	public void setConformation(){
		conformation = new IConformation(){
			public void initializePositions(IAtomList atomList) {
                // atoms are C, O and O, so we arrange them as 1-0-2
               
                atomList.getAtom(0).getPosition().E(0);
                atomList.getAtom(1).getPosition().E(0);
                atomList.getAtom(1).getPosition().setX(0, -bondL);
                atomList.getAtom(2).getPosition().E(0);
                atomList.getAtom(2).getPosition().setX(0, +bondL);
            }
		};
	}
		
        
        */
    
	
	 public Object clone(){
		 try{
			 MolecularModelEMP2_SpeciesCO2 cloned = (MolecularModelEMP2_SpeciesCO2)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
	 }

	
	//Creates the LJAtom Species
	public ISpecies createSpecies(){
		SpeciesFactory factory = new SpeciesFactory() {
	        public ISpecies makeSpecies(Space space) {
	            SpeciesSpheresHetero species = new SpeciesSpheresHetero(space, new IElement[]{Carbon.INSTANCE, Oxygen.INSTANCE});
	            species.setChildCount(new int[]{1,2});
	            
	            IConformation conformation = new IConformation(){
	    			public void initializePositions(IAtomList atomList) {
	                    // atoms are C, O and O, so we arrange them as 1-0-2
	                   
	                    atomList.getAtom(0).getPosition().E(0);
	                    atomList.getAtom(1).getPosition().E(0);
	                    atomList.getAtom(1).getPosition().setX(0, -bondL);
	                    atomList.getAtom(2).getPosition().E(0);
	                    atomList.getAtom(2).getPosition().setX(0, +bondL);
	                }
	    		};
	            
	            species.setConformation(conformation);
	            return species;
	        }
	    };
	    return factory.makeSpecies(this.space);
	}

	
	public SpeciesFactory createSpeciesFactory(){
		SpeciesFactory factory = new SpeciesFactory() {
	        public ISpecies makeSpecies(Space space) {
	            SpeciesSpheresHetero species = new SpeciesSpheresHetero(space, new IElement[]{Carbon.INSTANCE, Oxygen.INSTANCE});
	            species.setChildCount(new int[]{1,2});
	            
	            IConformation conformation = new IConformation(){
	    			public void initializePositions(IAtomList atomList) {
	                    // atoms are C, O and O, so we arrange them as 1-0-2
	                   
	                    atomList.getAtom(0).getPosition().E(0);
	                    atomList.getAtom(1).getPosition().E(0);
	                    atomList.getAtom(1).getPosition().setX(0, -bondL);
	                    atomList.getAtom(2).getPosition().E(0);
	                    atomList.getAtom(2).getPosition().setX(0, +bondL);
	                }
	    		};
	            
	            species.setConformation(conformation);
	            return species;
	        }
	    };
	    return factory;
	}
	
	public int getParameterCount() {
		return 2;
	}


	public void setParameter(String Parameter, String ParameterValue) {
		// TODO Auto-generated method stub
		
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMA.toString()+PotentialSites[i])){
				setSigma(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.EPSILON.toString()+PotentialSites[i])){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.CHARGE.toString()+getPotentialSiteAtIndex(i))){
				setCharge(Double.parseDouble(ParameterValue),i); 
			}
		}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.BONDL.toString())){
				setBondL(Double.parseDouble(ParameterValue)); 
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

			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.CHARGE.toString()+PotentialSites[i])){
				Description = EnumPotentialParamDescription.CHARGE.description();
			}
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.BONDL.toString())){
			Description = EnumPotentialParamDescription.BONDL.description();
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
		
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.CHARGE.toString()+PotentialSites[i])){
				parameterValue = getCharge(i);
			}
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.BONDL.toString())){
			parameterValue = getBondL();
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
		return "EPM2";
	}

	public String getPotentialSites(int index) {
		return PotentialSites[index];
	}

	@Override
	public String[][] getParamAndValues() {
		// TODO Auto-generated method stub
		return ParamAndValues;
	}

	@Override
	public String[] getPotentialSites() {
		// TODO Auto-generated method stub
		return PotentialSites;
	}

	@Override
	public String getMoleculeDisplayName() {
		// TODO Auto-generated method stub
		return MoleculeDisplayName;
	}

	@SuppressWarnings("rawtypes")
	
	public Class getPotential() {
		// TODO Auto-generated method stub
		return P2CO2EMP2.class;
	}

	@Override
	public Space getSpace() {
		// TODO Auto-generated method stub
		return this.space;
	}


	public boolean hasElectrostaticInteraction() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String getNonBondedInteractionModel() {
		// TODO Auto-generated method stub
		return "LennardJones";
	}

}
