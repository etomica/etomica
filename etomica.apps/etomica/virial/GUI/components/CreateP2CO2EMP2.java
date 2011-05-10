package etomica.virial.GUI.components;


import etomica.api.IAtomList;
import etomica.api.IElement;
import etomica.api.ISpecies;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.config.IConformation;
import etomica.potential.P2CO2EMP2;
import etomica.potential.P2LennardJones;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateP2CO2EMP2 implements ParameterMapping,Cloneable{
	
	private ISpace space;
	private double[] sigma;
	private double[] epsilon;
	private double[] charge;
	private double bondL;
	//private IConformation conformation;
	
	private int id;
	private static int numberOfInstances = 0;
	
	
	private String[] ParametersArray  = {"SIGMA","EPSILON","CHARGE","BONDL"};
	
	private String[][] PotentialSites = {{"C","2"},{"O","3"},{"CO","4"}};

	public String[] getPotentialSiteAtIndex(int index) {
		String[] tempReturn = new String[2];
		tempReturn[0]= PotentialSites[index][0];
		tempReturn[1]= PotentialSites[index][1];
		return tempReturn;
	}
	
	public boolean IsPotentialSiteMoreThanOne(){
		return true;
	}
	//Potentials references are created as Private members
	private P2CO2EMP2 p2CO2EPM2;
	

	//Constructors for different Instantiations
	
	public CreateP2CO2EMP2(){
		space = Space3D.getInstance();
		sigma = new double[3];
		epsilon = new double[3];
		charge = new double[2];
		sigma[0] = 2.757;
		sigma[1] = 3.033;
		sigma[2] = 2.8921;
		epsilon[0]= Kel;
		epsilon[1] = 1.0;
		epsilon[2] = 1.0;
		charge[0] = 1.0;
		charge[1] = 1.0;
		bondL = 1.149;
		//setConformation();
		id=++numberOfInstances;
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
			 CreateP2CO2EMP2 cloned = (CreateP2CO2EMP2)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
	 }

	
	//Creates the LJAtom Species
	public SpeciesFactory createSpeciesFactory(){
		SpeciesFactory factory = new SpeciesFactory() {
	        public ISpecies makeSpecies(ISpace space) {
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
		
		for(int i = 0;i< PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+getPotentialSiteAtIndex(i))){
				setSigma(Double.parseDouble(ParameterValue),i); 
			}
			
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+getPotentialSiteAtIndex(i))){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.CHARGE.toString()+getPotentialSiteAtIndex(i))){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
			}
		}
	}


	public String getDescription(String Parameter) {
		String Description = null;
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+"C") || 
				Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+"O") ||
				Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+"CO") ){
			Description = ParametersDouble.SIGMA.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+"C") || 
				Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+"O") ||
				Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+"CO") ){
			Description = ParametersDouble.EPSILON.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.CHARGE.toString()+"C") || 
				Parameter.toUpperCase().equals(ParametersDouble.CHARGE.toString()+"O")){
			Description = ParametersDouble.CHARGE.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
			Description = ParametersDouble.BONDL.Description();
		}
		return Description;
	}


	public Double getDoubleDefaultParameters(String Parameter) {
		// TODO Auto-generated method stub
		Double parameterValue = null;
		
		for(int i = 0;i< PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+getPotentialSiteAtIndex(i))){
				
				if(getPotentialSiteAtIndex(i).equals("C")){
					parameterValue = ParametersDouble.SIGMA.DefaultValue(2); 
				}
				if(getPotentialSiteAtIndex(i).equals("O")){
					parameterValue = ParametersDouble.SIGMA.DefaultValue(3); 
				}
				if(getPotentialSiteAtIndex(i).equals("CO")){
					parameterValue = ParametersDouble.SIGMA.DefaultValue(4); 
				}
			}
			
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+getPotentialSiteAtIndex(i))){
				
				if(getPotentialSiteAtIndex(i).equals("C")){
					parameterValue = ParametersDouble.EPSILON.DefaultValue(2); 
				}
				if(getPotentialSiteAtIndex(i).equals("O")){
					parameterValue = ParametersDouble.EPSILON.DefaultValue(3); 
				}
				if(getPotentialSiteAtIndex(i).equals("CO")){
					parameterValue = ParametersDouble.EPSILON.DefaultValue(4); 
				}
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.CHARGE.toString()+getPotentialSiteAtIndex(i))){
				
				if(getPotentialSiteAtIndex(i).equals("C")){
					parameterValue = ParametersDouble.CHARGE.DefaultValue(2); 
				}
				if(getPotentialSiteAtIndex(i).equals("O")){
					parameterValue = ParametersDouble.CHARGE.DefaultValue(3); 
				}
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
				parameterValue = ParametersDouble.BONDL.DefaultValue(2);
			}
		}
		return parameterValue;
	}
	
	public String[] getParametersArray() {
		return ParametersArray;
	}

	@Override
	public String getCustomName() {
		// TODO Auto-generated method stub
		return "EPM2";
	}

	@Override
	public String[][] getPotentialSites() {
		// TODO Auto-generated method stub
		return PotentialSites;
	}

	
	
}
