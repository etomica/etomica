package etomica.virial.GUI.components;


import etomica.potential.P2LennardJones;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateP2LJ implements ParameterMapping,Cloneable{
	
	private ISpace space;
	private double sigma;
	private double epsilon;
	private int id;
	private static int numberOfInstances = 0;
	
	
	private String[] ParametersArray  = {"SIGMA","EPSILON"};
	
	//Potentials references are created as Private members
	private P2LennardJones p2LJ;
	

	//Constructors for different Instantiations
	
	public CreateP2LJ(){
		space = Space3D.getInstance();
		sigma = 1.0;
		epsilon = 1.0;
		
		id=++numberOfInstances;
	}
	
	public int getId() {
		return id;
	}

	//Setter method for LJ atomic potentials
	public void setP2LJ(){
		this.p2LJ = new P2LennardJones(this.space,this.sigma,this.epsilon); 
	}

	//Getter for p2LJ
	public P2LennardJones getP2LJ() {
		return p2LJ;
	}
	
	 public Object clone(){
		 try{
			 CreateP2LJ cloned = (CreateP2LJ)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
	 }

	public double getSigma() {
		return sigma;
	}


	public void setSigma(double sigma) {
		this.sigma = sigma;
	}


	public double getEpsilon() {
		return epsilon;
	}

	public void setEpsilon(double epsilon) {
		this.epsilon = epsilon;
	}

	//Creates the LJAtom Species
	public SpeciesFactory createSpeciesFactory(){
		SpeciesFactory speciesFactory = new SpeciesFactorySpheres();
        return speciesFactory;
	}


	public int getParameterCount() {
		return 2;
	}


	public void setParameter(String Parameter, String ParameterValue) {
		// TODO Auto-generated method stub
		
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString())){
			setSigma(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString())){
			setEpsilon(Double.parseDouble(ParameterValue)); 
		}
	}


	public String getDescription(String Parameter) {
		String Description = null;
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString())){
			Description = ParametersDouble.SIGMA.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString())){
			Description = ParametersDouble.SIGMA.Description();
		}
		return Description;
	}


	public Double getDoubleDefaultParameters(String Parameter) {
		// TODO Auto-generated method stub
		Double parameterValue = null;
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString())){
			parameterValue = ParametersDouble.SIGMA.DefaultValue();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString())){
			parameterValue = ParametersDouble.EPSILON.DefaultValue();
		}
		
		return parameterValue;
	}
	
	public String[] getParametersArray() {
		return ParametersArray;
	}

	
	
}
