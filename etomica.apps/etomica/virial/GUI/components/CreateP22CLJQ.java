package etomica.virial.GUI.components;


import etomica.config.ConformationLinear;
import etomica.potential.P22CLJQ;

import etomica.space.ISpace;
import etomica.space3d.Space3D;

import etomica.virial.SpeciesFactory;

import etomica.virial.SpeciesFactoryTangentSpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateP22CLJQ implements ParameterMapping,Cloneable{
	
	private ISpace space;
	private double sigma;
	private double epsilon;
	private double moment;
	private double bondLength;
	private ConformationLinear conformation;
	
	private int id;
	private static int numberOfInstances = 0;
	
	private String[] ParametersArray  = {"SIGMA","EPSILON", "MOMENT","BONDL"};
	
	
	
	//Potentials references are created as Private members
	private P22CLJQ p22CLJQ;
	
	//Constructors for different Instantiations
	
	
	public CreateP22CLJQ(){
		space = Space3D.getInstance();
		this.sigma = 1.0;
		this.epsilon = 1.0;
		this.moment = 1.0;
		this.bondLength = 1.0;
		
		id=++numberOfInstances;
	}
	
	public int getId() {
		return id;
	}
	
	 public Object clone(){
		 try{
			 CreateP22CLJQ cloned = (CreateP22CLJQ)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
	 }

	//Sets the LJ Molecular Potential
	public void setP22CLJQ(){
		this.p22CLJQ = new P22CLJQ(this.space);
		
	}
	
	//Gets the LJ Molecular Potential
	public P22CLJQ getP22CLJQ(){
		return this.p22CLJQ;
	}
	
	
	//Getter for ConformationLinear class
	public ConformationLinear getConformation() {
		return this.conformation;
	}

	public void setConformation() {
		this.conformation = new ConformationLinear(this.space, this.bondLength);
	}
	
	//Creates the LJ Molecule Species
	public SpeciesFactory createP22CLJQSpecies(){
		SpeciesFactory speciesFactory = new SpeciesFactoryTangentSpheres(2,this.getConformation());
		return speciesFactory;
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

	public double getMoment() {
		return moment;
	}

	public void setMoment(double moment) {
		this.moment = moment;
	}

	public double getBondLength() {
		return bondLength;
	}
	
	@Override
	public int getParameterCount() {
		return 4;
	}

	@Override
	public void setParameter(String Parameter, String ParameterValue) {
		// TODO Auto-generated method stub
		
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString())){
			setSigma(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString())){
			setEpsilon(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.MOMENT.toString())){
			setEpsilon(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
			setEpsilon(Double.parseDouble(ParameterValue)); 
		}
		
	}

	@Override
	public String getDescription(String Parameter) {
		String Description = null;
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString())){
			Description = ParametersDouble.SIGMA.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString())){
			Description = ParametersDouble.EPSILON.Description();
		}
		
		if(Parameter.toUpperCase().equals(ParametersDouble.MOMENT.toString())){
			Description = ParametersDouble.MOMENT.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
			Description = ParametersDouble.BONDL.Description();
		
		}
		return Description;
	}

	//Testing Class
	public static void main(String[] args){
		CreateP22CLJQ lj = new CreateP22CLJQ();
		
	}

	@Override
	public Double getDoubleDefaultParameters(String Parameter) {
		// TODO Auto-generated method stub
		Double parameterValue = null;
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString())){
			parameterValue = ParametersDouble.SIGMA.DefaultValue();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString())){
			parameterValue = ParametersDouble.EPSILON.DefaultValue();
		}
		
		if(Parameter.toUpperCase().equals(ParametersDouble.MOMENT.toString())){
			parameterValue = ParametersDouble.MOMENT.DefaultValue();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
			parameterValue = ParametersDouble.BONDL.DefaultValue();
		}
		return parameterValue;
	}



	public String[] getParametersArray() {
		return ParametersArray;
	}
	
}
