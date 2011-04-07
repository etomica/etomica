package etomica.virial.GUI.components;

import etomica.potential.P2LJQ;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;

public class CreateP2LJQ implements ParameterMapping{
	
	private ISpace space;
	private double sigma;
	private double epsilon;
	private double momentSquare;
	
	
	//Potentials references are created as Private members
	private P2LJQ p2LJQ;
	
	//Constructors for different Instantiations
	
	CreateP2LJQ(){
		space = Space3D.getInstance();
		this.sigma = 1.0;
		this.epsilon = 1.0;
		this.momentSquare = 1.0;
	}
	
	//Setter method for LJ atomic potentials
	public void setP2LJQ(){
		this.p2LJQ = new P2LJQ(this.space); 
	}

	//Getter for p2LJQ
	public P2LJQ getP2LJQ() {
		return p2LJQ;
	}
	
	
	//Creates the LJAtom Species
	public SpeciesFactory createP2LJQSpecies(){
		SpeciesFactory speciesFactory;
		speciesFactory = new SpeciesFactorySpheres();
        return speciesFactory;
	}
	
	
	//Testing Class
	public static void main(String[] args){
		CreateP2LJQ lj = new CreateP2LJQ();
		System.out.println(lj.getDescription("epsilon"));
		System.out.println(lj.getDoubleDefaultParameters("epsilon"));
		lj.setParameter("epsilon", "1.5");
		System.out.println(lj.getEpsilon());
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

	public double getMomentSquare() {
		return momentSquare;
	}

	public void setMomentSquare(double momentSquare) {
		this.momentSquare = momentSquare;
	}
	
	
	@Override
	public int getParameterCount() {
		return 3;
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
		if(Parameter.toUpperCase().equals(ParametersDouble.MOMENTSQUARED.toString())){
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
		
		if(Parameter.toUpperCase().equals(ParametersDouble.MOMENTSQUARED.toString())){
			Description = ParametersDouble.MOMENTSQUARED.Description();
		}
		
		return Description;
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
		
		if(Parameter.toUpperCase().equals(ParametersDouble.MOMENTSQUARED.toString())){
			parameterValue = ParametersDouble.MOMENTSQUARED.DefaultValue();
		}
		
		return parameterValue;
	}
	
}
