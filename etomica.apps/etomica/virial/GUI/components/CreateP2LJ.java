package etomica.virial.GUI.components;


import etomica.potential.P2LennardJones;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;

public class CreateP2LJ implements ParameterMapping{
	
	private ISpace space;
	private double sigma;
	private double epsilon;
	
	//Potentials references are created as Private members
	private P2LennardJones p2LJ;
	

	//Constructors for different Instantiations
	
	CreateP2LJ(){
		space = Space3D.getInstance();
		sigma = 1.0;
		epsilon = 1.0;
	}
	
	//Setter method for LJ atomic potentials
	public void setP2LJ(){
		this.p2LJ = new P2LennardJones(this.space,this.sigma,this.epsilon); 
	}

	//Getter for p2LJ
	public P2LennardJones getP2LJ() {
		return p2LJ;
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
	public SpeciesFactory createP2LJSpeciesFactory(){
		SpeciesFactory speciesFactory = new SpeciesFactorySpheres();
        return speciesFactory;
	}

	@Override
	public int getParameterCount() {
		return 2;
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
	}

	@Override
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
		
		return parameterValue;
	}
	
	
	
}
