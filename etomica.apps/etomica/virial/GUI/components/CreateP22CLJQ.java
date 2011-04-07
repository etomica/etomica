package etomica.virial.GUI.components;

import etomica.api.ISpecies;
import etomica.config.ConformationLinear;
import etomica.potential.P22CLJQ;
import etomica.potential.P2LJQ;
import etomica.potential.P2LennardJones;
import etomica.potential.P2LennardJonesDreiding;
import etomica.potential.Potential;
import etomica.potential.PotentialMolecular;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactoryOrientedSpheres;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.SpeciesFactoryTangentSpheres;

public class CreateP22CLJQ implements ParameterMapping{
	
	private ISpace space;
	private double sigma;
	private double epsilon;
	private double moment;
	private double bondLength;
	private ConformationLinear conformation;
	
	//Potentials references are created as Private members
	private P22CLJQ p22CLJQ;
	
	//Constructors for different Instantiations
	
	
	CreateP22CLJQ(){
		space = Space3D.getInstance();
		this.sigma = 1.0;
		this.epsilon = 1.0;
		this.moment = 1.0;
		this.bondLength = 1.0;
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
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDLENGTH.toString())){
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
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDLENGTH.toString())){
			Description = ParametersDouble.BONDLENGTH.Description();
		
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
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDLENGTH.toString())){
			parameterValue = ParametersDouble.BONDLENGTH.DefaultValue();
		}
		return parameterValue;
	}
	
}
