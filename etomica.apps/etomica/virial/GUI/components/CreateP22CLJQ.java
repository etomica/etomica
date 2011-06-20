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
	private double[] sigma;
	private double[] epsilon;
	public void setPotentialSites(String[] potentialSites) {
		PotentialSites = potentialSites;
	}

	private double[] moment;
	private double bondLength;
	private ConformationLinear conformation;
	
	private int SpeciesID;
	private String[][] ParamAndValues; 

	

	private int id;
	private static int numberOfInstances = 0;
	
	private String[] ComponentParameters  = {"SIGMA","EPSILON","MOMENT"};
	
	private String[] SharedComponentParameters ={"BONDL"};
	
	private String[] PotentialSites = {"LJ"};
	
	private String[][] ComponentValues = {{"1.0","1.0","1.0"}};
	
	private String[] SharedComponentValues = {"1.0"};
	
	
	//Potentials references are created as Private members
	private P22CLJQ p22CLJQ;
	
	//Constructors for different Instantiations
	
	
	public CreateP22CLJQ(){
		space = Space3D.getInstance();
		sigma = new double[PotentialSites.length];
		epsilon = new double[PotentialSites.length];
		moment = new double[PotentialSites.length];
		ParamAndValues=setParameterValues();
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
				if(ComponentParameters[j]=="MOMENT"){
					setMoment(Double.parseDouble(ComponentValues[i][j]),i);
					
				}
				
				ReturnArray[index][0] = ComponentParameters[j]+PotentialSites[i];
				ReturnArray[index][1] = ComponentValues[i][j];
				index++;
			}
		}
		for(int k = 0;k<NoOfCommonParam;k++){
			if(SharedComponentParameters[k]=="BONDL"){
				setBondLength(Double.parseDouble(SharedComponentValues[k]));
			}
		}
		return ReturnArray;
		
		
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
	
	public int getSpeciesID() {
		return SpeciesID;
	}

	public void setSpeciesID(int speciesID) {
		SpeciesID = speciesID;
	}
	//Getter for ConformationLinear class
	public ConformationLinear getConformation() {
		return this.conformation;
	}

	public void setConformation() {
		this.conformation = new ConformationLinear(this.space, this.bondLength);
	}
	
	//Creates the LJ Molecule Species
	public SpeciesFactory createSpeciesFactory(){
		SpeciesFactory speciesFactory = new SpeciesFactoryTangentSpheres(2,this.getConformation());
		return speciesFactory;
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

	public double getMoment(int index) {
		return moment[index];
	}

	public void setMoment(double moment,int index) {
		this.moment[index] = moment;
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
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+PotentialSites[i])){
				setSigma(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+PotentialSites[i])){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.MOMENT.toString()+PotentialSites[i])){
				setMoment(Double.parseDouble(ParameterValue),i); 
			}
			
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
			setBondLength(Double.parseDouble(ParameterValue)); 
		}
	}

	public void setBondLength(double bondLength) {
		this.bondLength = bondLength;
	}

	@Override
	public String getDescription(String Parameter) {
		String Description = null;
		for(int i = 0;i <PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+PotentialSites[i])){
				Description = ParametersDouble.SIGMA.Description();
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+PotentialSites[i])){
				Description = ParametersDouble.EPSILON.Description();
			}
		
			if(Parameter.toUpperCase().equals(ParametersDouble.MOMENT.toString()+PotentialSites[i])){
				Description = ParametersDouble.MOMENT.Description();
			}
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
			Description = ParametersDouble.BONDL.Description();
		}
		return Description;
	}

	//Testing Class
	public static void main(String[] args){
		CreateP22CLJQ lj = new CreateP22CLJQ();
		for(int j=0;j<lj.ComponentParameters.length*lj.PotentialSites.length;j++){
			
				System.out.println(lj.ParamAndValues[j][0]+"\n");
				System.out.println(lj.ParamAndValues[j][1]+"\n");
		}
		
		
	}

	public String[][] getParamAndValues() {
		return ParamAndValues;
	}

	@Override
	public Double getDoubleDefaultParameters(String Parameter) {
		// TODO Auto-generated method stub
		Double parameterValue = null;
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+PotentialSites[i])){
				parameterValue = getSigma(i);
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+PotentialSites[i])){
				parameterValue = getEpsilon(i);
			}
		
			if(Parameter.toUpperCase().equals(ParametersDouble.MOMENT.toString()+PotentialSites[i])){
				parameterValue = getMoment(i);
			}
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
			parameterValue = getBondLength();
		}
		return parameterValue;
	}



	
	@Override
	public String getCustomName() {
		// TODO Auto-generated method stub
		return "2CLJQ";
	}


	public String[] getPotentialSites() {
		return PotentialSites;
	}

	
	public String getPotentialSiteAtIndex(int index) {
		return PotentialSites[index];
	}
	
	@Override
	public String[] getParametersArray() {
		// TODO Auto-generated method stub
		return ComponentParameters;
	}

	
	
}
