package etomica.virial.GUI.components;


import etomica.config.ConformationLinear;
import etomica.potential.P22CLJQ;

import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;

import etomica.units.Kelvin;
import etomica.virial.SpeciesFactory;

import etomica.virial.SpeciesFactoryTangentSpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateP2CO22CLJQ implements ParameterMapping,Cloneable{
	private static String MoleculeDisplayName = "CO2 - 2CCLJQ";
	private Space space;
	private double[] sigma;
	private double[] epsilon;
	private double[] moment;
	
	private double temperature;
	private int noOfSteps;
	private double sigmaHSRef;
	
	
	private double bondLength;
	
	
	private ConformationLinear conformation;
	
	private int SpeciesID;
	
	private String[] PotentialSites = {"CO2"};
	
	
	
	public String getPotentialSiteAtIndex(int index) {
		return PotentialSites[index];
	}

	private int id;
	private static int numberOfInstances = 0;
	
	private String[] ComponentParameters  =  {"SIGMA","EPSILON", "MOMENT"};
			
	private String[] SharedComponentParameters ={"BONDL"};
	
	private String[][] ParamAndValues; 

	
	private String[][] ComponentValues = {{"3.0354",Double.toString(Kelvin.UNIT.toSim(125.317)),Double.toString(3.0255*Kelvin.UNIT.toSim(125.317)*Math.pow(3.0354,5))}};
	
	private String[] SharedComponentValues = {"2.1347"};
	
	
	private String[] SimEnvParameters = {"TEMPERATURE","STEPS","SIGMAHSREF"};
	
	private String[] SimEnvValues = {"250.0","10000",Double.toString(1.5*3.0354)};
	
	
	
	
	//Constructors for different Instantiations
	
	
	public CreateP2CO22CLJQ(){
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
		
		int NoOfSimEnvParam = 3;
		for(int l = 0;l<NoOfSimEnvParam;l++){
			if(SimEnvParameters[l]=="TEMPERATURE"){
				setTemperature(Double.parseDouble(SimEnvValues[l]));
			}
			
			if(SimEnvParameters[l]=="STEPS"){
				setNoOfSteps(Integer.parseInt(SimEnvValues[l]));
			}
			
			if(SimEnvParameters[l]=="SIGMAHSREF"){
				setSigmaHSRef(Double.parseDouble(SimEnvValues[l]));
			}
		}
		return ReturnArray;
		
		
	}
	
	
	public double getTemperature() {
		return temperature;
	}

	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	public int getNoOfSteps() {
		return noOfSteps;
	}

	public void setNoOfSteps(int noOfSteps) {
		this.noOfSteps = noOfSteps;
	}

	public double getSigmaHSRef() {
		return sigmaHSRef;
	}

	public void setSigmaHSRef(double sigmaHSRef) {
		this.sigmaHSRef = sigmaHSRef;
	}

	public int getId() {
		return id;
	}
	
	 public Object clone(){
		 try{
			 CreateP2CO22CLJQ cloned = (CreateP2CO22CLJQ)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
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

	public void setSigma(double sigma, int index) {
		this.sigma[index] = sigma;
	}

	public double getEpsilon(int index) {
		return epsilon[index];
	}

	public void setEpsilon(double epsilon, int index) {
		this.epsilon[index] = epsilon;
	}

	public double getMoment(int index) {
		return moment[index];
	}

	public void setMoment(double moment, int index) {
		this.moment[index] = moment;
	}

	public double getBondLength() {
		return bondLength;
	}
	
	public void setBondLength(double bondLength) {
		this.bondLength = bondLength;
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
		if(Parameter.toUpperCase().equals(ParametersDouble.TEMPERATURE.toString())){
			setTemperature(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.STEPS.toString())){
			setNoOfSteps(Integer.parseInt(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMAHSREF.toString())){
			setSigmaHSRef(Double.parseDouble(ParameterValue)); 
		}
		
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
		
		if(Parameter.toUpperCase().equals(ParametersDouble.TEMPERATURE.toString())){
			Description = ParametersDouble.TEMPERATURE.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.STEPS.toString())){
			Description = ParametersDouble.STEPS.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMAHSREF.toString())){
			Description = ParametersDouble.SIGMAHSREF.Description();
		}
		return Description;
	}

	//Testing Class
	public static void main(String[] args){
		CreateP2CO22CLJQ lj = new CreateP2CO22CLJQ();
		
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
		if(Parameter.toUpperCase().equals(ParametersDouble.TEMPERATURE.toString())){
			parameterValue = getTemperature();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMAHSREF.toString())){
			parameterValue = getSigmaHSRef();
		}
		
		if(Parameter.toUpperCase().equals(ParametersDouble.STEPS.toString())){
			parameterValue = (double) getNoOfSteps();
		}
		return parameterValue;
	}



	public String[] getParametersArray() {
		return ComponentParameters;
	}

	@Override
	public String getCustomName() {
		// TODO Auto-generated method stub
		return "2CLJQ";
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
		return P22CLJQ.class;
	}
	
	@Override
	public Space getSpace() {
		// TODO Auto-generated method stub
		return this.space;
	}


}
