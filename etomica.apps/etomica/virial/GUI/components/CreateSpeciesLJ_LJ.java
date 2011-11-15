package etomica.virial.GUI.components;


import etomica.api.ISpecies;
import etomica.potential.P2LennardJones;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateSpeciesLJ_LJ implements MixtureBuilderSpeciesFactory,Cloneable{
	
	private static String MoleculeDisplayName = "P2LennardJones";
	private Space space;
	private double sigma[];
	private double epsilon[];
	

	private double sigmaHSRef;
	
	private int id;
	private static int numberOfInstances = 0;
	private String[][] ParamAndValues; 
	private String[] ComponentParameters  = {"SIGMA","EPSILON"};
	
	private String[] SharedComponentParameters = null;
	
	private String[] PotentialSites = {"LJ"};
	
	private String[][] ComponentValues = {{"1.0","1.0"}};
	
	private String[] SharedComponentValues = null;
	
	private String[] SimEnvParameters = {"SIGMAHSREF"};
	
	private String[] SimEnvValues = {"1.5"};
	
	

	//Constructors for different Instantiations
	
	public CreateSpeciesLJ_LJ(){
		space = Space3D.getInstance();
		sigma = new double[PotentialSites.length];
		epsilon = new double[PotentialSites.length];
		ParamAndValues=setParameterValues();
		//p2LJ = new P2LennardJones[PotentialSites.length];
		id=++numberOfInstances;
	}
	
private String[][] setParameterValues() {
		
		int NoOfParam = ComponentParameters.length;
		int NoOfSimEnvParam = 1;
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
				ReturnArray[index][0] = ComponentParameters[j]+PotentialSites[i];
				ReturnArray[index][1] = ComponentValues[i][j];
				index++;
			}
		}
		if (SharedComponentParameters != null){
			
		}
		
		for(int l = 0;l<NoOfSimEnvParam;l++){
			
			if(SimEnvParameters[l]=="SIGMAHSREF"){
				setSigmaHSRef(Double.parseDouble(SimEnvValues[l]));
			}
		}
		return ReturnArray;
		
	}

	

	
	 public Object clone(){
		 try{
			 CreateSpeciesLJ_LJ cloned = (CreateSpeciesLJ_LJ)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
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
	
	

	public double getSigmaHSRef() {
		return sigmaHSRef;
	}

	public void setSigmaHSRef(double sigmaHSRef) {
		this.sigmaHSRef = sigmaHSRef;
	}

	public String[][] getParamAndValues() {
		return ParamAndValues;
	}

	public String getPotentialSiteAtIndex(int index) {
		return PotentialSites[index];
	}
		
	public int getId() {
		return id;
	}
	
	//Creates the LJAtom Species
	public ISpecies createSpecies(){
		SpeciesFactory speciesFactory = new SpeciesFactorySpheres();
        return speciesFactory.makeSpecies(this.space);
	}

	public SpeciesFactory createSpeciesFactory(){
		SpeciesFactory speciesFactory = new SpeciesFactorySpheres();
        return speciesFactory;
	}

	public int getParameterCount() {
		return 2;
	}


	public void setParameter(String Parameter, String ParameterValue) {
		// TODO Auto-generated method stub
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+PotentialSites[i])){
				setSigma(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+PotentialSites[i])){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
			}
		}
		
		
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMAHSREF.toString())){
			setSigmaHSRef(Double.parseDouble(ParameterValue)); 
		}
		
	}


	public String getDescription(String Parameter) {
		String Description = null;
		for(int i = 0;i <PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+PotentialSites[i])){
				Description = ParametersDouble.SIGMA.Description();
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+PotentialSites[i])){
				Description = ParametersDouble.EPSILON.Description();
			}
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
		}
		
		
		if(Parameter.toUpperCase().equals(ParametersDouble.SIGMAHSREF.toString())){
			parameterValue = getSigmaHSRef();
		}
		
		
		
		return parameterValue;
	}
	
	public String[] getParametersArray() {
		return ComponentParameters;
	}

	@Override
	public String getCustomName() {
		// TODO Auto-generated method stub
		return "Spherical-2-Body";
	}


	public String[] getPotentialSites() {
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
		return P2LennardJones.class;
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
