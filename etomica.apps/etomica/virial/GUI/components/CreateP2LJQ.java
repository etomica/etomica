package etomica.virial.GUI.components;

import etomica.potential.P2LJQ;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateP2LJQ implements ParameterMapping,Cloneable{
	
	private ISpace space;
	private double sigma[];
	private double epsilon[];
	private double momentSquare[];
	private int id;
	private static int numberOfInstances = 0;
	
	
	
	private String CustomClassName = "Spherical-2-Body-With-Quad";
	
	private String[] ComponentParameters  = {"SIGMA","EPSILON","MOMENTSQR"};
	
	private String[] SharedComponentParameters =null;
	
	private String[] PotentialSites = {"LJ"};
	
	private String[][] ComponentValues = {{"1.0","1.0","1.0"}};
	
	private String[] SharedComponentValues = null;
	private String[][] ParamAndValues; 
	
	//Potentials references are created as Private members
	private P2LJQ p2LJQ;
	
	//Constructors for different Instantiations
	
	public CreateP2LJQ(){
		space = Space3D.getInstance();
		sigma = new double[PotentialSites.length];
		epsilon = new double[PotentialSites.length];
		momentSquare = new double[PotentialSites.length];
		ParamAndValues = setParameterValues();
		id = ++numberOfInstances;
	}
	
private String[][] setParameterValues() {
		
		int NoOfParam = ComponentParameters.length;
		
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
				if(ComponentParameters[j]=="MOMENTSQR"){
					setMomentSquare(Double.parseDouble(ComponentValues[i][j]),i);
				}
				ReturnArray[index][0] = ComponentParameters[j]+PotentialSites[i];
				ReturnArray[index][1] = ComponentValues[i][j];
				index++;
				
			}
		}
		return ReturnArray;
	}
	
	public String[][] getComponentValues() {
	return ComponentValues;
}

public void setComponentValues(String[][] componentValues) {
	ComponentValues = componentValues;
}

public String[][] getParamAndValues() {
	return ParamAndValues;
}

	public int getId() {
		return id;
	}

	//Setter method for LJ atomic potentials
	public void setP2LJQ(){
		this.p2LJQ = new P2LJQ(this.space); 
	}

	//Getter for p2LJQ
	public P2LJQ getP2LJQ() {
		return p2LJQ;
	}
	
	 public Object clone(){
		 try{
			 CreateP2LJQ cloned = (CreateP2LJQ)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
	 }
	
	//Creates the LJAtom Species
	public SpeciesFactory createSpeciesFactory(){
		SpeciesFactory speciesFactory;
		speciesFactory = new SpeciesFactorySpheres();
        return speciesFactory;
	}
	
	
	//Testing Class
	public static void main(String[] args){
		CreateP2LJQ lj = new CreateP2LJQ();
		
		System.out.println(lj.getDescription("EPSILONLJ"));
		System.out.println(lj.getDoubleDefaultParameters("EPSILONLJ"));
		for(int j=0;j<lj.ComponentParameters.length*lj.PotentialSites.length;j++){
			
			System.out.println(lj.ParamAndValues[j][0]+"\n");
			System.out.println(lj.ParamAndValues[j][1]+"\n");
	}
		//lj.setParameter("epsilon", "1.5");
		//System.out.println(lj.getEpsilon(i));
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

	public double getMomentSquare(int index) {
		return momentSquare[index];
	}

	public void setMomentSquare(double momentSquare,int index) {
		this.momentSquare[index] = momentSquare;
	}
	
	

	public int getParameterCount() {
		return 3;
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
			if(Parameter.toUpperCase().equals(ParametersDouble.MOMENTSQR.toString()+PotentialSites[i])){
				setMomentSquare(Double.parseDouble(ParameterValue),i); 
			}
			
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
		
			if(Parameter.toUpperCase().equals(ParametersDouble.MOMENTSQR.toString()+PotentialSites[i])){
				Description = ParametersDouble.MOMENTSQR.Description();
			}
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
		
			if(Parameter.toUpperCase().equals(ParametersDouble.MOMENTSQR.toString()+PotentialSites[i])){
				parameterValue = getMomentSquare(i);
			}
		}
		
		return parameterValue;
	}
	
	public String[] getParametersArray() {
		return ComponentParameters;
	}

	@Override
	public String getCustomName() {
		// TODO Auto-generated method stub
		return "Spherical-2-Body-With-Q";
	}

	

	@Override
	public String[] getPotentialSites() {
		// TODO Auto-generated method stub
		return PotentialSites;
	}

	@Override
	public String getPotentialSiteAtIndex(int index) {
		
		return PotentialSites[index];
	
	}

	
}
