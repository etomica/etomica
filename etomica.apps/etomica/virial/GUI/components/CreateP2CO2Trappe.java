package etomica.virial.GUI.components;


import etomica.api.IAtomList;
import etomica.api.IElement;
import etomica.api.ISpecies;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.config.IConformation;
import etomica.potential.P22CLJQ;
import etomica.potential.P2CO2EMP2;
import etomica.potential.P2CO2TraPPE;
import etomica.potential.P2LennardJones;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateP2CO2Trappe implements ParameterMapping,Cloneable{
	private static String MoleculeDisplayName = "CO2 - Trappe";
	private Space space;
	private double[] sigma;
	private double[] epsilon;
	private double[] charge;
	private double bondL;
	//private IConformation conformation;

	private double temperature;
	private int noOfSteps;
	private double sigmaHSRef;
	
	private int id;
	private static int numberOfInstances = 0;
	
	
	private String[] ComponentParameters  =  {"SIGMA","EPSILON","CHARGE"};
			
	private String[] SharedComponentParameters ={"BONDL"};
	
	private String[][] ParamAndValues; 
	
	private String[] PotentialSites = {"C","H"};
	
	private String[][] ComponentValues = {
			//{"2.8000",Double.toString(Kelvin.UNIT.toSim(27)),Double.toString(Electron.UNIT.toSim(0.70))},
			{"3.1670",Double.toString(Kelvin.UNIT.toSim(78.23)),Double.toString(Electron.UNIT.toSim(-0.82))},
			{"3.0500",Double.toString(Kelvin.UNIT.toSim(79)),Double.toString((-0.5)*Electron.UNIT.toSim(0.70))}
			//{"3.1670",Double.toString(Kelvin.UNIT.toSim(78.23)),Double.toString(Electron.UNIT.toSim(-0.82))}
	};


	private String[] SharedComponentValues = {"1.1491"};
	
	private String[] SimEnvParameters = {"TEMPERATURE","STEPS","SIGMAHSREF"};
	
	private String[] SimEnvValues = {"250.0","10000",Double.toString(5.0*1.1491)};
	
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
	
	
	public String getPotentialSiteAtIndex(int index) {
		
		return PotentialSites[index];
	}
	
	public boolean IsPotentialSiteMoreThanOne(){
		return true;
	}
	
	

	//Constructors for different Instantiations


	public CreateP2CO2Trappe(){
		space = Space3D.getInstance();
		sigma = new double[PotentialSites.length];
		epsilon = new double[PotentialSites.length];
		charge = new double[PotentialSites.length];
		ParamAndValues=setParameterValues();
		//setConformation();
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
				if(ComponentParameters[j]=="CHARGE"){
					setCharge(Double.parseDouble(ComponentValues[i][j]),i);
				}
				
				ReturnArray[index][0] = ComponentParameters[j]+PotentialSites[i];
				ReturnArray[index][1] = ComponentValues[i][j];
				index++;
			}
		}
		for(int k = 0;k<NoOfCommonParam;k++){
			if(SharedComponentParameters[k]=="BONDL"){
				setBondL(Double.parseDouble(SharedComponentValues[k]));
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
			 CreateP2CO2Trappe cloned = ( CreateP2CO2Trappe)super.clone();
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
		
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+PotentialSites[i])){
				setSigma(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+PotentialSites[i])){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.CHARGE.toString()+getPotentialSiteAtIndex(i))){
				setCharge(Double.parseDouble(ParameterValue),i); 
			}
		}
			if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
				setBondL(Double.parseDouble(ParameterValue)); 
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


	public String getDescription(String Parameter) {
		String Description = null;
		for(int i = 0;i <PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(ParametersDouble.SIGMA.toString()+PotentialSites[i])){
				Description = ParametersDouble.SIGMA.Description();
			}
			if(Parameter.toUpperCase().equals(ParametersDouble.EPSILON.toString()+PotentialSites[i])){
				Description = ParametersDouble.EPSILON.Description();
			}

			if(Parameter.toUpperCase().equals(ParametersDouble.CHARGE.toString()+PotentialSites[i])){
				Description = ParametersDouble.CHARGE.Description();
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
		
			if(Parameter.toUpperCase().equals(ParametersDouble.CHARGE.toString()+PotentialSites[i])){
				parameterValue = getCharge(i);
			}
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.BONDL.toString())){
			parameterValue = getBondL();
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
		return "TRAPPE";
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
		return P2CO2TraPPE.class;
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
