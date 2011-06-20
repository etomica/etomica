package etomica.virial.GUI.components;


import etomica.api.IAtomList;
import etomica.api.IElement;
import etomica.api.ISpecies;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.config.IConformation;
import etomica.potential.P22CLJQ;
import etomica.potential.P2CO2EMP2;
import etomica.potential.P2LennardJones;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateP2CO2EMP2 implements ParameterMapping,Cloneable{
	
	private ISpace space;
	private double[] sigma;
	private double[] epsilon;
	private double[] charge;
	private double bondL;
	//private IConformation conformation;
	
	private int id;
	private static int numberOfInstances = 0;
	
	
	private String[] ComponentParameters  =  {"SIGMA","EPSILON","CHARGE"};
			
	private String[] SharedComponentParameters ={"BONDL"};
	
	private String[][] ParamAndValues; 
	
	private String[] PotentialSites = {"C","O","CO"};
	
	private String[][] ComponentValues = {
			{"2.7570",Double.toString(Kelvin.UNIT.toSim(28.129)),Double.toString(Electron.UNIT.toSim(0.6512))},
			{"3.0330",Double.toString(Kelvin.UNIT.toSim(80.507)),Double.toString((-0.5)*Electron.UNIT.toSim(0.6512))},
			{"2.8921",Double.toString(Kelvin.UNIT.toSim(47.588)),Double.toString(Electron.UNIT.toSim(0.6512)*(-0.5)*Electron.UNIT.toSim(0.6512))}
	};

	private String[] SharedComponentValues = {"1.149"};
	
	public String getPotentialSiteAtIndex(int index) {
		
		return PotentialSites[index];
	}
	
	public boolean IsPotentialSiteMoreThanOne(){
		return true;
	}
	//Potentials references are created as Private members
	private P2CO2EMP2 p2CO2EPM2;
	

	//Constructors for different Instantiations
	
	public P2CO2EMP2 getP2CO2EPM2() {
		return p2CO2EPM2;
	}

	public void setP2CO2EPM2() {
		p2CO2EPM2 = new P2CO2EMP2(this.space);
	}

	public CreateP2CO2EMP2(){
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
			 CreateP2CO2EMP2 cloned = (CreateP2CO2EMP2)super.clone();
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
		return parameterValue;
	}
	
	public String[] getParametersArray() {
		return ComponentParameters;
	}

	@Override
	public String getCustomName() {
		// TODO Auto-generated method stub
		return "EPM2";
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

	
	
}
