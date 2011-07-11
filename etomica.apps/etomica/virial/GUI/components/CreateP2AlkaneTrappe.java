package etomica.virial.GUI.components;

import etomica.api.IAtomList;
import etomica.api.IElement;
import etomica.api.ISpecies;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.config.IConformation;
import etomica.potential.P2CO2TraPPE;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Electron;
import etomica.units.Kelvin;

import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySiepmannSpheres;
import etomica.virial.GUI.models.ParametersDouble;

public class CreateP2AlkaneTrappe implements ParameterMapping,Cloneable{
	private String MoleculeDisplayName;
	private ISpace space;
	private double[] sigma;
	private double[] epsilon;
	
	private double bondL;
	private double theta;
	private double forceconstant;
	private double torsionpotentiala0;
	private double torsionpotentiala1;
	private double torsionpotentiala2;
	private double torsionpotentiala3;
	//private IConformation conformation;
	
	private int id;
	private int AlkaneIndex;
	private static int numberOfInstances = 0;
	
	
	private String[] ComponentParameters  =  {"SIGMA","EPSILON"};
			
	private String[] SharedComponentParameters;
	
	
	private String[][] ParamAndValues; 
	
	private String[] PotentialSites; 
	
	private String[][] ComponentValues; 
	
	private String[] SharedComponentValues;
	
	public String getPotentialSiteAtIndex(int index) {
		
		return PotentialSites[index];
	}
	
	public boolean IsPotentialSiteMoreThanOne(){
		return true;
	}
	//Potentials references are created as Private members
	private P2CO2TraPPE p2CO2TraPPE;
	

	//Constructors for different Instantiations
	
	public P2CO2TraPPE getP2CO2TraPPE() {
		return p2CO2TraPPE;
	}

	public void setP2CO2TraPPE() {
		p2CO2TraPPE = new P2CO2TraPPE(this.space);
	}

	
	public CreateP2AlkaneTrappe(){}
	
	public CreateP2AlkaneTrappe(int index){
		space = Space3D.getInstance();
		switch(index){
		case 1:
			PotentialSites = new String[]{"CH4"};
			ComponentValues = new String[][]{{"3.73",Double.toString(Kelvin.UNIT.toSim(148.0))}};
			MoleculeDisplayName = "Methane-TRAPPE";
			AlkaneIndex = 1;
			break;
			
		case 2:
			PotentialSites = new String[]{"CH3"};
			ComponentValues = new String[][]{{"3.75",Double.toString(Kelvin.UNIT.toSim(98.0))}};
			SharedComponentParameters = new String[]{"NominalbondL"};
			SharedComponentValues = new String[]{"1.54"};
			MoleculeDisplayName = "Ethane-TRAPPE";
			AlkaneIndex = 2;
			break;
		case 3:
			PotentialSites = new String[]{"CH3","CH2"};
			ComponentValues = new String[][]{{"3.75",Double.toString(Kelvin.UNIT.toSim(98.0))},{"3.95",Double.toString(Kelvin.UNIT.toSim(46.0))}};
			SharedComponentParameters = new String[]{"NominalbondL","theta","forceconstant"};
			SharedComponentValues = new String[]{"1.54","114.0","62500"};
			MoleculeDisplayName = "Propane-TRAPPE";
			AlkaneIndex = 3;
			break;
		case 4:
			PotentialSites = new String[]{"CH3","CH2"};
			ComponentValues = new String[][]{{"3.75",Double.toString(Kelvin.UNIT.toSim(98.0))},{"3.95",Double.toString(Kelvin.UNIT.toSim(46.0))}};
			SharedComponentParameters = new String[]{"NominalbondL","theta","forceconstant","c1/kB",
			"c2/kB","c3/kB"};
			SharedComponentValues = new String[]{"1.54","114.0","62500","355.03","-68.19","791.32"};
			MoleculeDisplayName = "Higher n-Alkane TRAPPE";
			AlkaneIndex = 4;
			break;
		}
		sigma = new double[PotentialSites.length];
		epsilon = new double[PotentialSites.length];
		ParamAndValues=setParameterValues();
		//setConformation();
		id=++numberOfInstances;
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
				
				ReturnArray[index][0] = ComponentParameters[j]+PotentialSites[i];
				ReturnArray[index][1] = ComponentValues[i][j];
				index++;
			}
		}
		
		if(SharedComponentParameters != null){
			int NoOfCommonParam = SharedComponentParameters.length;
			for(int k = 0;k<NoOfCommonParam;k++){
				if(SharedComponentParameters[k]=="NominalbondL"){
					setBondL(Double.parseDouble(SharedComponentValues[k]));
				}
				if(SharedComponentParameters[k]=="theta"){
					setTheta(Double.parseDouble(SharedComponentValues[k]));
				}
				if(SharedComponentParameters[k]=="forceconstant"){
					setForceconstant(Double.parseDouble(SharedComponentValues[k]));
				}
			}
		}
		return ReturnArray;
		
		
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

	
    
	
	 public Object clone(){
		 try{
			 CreateP2AlkaneTrappe cloned = ( CreateP2AlkaneTrappe)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
	 }

	
	//Creates the LJAtom Species
	public SpeciesFactory createSpeciesFactory(){
		SpeciesFactory factory = new SpeciesFactorySiepmannSpheres(this.space,this.AlkaneIndex);
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
		}
		
		if(Parameter.toUpperCase().equals(ParametersDouble.NominalbondL.toString())){
			setBondL(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.theta.toString())){
			setTheta(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.forceconstant.toString())){
			setForceconstant(Double.parseDouble(ParameterValue)); 
		}
	}


	public double getTheta() {
		return theta;
	}

	public void setTheta(double theta) {
		this.theta = theta;
	}

	public double getForceconstant() {
		return forceconstant;
	}

	public void setForceconstant(double forceconstant) {
		this.forceconstant = forceconstant;
	}

	public double getTorsionpotentiala0() {
		return torsionpotentiala0;
	}

	public void setTorsionpotentiala0(double torsionpotentiala0) {
		this.torsionpotentiala0 = torsionpotentiala0;
	}

	public double getTorsionpotentiala1() {
		return torsionpotentiala1;
	}

	public void setTorsionpotentiala1(double torsionpotentiala1) {
		this.torsionpotentiala1 = torsionpotentiala1;
	}

	public double getTorsionpotentiala2() {
		return torsionpotentiala2;
	}

	public void setTorsionpotentiala2(double torsionpotentiala2) {
		this.torsionpotentiala2 = torsionpotentiala2;
	}

	public double getTorsionpotentiala3() {
		return torsionpotentiala3;
	}

	public void setTorsionpotentiala3(double torsionpotentiala3) {
		this.torsionpotentiala3 = torsionpotentiala3;
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
		if(Parameter.toUpperCase().equals(ParametersDouble.NominalbondL.toString())){
			Description = ParametersDouble.NominalbondL.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.theta.toString())){
			Description = ParametersDouble.theta.Description();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.forceconstant.toString())){
			Description = ParametersDouble.forceconstant.Description();
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
		if(Parameter.toUpperCase().equals(ParametersDouble.NominalbondL.toString())){
			parameterValue = getBondL();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.theta.toString())){
			parameterValue = getTheta();
		}
		if(Parameter.toUpperCase().equals(ParametersDouble.forceconstant.toString())){
			parameterValue = getTheta();
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

	
	public static void main(String[] args){
		CreateP2AlkaneTrappe lj = new CreateP2AlkaneTrappe(3);
		for(int j=0;j<lj.ComponentParameters.length*lj.PotentialSites.length;j++){
			
			System.out.println(lj.ParamAndValues[j][0]+"\n");
			System.out.println(lj.ParamAndValues[j][1]+"\n");
		}
	}
}
