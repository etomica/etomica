/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import etomica.atom.AtomType;
import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.virial.GUI.components.SimpleElementForSimilarSpecies;
import etomica.virial.SpeciesAlkane;
import etomica.virial.SpeciesFactory;

public class MolecularModelTRAPPE_SpeciesAlkane implements IMolecularModel_SpeciesFactory,Cloneable{
	private String MoleculeDisplayName;
	private Space space;
	private double[] sigma;
	private double[] epsilon;
	


	private double sigmaHSRef;
	private int NSpheres;
	
	
	public int getNSpheres() {
		return NSpheres;
	}

	public void setNSpheres(int nSpheres) {
		NSpheres = nSpheres;
	}

	private double bondL;
	private double theta;
	private double forceconstant;
	private double torsionpotentiala0;
	private double torsionpotentiala1;
	private double torsionpotentiala2;
	private double torsionpotentiala3;
	//private IConformation conformation;
	private SimpleElementForSimilarSpecies simpleElementCreator;
	
	private int id;
	private int AlkaneIndex;
	private static int numberOfInstances = 0;
		
	private String[] ComponentParameters  =  {"SIGMA","EPSILON"};
			
	private String[] SharedComponentParameters;
	
	private String[] SimEnvParameters = {"SIGMAHSREF","NUMBER"};
	
	private String[] SimEnvValues = {"1.5","4"};
	
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
	
	
	public MolecularModelTRAPPE_SpeciesAlkane(){}
	
	public MolecularModelTRAPPE_SpeciesAlkane(int index){
		space = Space3D.getInstance();
		switch(index){
		case 1:
			PotentialSites = new String[]{"CH3"};
			ComponentValues = new String[][]{{"3.73",Double.toString(Kelvin.UNIT.toSim(148.0))}};
			MoleculeDisplayName = "Methane-TRAPPE";
			AlkaneIndex = 1;
			setNSpheres(AlkaneIndex);
			break;
			
		case 2:
			PotentialSites = new String[]{"CH3"};
			ComponentValues = new String[][]{{"3.75",Double.toString(Kelvin.UNIT.toSim(98.0))}};
			SharedComponentParameters = new String[]{"NominalbondL"};
			SharedComponentValues = new String[]{"1.54"};
			MoleculeDisplayName = "Ethane-TRAPPE";
			AlkaneIndex = 2;
			setNSpheres(AlkaneIndex);
			break;
		case 3:
			PotentialSites = new String[]{"CH3","CH2"};
			ComponentValues = new String[][]{{"3.75",Double.toString(Kelvin.UNIT.toSim(98.0))},{"3.95",Double.toString(Kelvin.UNIT.toSim(46.0))}};
			SharedComponentParameters = new String[]{"NominalbondL","theta","forceconstant"};
			SharedComponentValues = new String[]{"1.54","114.0","62500"};
			MoleculeDisplayName = "Propane-TRAPPE";
			AlkaneIndex = 3;
			setNSpheres(AlkaneIndex);
			break;
		case 4:
			PotentialSites = new String[]{"CH3","CH2"};
			ComponentValues = new String[][]{{"3.75",Double.toString(Kelvin.UNIT.toSim(98.0))},{"3.95",Double.toString(Kelvin.UNIT.toSim(46.0))}};
			SharedComponentParameters = new String[]{"NominalbondL","theta","forceconstant","c1/kB",
			"c2/kB","c3/kB"};
			SharedComponentValues = new String[]{"1.54","114.0","62500","355.03","-68.19","791.32"};
			MoleculeDisplayName = "Higher n-Alkane TRAPPE";
			AlkaneIndex = 0;
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
		
		int NoOfSimEnvParam = 2;
		for(int l = 0;l<NoOfSimEnvParam;l++){
			
			
			if(SimEnvParameters[l]=="SIGMAHSREF"){
				setSigmaHSRef(Double.parseDouble(SimEnvValues[l]));
			}
			
			if(SimEnvParameters[l]=="NUMBER"){
				if(NSpheres > 4){
					setNSpheres(Integer.parseInt(SimEnvValues[l]));}
			}
		}
		return ReturnArray;
		
		
	}

	

public double getSigmaHSRef() {
	return sigmaHSRef;
}

public void setSigmaHSRef(double sigmaHSRef) {
	this.sigmaHSRef = sigmaHSRef;
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
			 MolecularModelTRAPPE_SpeciesAlkane cloned = ( MolecularModelTRAPPE_SpeciesAlkane)super.clone();
			 return cloned;
		  }
		  catch(CloneNotSupportedException e){
		     System.out.println(e);
		     return null;
		   }
	 }

	
	//Creates the LJAtom Species
	public ISpecies createSpecies(){
		SpeciesGeneral species;
		simpleElementCreator = SimpleElementForSimilarSpecies.getInstance();
		if(AlkaneIndex != 0){
			species = SpeciesAlkane.create(this.AlkaneIndex, AtomType.element(simpleElementCreator.getCH3element()), AtomType.element(simpleElementCreator.getCH2element()));
		}
		else{
			String number = Double.toString(this.getDoubleDefaultParameters("NUMBER"));
			String[] IntSteps= number.split("\\.");
			species = SpeciesAlkane.create(Integer.parseInt(IntSteps[0]), AtomType.element(simpleElementCreator.getCH3element()), AtomType.element(simpleElementCreator.getCH2element()));
		}
		return species;
	}

	public int getAlkaneIndex() {
		return AlkaneIndex;
	}

	public void setAlkaneIndex(int alkaneIndex) {
		AlkaneIndex = alkaneIndex;
	}

	public int getParameterCount() {
		return 2;
	}


	public void setParameter(String Parameter, String ParameterValue) {
		// TODO Auto-generated method stub
		
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMA.toString()+PotentialSites[i])){
				setSigma(Double.parseDouble(ParameterValue),i); 
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.EPSILON.toString()+PotentialSites[i])){
				setEpsilon(Double.parseDouble(ParameterValue),i); 
			}
		}
		
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.NominalbondL.toString())){
			setBondL(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.theta.toString())){
			setTheta(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.forceconstant.toString())){
			setForceconstant(Double.parseDouble(ParameterValue)); 
		}
		
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMAHSREF.toString())){
			setSigmaHSRef(Double.parseDouble(ParameterValue)); 
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.NUMBER.toString())){
			setNSpheres(Integer.parseInt(ParameterValue)); 
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
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMA.toString()+PotentialSites[i])){
				Description = EnumPotentialParamDescription.SIGMA.description();
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.EPSILON.toString()+PotentialSites[i])){
				Description = EnumPotentialParamDescription.EPSILON.description();
			}

			
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.NominalbondL.toString())){
			Description = EnumPotentialParamDescription.NominalbondL.description();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.theta.toString())){
			Description = EnumPotentialParamDescription.theta.description();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.forceconstant.toString())){
			Description = EnumPotentialParamDescription.forceconstant.description();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.TEMPERATURE.toString())){
			Description = EnumPotentialParamDescription.TEMPERATURE.description();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.STEPS.toString())){
			Description = EnumPotentialParamDescription.STEPS.description();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMAHSREF.toString())){
			Description = EnumPotentialParamDescription.SIGMAHSREF.description();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.NUMBER.toString())){
			Description = EnumPotentialParamDescription.NUMBER.description();
		}
		return Description;
	}


	public Double getDoubleDefaultParameters(String Parameter) {
		// TODO Auto-generated method stub
		Double parameterValue = null;
		
		for(int i=0;i<PotentialSites.length;i++){
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMA.toString()+PotentialSites[i])){
				parameterValue = getSigma(i);
			}
			if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.EPSILON.toString()+PotentialSites[i])){
				parameterValue = getEpsilon(i);
			}
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.NominalbondL.toString())){
			parameterValue = getBondL();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.theta.toString())){
			parameterValue = getTheta();
		}
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.forceconstant.toString())){
			parameterValue = getTheta();
		}
		
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.SIGMAHSREF.toString())){
			parameterValue = getSigmaHSRef();
		}
		
		
		if(Parameter.toUpperCase().equals(EnumPotentialParamDescription.NUMBER.toString())){
			parameterValue = (double) getNSpheres();
		}
		return parameterValue;
	}
	
	public String[] getParametersArray() {
		return ComponentParameters;
	}

	@Override
	public String getMolecularModelDisplayName() {
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
		MolecularModelTRAPPE_SpeciesAlkane lj = new MolecularModelTRAPPE_SpeciesAlkane(3);
		for(int j=0;j<lj.ComponentParameters.length*lj.PotentialSites.length;j++){
			
			System.out.println(lj.ParamAndValues[j][0]+"\n");
			System.out.println(lj.ParamAndValues[j][1]+"\n");
		}
	}

	@Override
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
		
		return false;
	}

	@Override
	public String getNonBondedInteractionModel() {
		// TODO Auto-generated method stub
		return "LennardJones";
	}
	

}
