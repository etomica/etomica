package etomica.virial.GUI.models;

public class SingleSpeciesModel {
	private static final int INTIAL_VALUE_ComponentDropDown = 0;
	private static final double INITIAL_VALUE_SigmaHSRefLJ = 1.5;
	private static final double INITIAL_VALUE_SigmaHSRefCO2LJ= 1.5;
	private static final double INITIAL_VALUE_SigmaHSRefCO2EPM2= 5*1.1491;
	
	private int ComponentDropDown;
	private double SigmaHSRefLJ;
	private double SigmaHSRefCO2LJ;
	private double SigmaHSRefCO2EPM2;
	
	SingleSpeciesModel(){
		reset();
	}
	
	
	public void reset(){
		ComponentDropDown = INTIAL_VALUE_ComponentDropDown;
		SigmaHSRefLJ = INITIAL_VALUE_SigmaHSRefLJ;
		SigmaHSRefCO2LJ = INITIAL_VALUE_SigmaHSRefCO2LJ;
		SigmaHSRefCO2EPM2 = INITIAL_VALUE_SigmaHSRefCO2EPM2;
	}
	
	
	public static int getIntialValueComponentdropdown() {
		return INTIAL_VALUE_ComponentDropDown;
	}

	public static double getInitialValueSigmahsreflj() {
		return INITIAL_VALUE_SigmaHSRefLJ;
	}
	
	public double getSigmaHSRefCO2LJ() {
		return SigmaHSRefCO2LJ;
	}

	public void setSigmaHSRefCO2LJ(double sigmaHSRefCO2LJ) {
		SigmaHSRefCO2LJ = sigmaHSRefCO2LJ;
	}

	public double getSigmaHSRefCO2EPM2() {
		return SigmaHSRefCO2EPM2;
	}

	public void setSigmaHSRefCO2EPM2(double sigmaHSRefCO2EPM2) {
		SigmaHSRefCO2EPM2 = sigmaHSRefCO2EPM2;
	}


	public int getComponentDropDown() {
		return ComponentDropDown;
	}

	public void setComponentDropDown(int componentDropDown) {
		ComponentDropDown = componentDropDown;
	}
	
	public double getSigmaHSRefLJ() {
		return SigmaHSRefLJ;
	}

	public void setSigmaHSRefLJ(double sigmaHSRefLJ) {
		SigmaHSRefLJ = sigmaHSRefLJ;
	}

}
