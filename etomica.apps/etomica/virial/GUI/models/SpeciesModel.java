package etomica.virial.GUI.models;

public class SpeciesModel {
	
	private static final int INTIAL_VALUE_PureOrBinaryRButtonValue = 0;
	
	
	private int PureOrBinaryRButtonValue;
	
	
	

	public static int getIntialValuePureorbinaryrbuttonvalue() {
		return INTIAL_VALUE_PureOrBinaryRButtonValue;
	}

	

	public SpeciesModel(){
		reset();
	}
	
	public void reset(){
		PureOrBinaryRButtonValue = INTIAL_VALUE_PureOrBinaryRButtonValue;
	}

	

	public int getPureBinaryButtonValue() {
		return PureOrBinaryRButtonValue;
	}

	public void setPureBinaryButtonValue(int componentDropDown) {
		PureOrBinaryRButtonValue = componentDropDown;
	}

	

	
}
