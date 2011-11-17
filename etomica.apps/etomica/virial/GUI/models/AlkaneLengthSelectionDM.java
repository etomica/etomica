package etomica.virial.GUI.models;

public class AlkaneLengthSelectionDM {
	
	private int alkane1Length;
	private int alkane2Length;
	
	AlkaneLengthSelectionDM(){
		reset();
	}

	private void reset() {
		// TODO Auto-generated method stub
		alkane1Length = 4;
		alkane2Length = 4;
	}

	public int getAlkane1Length() {
		return alkane1Length;
	}

	public void setAlkane1Length(int alkane1Length) {
		this.alkane1Length = alkane1Length;
	}

	public int getAlkane2Length() {
		return alkane2Length;
	}

	public void setAlkane2Length(int alkane2Length) {
		this.alkane2Length = alkane2Length;
	}

}
