package etomica.virial.GUI.containers;


public class SuperContainer {
	
	
	private SingleSpeciesParam SSP;
	
	private Species1 S1;
	public Species1 getS1() {
		return S1;
	}

	public void setS1(Species1 s1) {
		S1 = s1;
	}

	
	private DialogBoxPanel DBox;
	
	
	public SuperContainer(){
		
		SSP = new SingleSpeciesParam();
		
		S1 = new Species1();
		
	}

	
	public void setSSP(SingleSpeciesParam sSP) {
		SSP = sSP;
	}

	

	
	public SingleSpeciesParam getSSP() {
		return SSP;
	}

	public void setDBox(DialogBoxPanel dBox) {
		DBox = dBox;
	}

	public DialogBoxPanel getDBox() {
		return DBox;
	}

	
	

}
