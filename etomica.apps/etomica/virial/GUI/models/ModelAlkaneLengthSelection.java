package etomica.virial.GUI.models;

public class ModelAlkaneLengthSelection {
	
	private int[] alkaneLength;
	
	
	public ModelAlkaneLengthSelection(){
		reset();
	}

	public void reset() {
		// TODO Auto-generated method stub
		alkaneLength = new int[10];
		for(int i=0;i<10;i++){
			alkaneLength[i]=4;
		}
		
		
	}

	public int getAlkaneLength(int index) {
		return alkaneLength[index];
	}

	public void setAlkaneLength(int alkane1Length,int index) {
		this.alkaneLength[index] = alkane1Length;
	}

	

}
