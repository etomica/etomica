package etomica.virial.GUI.components;

public class ExceptionDoNotRemoveSpecies extends Exception{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public int indexSelectedSpecies;
	
	public ExceptionDoNotRemoveSpecies(int index){
		indexSelectedSpecies = index;
	}

	public int getIndexSelectedSpecies() {
		return indexSelectedSpecies;
	}

	public void setIndexSelectedSpecies(int indexSelectedSpecies) {
		this.indexSelectedSpecies = indexSelectedSpecies;
	}
	
	
}
