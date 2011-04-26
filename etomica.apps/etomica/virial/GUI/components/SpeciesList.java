package etomica.virial.GUI.components;




public class SpeciesList {
	
	


	private static int Index=0;
	private int id;
	private Object[] SpeciesL;
	
	public SpeciesList(){
		id = Index;
		SpeciesL = new Object[8];
	}
	
	public void addSpecies(Object object){
		SpeciesL[id] = object;
		id++;
	}
	
	
	public void removeSpecies(){
		SpeciesL[id] = null;
		id--;
	}
	
	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}
	
	public Object getObject(int Index){
		return SpeciesL[Index];
	}
	
	public void displayList(){
		for (int i=0;i<SpeciesL.length;i++){
			if(SpeciesL[i] != null){
				System.out.println(SpeciesL[i]);
			}
		}
	}
}
