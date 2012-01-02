package etomica.virial.GUI.models;

public enum EnumSiteName {
	
	C("Carbon"),
	O("Oxygen"),
	H("Hydrogen"),
	CH3("Methyl"),
	CH2("Methylene"),
	CH("Methyne"),
	OH("Hydroxyl"),
	aC("AlphaCarbon"),
	aH("AlphaHydrogen"),
	X("Satellite"),
	LJ("Leaf");
	
	
	
	
	private String Site;
	
	public String getSite() {
		return Site;
	}

	public void setSite(String site) {
		Site = site;
	}

	EnumSiteName(String sitename){
		Site = sitename;
	}
	

}
