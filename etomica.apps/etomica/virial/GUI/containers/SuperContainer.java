package etomica.virial.GUI.containers;
public class SuperContainer {
	
	private BinaryMixtureParam BMP;
	private SingleSpeciesParam SSP;
	private Species S;
	private Species1 S1;
	public Species1 getS1() {
		return S1;
	}

	public void setS1(Species1 s1) {
		S1 = s1;
	}

	private RunParam RP;
	private DialogBoxPanel DBox;
	private LJDefaultViewPanel LJD;
	
	public SuperContainer(){
		BMP = new BinaryMixtureParam();
		SSP = new SingleSpeciesParam();
		S = new Species();
		S1 = new Species1();
		RP = new RunParam();
		DBox = new DialogBoxPanel("");
		LJD = new LJDefaultViewPanel();
	}

	public LJDefaultViewPanel getLJD() {
		return LJD;
	}

	public void setLJD(LJDefaultViewPanel lJD) {
		LJD = lJD;
	}

	public void setBMP(BinaryMixtureParam bMP) {
		BMP = bMP;
	}

	public void setSSP(SingleSpeciesParam sSP) {
		SSP = sSP;
	}

	public void setS(Species s) {
		S = s;
	}

	public void setRP(RunParam rP) {
		RP = rP;
	}

	public BinaryMixtureParam getBMP() {
		return BMP;
	}

	public SingleSpeciesParam getSSP() {
		return SSP;
	}

	public Species getS() {
		return S;
	}

	public RunParam getRP() {
		return RP;
	}

	public void setDBox(DialogBoxPanel dBox) {
		DBox = dBox;
	}

	public DialogBoxPanel getDBox() {
		return DBox;
	}

	
	

}
