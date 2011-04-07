package etomica.virial.GUI.models;

public class SuperModel {

	private SpeciesModel sm;
	private SingleSpeciesModel ssm;
	private RunParametersModel rpm;
	private LJDefaultParametersModel LJD;
	
	public SuperModel(){
		this.sm = new SpeciesModel();
		this.ssm = new SingleSpeciesModel();
		this.rpm = new RunParametersModel();
		this.LJD = new LJDefaultParametersModel();
	}

	public LJDefaultParametersModel getLJD() {
		return LJD;
	}

	public void setLJD(LJDefaultParametersModel lJD) {
		LJD = lJD;
	}

	public SpeciesModel getSm() {
		return sm;
	}

	

	public SingleSpeciesModel getSsm() {
		return ssm;
	}


	public RunParametersModel getRpm() {
		return rpm;
	}

	
	
}
