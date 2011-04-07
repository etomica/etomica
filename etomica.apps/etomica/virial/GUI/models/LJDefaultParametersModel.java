package etomica.virial.GUI.models;

import etomica.util.Function.Sqrt;

public class LJDefaultParametersModel {
	
	private static final double INITIAL_SigmaValue = 1.0;
	private static final double INITIAL_EpsilonValue = 1.0;
	private static final double INITIAL_MomentValue = 1.0;
	private static final double INITIAL_BondLenthValue = 1.0;
	
	private static final double INITIAL_SigmaAAValueMixture1 = 1.0;
	private static final double INITIAL_SigmaBBValueMixture1= 1.0;
	private static final double INITIAL_SigmaABValueMixture1 = 1.0;
	
	private static final double INITIAL_EpsilonValueAAMixture1 = 1.0;
	private static final double INITIAL_EpsilonValueABMixture1 = 1.0;
	private static final double INITIAL_EpsilonValueBBMixture1 = 1.0;
	
	private static final double INITIAL_EpsilonValueABMixture2 = 0.75;
	
	private static final double INITIAL_SigmaABValueMixture3= 0.864;
	private static final double INITIAL_EpsilonValueABMixture3 = 0.773;
	private static final double INITIAL_SigmaBBValueMixture3= 0.768;
	private static final double INITIAL_EpsilonValueBBMixture3 = 0.597;
	
	private static final double INITIAL_EpsilonValueABMixture4 = Sqrt.INSTANCE(0.5);
	private static final double INITIAL_EpsilonValueBBMixture4 = Sqrt.INSTANCE(0.5);
	
	private double SigmaValue;
	private double EpsilonValue;
	private double MomentValue;
	private double BondLenthValue;
	
	private double SigmaAAValueMixture1;
	private double SigmaBBValueMixture1;
	private double SigmaABValueMixture1;
	
	private double  EpsilonValueAAMixture1;
	private double EpsilonValueABMixture1;
	private double EpsilonValueBBMixture1;
	
	private double EpsilonValueABMixture2;
	
	private double SigmaABValueMixture3;
	private double EpsilonValueABMixture3;
	private double SigmaBBValueMixture3;
	private double EpsilonValueBBMixture3;
	
	private double EpsilonValueABMixture4;
	private double EpsilonValueBBMixture4;
	
	
	LJDefaultParametersModel(){
		reset();
	}

	public void reset() {
		// TODO Auto-generated method stub
		this.SigmaValue = INITIAL_SigmaValue;
		this.EpsilonValue = INITIAL_EpsilonValue;
		this.MomentValue = INITIAL_MomentValue;
		this.BondLenthValue = INITIAL_BondLenthValue;
		
		this.SigmaAAValueMixture1 = INITIAL_SigmaAAValueMixture1;
		this.SigmaABValueMixture1 = INITIAL_SigmaABValueMixture1;
		this.SigmaBBValueMixture1 = INITIAL_SigmaBBValueMixture1;
		
		this.EpsilonValueAAMixture1 = INITIAL_EpsilonValueAAMixture1;
		this.EpsilonValueABMixture1 = INITIAL_EpsilonValueABMixture1;
		this.EpsilonValueBBMixture1 = INITIAL_EpsilonValueBBMixture1;
		
		this.EpsilonValueABMixture2 = INITIAL_EpsilonValueABMixture2;
		
		this.EpsilonValueABMixture3 = INITIAL_EpsilonValueABMixture3;
		this.EpsilonValueBBMixture3 = INITIAL_EpsilonValueBBMixture3;
		this.SigmaABValueMixture3 = INITIAL_SigmaABValueMixture3;
		this.SigmaBBValueMixture3 = INITIAL_SigmaBBValueMixture3;
		
		this.EpsilonValueABMixture4 = INITIAL_EpsilonValueABMixture4;
		this.EpsilonValueBBMixture4 = INITIAL_EpsilonValueBBMixture4;
		
		
		
	}

	public double getSigmaValue() {
		return SigmaValue;
	}

	public void setSigmaValue(double sigmaValue) {
		SigmaValue = sigmaValue;
	}

	public double getEpsilonValue() {
		return EpsilonValue;
	}

	public void setEpsilonValue(double epsilonValue) {
		EpsilonValue = epsilonValue;
	}

	public double getMomentValue() {
		return MomentValue;
	}

	public void setMomentValue(double momentValue) {
		MomentValue = momentValue;
	}

	public double getBondLenthValue() {
		return BondLenthValue;
	}

	public void setBondLenthValue(double bondLenthValue) {
		BondLenthValue = bondLenthValue;
	}

	public double getSigmaAAValueMixture1() {
		return SigmaAAValueMixture1;
	}

	public void setSigmaAAValueMixture1(double sigmaAAValueMixture1) {
		SigmaAAValueMixture1 = sigmaAAValueMixture1;
	}

	public double getSigmaBBValueMixture1() {
		return SigmaBBValueMixture1;
	}

	public void setSigmaBBValueMixture1(double sigmaBBValueMixture1) {
		SigmaBBValueMixture1 = sigmaBBValueMixture1;
	}

	public double getSigmaABValueMixture1() {
		return SigmaABValueMixture1;
	}

	public void setSigmaABValueMixture1(double sigmaABValueMixture1) {
		SigmaABValueMixture1 = sigmaABValueMixture1;
	}

	public double getEpsilonValueAAMixture1() {
		return EpsilonValueAAMixture1;
	}

	public void setEpsilonValueAAMixture1(double epsilonValueAAMixture1) {
		EpsilonValueAAMixture1 = epsilonValueAAMixture1;
	}

	public double getEpsilonValueABMixture1() {
		return EpsilonValueABMixture1;
	}

	public void setEpsilonValueABMixture1(double epsilonValueABMixture1) {
		EpsilonValueABMixture1 = epsilonValueABMixture1;
	}

	public double getEpsilonValueBBMixture1() {
		return EpsilonValueBBMixture1;
	}

	public void setEpsilonValueBBMixture1(double epsilonValueBBMixture1) {
		EpsilonValueBBMixture1 = epsilonValueBBMixture1;
	}

	public double getEpsilonValueABMixture2() {
		return EpsilonValueABMixture2;
	}

	public void setEpsilonValueABMixture2(double epsilonValueABMixture2) {
		EpsilonValueABMixture2 = epsilonValueABMixture2;
	}

	public double getSigmaABValueMixture3() {
		return SigmaABValueMixture3;
	}

	public void setSigmaABValueMixture3(double sigmaABValueMixture3) {
		SigmaABValueMixture3 = sigmaABValueMixture3;
	}

	public double getEpsilonValueABMixture3() {
		return EpsilonValueABMixture3;
	}

	public void setEpsilonValueABMixture3(double epsilonValueABMixture3) {
		EpsilonValueABMixture3 = epsilonValueABMixture3;
	}

	public double getSigmaBBValueMixture3() {
		return SigmaBBValueMixture3;
	}

	public void setSigmaBBValueMixture3(double sigmaBBValueMixture3) {
		SigmaBBValueMixture3 = sigmaBBValueMixture3;
	}

	public double getEpsilonValueBBMixture3() {
		return EpsilonValueBBMixture3;
	}

	public void setEpsilonValueBBMixture3(double epsilonValueBBMixture3) {
		EpsilonValueBBMixture3 = epsilonValueBBMixture3;
	}

	public double getEpsilonValueABMixture4() {
		return EpsilonValueABMixture4;
	}

	public void setEpsilonValueABMixture4(double epsilonValueABMixture4) {
		EpsilonValueABMixture4 = epsilonValueABMixture4;
	}

	public double getEpsilonValueBBMixture4() {
		return EpsilonValueBBMixture4;
	}

	public void setEpsilonValueBBMixture4(double epsilonValueBBMixture4) {
		EpsilonValueBBMixture4 = epsilonValueBBMixture4;
	}
	
}
