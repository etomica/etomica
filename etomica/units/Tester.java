package etomica.units;

public class Tester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		UnitGraphics gUI = new UnitGraphics();
		Dimension targetDimension = Force.DIMENSION;
		gUI.startWithDim(targetDimension);
	}
}