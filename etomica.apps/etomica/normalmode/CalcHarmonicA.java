package etomica.normalmode;


/**
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * 
 * @author Andrew Schultz
 */
public class CalcHarmonicA {
    /**
     * @param args filename containing simulation parameters
     * @see CalcHarmonicA.SimOverlapParam
     */
    public static void main(String[] args) {
        
        //set up simulation parameters
        String filename = args[0];
        int D = Integer.parseInt(args[1]);
        int totalCells = Integer.parseInt(args[2]);
        double temperature = 1.0;
        double harmonicFudge = 1.0;
        int basisSize = 1;
        
        if (args.length > 3) {
            harmonicFudge = Double.parseDouble(args[3]);
        }
        if (args.length > 4) {
            temperature = Double.parseDouble(args[4]);
        }
        if (args.length > 5) {
            basisSize = Integer.parseInt(args[5]);
        }

        NormalModesFromFile normalModes = new NormalModesFromFile(filename, D);
        normalModes.setHarmonicFudge(harmonicFudge);
        normalModes.setTemperature(temperature);
        double[][] omega2 = normalModes.getOmegaSquared(null);
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        double AHarmonic = 0;
        for(int i=0; i<omega2.length; i++) {
            for(int j=0; j<omega2[0].length; j++) {
                if (Math.abs(omega2[i][j]) < 1.e15) {
                    AHarmonic += coeffs[i]*Math.log(omega2[i][j]*coeffs[i]/(temperature*Math.PI));
                }
            }
        }
        
        if (totalCells % 2 == 0) {
            AHarmonic += Math.log(Math.pow(2.0, basisSize*D*(totalCells - Math.pow(2,D))/2.0) / Math.pow(totalCells,0.5*D));
        }
        else {
            AHarmonic += Math.log(Math.pow(2.0, basisSize*D*(totalCells - 1)/2.0) / Math.pow(totalCells,0.5*D));
        }
        
        System.out.println("Harmonic-reference free energy: "+AHarmonic*temperature);
    }
}
