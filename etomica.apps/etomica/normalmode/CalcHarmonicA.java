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
        String filename = "normal_modes_LJ_3D_108";
        int D = 3;
        int totalCells = 27;
        double temperature = 1;
        double harmonicFudge = 1.0;
        int basisSize = 4;

        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            D = Integer.parseInt(args[1]);
        }
        if (args.length > 2) {
            totalCells = Integer.parseInt(args[2]);
        }
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
        normalModes.setTemperature(temperature);
        normalModes.setHarmonicFudge(harmonicFudge);
        doit(normalModes, D, temperature, basisSize, totalCells);
    }
    
    public static double doit(NormalModes normalModes, int D, double temperature, int basisSize, int totalCells) {
        double AHarmonic = 0;
        double[][] omega2 = normalModes.getOmegaSquared();
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        for(int i=0; i<omega2.length; i++) {
            for(int j=0; j<omega2[0].length; j++) {
                if (!Double.isInfinite(omega2[i][j])) {
                    AHarmonic += coeffs[i]*Math.log(omega2[i][j]/(2*temperature*Math.PI));
                }
            }
        }

        // and COM correction
        AHarmonic += 0.5*D*Math.log(totalCells*basisSize);

        System.out.println("Harmonic-reference free energy: "+AHarmonic*temperature);
        return AHarmonic*temperature;
    }
}
