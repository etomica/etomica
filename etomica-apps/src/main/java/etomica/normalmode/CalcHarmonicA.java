/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


/**
 * Convenience class to calculate the Helmholtz free energy based on the
 * normal mode description of the system with harmonic springs governing the
 * system's exploration of the modes.
 * 
 * For derivation of harmonic free energy, see
 *   http://rheneas.eng.buffalo.edu/~andrew/harmonicF.pdf
 * 
 * @author Andrew Schultz
 */
public class CalcHarmonicA {

    public static void main(String[] args) {
        
        String filename = "normal_modes_LJ_3D_108";
        int D = 3;
        int numMolecules = 108;
        double temperature = 1;
        double harmonicFudge = 1.0;

        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            D = Integer.parseInt(args[1]);
        }
        if (args.length > 2) {
            numMolecules = Integer.parseInt(args[2]);
        }
        if (args.length > 3) {
            harmonicFudge = Double.parseDouble(args[3]);
        }
        if (args.length > 4) {
            temperature = Double.parseDouble(args[4]);
        }

        NormalModesFromFile normalModes = new NormalModesFromFile(filename, D);
        normalModes.setTemperature(temperature);
        normalModes.setHarmonicFudge(harmonicFudge);
        doit(normalModes, D, temperature, numMolecules);
    }
    
    public static double doit(NormalModes normalModes, int D, double temperature, int numMolecules) {
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

        // FE contribution from COM normal modes
        // see http://rheneas.eng.buffalo.edu/~andrew/harmonicF.pdf
        AHarmonic -= 0.5*D*Math.log(numMolecules);

        System.out.println("Harmonic-reference free energy: "+AHarmonic*temperature);
        return AHarmonic*temperature;
    }
}
