package etomica.virial.simulations.mbnrg;

import java.util.ArrayList;

/**
 * @file energy2b
 * @brief Contains the implementation of the 2b energy calls
 * unittest-energy2b.cpp has sample xyz1 and xyz2
 */

public class energy2b {
    Double get_2b_energy(String mon1, String mon2, Integer nm, Double[] xyz1, Double[] xyz2) {
        Double energy = 0.0;

        // Order the two monomer names && corresponding xyz
        if (mon2.compareTo(mon1) < 0) {
            String tmp = mon1;
            mon1 = mon2;
            mon2 = tmp;
            Double[] tmp2 = xyz1;
            xyz1 = xyz2;
            xyz2 = tmp2;
        }
        if (mon1.equals("co2_archive") && mon2.equals("co2_archive")) {
            A1B2_2b pot = new A1B2_2b(mon1, mon2);
            energy = pot.eval(xyz1, xyz2, nm);
        } else if (mon1.equals("co2") && mon2.equals("co2")) {
            A1B2_2b pot = new A1B2_2b(mon1, mon2);
            energy = pot.eval(xyz1, xyz2, nm);
        } else if (mon1.equals("co2cm5100") && mon2.equals("co2cm5100")) {
            A1B2_2b pot = new A1B2_2b(mon1, mon2);
            energy = pot.eval(xyz1, xyz2, nm);
        } else if (mon1.equals("co2cm595") && mon2.equals("co2cm595")) {
            A1B2_2b pot = new A1B2_2b(mon1, mon2);
            energy = pot.eval(xyz1, xyz2, nm);
        } else if (mon1.equals("co2cm5875") && mon2.equals("co2cm5875")) {
            A1B2_2b pot = new A1B2_2b(mon1, mon2);
            energy = pot.eval(xyz1, xyz2, nm);
        } else if (mon1.equals("co2cm585") && mon2.equals("co2cm585")) {
            A1B2_2b pot = new A1B2_2b(mon1, mon2);
            energy = pot.eval(xyz1, xyz2, nm);
        } else if (mon1.equals("co2cm580") && mon2.equals("co2cm580")) {
            A1B2_2b pot = new A1B2_2b(mon1, mon2);
            energy = pot.eval(xyz1, xyz2, nm);
        } else if (mon1.equals("co2cm590") && mon2.equals("co2cm590")) {
            A1B2_2b pot = new A1B2_2b(mon1, mon2);
            energy = pot.eval(xyz1, xyz2, nm);
        }

        return energy;

    }

}