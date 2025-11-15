package etomica.virial.simulations.mbnrg;

import java.util.ArrayList;

/**
 * @file energy2b
 * @brief Contains the implementation of the 2b energy calls
 */

public class energy2b {
    double get_2b_energy(String mon1, String mon2, Integer nm, ArrayList<Double> xyz1, ArrayList<Double> xyz2) {
        double energy = 0.0;

        // Order the two monomer names && corresponding xyz
        if (mon2 < mon1) {
            String tmp = mon1;
            mon1 = mon2;
            mon2 = tmp;
            ArrayList<Double> tmp2 = std::move (xyz1);
            xyz1 = std::move (xyz2);
            xyz2 = std::move (tmp2);
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
        } else if (mon1.equals("co2_archive") && mon2.equals("h2o")) {
            x2b_A1B2Z2_C1D2_deg4::x2b_A1B2Z2_C1D2_v1x pot(mon2, mon1);
            energy = pot.eval(xyz2, xyz1, nm);
        } else if ((mon1.equals("co2") || mon1.equals("co2cm5100") || mon1.equals("co2cm595") || mon1.equals("co2cm590") ||
                mon1.equals("co2cm5875") || mon1.equals("co2cm585") || mon1.equals("co2cm580")) && mon2.equals("h2o")) {
            x2b_A1B2Z2_C1D2_deg4::x2b_A1B2Z2_C1D2_v1x pot(mon2, "co2");
            energy = pot.eval(xyz2, xyz1, nm);

        }
        return energy;

    }

}