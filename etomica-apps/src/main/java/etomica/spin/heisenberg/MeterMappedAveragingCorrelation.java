package etomica.spin.heisenberg;

import etomica.space.Vector;

public class MeterMappedAveragingCorrelation {

    //Computes derivatives of thetaDot_k with respect to eta_k and theta0_k for all spins
    //The meterMeanField passed as an argument should be set up by invoking its getData() method before calling this
    //dtdotkdetak and dtdotkdt0k should be instantiated and filled with Vector2D instances prior to calling
    public static void computeTdotDerivs(double beta, MeterMeanFieldFasterer meterMeanField, Vector[] dtdotkdetak, Vector[] dtdotkdt0k) {
        int N = dtdotkdetak.length;
        for(int k=0; k < N; k++) {
            double dtheta = meterMeanField.getDtheta()[k];
            double cosdtheta = meterMeanField.getCosdtheta()[k];
            double sindtheta = meterMeanField.getSindtheta()[k];
            double t0 = meterMeanField.getTheta0()[k];
            double cost0 = meterMeanField.getCost0()[k];
            double sint0 = meterMeanField.getSint0()[k];
            double I1I0 = meterMeanField.getI1I0()[k];
            double I2etc = meterMeanField.getI2etc()[k];
            double bh = beta * meterMeanField.getEta()[k];
            double pInv = Math.exp(-bh * cosdtheta);

            // Cmn is integral of (cos(t-t0))^m (sin(t-t0))^n Exp[bh cos(t-t0)] for t = t0 to theta
            // Cmnc is integral from 0 to theta
            // Cmns is integral from Pi/2 to theta
            // integrals are computed using functions that return integral of cosx^m sinx^n Exp[bh cosx] for x = 0 to [first argument]
            double C00 = FunctionCosIntegral.cosInt(dtheta, bh);
            double C00c = C00 - FunctionCosIntegral.cosInt(0 - t0, bh);
            double C00s = C00 - FunctionCosIntegral.cosInt(Math.PI/2. - t0, bh);

            double C10 = FunctionCosIntegral.coscosInt(dtheta, bh);
            double C10c = C10 - FunctionCosIntegral.coscosInt(0 - t0, bh);
            double C10s = C10 - FunctionCosIntegral.coscosInt(Math.PI/2. - t0, bh);

            double C20 = FunctionCosIntegral.cos2cosInt(dtheta, bh);
            double C20c = C20 - FunctionCosIntegral.cos2cosInt(0 - t0, bh);
            double C20s = C20 - FunctionCosIntegral.cos2cosInt(Math.PI/2. - t0, bh);

            double C11 = FunctionCosIntegral.sincoscosInt(dtheta, bh);
            double C11c = C11 - FunctionCosIntegral.sincoscosInt(0 - t0, bh);
            double C11s = C11 - FunctionCosIntegral.sincoscosInt(Math.PI/2. - t0, bh);

            double C01 = FunctionCosIntegral.sincosInt(dtheta, bh);
            double C01c = C01 - FunctionCosIntegral.sincosInt(0 - t0, bh);
            double C01s = C01 - FunctionCosIntegral.sincosInt(Math.PI/2. - t0, bh);

            double C02c = C00c - C20c;//integral of sin^2 is integral of 1 - cos^2
            double C02s = C00s - C20s;

            //compute d(thetaDot_k)/d(eta_k)
            double c = -I2etc * C00c * cost0;
            double s = -I2etc * C00s * sint0;

            c += cost0 * C20c - sint0 * C11c;
            s += cost0 * C11s + sint0 * C20s;

            c -= cost0 * I1I0 * C10c;
            s -= sint0 * I1I0 * C10s;

            c -= cosdtheta * (cost0 * C10c - sint0 * C01c);
            s -= cosdtheta * (cost0 * C01s + sint0 * C10s);

            c += cosdtheta * I1I0 * cost0 * C00c;
            s += cosdtheta * I1I0 * sint0 * C00s;

            dtdotkdetak[k].setX(0, c);
            dtdotkdetak[k].setX(1, s);
            dtdotkdetak[k].TE(-beta * pInv);

            //compute d(thetaDot_k)/d(theta0_k)
            c = + I1I0 * sint0 * C00c;
            s = - I1I0 * cost0 * C00s;

            c += bh * (cost0 * C11c - sint0 * C02c);
            s += bh * (cost0 * C02s + sint0 * C11s);

            c -= bh * cost0 * I1I0 * C01c;
            s -= bh * sint0 * I1I0 * C01s;

            c -= bh * sindtheta * (cost0 * C10c - sint0 * C01c);
            s -= bh * sindtheta * (cost0 * C01s + sint0 * C10s);

            c += bh * sindtheta * I1I0 * cost0 * C00c;
            s += bh * sindtheta * I1I0 * sint0 * C00s;

            dtdotkdt0k[k].setX(0, c);
            dtdotkdt0k[k].setX(1, s);
            dtdotkdt0k[k].TE(-pInv);
        }
    }
}
