package etomica.potential.UFF;

import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.potential.TruncationFactory;

public class LJUFF extends P2LennardJones {

    public static LJUFF maketruncated(double sigma, double epsilon, double sci, TruncationFactory tf){
        return (LJUFF) tf.make(new LJUFF(sigma, epsilon, sci));
    }

    public LJUFF(){
        this(1.0, 1.0, 12);
    }
    public LJUFF(double sigma, double epsilon, double sci){
        setSigmaNew(sigma);
        setEpsilonNew(epsilon);
        setSciNew(sci);
    }

    public double u(double r2) {
        double r = Math.sqrt(r2);
        double r6 = r2*r2*r2;
        double epsilon6 = Math.pow(epsilon,6);
        double u = epsilon * (6/(scale-6)) * ((Math.exp(scale-(scale*r/sigma)))-(epsilon6/r6));
        return u;
    }

    public double udu(double r2) {
        double r = Math.sqrt(r2);
        double r6 = r2*r2*r2;
        double r7 = r6*r;
        double epsilon6 = Math.pow(epsilon,6);
        u = epsilon * (6/(scale-6)) * ((Math.exp(scale-(scale*r/sigma)))-(epsilon6/r6));
        double du = (epsilon * (6/scale-6))*((6*epsilon6/(r7))-(scale*Math.exp(scale-(scale*r/sigma))));
        return u*du;
    }

   /* public double u2du2 (double r2){
        double r = Math.sqrt(r2);
        double r6 = r2*r2*r2;
        double r7 = r6*r;
        double scale2 = scale * scale;
        double epsilon6 = Math.pow(epsilon,6);
        u = u(r2);
        double u2 = epsilon * (6/scale - 6) * ((-scale2* Math.exp(scale - (scale*r/epsilon))/epsilon)-(6*r6/(r*epsilon6)));
        return u2;
    }*/
    public double getSigma(){ return sigma;}

    public final void setSigmaNew(double s){
        sigma = s;
    }

    public double getEpsilon(){ return epsilon;}

    public final void setEpsilonNew(double eps){
        epsilon = eps;

    }
    public double getSci(){ return scale;}

    public final void setSciNew(double sci){
        scale = sci;
    }

    private double sigma, epsilon, scale;
    private double u, u2;

}
