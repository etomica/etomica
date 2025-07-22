package etomica.potential.COMPASS;

import etomica.potential.IPotential2;
import etomica.potential.TruncationFactory;

public class LJCOMPASS implements IPotential2 {
    public static IPotential2 maketruncated(double sigma, double epsilon, TruncationFactory tf){
        return tf.make(new LJCOMPASS(sigma, epsilon));
    }
    public LJCOMPASS(){
        this(1.0, 1.0);
    }
    public LJCOMPASS(double sigma, double epsilon){
        super();
        setSigmaNew(sigma);
        setEpsilonNew(epsilon);
        // setSciNew(sci);
    }
    public double u(double r2) {
        double sigmaSquared = sigma*sigma;
        double s2 = sigmaSquared/r2;
        double sr = Math.sqrt(sigmaSquared);
        if( s2 > 150){
            return Double.POSITIVE_INFINITY;
        }
        double s6 = s2*s2*s2;
        double s9 = s6*s2*sr;
        double u = epsilon*(2*s9-3*s6);
        return u;
        /*System.out.println(" sigma : " + sigma + " epsilon: " +epsilon + " sci " + scale);
        //System.out.println( r + " " + r6 + " " + sigma6 + " " + scale_6 + " " + expo + " " + epsilon +  " " + scale + " ");
       // double u =  epsilon * (6 / scale_6) * expo - epsilon * ( scale / scale_6) * sigma6 / r6;
       //System.out.println(" u " +  u  );

       // return  epsilon * (6 / scale_6) * expo - epsilon * ( scale / scale_6) * sigma6 / r6;*/

    }

    public double getSigma(){ return sigma;}

    public final void setSigmaNew(double s){
        sigma = s;

    }

    public double getEpsilon(){ return epsilon;}

    public final void setEpsilonNew(double eps){
        epsilon = eps;

    }
    /* public double getSci(){ return scale;}

     public final void setSciNew(double sci){
         scale = sci;
     }*/
    private double sigma, epsilon;
}
