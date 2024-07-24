package etomica.potential.UFF;

import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.potential.TruncationFactory;
import etomica.units.*;

public class LJUFF implements IPotential2 {

    public static IPotential2 maketruncated(double sigma, double epsilon, TruncationFactory tf){
        return tf.make(new LJUFF(sigma, epsilon));
    }

    public LJUFF(){
        this(1.0, 1.0);
    }
    public LJUFF(double sigma, double epsilon){
        super();
        setSigma(sigma);
        setEpsilon(epsilon);
       // setSciNew(sci);
    }

    public double u(double r2) {
        sigmaSquared = sigma*sigma;
        s2 = sigmaSquared/r2;
        if( s2 > 150){
            return Double.POSITIVE_INFINITY;
        }
        s6 = s2*s2*s2;
       // System.out.println("in uff" +epsilon*(s6*s6 - (2*s6)));
        return epsilon*(s6*s6 - (2*s6));
        /*System.out.println(" sigma : " + sigma + " epsilon: " +epsilon + " sci " + scale);
        //System.out.println( r + " " + r6 + " " + sigma6 + " " + scale_6 + " " + expo + " " + epsilon +  " " + scale + " ");
       // double u =  epsilon * (6 / scale_6) * expo - epsilon * ( scale / scale_6) * sigma6 / r6;
       //System.out.println(" u " +  u  );

       // return  epsilon * (6 / scale_6) * expo - epsilon * ( scale / scale_6) * sigma6 / r6;*/

    }
    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        sigmaSquared = sigma*sigma;
        s2 = sigmaSquared/r2;
        s6 = s2*s2*s2;
   //     double rdudr = -12*epsilon*((s6*s6)-s6);
       // System.out.println("in uff" + -12*epsilon*((s6*s6)));
        return  -12*epsilon*((s6*s6)-s6);
       // System.out.println(epsilon + " " + scale + " " + sigma6 + " " + sigma6 + " " + r7 + " "+ expo);
       // double dudr = 6 * epsilon * scale * sigma6 / (scale_6 * r7) - ( 6 * epsilon * sigma * expo / (scale_6 * sigma) );
        //System.out.println( " rdudr "+  r *dudr);
       // double dudr = epsilon * ((12*sigma6/r7)  - (12*sigma6*sigma6/(r6*r7)));
       // return r * dudr;
    }

    public void u012add(double r2, double[] u012) {
        sigmaSquared = sigma*sigma;
        s2 = sigmaSquared/r2;
        if(s2 > 300){
            u012[0] = Double.POSITIVE_INFINITY;
            return;
        }
        s6 = s2*s2*s2;
        //u012[0] += epsilon * (s6*s6);
        //u012[1] += -12*epsilon*((s6*s6));
        //u012[2] += -156*epsilon*s6*s6;
        u012[0] += epsilon * (s6*s6 - 2*s6);
        u012[1] += -12*epsilon*((s6*s6)-s6);
        u012[2] += -84*epsilon*s6 -156*epsilon*s6*s6;
        /*u012[0] += epsilon * (6 / (scale - 6)) * expo - epsilon * ( scale / ( scale - 6)) * sigma6 / r6;
        u012[1] += 6 * epsilon * scale * sigma6 / (scale_6 * r7) - ( 6 * epsilon * sigma * expo / (scale_6 * sigma) );
        u012[2] += 6 * epsilon * scale2 * expo / ( scale_6 * sigma2) - 42 * epsilon * scale * sigma6 / ( scale_6 * r8);*/
        /*u012[0] +=  epsilon*((sigma6*sigma6/(r6*r6)) - (2*sigma6/r6));
        u012[0] += epsilon * ((12*sigma6/r7)  - (12*sigma6*sigma6/(r6*r7)));
        u012[0] += epsilon * (12*sigma6/r7) - (12*sigma6*sigma6/(r6*r7));*/
    }

    public double getSigma(){ return sigma;}

    public final void setSigma(double s){
        sigma = s;

    }

    public double getEpsilon(){ return epsilon;}

    public final void setEpsilon(double eps){
        epsilon = eps;

    }
   /* public double getSci(){ return scale;}

    public final void setSciNew(double sci){
        scale = sci;
    }*/
    private double sigma, epsilon;
    private double s2, s6, sigmaSquared;

}
