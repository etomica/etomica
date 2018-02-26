package etomica.rotation;


public class KappaIntergration {
    public static double theta;
    public double k1, k2, k3;
    public double x, y, z;
    public double q0, q1, q2, q3;
    public static double dt = 1E-2;
    //	public static int numSteps = (int) Math.round(10/dt);
    public static int numSteps = 200;


    public static void main(String[] args) {
        //x y z is the component of angular velocity
        double x = 1.0;
        double y = 0;
        double z = 0;
        double theta = 0;//  0< theta < Pi/2
        //ux uy uz is the ratio of the component of the rotation axis vector
        double ux = 0;
        double uy = 0;
        double uz = 1;
        double norm = Math.sqrt(ux * ux + uy * uy + uz * uz);
        double q1 = Math.sin(theta) * ux / norm;
        double q2 = Math.sin(theta) * uy / norm;
        double q3 = Math.sin(theta) * uz / norm;
        double k1 = k(theta) * q1;
        double k2 = k(theta) * q2;
        double k3 = k(theta) * q3;// 0 <= k_i <= (0.75*Pi)^(1/3)
//	System.out.println("0 " + k1);
//	System.out.println(k2);
//	System.out.println(k3);

        for (int i = 1; i < numSteps; i++) {

            double dk1dt = -k3 * y + k2 * z + B(theta) * x + A(theta) * k1 * (k1 * x + k2 * y + k3 * z);
            double dk2dt = k3 * x - k1 * z + B(theta) * y + A(theta) * k2 * (k1 * x + k2 * y + k3 * z);
            double dk3dt = -k2 * x + k1 * y + B(theta) * z + A(theta) * k3 * (k1 * x + k2 * y + k3 * z);


            k1 += dk1dt * dt;
            k2 += dk2dt * dt;
            k3 += dk3dt * dt;

            if (Math.sqrt(k1 * k1 + k2 * k2 + k3 * k3) > Math.cbrt((0.75 * Math.PI))) {
                double LH = Math.sqrt(k1 * k1 + k2 * k2 + k3 * k3);
                double RH = Math.cbrt((0.75 * Math.PI));
                double ratio = -1 * (2 * RH - LH) / LH;
//    			double ratio = (2*RH - LH)/LH;
//    			x *= -1;
//    			y *= -1;
//    			z *= -1;
                k1 *= ratio;
                k2 *= ratio;
                k3 *= ratio;
            }

            double theta_plus = Math.PI / 2;
            double theta_minus = 0;
            while (true) {
                theta = (theta_minus + theta_plus) / 2;
                if (theta == theta_minus || theta == theta_plus) {
                    break;
                }
                double f = Math.cbrt(0.75 * (2 * theta - Math.sin(2 * theta))) * Math.cbrt(0.75 * (2 * theta - Math.sin(2 * theta))) - (k1 * k1 + k2 * k2 + k3 * k3);
                if (f == 0) break;
                if (f > 0) {
                    theta_plus = theta;
                } else {
                    theta_minus = theta;
                }
            }


            double kNorm = Math.sqrt(k1 * k1 + k2 * k2 + k3 * k3);
            double uxt = theta * k1 / kNorm;
            double uyt = theta * k2 / kNorm;
            double uzt = theta * k3 / kNorm;
            if (i % (numSteps / 100) == 0) {
//					System.out.println(i/(numSteps/100) + " " + uxt + " " + uyt + " " + uzt);
//    				System.out.println(i*dt + " " + theta + " " + uxt);
//        			System.out.println(i/(numSteps/100) + " " +  theta + " " + k1);
//	        		System.out.println(i/(numSteps/100) + " " + theta );

            }

            System.out.println(i * dt + " " + theta);
        }

    }


    public static double k(double theta) {
        if (theta == 0) {
            return 1;
        }
        if (theta == Math.PI / 2) {
            return Math.cbrt((3 * Math.PI)) / Math.cbrt(4);
        }
        double k = Math.cbrt(0.75) * Math.cbrt((2 * theta - Math.sin(2 * theta))) / Math.sin(theta);
        return k;
    }

    public static double g(double theta) {
        if (theta == 0) {
            return 0.2;
        }
        if (theta == Math.PI / 2) {
            return 2 * Math.cbrt(2) / Math.cbrt(9 * Math.PI * Math.PI);
        }
        double g = (-12.0 * theta * Math.cos(theta) + 9 * Math.sin(theta) + Math.sin(3 * theta)) / Math.sin(theta) / Math.sin(theta) / Math.sin(theta) / 2 / Math.cbrt(36) / Math.cbrt((2 * theta - Math.sin(theta)) * (2 * theta - Math.sin(theta)));
        return g;
    }

    public static double A(double theta) {
        if (theta == 0) {
            return 0.2;
        }
        if (theta == Math.PI / 2) {
            return 4.0 * Math.cbrt(4) / 3 / Math.cbrt(3) / Math.PI / Math.cbrt(Math.PI);
        }

        double A = g(theta) / k(theta) / k(theta);
        return A;
    }

    public static double B(double theta) {
        if (theta == 0) {
            return 1;
        }
        if (theta == Math.PI / 2) {
            return 0;
        }
        double B = 1.0 * Math.cos(theta) * k(theta);
        return B;
    }


//  this is for dkdt expression with generized momenta	future usage

//	public static double C(double theta){
//		double C = (Math.sin(theta))*(Math.sin(theta)) * k(theta) * k(theta);
//		return C;
//		}
//
//	public double H(double theta){
//		double H = A(theta) * C(theta) + B(theta);
//		return H;
//	}
//	
//	public double D(double theta){
//		double D = -1*(B(theta)*B(theta)+C(theta))*H(theta);
//		return D;
//	}
//	
//	public double F(double theta){
//		double F = A(theta)*B(theta) -1 ;
//		return F;
//	}

}
