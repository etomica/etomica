package etomica.rotation;


public class QuaternionIntergration {
    public static double theta;
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
        double q0 = Math.cos(theta);
        double q1 = Math.sin(theta) * ux / norm;
        double q2 = Math.sin(theta) * uy / norm;
        double q3 = Math.sin(theta) * uz / norm;

        for (int i = 1; i < numSteps; i++) {


            double dq0dt = -q1 * x - q2 * y - q3 * z;
            double dq1dt = q0 * x + q2 * z - q3 * y;
            double dq2dt = q0 * y - q1 * z + q3 * x;
            double dq3dt = q0 * z + q1 * y - q2 * x;


            q0 += dq0dt * dt;
            q1 += dq1dt * dt;
            q2 += dq2dt * dt;
            q3 += dq3dt * dt;

            if (q0 < 0) {
                q0 *= -1;
                q1 *= -1;
                q2 *= -1;
                q3 *= -1;
            }

            double qNorm = Math.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
            q0 = q0 / qNorm;
            q1 = q1 / qNorm;
            q2 = q2 / qNorm;
            q3 = q3 / qNorm;
//    		theta = Math.asin(q1);
            theta = Math.acos(q0);


            qNorm = Math.sqrt(q1 * q1 + q2 * q2 + q3 * q3);
            double uxt = theta * q1 / qNorm;
            double uyt = theta * q2 / qNorm;
            double uzt = theta * q3 / qNorm;
            if (i % (numSteps / 100) == 0) {
//				System.out.println(i/(numSteps/100) + " " + uxt + " " + uyt + " " + uzt);
//				System.out.println(i*dt + " " + 2*theta + " " + (1-q0) +" " + q1 + " " + q2 + " " + q3 );
//				System.out.println(i*dt + " " + 2*theta + " " + 2*Math.acos(q0)+ " " +(1-q0)  );
//				System.out.println(i*dt + " " + theta+ " " + uxt);
            }


            System.out.println(i * dt + " " + theta + " " + 2 * Math.acos(q0) + " " + (1 - q0));


        }


    }


}
