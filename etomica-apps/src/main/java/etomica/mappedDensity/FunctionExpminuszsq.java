package etomica.mappedDensity;

import etomica.math.function.FunctionDifferentiable;

public class FunctionExpminuszsq implements FunctionDifferentiable {

    private final double L;
    private final double w;
    private final double aa0;
    private final double aa1;
    private final double aa2;
    private final double aa3;
    private final double aa4;
    private final double aa5;
    private final double aa6;
    private final double aa7;
    private final double aa8;
    private final double bb1;
    private final double bb2;
    private final double bb3;
    private final double bb4;
    private final double bb5;
    private final double bb6;
    private final double bb7;
    private final double bb8;

    public FunctionExpminuszsq(double L,double aa0,double aa1,double aa2,double aa3,double aa4,double aa5,double aa6,double aa7,double aa8,double bb1,double bb2,double bb3,double bb4,double bb5,double bb6,double bb7,double bb8,double w) {
        this.L = L;
        this.w = w;
        this.aa0 = aa0;
        this.aa1 = aa1;
        this.aa2 = aa2;
        this.aa3 = aa3;
        this.aa4 = aa4;
        this.aa5 = aa5;
        this.aa6 = aa6;
        this.aa7 = aa7;
        this.aa8 = aa8;
        this.bb1 = bb1;
        this.bb2 = bb2;
        this.bb3 = bb3;
        this.bb4 = bb4;
        this.bb5 = bb5;
        this.bb6 = bb6;
        this.bb7 = bb7;
        this.bb8 = bb8; }

    @Override
    public double df(int nder, double x) {
        switch (nder) {
            //c
            case 0:
                return (aa0*L/2)+aa0*x+((Math.cos(L*w/2)-Math.cos(w*x))*bb1/w)+((Math.cos(L*w)-Math.cos(2*w*x))*bb2/(2*w))+((Math.cos(3*L*w/2)-Math.cos(3*w*x))*bb3/(3*w))+((Math.cos(2*L*w)-Math.cos(4*w*x))*bb4/(4*w))+((Math.cos(5*L*w/2)-Math.cos(5*w*x))*bb5/(5*w))+((Math.cos(6*L*w/2)-Math.cos(6*w*x))*bb6/(6*w))+((Math.cos(7*L*w/2)-Math.cos(7*w*x))*bb7/(7*w))+((Math.cos(8*L*w/2)-Math.cos(8*w*x))*bb8/(8*w))+ ((Math.sin(L*w/2)+Math.sin(w*x))*aa1/w)+((Math.sin(2*L*w/2)+Math.sin(2*w*x))*aa2/(2*w))+((Math.sin(3*L*w/2)+Math.sin(3*w*x))*aa3/(3*w))+((Math.sin(4*L*w/2)+Math.sin(4*w*x))*aa4/(4*w))+((Math.sin(5*L*w/2)+Math.sin(5*w*x))*aa5/(5*w))+((Math.sin(6*L*w/2)+Math.sin(6*w*x))*aa6/(6*w))+((Math.sin(7*L*w/2)+Math.sin(7*w*x))*aa7/(7*w))+((Math.sin(8*L*w/2)+Math.sin(8*w*x))*aa8/(8*w));
            // p
            case 1:
                return (aa0+aa1*Math.cos(x*w)+bb1*Math.sin(x*w)+aa2*Math.cos(2*x*w)+bb2*Math.sin(2*x*w)+aa3*Math.cos(3*x*w)+bb3*Math.sin(3*x*w)+aa4*Math.cos(4*x*w)+bb4*Math.sin(4*x*w)+aa5*Math.cos(5*x*w)+bb5*Math.sin(5*x*w)+aa6*Math.cos(6*x*w)+bb6*Math.sin(6*x*w)+aa7*Math.cos(7*x*w)+bb7*Math.sin(7*x*w)+aa8*Math.cos(8*x*w)+bb8*Math.sin(8*x*w));
            // dp/dz
            case 2:
                return (bb1*w*Math.cos(x*w)+2*bb2*w*Math.cos(2*x*w)+3*bb3*w*Math.cos(3*x*w)+4*bb4*w*Math.cos(4*x*w)+5*bb5*w*Math.cos(5*x*w)+6*bb6*w*Math.cos(6*x*w)+7*bb7*w*Math.cos(7*x*w)+8*bb8*w*Math.cos(8*x*w)-aa1*w*Math.sin(x*w)-2*aa2*w*Math.sin(2*x*w)-3*aa3*w*Math.sin(3*x*w)-4*aa4*w*Math.sin(4*x*w)-5*aa5*w*Math.sin(5*x*w)-6*aa6*w*Math.sin(6*x*w)-7*aa7*w*Math.sin(7*x*w)-8*aa8*w*Math.sin(8*x*w));
            default:
                throw new RuntimeException("can't do that");
        }
    }

    @Override
    public double f(double x) {
        return df(0, x);
    }

}
