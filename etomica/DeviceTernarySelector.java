package etomica;

/**
 * Selector for fractions in a ternary system.  Presents a triangle that can
 * be clicked upon to select (for example) the three mole fractions of a ternary
 * mixture.
 *
 * @author David Kofke
 */
public class DeviceTernarySelector extends Device implements EtomicaElement {
    
    public String getVersion() {return "DeviceTernarySelector:01.11.17/"+Device.VERSION;}

    private Triangle triangle;
    private int sideLength;
    private javax.swing.JPanel panel = new TrianglePanel();
    
    public DeviceTernarySelector() {
        this(Simulation.instance);
    }
    
    public DeviceTernarySelector(Simulation sim) {
        super(sim);
        sideLength = 200;
        triangle = new Triangle(sideLength);
        panel = new TrianglePanel();
        panel.addMouseListener(new ClickListener());
        panel.setPreferredSize(new java.awt.Dimension(sideLength+10, sideLength+10));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Select fractions in a ternary system");
        return info;
    }
    
    public java.awt.Component graphic(Object obj) {
        return panel;
    }
    
    /**
     * Simple panel used to display the triangle and receive mouse events.
     */
    private class TrianglePanel extends javax.swing.JPanel {
        public void paintComponent(java.awt.Graphics g) {
            super.paintComponent(g);
            g.drawPolygon(triangle);
        }
    }
    
    //listener for mouse events.  Reads mouse position and converts
    //to fractions
    private class ClickListener extends java.awt.event.MouseAdapter {
        public void mousePressed(java.awt.event.MouseEvent evt) {
            int x = evt.getX();
            int y = evt.getY();
//            System.out.println(x + "  " + y + "  " + triangle.contains(x, y));
            double[] frac = triangle.fractions(x,y);
            if(triangle.contains(x,y)) System.out.println(frac[0]+" "+frac[1]+" "+frac[2]);
        }
    }
    
    /**
     * Class used to draw and perform calculations for the triangle.
     */
    private static class Triangle extends java.awt.Polygon {
        
        public int height;
        private double[] x = new double[3];
        public Triangle(int L) {
            super(new int[] {   0,    L, L/2}, //x points       2
                  new int[] {h(L), h(L),   0}, //y points      0 1
                  3);
            height = h(L);
        }
        
        /**
         * The height of the triangle, equal to sqrt(3)/2 times the length of a side.
         */
        private static int h(int L) {return (int)(0.5*Math.sqrt(3.)*L);}
        
        /**
         * Returns the ternary fractions corresponding to the given pixel-point 
         * in the triangle.
         */
        public double[] fractions(int xP, int yP) {
            x[0] = 1.0 - (0.5*Math.sqrt(3.)*xP + 0.5*(height-yP))/height;
            x[1] = 1.0 - (double)yP/(double)height;
            x[2] = 1.0 - x[0] - x[1];
            return x;
        }
        
    }
    
    public static void main(String[] args) {
        Simulation sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        DeviceTernarySelector selector = new DeviceTernarySelector(sim);
        sim.elementCoordinator.go();
        sim.makeAndDisplayFrame();
    }
}