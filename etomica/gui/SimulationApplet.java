package etomica.gui;

import etomica.*;
import javax.swing.JApplet;
import java.awt.GridLayout;
import java.awt.Frame;
import java.beans.*;
import javax.swing.JPanel;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class SimulationApplet extends JApplet {
    transient String space;
    javax.swing.JApplet instance;
    
    public void init(){
        getContentPane().setLayout(new GridLayout(1,1)); 
        System.out.println("input stream"); 
        try { 
            //java.io.InputStream ins = ClassLoader.getSystemResourceAsStream("SimulationApplet.ser"); 
            java.io.FileInputStream ins = new java.io.FileInputStream("SimulationApplet.ser"); 
            System.out.println("ins: " + ins); 
            ObjectInputStream ois = new ObjectInputStream(ins);
            System.out.println(ois);
            instance = ((javax.swing.JApplet) ois.readObject());
            System.out.println(instance);
            getContentPane().add(instance);
            System.out.println("yeah"); 
        } 
        catch (Exception ex) { ex.printStackTrace(); } 
    } 

//    public void paint(java.awt.Graphics g){
//        g.setColor(java.awt.Color.blue);
//        g.drawRect(20, 20, 100, 100);
//    }

    public static void main(String args[]) { 
        Frame f = new Frame("SimulationApplet"); 
        SimulationApplet applet = new etomica.gui.SimulationApplet(); 
        applet.init(); 
        f.setLayout(new GridLayout(1,1)); 
        f.add(applet); 
        f.pack(); 
        f.show(); 
        f.setSize(450, 350); 
        // Handle window close requests 
        f.addWindowListener(new java.awt.event.WindowAdapter() { 
        public void windowClosing(java.awt.event.WindowEvent e) { System.exit(0); } 
        }); 
    } 
}