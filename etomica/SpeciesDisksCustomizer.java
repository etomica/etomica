package simulate;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;

public class SpeciesDisksCustomizer extends Panel implements ItemListener, Customizer, ActionListener {

    public SpeciesDisksCustomizer(){
        setBackground(Color.yellow);
        setBounds(20,20,200,200);
    }
    
    public void setObject(Object o){
        species = (Species)o;
        nMolecules.setText(Double.toString(species.getNMolecules()));
        nMolecules.addActionListener(this);
        nMolecules.setBounds(30, 10, 40, 20);
        add(nMolecules);
        color.addItem("black");
        color.addItem("blue");
        color.addItem("cyan");
        color.addItem("red");
        color.addItem("gray");
        color.addItem("green");
        color.addItem("lightGray");
        color.addItem("magenta");
        color.setBounds(50, 40, 100, 100);
        color.addItemListener(this);
        add(color);
    }
    
    public void itemStateChanged(ItemEvent evt){
        int i = color.getSelectedIndex();
        System.out.println("inside itemStateChanged; index is " + i);
        //a.valuesChanged();
        repaint();
        pcs.firePropertyChange(null, null, null);
    }
    
    public void actionPerformed(ActionEvent evt){
        int d = 10;
        String s = "";
        if(evt.getSource() == nMolecules){
            s=nMolecules.getText();
            if(s.length()>0){
                d = Integer.valueOf(s.trim()).intValue();
                species.setNMolecules(d);
            }
        }
        //a.valuesChanged();

        pcs.firePropertyChange(null, null, null);
    }
    
    public void paint(Graphics g){
//        g.setColor(a.color);
        g.fillOval(100, 170, 15, 15);
    }
    
    public void addPropertyChangeListener(PropertyChangeListener l){
        pcs.addPropertyChangeListener(l);
    }
    
    public void removePropertyChangeListener(PropertyChangeListener l){
        pcs.removePropertyChangeListener(l);
    }
    
    private TextField nMolecules = new TextField(5);
    private List color = new List();
    Species species;
    PropertyChangeSupport pcs = new PropertyChangeSupport(this);
}