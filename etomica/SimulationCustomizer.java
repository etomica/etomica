package simulate;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;

public class SimulationCustomizer extends Panel implements Customizer,ActionListener {
    public SimulationCustomizer(){
        setBackground(Color.gray);
        setBounds(20,20,200,200);
    }
    
    public void setObject(Object o){
        species = (Species)o;
        nMolecules.setText(Double.toString(species.getNMolecules()));
        nMolecules.addActionListener(this);
        nMolecules.setBounds(30, 10, 40, 20);
        add(species);
    }
    
/*    public void itemStateChanged(ItemEvent evt){
        int i = color.getSelectedIndex();
        System.out.println("inside itemStateChanged; index is " + i);
        switch(i){
            case 0: a.setColor(Color.black);
                    break;
            case 1: a.setColor(Color.blue);
                    break;
            case 2: a.setColor(Color.cyan);
                    break;
            case 3: a.setColor(Color.red);
                    break;
            case 4: a.setColor(Color.gray);
                    break;
            case 5: a.setColor(Color.green);
                    break;
            case 6: a.setColor(Color.lightGray);
                    break;
            case 7: a.setColor(Color.magenta);
                    break;
            default: a.setColor(Color.black);
        }
        //a.valuesChanged();
        repaint();
        //pcs.firePropertyChange(null, null, null);
    }*/
    
    public void actionPerformed(ActionEvent evt){
        int d = 0;
        String s = "";
        if(evt.getSource() == nMolecules){
            s=nMolecules.getText();
            if(s.length()>0){
                d = Integer.valueOf(s.trim()).intValue();
                species.setNMolecules(d);
            }
        }
        //a.valuesChanged();

        //pcs.firePropertyChange(null, null, null);
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
    Species species;
    PropertyChangeSupport pcs = new PropertyChangeSupport(this);
}