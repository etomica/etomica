package etomica.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.Dimension;
import javax.swing.JButton;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class PeriodicTable extends JInternalFrame {
    java.awt.GridBagConstraints gbc = new java.awt.GridBagConstraints();
    java.awt.GridBagLayout gbl = new java.awt.GridBagLayout();
    JPanel mainPanel = new JPanel();
    JButton[] buttons = new JButton[162];
    JButton button;
    JLabel[] labels = new JLabel[162];
    JLabel label;
    String[] buttonLabels = new String[162];
    
    PeriodicTable(){
        setBounds(500,200, 600,600);
        setMaximizable(true);
        setIconifiable(true);
        setClosable(true);
        setResizable(true);
        setVisible(true);
        mainPanel.setLayout(gbl);
        mainPanel.setSize(600,600);
        gbc.gridx = -1;
        gbc.gridy = 0;
        gbc.gridheight = 1;
        gbc.gridwidth = 1;
        
        buttonLabels[0] = "1  H";
        buttonLabels[17] = "2  He";
        buttonLabels[18] = "3  Li";
        buttonLabels[19] = "4  Be";
        buttonLabels[30] = "5  B";
        buttonLabels[33] = "8  O";
        
        for (int i = 0; i < 162; i++){
            gbc.gridx++;
            if (gbc.gridx == 18) {
                gbc.gridx = 0;
                gbc.gridy++;
            }
            if (buttonLabels[i] != null) {
                button = new JButton(buttonLabels[i]);
                button.setMinimumSize(new Dimension(50,40));
                button.setMaximumSize(new Dimension(50,40));
                button.setPreferredSize(new Dimension(50,40));
                buttons[i] = button;
                button.addActionListener(new ActionListener(){
                    public void actionPerformed(ActionEvent e){
                        System.out.println(((JButton)e.getSource()).getText());
                    }});
                gbl.addLayoutComponent(button,gbc);
                mainPanel.add(button);
            }
            else {
                label = new JLabel();
                label.setMinimumSize(new Dimension(50,40));
                label.setMaximumSize(new Dimension(50,40));
                label.setPreferredSize(new Dimension(50,40));
                labels[i] = label;
                gbl.addLayoutComponent(label,gbc);
                mainPanel.add(label);
            }                
        }
        getContentPane().add(mainPanel);
    }
}