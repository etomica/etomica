package etomica.modules;

import etomica.graphics.SimulationGraphic;

import javax.swing.*;
import java.awt.*;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

public class ModuleLauncher extends JPanel {
    private EtomicaModuleInfo currentModule;
    private final JTextPane descriptionPane;
    private JTextField argsField;

    public ModuleLauncher() {
        super(new BorderLayout(10, 10));
        this.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        descriptionPane = new JTextPane();
        descriptionPane.setEditable(false);
        descriptionPane.setText("Module description");
        this.add(descriptionPane, BorderLayout.CENTER);

        JComboBox<EtomicaModuleInfo> selector = new JComboBox<>(EtomicaModuleInfo.ETOMICA_MODULES);
        this.currentModule = (EtomicaModuleInfo) selector.getSelectedItem();
        selector.addActionListener(e -> {
            this.currentModule = (EtomicaModuleInfo) ((JComboBox) e.getSource()).getSelectedItem();
            this.descriptionPane.setText(this.currentModule.description);
            this.argsField.setText(this.currentModule.args);
        });
        this.add(selector, BorderLayout.NORTH);

        JButton launchButton = new JButton("Launch");
        launchButton.addActionListener(e -> {
            Class<?> cls = this.currentModule.moduleClass;
            try {
                Method method = cls.getMethod("main", String[].class);
                String[] args = this.argsField.getText().trim().isEmpty() ? new String[0] : this.argsField.getText().split(" ");
                SwingUtilities.invokeLater(() -> {
                    try {
                        method.invoke(null, (Object) args);
                    } catch (IllegalAccessException | InvocationTargetException e1) {
                        e1.printStackTrace();
                    }
                });

            } catch (NoSuchMethodException e1) {
                e1.printStackTrace();
            }
        });


        argsField = new JTextField(10);
        argsField.setToolTipText("Command-line arguments for the module main method. This will usually " +
                "be empty unless you know otherwise");
        JPanel launchPanel = new JPanel();
        launchPanel.setLayout(new BorderLayout(5, 5));
        JLabel argsLabel = new JLabel("Module args:");
        argsLabel.setOpaque(true);
        argsLabel.setLabelFor(argsField);
        launchPanel.add(argsLabel, BorderLayout.WEST);
        launchPanel.add(argsField, BorderLayout.CENTER);
        launchPanel.add(launchButton, BorderLayout.EAST);
        this.add(launchPanel, BorderLayout.SOUTH);

    }

    private static void createGUI() {
        SimulationGraphic.initGraphics();
        JFrame frame = new JFrame("Etomica Module Launcher");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        JPanel panel = new ModuleLauncher();
        panel.setOpaque(true);

        frame.setContentPane(panel);
        frame.pack();
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(ModuleLauncher::createGUI);
    }

}
