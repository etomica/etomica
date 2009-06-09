package etomica.virial.cluster2.ui;

import java.awt.Dimension;
import java.net.URL;

import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JRadioButton;
import javax.swing.LookAndFeel;
import javax.swing.UIManager;
import javax.swing.plaf.metal.DefaultMetalTheme;
import javax.swing.plaf.metal.MetalLookAndFeel;

import com.jgoodies.looks.Options;
import com.jgoodies.looks.plastic.PlasticLookAndFeel;
import com.jgoodies.looks.plastic.PlasticTheme;
import com.jgoodies.looks.plastic.theme.SkyKrupp;
import com.jgoodies.looks.plastic.theme.DesertBluer;
import com.jgoodies.looks.plastic.theme.ExperienceBlue;
import com.jgoodies.looks.plastic.theme.ExperienceRoyale;
import com.jgoodies.looks.plastic.theme.LightGray;
import com.jgoodies.looks.plastic.theme.Silver;
import com.jgoodies.looks.plastic.theme.SkyBlue;

public class ApplicationUI {

  /**
   * The constants below provide a number of choices for each of the
   * configurable options of the application.
   */
  // Icon Size options
  public static final Dimension        IS_DIM1         = new Dimension(18, 18);
  // Look and Feel options
  public static final String           LF_WINDOWS      = "Windows";
  public static final String           LF_PLASTIC      = "Plastic";
  public static final String           LF_PLASTIC3D    = "Plastic3D";
  public static final String           LF_PLASTICXP    = "PlasticXP";
  // Theme options
  public static final String           THM_DESERTBR    = "DesertBluer";
  public static final String           THM_EXPERIENCEB = "ExperienceBlue";
  public static final String           THM_EXPERIENCER = "ExperienceRoyale";
  public static final String           THM_LIGHTGRAY   = "LightGray";
  public static final String           THM_SILVER      = "Silver";
  public static final String           THM_SKYBLUE     = "SkyBlue";
  public static final String           THM_SKYKRUPP    = "SkyKrupp";
  // Main Window Dimensions
  public static final Dimension        MWP_DIM1        = new Dimension(640, 480);
  public static final Dimension        MWP_DIM2        = new Dimension(720, 540);
  public static final Dimension        MWP_DIM3        = new Dimension(800, 600);
  /**
   * The constants below define the actual values for each of the configurable
   * options above.
   */
  // Default Icon Size
  public static final Dimension        DF_ICON_SIZE    = IS_DIM1;
  // Look and Feel
  public static final String           LF_CHOICE       = LF_PLASTIC3D;
  // Look and Feel
  public static final String           THM_CHOICE      = THM_EXPERIENCER;
  // UI standard settings
  public static final JGoodiesSettings uiSettings      = defaultUISettings();

  /**
   * Configures the user interface; requests Swing settings and JGoodies Looks
   * options from the launcher.
   */
  public static void configure() {

    // UIManager.put("ToolTip.hideAccelerator", Boolean.FALSE);
    Options.setDefaultIconSize(DF_ICON_SIZE);
    Options.setUseNarrowButtons(uiSettings.isUseNarrowButtons());
// Global options
    Options.setTabIconsEnabled(uiSettings.isTabIconsEnabled());
    UIManager.put(Options.POPUP_DROP_SHADOW_ENABLED_KEY, uiSettings
        .isPopupDropShadowEnabled());
// Swing Settings
    LookAndFeel selectedLaf = uiSettings.getSelectedLookAndFeel();
    if (selectedLaf instanceof PlasticLookAndFeel) {
      PlasticLookAndFeel.setPlasticTheme(uiSettings.getSelectedTheme());
      PlasticLookAndFeel.setTabStyle(uiSettings.getPlasticTabStyle());
      PlasticLookAndFeel.setHighContrastFocusColorsEnabled(uiSettings
          .isPlasticHighContrastFocusEnabled());
    }
    else if (selectedLaf.getClass() == MetalLookAndFeel.class) {
      MetalLookAndFeel.setCurrentTheme(new DefaultMetalTheme());
    }
// Work around caching in MetalRadioButtonUI
    JRadioButton radio = new JRadioButton();
    radio.getUI().uninstallUI(radio);
    JCheckBox checkBox = new JCheckBox();
    checkBox.getUI().uninstallUI(checkBox);
// try setting the look and feel
    try {
      UIManager.setLookAndFeel(selectedLaf);
    }
    catch (Exception e) {
      System.out.println("Can't change the look and feel: " + e);
    }
  }

  /**
   * Configures the default UI settings, including look and feel and theme.
   */
  protected static JGoodiesSettings defaultUISettings() {

    JGoodiesSettings settings = JGoodiesSettings.createDefault();
    settings.setSelectedLookAndFeel(getLookAndFeel());
    settings.setTabIconsEnabled(true);
    settings.setSelectedTheme(getTheme());
    // Configure more settings here.
    return settings;
  }

  /**
   * Looks up and returns an icon for the specified filename suffix.
   */
  public static ImageIcon readImageIcon(String filename) {

    URL url = ApplicationUI.class.getResource("resources/images/" + filename);
    return new ImageIcon(url);
  }

  /**
   * Translates the THM_CHOICE field into a PlasticThem class instance.
   */
  public static PlasticTheme getTheme() {

    if (THM_DESERTBR.equalsIgnoreCase(THM_CHOICE)) {
      return new DesertBluer();
    }
    else if (THM_EXPERIENCEB.equalsIgnoreCase(THM_CHOICE)) {
      return new ExperienceBlue();
    }
    else if (THM_EXPERIENCER.equalsIgnoreCase(THM_CHOICE)) {
      return new ExperienceRoyale();
    }
    else if (THM_LIGHTGRAY.equalsIgnoreCase(THM_CHOICE)) {
      return new LightGray();
    }
    else if (THM_SILVER.equalsIgnoreCase(THM_CHOICE)) {
      return new Silver();
    }
    else if (THM_SKYBLUE.equalsIgnoreCase(THM_CHOICE)) {
      return new SkyBlue();
    }
    else {
      return new SkyKrupp();
    }
  }

  /**
   * Translates the LF_CHOICE field into a LookAndFeel class name.
   */
  private static String getLookAndFeel() {

    if (LF_WINDOWS.equalsIgnoreCase(LF_CHOICE)) {
      return Options.JGOODIES_WINDOWS_NAME;
    }
    else if (LF_PLASTIC.equalsIgnoreCase(LF_CHOICE)) {
      return Options.PLASTIC_NAME;
    }
    else if (LF_PLASTIC3D.equalsIgnoreCase(LF_CHOICE)) {
      return Options.PLASTIC3D_NAME;
    }
    else if (LF_PLASTICXP.equalsIgnoreCase(LF_CHOICE)) {
      return Options.PLASTICXP_NAME;
    }
    else {
      return LF_CHOICE;
    }
  }
}