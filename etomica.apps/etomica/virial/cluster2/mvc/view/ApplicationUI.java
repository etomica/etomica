/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.AlphaComposite;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.URL;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.LookAndFeel;
import javax.swing.UIManager;
import javax.swing.plaf.metal.MetalLookAndFeel;
import javax.swing.plaf.metal.DefaultMetalTheme;

import com.jgoodies.looks.Options;
import com.jgoodies.looks.plastic.PlasticLookAndFeel;
import com.jgoodies.looks.plastic.PlasticTheme;
import com.jgoodies.looks.plastic.theme.DesertBluer;
import com.jgoodies.looks.plastic.theme.ExperienceBlue;
import com.jgoodies.looks.plastic.theme.ExperienceRoyale;
import com.jgoodies.looks.plastic.theme.LightGray;
import com.jgoodies.looks.plastic.theme.Silver;
import com.jgoodies.looks.plastic.theme.SkyBlue;
import com.jgoodies.looks.plastic.theme.SkyKrupp;

public class ApplicationUI {

  /**
   * The constants below provide a number of choices for each of the configurable options
   * of the application.
   */
  // Icon Size options
  public static final Dimension IS_DIM1 = new Dimension(18, 18);
  // Look and Feel options
  public static final String LF_WINDOWS = "Windows";
  public static final String LF_PLASTIC = "Plastic";
  public static final String LF_PLASTIC3D = "Plastic3D";
  public static final String LF_PLASTICXP = "PlasticXP";
  // Theme options
  public static final String THM_DESERTBR = "DesertBluer";
  public static final String THM_EXPERIENCEB = "ExperienceBlue";
  public static final String THM_EXPERIENCER = "ExperienceRoyale";
  public static final String THM_LIGHTGRAY = "LightGray";
  public static final String THM_SILVER = "Silver";
  public static final String THM_SKYBLUE = "SkyBlue";
  public static final String THM_SKYKRUPP = "SkyKrupp";
  // Resource relative base folder
  public static final String RESOURCE_BASE_FOLDER = "images/";
  /**
   * The constants below define the actual values for each of the configurable options
   * above.
   */
  // Default Icon Size
  public static final Dimension DF_ICON_SIZE = IS_DIM1;
  // Look and Feel
  public static final String LF_CHOICE = LF_PLASTIC3D;
  // Look and Feel
  public static final String THM_CHOICE = THM_EXPERIENCER;
  // use the theme's control background color for our custom dialogs instead of the
  // theme's dialog background
  public static boolean overrideDialogBackground = false;
  // UI standard settings
  public static final JGoodiesSettings uiSettings = defaultUISettings();

  /**
   * Configures the user interface; requests Swing settings and JGoodies Looks options
   * from the launcher.
   */
  public static void configure() {

    // UIManager.put("ToolTip.hideAccelerator", Boolean.FALSE);
    Options.setDefaultIconSize(DF_ICON_SIZE);
    Options.setUseNarrowButtons(uiSettings.isUseNarrowButtons());
    // Global options
    Options.setTabIconsEnabled(uiSettings.isTabIconsEnabled());
    UIManager.put(Options.POPUP_DROP_SHADOW_ENABLED_KEY, uiSettings.isPopupDropShadowEnabled());
    // Swing Settings
    LookAndFeel selectedLaf = uiSettings.getSelectedLookAndFeel();
    if (selectedLaf instanceof PlasticLookAndFeel) {
      PlasticLookAndFeel.setPlasticTheme(uiSettings.getSelectedTheme());
      PlasticLookAndFeel.setTabStyle(uiSettings.getPlasticTabStyle());
      PlasticLookAndFeel.setHighContrastFocusColorsEnabled(uiSettings.isPlasticHighContrastFocusEnabled());
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
    settings.setSelectedTheme(getTheme()); // PlasticLookAndFeel.createMyDefaultTheme()
    // Configure more settings here.
    UIManager.put("Application.useSystemFontSettings", Boolean.TRUE);
    UIManager.put(Options.USE_SYSTEM_FONTS_APP_KEY, Boolean.TRUE);
    return settings;
  }

  /**
   * Looks up and returns an icon for the specified filename suffix.
   */
  public static ImageIcon readImageIcon(String filename) {

    URL url = ApplicationUI.class.getResource(RESOURCE_BASE_FOLDER + filename);
    return new ImageIcon(url);
  }

  public static class ImagePanel extends JPanel {

    private static final long serialVersionUID = -7302863849251794527L;
    private BufferedImage image;
    private static int arNum = 10;
    private static int arDen = 20;
    private static final int delta = 10;
    private static final float transparency = 0.4f;

    public ImagePanel(final URL resourceURL) {

      try {
        BufferedImage memImage = ImageIO.read(resourceURL);
        image = new BufferedImage(memImage.getWidth(), memImage.getHeight(), BufferedImage.TRANSLUCENT);
        // Get the image's graphics
        Graphics2D g = image.createGraphics();
        // Set the Graphics composite to Alpha
        g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, transparency));
        // Draw the LOADED image onto the actual image
        g.drawImage(memImage, null, 0, 0);
        // let go of all system resources in this Graphics
        g.dispose();
      }
      catch (IOException ex) {
        // handle exception...
      }
    }

    public ImagePanel(final String resourceName) {

      this(ApplicationUI.class.getResource(RESOURCE_BASE_FOLDER + resourceName));
    }

    @Override
    public void paintComponent(Graphics g) {

      super.paintComponent(g);
      int arWidth = image.getWidth();
      int arHeight = image.getHeight();
      Dimension prefSize = new Dimension(arWidth, arHeight);
      if (image.getWidth() > (arNum * getWidth() / arDen)) {
        prefSize.width = arNum * getWidth() / arDen;
        prefSize.height = prefSize.width * arHeight / arWidth;
      }
      else if (image.getHeight() > (arNum * getHeight() / arDen)) {
        prefSize.height = arNum * getHeight() / arDen;
        prefSize.width = prefSize.height * arWidth / arHeight;
      }
      // System.out.println(getWidth() + ":" + getHeight() + ":"
      // + image.getWidth() + ":" + image.getHeight() + ":"
      // + prefSize.getWidth() + ":" + prefSize.getHeight());
      g.drawImage(image, getWidth() - prefSize.width - delta, getHeight() - prefSize.height - delta,
          prefSize.width, prefSize.height, null);
      // System.err.println("components: " + getComponentCount());
      // for (int i = 0; i < getComponentCount(); i++) {
      // System.err.println("component #" + i + ": " +
      // getComponent(i).getClass().getName());
      // }
    }
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