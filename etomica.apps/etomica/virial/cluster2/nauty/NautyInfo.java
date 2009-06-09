package etomica.virial.cluster2.nauty;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.virial.cluster2.graph.Tagged;

public class NautyInfo implements ProcessInfo, Tagged {

  private static final String SHELL = "/bin/bash";
  private static final String SHELL_FLAG = "-c";

  private static final String CMD_GENG = "./geng";
  private static final String CMD_GENBG = "./genbg";

  private static final String GEN_FLAG_CONNECTED = "-c";
  private static final String GEN_FLAG_BICONNECTED = "-C";

  private String basePath = "";
  private int numBlackNodes = 0;
  private int numWhiteNodes = 0;

  private boolean isMono = true;
  private boolean isConnected = false;
  private boolean isBiconnected = false;
  private Set<String> tags = new HashSet<String>();

  public NautyInfo(final String runFrom, int nodeCount) {

    isMono = true;
    basePath = runFrom;
    numBlackNodes = nodeCount;
    computeTags();
  }

  public NautyInfo(final String runFrom, int nodeCount1, int nodeCount2) {

    isMono = false;
    basePath = runFrom;
    numBlackNodes = nodeCount1;
    numWhiteNodes = nodeCount2;
    computeTags();
  }

  protected void computeTags() {
    
    if (isMono()) {
      tags.add("Monochromatic");
    } else {
      tags.add("Bichromatic");
    }
    tags.add("Isomorph-Free");
    if (isBiconnected()) {
      tags.add("Biconnected");
    } else if (isConnected()) {
      tags.add("Connected");
    }
  }
  
  protected String getGenCommand() {

    String result = "";
    if (isMono()) {
      result += NautyInfo.CMD_GENG;
    } else {
      result += NautyInfo.CMD_GENBG;
    }
    if (isBiconnected()) {
      result += " " + NautyInfo.GEN_FLAG_BICONNECTED;
    } else if (isConnected()) {
      result += " " + NautyInfo.GEN_FLAG_CONNECTED;
    }
    result += " " + getNumBlackNodes();
    if (!isMono()) {
      result += " " + getNumWhiteNodes();
    }
    return result;
  }

  @Override
  public List<String> getCommand() {

    List<String> result = new ArrayList<String>();
    result.add(NautyInfo.SHELL);
    result.add(NautyInfo.SHELL_FLAG);
    result.add(getGenCommand());
    return result;
  }

  public boolean isMono() {

    return isMono;
  }

  public void setConnected(final boolean value) {

    isConnected = value;
  }

  public boolean isConnected() {

    return isConnected;
  }

  public void setBiconnected(final boolean value) {

    isBiconnected = value;
  }

  public boolean isBiconnected() {

    return isBiconnected;
  }

  public int getNumBlackNodes() {

    return numBlackNodes;
  }

  public int getNumWhiteNodes() {

    return numWhiteNodes;
  }

  @Override
  public String getPath() {

    return basePath;
  }

  @Override
  public Set<String> getTags() {

    return Collections.unmodifiableSet(tags);
  }

  @Override
  public String toString() {

    return getCommand() + "\n" + getTags();
  }

  @Override
  public boolean redirectErrorStream() {

    return true;
  }
}