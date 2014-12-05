/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.nauty;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import etomica.virial.cluster2.graph.Tagged;
import etomica.virial.cluster2.util.TagsList;

public class NautyInfo implements ProcessInfo, Tagged {

  // public tags used by NautyInfo
  public static final String TAG_MONOCHROMATIC = "Monochromatic";
  public static final String TAG_BICHROMATIC = "Bichromatic";
  public static final String TAG_ISOMORPH_FREE = "Isomorph-Free";
  public static final String TAG_CONNECTED = "Connected";
  public static final String TAG_BICONNECTED = "Biconnected";
  public static final String TAG_UPPER_TRIANGLE = "Upper Triangle";
  public static final String TAG_ADJACENCY_MATRIX = "Adjacency Matrix";
  // private fields to parametrize calls to nauty
  private static final String SHELL = "/bin/bash";
  private static final String SHELL_FLAG = "-c";
  private static final String CMD_GENG = "./geng";
  private static final String CMD_GENBG = "./genbg";
  private static final String GEN_FLAG_CONNECTED = "-c";
  private static final String GEN_FLAG_BICONNECTED = "-C";
  private static final String GEN_FLAG_UPPER_TRIANGLE = "-T";
  // implementation fields
  private String basePath = "";
  private int numBlackNodes = 0;
  private int numWhiteNodes = 0;
  private boolean isMono = true;
  private boolean isConnected = false;
  private boolean isBiconnected = false;
  private boolean isUpperTriangle = true;
  private List<String> tags = new TagsList();

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

    tags.clear();
    if (isMono()) {
      tags.add(TAG_MONOCHROMATIC);
    }
    else {
      tags.add(TAG_BICHROMATIC);
    }
    tags.add(TAG_ISOMORPH_FREE);
    if (isBiconnected()) {
      tags.add(TAG_BICONNECTED);
    }
    else if (isConnected()) {
      tags.add(TAG_CONNECTED);
    }
    if (isUpperTriangle()) {
      tags.add(TAG_UPPER_TRIANGLE);
    }
    else {
      tags.add(TAG_ADJACENCY_MATRIX);
    }
  }

  protected String getGenCommand() {

    String result = "";
    if (isMono()) {
      result += CMD_GENG;
    }
    else {
      result += CMD_GENBG;
    }
    if (isBiconnected()) {
      result += " " + GEN_FLAG_BICONNECTED;
    }
    else if (isConnected()) {
      result += " " + GEN_FLAG_CONNECTED;
    }
    if (isUpperTriangle()) {
      result += " " + GEN_FLAG_UPPER_TRIANGLE;
    }
    result += " " + getNumBlackNodes();
    if (!isMono()) {
      result += " " + getNumWhiteNodes();
    }
    return result;
  }

  public List<String> getCommand() {

    List<String> result = new ArrayList<String>();
    result.add(SHELL);
    result.add(SHELL_FLAG);
    result.add(getGenCommand());
    return result;
  }

  public boolean isMono() {

    return isMono;
  }

  public void setConnected(boolean value) {

    isConnected = value;
    computeTags();
  }

  public boolean isConnected() {

    return isConnected;
  }

  public boolean isUpperTriangle() {

    return isUpperTriangle;
  }

  public void setBiconnected(boolean value) {

    isBiconnected = value;
    computeTags();
  }

  public boolean isBiconnected() {

    return isBiconnected;
  }

  public void setUpperTriangle(boolean value) {

    isUpperTriangle = value;
    computeTags();
  }

  public int getNumBlackNodes() {

    return numBlackNodes;
  }

  public int getNumWhiteNodes() {

    return numWhiteNodes;
  }

  public String getPath() {

    return basePath;
  }

  public List<String> getTags() {

    return Collections.unmodifiableList(tags);
  }

  @Override
  public String toString() {

    return getCommand().toString();
  }

  public boolean redirectErrorStream() {

    return true;
  }
}