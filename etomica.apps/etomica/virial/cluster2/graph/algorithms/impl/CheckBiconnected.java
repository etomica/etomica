package etomica.virial.cluster2.graph.algorithms.impl;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.graph.algorithms.GraphProperty;

public class CheckBiconnected implements GraphProperty {

//  static boolean
//  isbiconnected(graph *g, int n)
//  /* test if g is biconnected */
//  {
//          int sp,v,w;
//          setword sw;
//          setword visited;
//          int numvis,num[MAXN],lp[MAXN],stack[MAXN];
//   
//          if (n <= 2) return FALSE;
//   
//          visited = bit[0];
//          stack[0] = 0;
//          num[0] = 0;
//          lp[0] = 0;
//          numvis = 1;
//          sp = 0;
//          v = 0;
//   
//          for (;;)
//          {
//              if ((sw = g[v] & ~visited))           /* not "==" */
//              {
//                  w = v;
//                  v = FIRSTBIT(sw);       /* visit next child */
//                  stack[++sp] = v;
//                  visited |= bit[v];
//                  lp[v] = num[v] = numvis++;
//                  sw = g[v] & visited & ~bit[w];
//                  while (sw)
//                  {
//                      w = FIRSTBIT(sw);
//                      sw &= ~bit[w];
//                      if (num[w] < lp[v])  lp[v] = num[w];
//                  }
//              }
//              else
//              {
//                  w = v;                  /* back up to parent */
//                  if (sp <= 1)          return numvis == n;
//                  v = stack[--sp];
//                  if (lp[w] >= num[v])  return FALSE;
//                  if (lp[w] < lp[v])    lp[v] = lp[w];
//              }
//          }
//  }
//

  @Override
  public boolean check(Nodes nodes, Edges edges) {

    if ((nodes.count() <= 2) || (edges.count() < (nodes.count() - 2))) {
      return false;
    }
    // one bit for each node we have to see
    int goal = BitmapUtils.bitMask(nodes.count());
    // one bit for each node seen
    int seen = 0x00000001;
    // one bit for each node explored
    int explored = 0x00000000;
    // nodes we still have to explore
    int toExplore = 0x00000001;
    // done when: (1) all nodes seen OR (2) no new nodes to explore
    while ((seen != goal) && (toExplore != 0)) {
      // pop the next node to explore
      int explore = BitmapUtils.leftmostBit(toExplore);
      for (int node = 0; node < nodes.count(); node++) {
        if ((node != explore) && edges.hasEdge(explore, node)) {
          // mark the node as seen
          seen = seen | BitmapUtils.bitOnMask(node);
          // if the node has not been explored yet, add it to the toExplore list
          if ((explored & BitmapUtils.bitOnMask(node)) == 0) {
            toExplore = toExplore | BitmapUtils.bitOnMask(node);
          }
        }
      }
      // remove the explored node from the bitmap
      toExplore = toExplore & BitmapUtils.bitOffMask(explore);
      // add the explored node to the bitmap
      explored = explored | BitmapUtils.bitOnMask(explore);
    }
    return (seen == goal);
  }
}