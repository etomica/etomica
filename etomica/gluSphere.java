package etomica;

import gl4java.GLFunc;

public final class gluSphere {
  /* Make it not a power of two to avoid cache thrashing on the chip */
  private static final int CACHE_SIZE = 240;
  private static final boolean lines = false;
  private static gl4java.GLEnum glt;

  public static /*float[][]*/void gluSmoothSphere(gl4java.GLFunc glf, double radius, int divisions) {
    float sinCache1[] = new float[CACHE_SIZE];
    float cosCache1[] = new float[CACHE_SIZE];
    float sinCache2[] = new float[CACHE_SIZE];
    float cosCache2[] = new float[CACHE_SIZE];
    float sintemp1=0, sintemp2=0;
    float costemp1=0, costemp2=0;
    float angle;
    int i,j,k;

    int stacks = divisions, slices = divisions*2;

    if (slices >= CACHE_SIZE) slices = CACHE_SIZE-1;
    if (stacks >= CACHE_SIZE) stacks = CACHE_SIZE-1;

    if (divisions < 2 || radius < 0.0) {
      //gluQuadricError(qobj, w.GLU_INVALID_VALUE);
      return;
    }

    /* Cache1 is the vertex locations cache */
    /* Cache2 is the various normals at the vertices themselves */
    for (j = 0; j <= stacks; j++) {
      angle = (float)(Math.PI * j / stacks);
      sinCache2[j] = (float)Math.sin(angle);
      cosCache2[j] = (float)Math.cos(angle);
      sinCache1[j] = (float)(radius * Math.sin(angle));
      cosCache1[j] = (float)(radius * Math.cos(angle));
    }
    /* Make sure it comes to a point */
    sinCache1[0] = 0;
    sinCache1[stacks] = 0;

    float normalArray[][] = new float[stacks+1][((stacks*4)+2)*3];
    float vertexArray[][] = new float[stacks+1][((stacks*4)+2)*3];
    //int vertexCount=0;
    int tempInt;
    normalArray[0][1] = sinCache2[0];
    normalArray[0][2] = cosCache2[0];
    vertexArray[0][2] = (float)radius;
    normalArray[stacks][0] = (float)Math.sin(2*(Math.PI*stacks/slices)) * sinCache2[stacks];
    normalArray[stacks][1] = (float)Math.cos(2*(Math.PI*stacks/slices)) * sinCache2[stacks];
    normalArray[stacks][2] = cosCache2[stacks];
    vertexArray[stacks][2] = (float)-radius;
    for (j = 1; j < stacks; j++) {
      sintemp1 = sinCache1[j];
      costemp1 = cosCache1[j];
      sintemp2 = sinCache2[j];
      costemp2 = cosCache2[j];

      if(j > stacks/2) {
        tempInt = ((stacks-j)*4);
      } else
        tempInt=j*4;
      if(lines)
        glf.glBegin(glt.GL_LINE_STRIP);
      for (i = 0; i <= tempInt; i++) {
        angle = 2 * (float)(Math.PI * i / tempInt);
        normalArray[j][(i*3)] = sintemp2*(float)Math.sin(angle);
        normalArray[j][(i*3)+1] = sintemp2*(float)Math.cos(angle);
        normalArray[j][(i*3)+2] = costemp2;
        vertexArray[j][(i*3)] = sintemp1*(float)Math.sin(angle);
        vertexArray[j][(i*3)+1] = sintemp1*(float)Math.cos(angle);
        vertexArray[j][(i*3)+2] = costemp1;
        if(lines) {
          glf.glNormal3f(normalArray[j][(i*3)], normalArray[j][(i*3)+1], normalArray[j][(i*3)+2]);
          glf.glVertex3f(vertexArray[j][(i*3)], vertexArray[j][(i*3)+1], vertexArray[j][(i*3)+2]);
        }
      }
      if(lines)
        glf.glEnd();
    }
    for (j = 1; j <= (stacks/2); j++) {
      if(lines)
        glf.glBegin(glt.GL_LINE_STRIP);
      else
        glf.glBegin(glt.GL_TRIANGLE_STRIP);
      tempInt=j*4;
      //vertexCount+=3;
      for (i = 0, k = 0; i < tempInt; i++) {
        if(!((i==0)||(i==j)||(i==(2*j))||(i==(3*j)))) {
          k++;
        }
        glf.glNormal3f(normalArray[j][(i*3)], normalArray[j][(i*3)+1], normalArray[j][(i*3)+2]);
        glf.glVertex3f(vertexArray[j][(i*3)], vertexArray[j][(i*3)+1], vertexArray[j][(i*3)+2]);
        glf.glNormal3f(normalArray[j-1][(k*3)], normalArray[j-1][(k*3)+1], normalArray[j-1][(k*3)+2]);
        glf.glVertex3f(vertexArray[j-1][(k*3)], vertexArray[j-1][(k*3)+1], vertexArray[j-1][(k*3)+2]);
        //vertexCount+=2;
      }
      glf.glNormal3f(normalArray[j][0], normalArray[j][1], normalArray[j][2]);
      glf.glVertex3f(vertexArray[j][0], vertexArray[j][1], vertexArray[j][2]);
      //vertexCount++;
      glf.glEnd();
    }
    for (j = stacks-1; j >= (stacks/2); j--) {
      if(lines)
        glf.glBegin(glt.GL_LINE_STRIP);
      else
        glf.glBegin(glt.GL_TRIANGLE_STRIP);
      tempInt = ((stacks-j)*4);
      //vertexCount+=3;
      if (!((j==(stacks/2))&&((divisions%2)==1))) {
        for (i = tempInt, k = ((stacks-j-1)*4); (i >= 0) && (k >= 0); i--) {
          if(!((i==0)||(i==(tempInt/4))||(i==(2*(tempInt/4)))||(i==(3*(tempInt/4)))||(i==tempInt))) {
            k--;
          }
          if(!((vertexArray[j][(i*3)]==0)&&(vertexArray[j][(i*3)+1]==0)&&(vertexArray[j][(i*3)+2]==0))) {
            glf.glNormal3f(normalArray[j][(i*3)], normalArray[j][(i*3)+1], normalArray[j][(i*3)+2]);
            glf.glVertex3f(vertexArray[j][(i*3)], vertexArray[j][(i*3)+1], vertexArray[j][(i*3)+2]);
            //vertexCount++;
          }
          if(!((vertexArray[j+1][(k*3)]==0)&&(vertexArray[j+1][(k*3)+1]==0)&&(vertexArray[j+1][(k*3)+2]==0))) {
            glf.glNormal3f(normalArray[j+1][(k*3)], normalArray[j+1][(k*3)+1], normalArray[j+1][(k*3)+2]);
            glf.glVertex3f(vertexArray[j+1][(k*3)], vertexArray[j+1][(k*3)+1], vertexArray[j+1][(k*3)+2]);
            //vertexCount++;
          }
        }
      } else {
        for (i = tempInt, k = tempInt; (i >= 0) && (k >= 0); i--) {
          if(!((vertexArray[j][(i*3)]==0)&&(vertexArray[j][(i*3)+1]==0)&&(vertexArray[j][(i*3)+2]==0))) {
            glf.glNormal3f(normalArray[j][(i*3)], normalArray[j][(i*3)+1], normalArray[j][(i*3)+2]);
            glf.glVertex3f(vertexArray[j][(i*3)], vertexArray[j][(i*3)+1], vertexArray[j][(i*3)+2]);
            //vertexCount++;
          }
          if(!((vertexArray[j+1][(k*3)]==0)&&(vertexArray[j+1][(k*3)+1]==0)&&(vertexArray[j+1][(k*3)+2]==0))) {
            glf.glNormal3f(normalArray[j+1][(k*3)], normalArray[j+1][(k*3)+1], normalArray[j+1][(k*3)+2]);
            glf.glVertex3f(vertexArray[j+1][(k*3)], vertexArray[j+1][(k*3)+1], vertexArray[j+1][(k*3)+2]);
            //vertexCount++;
          }
          k--;
        }
      }
      glf.glEnd();
    }
    /*float finalArray[][] = new float[2][vertexCount*3];
    vertexCount=0;
    for (j = 1; j <= (stacks/2); j++) {
      tempInt=j*4;
      finalArray[0][(vertexCount*3)]=normalArray[j][0];
      finalArray[0][(vertexCount*3)+1]=normalArray[j][1];
      finalArray[0][(vertexCount*3)+2]=normalArray[j][2];
      finalArray[1][(vertexCount*3)]=vertexArray[j][0];
      finalArray[1][(vertexCount*3)+1]=vertexArray[j][1];
      finalArray[1][(vertexCount*3)+2]=vertexArray[j][2];
      vertexCount++;
      finalArray[0][(vertexCount*3)]=normalArray[j][0];
      finalArray[0][(vertexCount*3)+1]=normalArray[j][1];
      finalArray[0][(vertexCount*3)+2]=normalArray[j][2];
      finalArray[1][(vertexCount*3)]=vertexArray[j][0];
      finalArray[1][(vertexCount*3)+1]=vertexArray[j][1];
      finalArray[1][(vertexCount*3)+2]=vertexArray[j][2];
      vertexCount++;
      finalArray[0][(vertexCount*3)]=normalArray[j][0];
      finalArray[0][(vertexCount*3)+1]=normalArray[j][1];
      finalArray[0][(vertexCount*3)+2]=normalArray[j][2];
      finalArray[1][(vertexCount*3)]=vertexArray[j][0];
      finalArray[1][(vertexCount*3)+1]=vertexArray[j][1];
      finalArray[1][(vertexCount*3)+2]=vertexArray[j][2];
      vertexCount++;
      for (i = 0, k = 0; i < tempInt; i++) {
        if(!((i==0)||(i==j)||(i==(2*j))||(i==(3*j)))) {
          k++;
        }
        finalArray[0][(vertexCount*3)]=normalArray[j][(i*3)];
        finalArray[0][(vertexCount*3)+1]=normalArray[j][(i*3)+1];
        finalArray[0][(vertexCount*3)+2]=normalArray[j][(i*3)+2];
        finalArray[1][(vertexCount*3)]=vertexArray[j][(i*3)];
        finalArray[1][(vertexCount*3)+1]=vertexArray[j][(i*3)+1];
        finalArray[1][(vertexCount*3)+2]=vertexArray[j][(i*3)+2];
        vertexCount++;
        finalArray[0][(vertexCount*3)]=normalArray[j-1][(k*3)];
        finalArray[0][(vertexCount*3)+1]=normalArray[j-1][(k*3)+1];
        finalArray[0][(vertexCount*3)+2]=normalArray[j-1][(k*3)+2];
        finalArray[1][(vertexCount*3)]=vertexArray[j-1][(k*3)];
        finalArray[1][(vertexCount*3)+1]=vertexArray[j-1][(k*3)+1];
        finalArray[1][(vertexCount*3)+2]=vertexArray[j-1][(k*3)+2];
        vertexCount++;
      }
      finalArray[0][(vertexCount*3)]=normalArray[j][0];
      finalArray[0][(vertexCount*3)+1]=normalArray[j][1];
      finalArray[0][(vertexCount*3)+2]=normalArray[j][2];
      finalArray[1][(vertexCount*3)]=vertexArray[j][0];
      finalArray[1][(vertexCount*3)+1]=vertexArray[j][1];
      finalArray[1][(vertexCount*3)+2]=vertexArray[j][2];
      vertexCount++;
    }
    for (j = stacks-1; j >= (stacks/2); j--) {
      tempInt = ((stacks-j)*4);
      if(!((vertexArray[j][0]==0)&&(vertexArray[j][1]==0)&&(vertexArray[j][2]==0))) {
        finalArray[0][(vertexCount*3)]=normalArray[j][0];
        finalArray[0][(vertexCount*3)+1]=normalArray[j][1];
        finalArray[0][(vertexCount*3)+2]=normalArray[j][2];
        finalArray[1][(vertexCount*3)]=vertexArray[j][0];
        finalArray[1][(vertexCount*3)+1]=vertexArray[j][1];
        finalArray[1][(vertexCount*3)+2]=vertexArray[j][2];
        vertexCount++;
        finalArray[0][(vertexCount*3)]=normalArray[j][0];
        finalArray[0][(vertexCount*3)+1]=normalArray[j][1];
        finalArray[0][(vertexCount*3)+2]=normalArray[j][2];
        finalArray[1][(vertexCount*3)]=vertexArray[j][0];
        finalArray[1][(vertexCount*3)+1]=vertexArray[j][1];
        finalArray[1][(vertexCount*3)+2]=vertexArray[j][2];
        vertexCount++;
        finalArray[0][(vertexCount*3)]=normalArray[j][0];
        finalArray[0][(vertexCount*3)+1]=normalArray[j][1];
        finalArray[0][(vertexCount*3)+2]=normalArray[j][2];
        finalArray[1][(vertexCount*3)]=vertexArray[j][0];
        finalArray[1][(vertexCount*3)+1]=vertexArray[j][1];
        finalArray[1][(vertexCount*3)+2]=vertexArray[j][2];
        vertexCount++;
      }
      if (!((j==(stacks/2))&&((divisions%2)==1))) {
        for (i = tempInt, k = ((stacks-j-1)*4); (i >= 0) && (k >= 0); i--) {
          if(!((i==0)||(i==(tempInt/4))||(i==(2*(tempInt/4)))||(i==(3*(tempInt/4)))||(i==tempInt))) {
            k--;
          }
          if(!((vertexArray[j][(i*3)]==0)&&(vertexArray[j][(i*3)+1]==0)&&(vertexArray[j][(i*3)+2]==0))) {
            finalArray[0][(vertexCount*3)]=normalArray[j][(i*3)];
            finalArray[0][(vertexCount*3)+1]=normalArray[j][(i*3)+1];
            finalArray[0][(vertexCount*3)+2]=normalArray[j][(i*3)+2];
            finalArray[1][(vertexCount*3)]=vertexArray[j][(i*3)];
            finalArray[1][(vertexCount*3)+1]=vertexArray[j][(i*3)+1];
            finalArray[1][(vertexCount*3)+2]=vertexArray[j][(i*3)+2];
            vertexCount++;
          }
          if(!((vertexArray[j+1][(k*3)]==0)&&(vertexArray[j+1][(k*3)+1]==0)&&(vertexArray[j+1][(k*3)+2]==0))) {
            finalArray[0][(vertexCount*3)]=normalArray[j-1][(k*3)];
            finalArray[0][(vertexCount*3)+1]=normalArray[j-1][(k*3)+1];
            finalArray[0][(vertexCount*3)+2]=normalArray[j-1][(k*3)+2];
            finalArray[1][(vertexCount*3)]=vertexArray[j-1][(k*3)];
            finalArray[1][(vertexCount*3)+1]=vertexArray[j-1][(k*3)+1];
            finalArray[1][(vertexCount*3)+2]=vertexArray[j-1][(k*3)+2];
            vertexCount++;
          }
        }
      } else {
        for (i = tempInt, k = tempInt; (i >= 0) && (k >= 0); i--, k--) {
          if(!((vertexArray[j][(i*3)]==0)&&(vertexArray[j][(i*3)+1]==0)&&(vertexArray[j][(i*3)+2]==0))) {
            finalArray[0][(vertexCount*3)]=normalArray[j][(i*3)];
            finalArray[0][(vertexCount*3)+1]=normalArray[j][(i*3)+1];
            finalArray[0][(vertexCount*3)+2]=normalArray[j][(i*3)+2];
            finalArray[1][(vertexCount*3)]=vertexArray[j][(i*3)];
            finalArray[1][(vertexCount*3)+1]=vertexArray[j][(i*3)+1];
            finalArray[1][(vertexCount*3)+2]=vertexArray[j][(i*3)+2];
            vertexCount++;
          }
          if(!((vertexArray[j+1][(k*3)]==0)&&(vertexArray[j+1][(k*3)+1]==0)&&(vertexArray[j+1][(k*3)+2]==0))) {
            finalArray[0][(vertexCount*3)]=normalArray[j-1][(k*3)];
            finalArray[0][(vertexCount*3)+1]=normalArray[j-1][(k*3)+1];
            finalArray[0][(vertexCount*3)+2]=normalArray[j-1][(k*3)+2];
            finalArray[1][(vertexCount*3)]=vertexArray[j-1][(k*3)];
            finalArray[1][(vertexCount*3)+1]=vertexArray[j-1][(k*3)+1];
            finalArray[1][(vertexCount*3)+2]=vertexArray[j-1][(k*3)+2];
            vertexCount++;
          }
          //k--;
        }
      }
    }
    return(finalArray);*/
  }
}