package etomica.graphics.vis;

import com.jogamp.opengl.*;
import com.jogamp.opengl.util.GLBuffers;
import com.jogamp.opengl.util.glsl.ShaderCode;
import com.jogamp.opengl.util.glsl.ShaderProgram;
import org.joml.Matrix4f;
import org.joml.Vector3f;

import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.util.Random;

public class BoxRenderer implements GLEventListener, MouseMotionListener, MouseWheelListener {

    private final IntBuffer intBuf = GLBuffers.newDirectIntBuffer(1);
    private final FloatBuffer matBuffer = GLBuffers.newDirectFloatBuffer(16);

    private int width;
    private int height;
    private float zoom = 20;
    private float pitch = 0.3f;
    private float yaw = 0.2f;
    private int mouseX;
    private int mouseY;

    private FBObject gBufferFBO;
    private FBObject.TextureAttachment gPosition;
    private FBObject.TextureAttachment gNormal;
    private FBObject.TextureAttachment gAlbedoSpec;
    private FBObject ssaoFBO;
    private FBObject.TextureAttachment ssaoColorBuffer;
    private FBObject ssaoBlurFBO;
    private FBObject.TextureAttachment ssaoColorBufferBlur;

    private ShaderProgram shaderGeomPass;
    private ShaderProgram shaderLightingPass;
    private ShaderProgram shaderSSAO;
    private ShaderProgram shaderSSAOBlur;

    private int spheresVAO;
    private int quadVAO;
    private int positionsVBO;
    FloatBuffer positions;
    int atomCount;

    private FloatBuffer ssaoKernel;
    private int noiseTexture;

    public BoxRenderer(int initialAtomCount) {
        atomCount = initialAtomCount;
        positions = GLBuffers.newDirectFloatBuffer(initialAtomCount * 3);
    }

    private void makeFramebuffers(GL3 gl) {
        gBufferFBO = new FBObject();
        gBufferFBO.init(gl, width, height, 0);
        gPosition = gBufferFBO.attachTexture2D(gl,
                0,
                gl.GL_RGBA16F,
                gl.GL_RGBA,
                gl.GL_FLOAT,
                gl.GL_NEAREST,
                gl.GL_NEAREST,
                gl.GL_CLAMP_TO_EDGE,
                gl.GL_CLAMP_TO_EDGE);

        gNormal = gBufferFBO.attachTexture2D(gl,
                1,
                gl.GL_RGBA16F,
                gl.GL_RGBA,
                gl.GL_FLOAT,
                gl.GL_NEAREST,
                gl.GL_NEAREST,
                gl.GL_CLAMP_TO_EDGE,
                gl.GL_CLAMP_TO_EDGE);

        gAlbedoSpec = gBufferFBO.attachTexture2D(gl,
                2,
                gl.GL_RGBA,
                gl.GL_RGBA,
                gl.GL_UNSIGNED_BYTE,
                gl.GL_NEAREST,
                gl.GL_NEAREST,
                gl.GL_REPEAT,
                gl.GL_REPEAT);
        gBufferFBO.bind(gl);
        gl.glDrawBuffers(3, new int[]{ gl.GL_COLOR_ATTACHMENT0, gl.GL_COLOR_ATTACHMENT1, gl.GL_COLOR_ATTACHMENT2}, 0);

        gBufferFBO.attachRenderbuffer(gl, FBObject.Attachment.Type.DEPTH, FBObject.DEFAULT_BITS);
        System.out.println(gBufferFBO.getStatusString());
        if (!gBufferFBO.isStatusValid()) {
            throw new RuntimeException();
        }
        gBufferFBO.unbind(gl);

        ssaoFBO = new FBObject();
        ssaoFBO.init(gl, width, height, 0);
        ssaoFBO.bind(gl);
        ssaoColorBuffer = ssaoFBO.attachTexture2D(gl,
                0,
                gl.GL_RED,
                gl.GL_RED,
                gl.GL_FLOAT,
                gl.GL_NEAREST,
                gl.GL_NEAREST,
                gl.GL_CLAMP_TO_EDGE,
                gl.GL_CLAMP_TO_EDGE);

        ssaoBlurFBO = new FBObject();
        ssaoBlurFBO.init(gl, width, height, 0);
        ssaoBlurFBO.bind(gl);
        ssaoColorBufferBlur = ssaoBlurFBO.attachTexture2D(gl,
                0,
                gl.GL_RED,
                gl.GL_RED,
                gl.GL_FLOAT,
                gl.GL_NEAREST,
                gl.GL_NEAREST,
                gl.GL_CLAMP_TO_EDGE,
                gl.GL_CLAMP_TO_EDGE);
    }

    private void resizeFBOs(GL3 gl) {
        gBufferFBO.reset(gl, width, height, 0);
        ssaoFBO.reset(gl, width, height, 0);
        ssaoBlurFBO.reset(gl, width, height, 0);
    }

    private void makeShaders(GL3 gl) {
        ShaderCode spheresVert = ShaderCode.create(gl, gl.GL_VERTEX_SHADER, 1, this.getClass(), new String[]{"spheres.vert"}, false);
        ShaderCode spheresGeom = ShaderCode.create(gl, gl.GL_GEOMETRY_SHADER, 1, this.getClass(), new String[]{"spheres.geom"}, false);
        ShaderCode spheresFrag = ShaderCode.create(gl, gl.GL_FRAGMENT_SHADER, 1, this.getClass(), new String[]{"spheres.frag"}, false);
        this.shaderGeomPass = new ShaderProgram();
        shaderGeomPass.add(gl, spheresVert, System.err);
        shaderGeomPass.add(gl, spheresGeom, System.err);
        shaderGeomPass.add(gl, spheresFrag, System.err);
        shaderGeomPass.link(gl, System.err);

        ShaderProgram lightingShader = new ShaderProgram();
        ShaderCode vertLighting = ShaderCode.create(gl, gl.GL_VERTEX_SHADER, 1, this.getClass(), new String[]{"screenQuad.vert"}, false);
        ShaderCode fragLighting = ShaderCode.create(gl, gl.GL_FRAGMENT_SHADER, 1, this.getClass(), new String[]{"sphereLightingPass.frag"}, false);
        lightingShader.add(gl, vertLighting, System.err);
        lightingShader.add(gl, fragLighting, System.err);
        lightingShader.link(gl, System.err);
        this.shaderLightingPass = lightingShader;
        this.shaderLightingPass.useProgram(gl, true);
        setUniform(gl, shaderLightingPass, new GLUniformData("gPosition", 0));
        setUniform(gl, shaderLightingPass, new GLUniformData("gNormal", 1));
        setUniform(gl, shaderLightingPass, new GLUniformData("gAlbedoSpec", 2));
        setUniform(gl, shaderLightingPass, new GLUniformData("ssao", 3));

        this.shaderSSAO = new ShaderProgram();
        ShaderCode vertSSAO = ShaderCode.create(gl, gl.GL_VERTEX_SHADER, 1, this.getClass(), new String[]{"screenQuad.vert"}, false);
        ShaderCode fragSSAO = ShaderCode.create(gl, gl.GL_FRAGMENT_SHADER, 1, this.getClass(), new String[]{"ssao.frag"}, false);
        shaderSSAO.add(gl, vertSSAO, System.err);
        shaderSSAO.add(gl, fragSSAO, System.err);
        shaderSSAO.link(gl, System.err);
        this.shaderSSAO.useProgram(gl, true);
        setUniform(gl, shaderSSAO, new GLUniformData("gPosition", 0));
        setUniform(gl, shaderSSAO, new GLUniformData("gNormal", 1));
        setUniform(gl, shaderSSAO, new GLUniformData("texNoise", 2));

        this.shaderSSAOBlur = new ShaderProgram();
        ShaderCode vertSSAOBlur = ShaderCode.create(gl, gl.GL_VERTEX_SHADER, 1, this.getClass(), new String[]{"screenQuad.vert"}, false);
        ShaderCode fragSSAOBlur = ShaderCode.create(gl, gl.GL_FRAGMENT_SHADER, 1, this.getClass(), new String[]{"ssaoBlur.frag"}, false);
        shaderSSAOBlur.add(gl, vertSSAOBlur, System.err);
        shaderSSAOBlur.add(gl, fragSSAOBlur, System.err);
        shaderSSAOBlur.link(gl, System.err);
        this.shaderSSAOBlur.useProgram(gl, true);
        setUniform(gl, shaderSSAOBlur, new GLUniformData("ssaoInput", 0));
    }

    private void makeBuffers(GL3 gl) {
        gl.glGenVertexArrays(1, intBuf);
        this.spheresVAO = intBuf.get(0);
        gl.glBindVertexArray(spheresVAO);

        gl.glGenBuffers(1, intBuf);
        this.positionsVBO = intBuf.get(0);
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, positionsVBO);
        gl.glBufferData(gl.GL_ARRAY_BUFFER, atomCount * 3 * Float.BYTES, positions, gl.GL_STREAM_DRAW);
        gl.glVertexAttribPointer(0, 3, gl.GL_FLOAT, false, 3 * Float.BYTES, 0);
        gl.glEnableVertexAttribArray(0);

        gl.glBindVertexArray(0);

        gl.glGenVertexArrays(1, intBuf);
        this.quadVAO = intBuf.get(0);


    }

    private void setupSSAO(GL3 gl) {
        ssaoKernel = GLBuffers.newDirectFloatBuffer(64 * 3);
        Random rand = new Random();
        for (int i = 0; i < 64; i++) {
            Vector3f vec3 = new Vector3f(
                    rand.nextFloat() * 2 - 1,
                    rand.nextFloat() * 2 - 1,
                    rand.nextFloat()
            );
            vec3.normalize();
            vec3.mul(rand.nextFloat());
            float scale = (float) i / 64f;
            scale = 0.1f + (scale * scale) * (1.0f - 0.1f);
            vec3.mul(scale);
            vec3.get(i * 3, ssaoKernel);
        }

        FloatBuffer ssaoNoise = GLBuffers.newDirectFloatBuffer(16 * 3);
        for (int i = 0; i < 16; i++) {
            Vector3f vec3 = new Vector3f(
                    rand.nextFloat() * 2 - 1,
                    rand.nextFloat() * 2 - 1,
                    0.0f
            );
            vec3.get(i * 3, ssaoNoise);
        }

        gl.glGenTextures(1, intBuf);
        this.noiseTexture = intBuf.get(0);
        gl.glBindTexture(gl.GL_TEXTURE_2D, noiseTexture);
        gl.glTexImage2D(gl.GL_TEXTURE_2D, 0, gl.GL_RGBA16F, 4, 4, 0, gl.GL_RGB, gl.GL_FLOAT, ssaoNoise);
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST);
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST);
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_WRAP_S, gl.GL_REPEAT);
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_WRAP_T, gl.GL_REPEAT);
    }

    @Override
    public void init(GLAutoDrawable glAutoDrawable) {
        final GL3 gl = glAutoDrawable.getGL().getGL3();
        gl.glEnable(gl.GL_DEPTH_TEST);

        makeShaders(gl);
        makeFramebuffers(gl);
        makeBuffers(gl);
        setupSSAO(gl);

    }

    @Override
    public void dispose(GLAutoDrawable glAutoDrawable) {

    }

    @Override
    public void display(GLAutoDrawable glAutoDrawable) {
        final GL3 gl = glAutoDrawable.getGL().getGL3();

        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, this.positionsVBO);
        synchronized (this.positions) {
            gl.glBufferSubData(gl.GL_ARRAY_BUFFER, 0, positions.capacity() * Float.BYTES, positions);
        }

        // Geometry pass
        this.gBufferFBO.bind(gl);
        gl.glClearColor( 0.0f, 0.0f, 0.0f, 1f );
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT);

        Matrix4f view = new Matrix4f()
                .translation(0, 0, -zoom)
                .rotateX(pitch)
                .rotateY(yaw);
        Matrix4f projection = new Matrix4f()
                .perspective(((float) Math.toRadians(45.0)), ((float) width) / height, 0.1f, 100.0f);

        gl.glUseProgram(this.shaderGeomPass.program());
        setUniform(gl, this.shaderGeomPass, new GLUniformData("view", 4, 4, view.get(matBuffer)));
        setUniform(gl, this.shaderGeomPass, new GLUniformData("projection", 4, 4, projection.get(matBuffer)));
        setUniform(gl, this.shaderGeomPass, new GLUniformData("nearZ", 0.1f));

        setUniform(gl, this.shaderGeomPass, new GLUniformData("material.diffuse", 3, GLBuffers.newDirectFloatBuffer(new float[]{1.0f, 0.5f, 0.31f})));
        setUniform(gl, this.shaderGeomPass, new GLUniformData("material.specular", 3, GLBuffers.newDirectFloatBuffer(new float[]{0.1f, 0.1f, 0.1f})));

        gl.glBindVertexArray(this.spheresVAO);
        gl.glDrawArrays(gl.GL_POINTS, 0, this.atomCount);
        gl.glBindFramebuffer(gl.GL_FRAMEBUFFER, 0);

        // SSAO pass
        ssaoFBO.bind(gl);
        gl.glClear(gl.GL_COLOR_BUFFER_BIT);
        gl.glUseProgram(shaderSSAO.program());
        setUniform(gl, shaderSSAO, new GLUniformData("ssaoKernel", 3, this.ssaoKernel));
        setUniform(gl, shaderSSAO, new GLUniformData("projection", 4, 4, projection.get(matBuffer)));
        setUniform(gl, shaderSSAO, new GLUniformData("width", (float)width));
        setUniform(gl, shaderSSAO, new GLUniformData("height", (float)height));
        gl.glActiveTexture(gl.GL_TEXTURE0);
        gl.glBindTexture(gl.GL_TEXTURE_2D, gPosition.getName());
        gl.glActiveTexture(gl.GL_TEXTURE1);
        gl.glBindTexture(gl.GL_TEXTURE_2D, gNormal.getName());
        gl.glActiveTexture(gl.GL_TEXTURE2);
        gl.glBindTexture(gl.GL_TEXTURE_2D, noiseTexture);
        gl.glBindVertexArray(this.quadVAO);
        gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, 0, 4);

        ssaoBlurFBO.bind(gl);
        gl.glClear(gl.GL_COLOR_BUFFER_BIT);
        gl.glUseProgram(shaderSSAOBlur.program());
        gl.glActiveTexture(gl.GL_TEXTURE0);
        gl.glBindTexture(gl.GL_TEXTURE_2D, ssaoColorBuffer.getName());
        gl.glBindVertexArray(this.quadVAO);
        gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, 0, 4);

        // Lighting pass
        // -------------

        gl.glBindFramebuffer(gl.GL_FRAMEBUFFER, 0);
        gl.glClearColor( 0.0f, 0.0f, 0.0f, 1f );
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT);

        gl.glUseProgram(shaderLightingPass.program());

        gl.glActiveTexture(gl.GL_TEXTURE0);
        gl.glBindTexture(gl.GL_TEXTURE_2D, gPosition.getName());
        gl.glActiveTexture(gl.GL_TEXTURE1);
        gl.glBindTexture(gl.GL_TEXTURE_2D, gNormal.getName());
        gl.glActiveTexture(gl.GL_TEXTURE2);
        gl.glBindTexture(gl.GL_TEXTURE_2D, gAlbedoSpec.getName());
        gl.glActiveTexture(gl.GL_TEXTURE3);
        gl.glBindTexture(gl.GL_TEXTURE_2D, ssaoColorBufferBlur.getName());

        Vector3f lightDir = new Vector3f(.5f, .5f, -1);
        lightDir.mulDirection(view);
        setUniform(gl, this.shaderLightingPass, new GLUniformData("light.direction", 3, lightDir.get(GLBuffers.newDirectFloatBuffer(3))));
        setUniform(gl, this.shaderLightingPass, new GLUniformData("light.ambient", 3, GLBuffers.newDirectFloatBuffer(new float[]{0.9f, 0.9f, 0.9f})));
        setUniform(gl, this.shaderLightingPass, new GLUniformData("light.diffuse", 3, GLBuffers.newDirectFloatBuffer(new float[]{0.4f, 0.4f, 0.4f})));
        setUniform(gl, this.shaderLightingPass, new GLUniformData("light.specular", 3, GLBuffers.newDirectFloatBuffer(new float[]{0.9f, 0.9f, 0.9f})));

        // render quad
        gl.glBindVertexArray(this.quadVAO);
        gl.glDrawArrays(gl.GL_TRIANGLE_STRIP, 0, 4);
        gl.glBindVertexArray(0);
    }

    @Override
    public void reshape(GLAutoDrawable glAutoDrawable, int x, int y, int width, int height) {
        final GL3 gl = glAutoDrawable.getGL().getGL3();
        System.out.println(width);
        this.width = width;
        this.height = height;

        this.resizeFBOs(gl);

    }

    public static void setUniform(GL3 gl, ShaderProgram program, GLUniformData uniform) {
        int loc = gl.glGetUniformLocation(program.program(), uniform.getName());
        if (loc < 0) {
            throw new RuntimeException(uniform.toString());
        }
        uniform.setLocation(loc);
        gl.glUniform(uniform);
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        yaw += Math.toRadians((e.getX() - mouseX) * 1.01f);
        pitch += Math.toRadians((e.getY() - mouseY) * 1.01f);
        mouseX = e.getX();
        mouseY = e.getY();
    }

    @Override
    public void mouseMoved(MouseEvent e) {
        this.mouseX = e.getX();
        this.mouseY = e.getY();
    }

    @Override
    public void mouseWheelMoved(MouseWheelEvent e) {
        if (e.getPreciseWheelRotation() < 0) {
            zoom /= 1.05f;
        } else {
            zoom *= 1.05f;
        }
    }
}
