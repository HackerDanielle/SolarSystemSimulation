#version 120
attribute vec4 aPos;
attribute vec2 aTex;
uniform mat4 P;
uniform mat4 MV;
uniform float radius;
varying vec2 vTex;

void main()
{
	// Billboarding: set the upper 3x3 to be the identity matrix
	mat4 MV0 = MV;
	MV0[0] = vec4(1.0, 0.0, 0.0, 0.0);
	MV0[1] = vec4(0.0, 1.0, 0.0, 0.0);
	MV0[2] = vec4(0.0, 0.0, 1.0, 0.0);
	gl_Position = P * MV0 * vec4(radius * aPos.xy, 0.0, 1.0);
	vTex = aTex;
}
