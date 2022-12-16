#version 120
varying vec2 vTex;
uniform sampler2D alphaTexture;
uniform vec4 color;

void main()
{
	float alpha = texture2D(alphaTexture, vTex).r;
	gl_FragColor = vec4(color.rgb, color.a*alpha);
}
