#version 150

in vec3 position;
in vec3 positionLeft;
in vec3 positionRight;
in vec3 positionUp;
in vec3 positionDown;
in vec4 color;
out vec4 col;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

uniform int mode = 0; // 0 default mode, 1 smooth mode
uniform float scale;
uniform float exponent;
uniform int jetColor = 0;



vec4 JetColorMap(float x)
{
    float a; // alpha
    float r, g, b;
    if (x < 0)
    {
        r = 0;
        g = 0;
        b = 0;
    }
    else if (x < 0.125)
    {
        a = x / 0.125;
        r = 0;
        g = 0;
        b = 0.5 + 0.5 * a;
    }
    else if (x < 0.375)
    {
        a = (x - 0.125) / 0.25;
        r = 0;
        g = a;
        b = 1;
    }
    else if (x < 0.625)
    {
        a = (x - 0.375) / 0.25;
        r = a;
        g = 1;
        b = 1 - a;
    }
    else if (x < 0.875)
    {
        a = (x - 0.625) / 0.25;
        r = 1;
        g = 1 - a;
        b = 0;
    }
    else if (x <= 1.0)
    {
        a = (x - 0.875) / 0.125;
        r = 1 - 0.5 * a;
        g = 0;
        b = 0;
    }
    else
    {
        r = 1;
        g = 1;
        b = 1;
    }
    return vec4(r, g, b, 1.0);
}


void main()
{

	if (mode == 0) {
		gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0f);
		col = color;
	}
	else {
		vec3 newPosition = (position + positionLeft + positionRight + positionUp + positionDown) / 5.0;
		col = vec4(pow(newPosition.y, exponent), pow(newPosition.y, exponent), pow(newPosition.y, exponent), 1.0);
        if (jetColor == 1) { // to make the terrian colourful
            col = JetColorMap(newPosition.y);
        }
		newPosition.y = scale * pow(newPosition.y, exponent);
		gl_Position = projectionMatrix * modelViewMatrix * vec4(newPosition, 1.0f);
	}	
}


