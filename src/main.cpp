#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#ifndef _GLIBCXX_USE_NANOSLEEP
#define _GLIBCXX_USE_NANOSLEEP
#endif
#include <thread>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Camera.h"
#include "GLSL.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Particle.h"
#include "Texture.h"

using namespace std;
using namespace Eigen;
using namespace glm;

bool keyToggles[256] = {false}; // only for English keyboards!

GLFWwindow *window; // Main application window
string RESOURCE_DIR = ""; // Where the resources are loaded from

enum types { none, low, high };
// white - none,none 
// red - none,low 
// lime - none,high 
// blue - low,none
// yellow - low,low
// aqua - low,high 
// magenta - high,none
// ??? - high,low
// ??? - high,high

shared_ptr<Program> progSimple;
shared_ptr<Program> prog;
shared_ptr<Camera> camera;
vector< shared_ptr<Particle> > particles;
shared_ptr<Texture> texture;
double t, h, e2, G;
int currInitV = none, currCollEff = none;
shared_ptr<Particle> sampleMeteor1;

Vector3f getColor() {
	Vector3f color;
	if (currInitV == none && currCollEff == none) {
		color << 255, 255, 255;
	} else if (currInitV == none && currCollEff == low) {
		color << 255, 0, 0;
	}
	else if (currInitV == none && currCollEff == high) {
		color << 0, 255, 0;
	}
	else if (currInitV == low && currCollEff == none) {
		color << 0, 0, 255;
	}
	else if (currInitV == low && currCollEff == low) {
		color << 255, 255, 0;
	}
	else if (currInitV == low && currCollEff == high) {
		color << 0, 255, 255;
	}
	else if (currInitV == high && currCollEff == none) {
		color << 255, 0, 255;
	}
	else if (currInitV == high && currCollEff == low) {
		color << 255, 128, 0;
	}
	else if (currInitV == high && currCollEff == high) {
		color << 127, 0, 255;
	}
	else {
		color << 255, 255, 255;
	}
	return color;
}

class Meteor {
public:
	Meteor(Vector3d p, Vector3d d) : meteorPart(), collisionEffect(currCollEff), active(true) {
		meteorPart = make_shared<Particle>();
		Vector3f color = getColor();
		//color << 255, 255, 255;
		
		meteorPart->setColor(color);
		meteorPart->setPosition(p);
		switch (currInitV) {
			case none:
				meteorPart->setVelocity(Vector3d::Zero());
				break;
			case low:
				meteorPart->setVelocity(d*0.1);
				break;
			case high:
				meteorPart->setVelocity(d*0.5);
				break;
		}
	}
	shared_ptr<Particle> meteorPart;
	int collisionEffect;
	bool active;
};
vector<shared_ptr<Meteor>> meteors;

// https://stackoverflow.com/questions/41470942/stop-infinite-loop-in-different-thread
std::atomic<bool> stop_flag;

static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
}

static void char_callback(GLFWwindow *window, unsigned int key)
{
	double xmouse, ymouse;
	glfwGetCursorPos(window, &xmouse, &ymouse);
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	keyToggles[key] = !keyToggles[key];
	if (key == 'e') {

		// --- Source: csce 441 lab 13 --- //
		auto MV = make_shared<MatrixStack>();
		MV->loadIdentity();
		camera->applyProjectionMatrix(MV);
		mat4 P = MV->topMatrix();
		MV->loadIdentity();
		camera->applyViewMatrix(MV);
		mat4 V = MV->topMatrix();
		mat4 C = inverse(V);

		vec3 glmP = vec3(C[3][0], C[3][1], C[3][2]);
		vec3 glmD = vec3((float)xmouse, (float)ymouse, 0.0);
		glmD.x = (2.0 * glmD.x / width) - 1.0;
		glmD.y = 1.0 - (2.0 * glmD.y / height);

		vec4 tempCoord = vec4(glmD.x, glmD.y, -1.0, 1.0);
		tempCoord = inverse(P) * tempCoord;
		tempCoord.w = 1.0;
		tempCoord = C * tempCoord;
		glmD = normalize(vec3(tempCoord) - glmP);
		// --- end --- //

		Vector3d pos;
		pos << glmP[0], glmP[1], glmP[2];
		Vector3d dir;
		dir << glmD[0], glmD[1], glmD[2];

		auto newMeteor = make_shared<Meteor>(pos, dir);
		meteors.push_back(newMeteor);
		newMeteor->meteorPart->init();
		
		/*cout << "New meteor created:" << endl;
		cout << pos << endl;*/
	}
	else if (key == 'v') {
		currInitV++;
		if (currInitV > high) {
			currInitV = none;
		}
		sampleMeteor1->setColor(getColor());
	}
	else if (key == 'b') {
		currCollEff++;
		if (currCollEff > high) {
			currCollEff = none;
		}
		sampleMeteor1->setColor(getColor());
	}
}

static void cursor_position_callback(GLFWwindow* window, double xmouse, double ymouse)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if(state == GLFW_PRESS) {
		camera->mouseMoved(xmouse, ymouse);
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the current mouse position.
	double xmouse, ymouse;
	glfwGetCursorPos(window, &xmouse, &ymouse);
	// Get current window size.
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	if(action == GLFW_PRESS) {
		bool shift = mods & GLFW_MOD_SHIFT;
		bool ctrl  = mods & GLFW_MOD_CONTROL;
		bool alt   = mods & GLFW_MOD_ALT;
		camera->mouseClicked(xmouse, ymouse, shift, ctrl, alt);
	}
}

static void initGL()
{
	GLSL::checkVersion();
	
	// Set background color
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	// Enable z-buffer test
	glEnable(GL_DEPTH_TEST);
	// Enable alpha blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	progSimple = make_shared<Program>();
	progSimple->setShaderNames(RESOURCE_DIR + "simple_vert.glsl", RESOURCE_DIR + "simple_frag.glsl");
	progSimple->setVerbose(false); // Set this to true when debugging.
	progSimple->init();
	progSimple->addUniform("P");
	progSimple->addUniform("MV");
	
	prog = make_shared<Program>();
	prog->setVerbose(true); // Set this to true when debugging.
	prog->setShaderNames(RESOURCE_DIR + "particle_vert.glsl", RESOURCE_DIR + "particle_frag.glsl");
	prog->init();
	prog->addUniform("P");
	prog->addUniform("MV");
	prog->addAttribute("aPos");
	prog->addAttribute("aTex");
	prog->addUniform("radius");
	prog->addUniform("alphaTexture");
	prog->addUniform("color");
	
	texture = make_shared<Texture>();
	texture->setFilename(RESOURCE_DIR + "alpha.jpg");
	texture->init();
	texture->setUnit(0);
	
	camera = make_shared<Camera>();
	
	// Initialize OpenGL for particles.
	for(auto p : particles) {
		p->init();
	}

	sampleMeteor1->init();
	/*sampleMeteor2->init(); */
	
	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

// Sort particles by their z values in camera space
class ParticleSorter {
public:
	bool operator()(size_t i0, size_t i1) const
	{
		// Particle positions in world space
		const Vector3d &x0 = particles[i0]->getPosition();
		const Vector3d &x1 = particles[i1]->getPosition();
		// Particle positions in camera space
		float z0 = V.row(2) * Vector4f(x0(0), x0(1), x0(2), 1.0f);
		float z1 = V.row(2) * Vector4f(x1(0), x1(1), x1(2), 1.0f);
		return z0 < z1;
	}
	
	void setViewMatrix(glm::mat4 V2)
	{
		for(int i = 0; i < 4; ++i) {
			for(int j = 0; j < 4; ++j) {
				V(i,j) = V2[j][i]; // indexing is different in Eigen and glm
			}
		}
	}
	
	Matrix4f V; // current view matrix
};
ParticleSorter sorter;

// http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
vector<size_t> sortIndices(const vector<T> &v) {
	// initialize original index locations
	vector<size_t> idx(v.size());
	for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(), sorter);
	return idx;
}

void renderGL()
{
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);
	
	// Use the window size for camera.
	glfwGetWindowSize(window, &width, &height);
	camera->setAspect((float)width/(float)height);
	
	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(keyToggles[(unsigned)'c']) {
		glEnable(GL_CULL_FACE);
	} else {
		glDisable(GL_CULL_FACE);
	}
	if(keyToggles[(unsigned)'l']) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	
	auto P = make_shared<MatrixStack>();
	auto MV = make_shared<MatrixStack>();
	
	// draw HUD meteor
	prog->bind();
	texture->bind(prog->getUniform("alphaTexture"));
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	if (keyToggles['h']) {
		sampleMeteor1->draw(prog, MV);
	}
	
	texture->unbind();
	prog->unbind();

	// Apply camera transforms
	P->pushMatrix();
	camera->applyProjectionMatrix(P);
	MV->pushMatrix();
	camera->applyViewMatrix(MV);
	// Set view matrix for the sorter
	sorter.setViewMatrix(MV->topMatrix());
	
	// Draw particles
	prog->bind();
	texture->bind(prog->getUniform("alphaTexture"));
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	// Sort particles by Z for transparency rendering.
	// Since we don't want to modify the contents of the vector, we compute the
	// sorted indices and traverse the particles in this sorted order.
	for(auto i : sortIndices(particles)) {
		particles[i]->draw(prog, MV);
	}
	texture->unbind();
	prog->unbind();

	// Draw meteors
	prog->bind();
	texture->bind(prog->getUniform("alphaTexture"));
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	// Sort particles by Z for transparency rendering.
	// Since we don't want to modify the contents of the vector, we compute the
	// sorted indices and traverse the particles in this sorted order.
	for (auto i : sortIndices(meteors)) {
		meteors[i]->meteorPart->draw(prog, MV);
	}
	texture->unbind();
	prog->unbind();
	
	//////////////////////////////////////////////////////
	// Cleanup
	//////////////////////////////////////////////////////
	
	// Pop stacks
	MV->popMatrix();
	P->popMatrix();
	
	GLSL::checkError(GET_FILE_LINE);
}

void saveParticles(const char *filename)
{
	ofstream out(filename);
	if(!out.good()) {
		cout << "Could not open " << filename << endl;
		return;
	}
	
	// 1st line:
	// <n> <h> <e2>
	out << particles.size() << " " << h << " " << " " << e2 << endl;

	// Rest of the lines:
	// <mass> <position> <velocity> <color> <radius>
	double m, r;
	Vector3d p, v;
	Vector3f c;

	for (int i = 0; i < particles.size(); i++) {
		auto temp = particles[i];
		m = temp->getMass();
		p = temp->getPosition();
		v = temp->getVelocity();
		c = temp->getColor();
		r = temp->getRadius();
		out << m << " ";
		out << p[0] << " " << p[1] << " " << p[2] << " ";
		out << v[0] << " " << v[1] << " " << v[2] << " ";
		out << c[0] << " " << c[1] << " " << c[2] << " ";
		out	<< r << endl;
	}
	
	out.close();
	cout << "Wrote galaxy to " << filename << endl;
}

void loadParticles(const char *filename)
{
	ifstream in;
	in.open(filename);
	if(!in.good()) {
		cout << "Cannot read " << filename << endl;
		return;
	}

	// 1st line:
	// <n> <h> <e2>
	int n;
	in >> n;
	in >> h;
	in >> e2;

	// Rest of the lines:
	// <mass> <position> <velocity> <color> <radius>
	double mass,p1,p2,p3,v1,v2,v3,r;
	float c1, c2, c3;

	for (int i = 0; i < n; i++) {
		in >> mass >> p1 >> p2 >> p3 >> v1 >> v2 >> v3 >> c1 >> c2 >> c3 >> r;
		auto temp = make_shared<Particle>();
		Vector3d p;
		p << p1, p2, p3;
		Vector3d v;
		v << v1, v2, v3;
		Vector3f color;
		color << c1, c2, c3;
		temp->setMass(mass);
		temp->setPosition(p);
		temp->setVelocity(v);
		temp->setColor(color);
		temp->setRadius(r);
		particles.push_back(temp);
	}

	in.close();
	cout << "Loaded galaxy from " << filename << endl;
}

void createParticles()
{
	srand(0);
	t = 0.0;
	h = 1.0;
	e2 = 1e-4;
	G = 1.0;
	
	Vector3d temp;
	Vector3f color;
	temp << 0.0, 0.0, 0.0;
	double setR = 3.0;
	double r = -setR;
	double a = 1.0;
	double y;

	auto tempSun = make_shared<Particle>();
	color << 255, 255, 0;
	tempSun->setMass(1e-3);
	tempSun->setPosition(temp);
	tempSun->setVelocity(temp);
	tempSun->setColor(color);
	tempSun->setRadius(0.5);
	particles.push_back(tempSun);

	int n = 1000;
	for (int i = 0; i < n; i++) {
		double sqr = G * 1e-3 * (2 / abs(r) - 1 / a);
		y = sqrt(abs(sqr));
		if (r < 0) {
			y *= -1.0;
		}

		auto tempStar = make_shared<Particle>(); 
		tempStar->setMass(1e-6);
		temp.Zero();
		temp << r, 0.0, 0.0;
		tempStar->setPosition(temp);
		temp.Zero();
		temp << 0.0, y, 0.0;
		tempStar->setVelocity(temp);
		particles.push_back(tempStar);

		if (i == n / 2) {
			r = 1.5;
		} else {
			r += setR / (double)n;
		}
		a = abs(r);
	}

	Vector3d HUDP;
	HUDP << 0.8, 0.8, 0.0;
	sampleMeteor1 = make_shared<Particle>();
	sampleMeteor1->setRadius(0.1); 
	sampleMeteor1->setPosition(HUDP);
	sampleMeteor1->setVelocity(Vector3d::Zero());
	sampleMeteor1->setColor(getColor());
}

void stepParticles()
{
	auto sun = particles[0];
	vector<Vector3d> forces;
	for (int i = 0; i < particles.size() + meteors.size(); i++) {
		forces.push_back(Vector3d::Zero());
	}

	bool isMeteor = false;
	for (int i = 0; i < particles.size() + meteors.size(); i++) {
		shared_ptr<Particle> obj;
		if (i < particles.size()) {
			obj = particles[i];
			isMeteor = false;
		}
		else {
			if (!meteors[i - particles.size()]->active) {
				continue;
			}
			obj = meteors[i - particles.size()]->meteorPart;
			isMeteor = true;
		}

		Vector3d fi = Vector3d::Zero();
		if (i != 0) {
			Eigen::Vector3d rij = sun->getPosition() - obj->getPosition();
			double fract = G * obj->getMass() * sun->getMass() / pow(rij.norm() * rij.norm() + e2, 3.0 / 2.0);
			fi = fi + fract * rij;

			// c = 1e1 "deletes" particles
			// c = 1e-5 nice explosion effect on impact
			// c = 1e-8 slower impact and less spread on impact
			// --- Source: Assignment 5 ---
			if (!isMeteor) {
				for (int j = 0; j < meteors.size(); j++) {
					auto meteor = meteors[j];
					if (!meteor->active) {
						continue;
					}
					auto meteorPart = meteor->meteorPart;
					Vector3d partP = obj->getPosition();
					Vector3d metP = meteorPart->getPosition();
					Vector3d deltaX = partP - metP;
					double l = deltaX.norm();
					double partR = obj->getRadius();
					double metR = meteorPart->getRadius();
					double d = partR + metR - l;
					if (d > 0) {
						Vector3d n = deltaX / l;
						double c;
						switch (meteor->collisionEffect) {
							case none:
								c = 1e1;
								break;
							case low:
								c = 1e-8;
								break;
							case high:
								c = 1e-5;
								break;
						}
						Vector3d fc = c * d * n;
						fi += fc;
						forces[i + j] -= fc;
					}
				}
			}
			// --- end ---
		}
		forces[i] = fi;
	}

	for (unsigned i = 1; i < particles.size() + meteors.size(); i++) {
		Vector3d fi = forces[i];
		shared_ptr<Particle> obj;
		if (i < particles.size()) {
			obj = particles[i];
		}
		else {
			if (!meteors[i - particles.size()]->active) {
				continue;
			}
			obj = meteors[i - particles.size()]->meteorPart;
		}

		Vector3d v = obj->getVelocity();
		double m = obj->getMass();
		Vector3d x = obj->getPosition();
		Vector3d newV = v + h / m * fi;
		Vector3d newX = x + h * newV;
		obj->setVelocity(newV);
		obj->setPosition(newX);

		double dist = (newX - Vector3d::Zero()).norm();
		if (dist > 50.0 && i >= particles.size()) {
			meteors[i - particles.size()]->active = false;
			meteors[i - particles.size()]->meteorPart->setColor(Vector3f::Zero());

		}
	}

	t += h;
}

void stepperFunc()
{
	while(!stop_flag) {
		if(keyToggles[(unsigned)' ']) {
			stepParticles();
		}
		this_thread::sleep_for(chrono::microseconds(1));
	}
}

int main(int argc, char **argv)
{
	if(argc != 2 && argc != 3) {
		// Wrong number of arguments
		cout << "Usage: Lab18 <RESOURCE_DIR> <(OPTIONAL) INPUT FILE>" << endl;
		cout << "   or: Lab18 <#steps>       <(OPTIONAL) INPUT FILE>" << endl;
		exit(0);
	}
	// Create the particles...
	if(argc == 2) {
		// ... without input file
		createParticles();
	} else {
		// ... with input file
		loadParticles(argv[2]);
	}
	// Try parsing `steps`
	int steps;
	if(sscanf(argv[1], "%i", &steps)) {
		// Success!
		cout << "Running without OpenGL for " << steps << " steps" << endl;
		// Run without OpenGL
		for(int k = 0; k < steps; ++k) {
			stepParticles();
		}
	} else {
		// `steps` could not be parsed
		cout << "Running with OpenGL" << endl;
		// Run with OpenGL until the window is closed
		RESOURCE_DIR = argv[1] + string("/");
		// Set error callback.
		glfwSetErrorCallback(error_callback);
		// Initialize the library.
		if(!glfwInit()) {
			return -1;
		}
		// Create a windowed mode window and its OpenGL context.
		window = glfwCreateWindow(640, 480, "Space Simulator", NULL, NULL);
		if(!window) {
			glfwTerminate();
			return -1;
		}
		// Make the window's context current.
		glfwMakeContextCurrent(window);
		// Initialize GLEW.
		glewExperimental = true;
		if(glewInit() != GLEW_OK) {
			cerr << "Failed to initialize GLEW" << endl;
			return -1;
		}
		glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
		cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
		cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
		// Set vsync.
		glfwSwapInterval(1);
		// Set keyboard callback.
		glfwSetKeyCallback(window, key_callback);
		// Set char callback.
		glfwSetCharCallback(window, char_callback);
		// Set cursor position callback.
		glfwSetCursorPosCallback(window, cursor_position_callback);
		// Set mouse button callback.
		glfwSetMouseButtonCallback(window, mouse_button_callback);
		// Initialize scene.
		initGL();
		// Start simulation thread.
		stop_flag = false;
		thread stepperThread(stepperFunc);
		// Loop until the user closes the window.
		while(!glfwWindowShouldClose(window)) {
			if(!glfwGetWindowAttrib(window, GLFW_ICONIFIED)) {
				// Render scene.
				renderGL();
				// Swap front and back buffers.
				glfwSwapBuffers(window);
			}
			// Poll for and process events.
			glfwPollEvents();
		}
		// Quit program.
		stop_flag = true;
		stepperThread.join();
		glfwDestroyWindow(window);
		glfwTerminate();
	}
	cout << "Elapsed time: " << (t*3.261539827498732e6) << " years" << endl;
	cout << "Total meteors created: " << meteors.size() << endl;
	//saveParticles("../resources/Output.txt");
	return 0;
}
