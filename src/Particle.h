#pragma once
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <memory>
#include <vector>
#include <list>
#include <cstdlib>
#include <cmath>
#include <limits>

#define GLEW_STATIC
#include <GL/glew.h>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

double randRange(double l, double h);
double generateGaussianNoise(double mu, double sigma);

class MatrixStack;
class Program;
class Texture;

class Particle
{
public:
	Particle();
	virtual ~Particle();

	// OpenGL methods
	void init();
	void draw(std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> MV) const;
	
	// Getters
	double getMass() const { return m; }
	const Eigen::Vector3d& getPosition() const { return x; }
	const Eigen::Vector3d& getVelocity() const { return v; }
	Eigen::Vector3f getColor() const { return color.segment<3>(0); }
	float getRadius() const { return radius; }
	
	// Setters
	void setMass(double m) { this->m = m; }
	void setPosition(const Eigen::Vector3d x) { this->x = x; }
	void setVelocity(const Eigen::Vector3d v) { this->v = v; }
	void setColor(const Eigen::Vector3f &c) { color.segment<3>(0) = c; }
	void setRadius(float r) { radius = r; }
	
private:
	// For physics
	double m; // mass
	Eigen::Vector3d x; // position
	Eigen::Vector3d v; // velocity
	Eigen::Vector3d a; // aceleration (force/mass)
	
	// For display only
	float radius;
	Eigen::Vector4f color;
	std::vector<float> posBuf;
	std::vector<float> texBuf;
	std::vector<unsigned int> indBuf;
	GLuint posBufID;
	GLuint texBufID;
	GLuint indBufID;
};

#endif
