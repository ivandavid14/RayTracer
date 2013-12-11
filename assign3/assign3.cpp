/*
CSCI 480
Assignment 3 Raytracer

Name: <Your name here>
*/

/*
	TODO:	
		-use tangent to calculate fov rays rather than my way, just to be safe
			even though my way works, however its hard to change
		-
*/

#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>

#include <stdio.h>
#include <string>
#include <iostream>
#include <random>
#include <vector>
#include <time.h>
#include <random>
#include <exception>

#define MAX_TRIANGLES 100000
#define MAX_SPHERES 1000
#define MAX_LIGHTS 100

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

int maxBounces = 5;
bool antiAliasing = false;
double shadowRayCount = 100;

//the field of view of the camera
#define fov 50.0

#define PI 3.14159265

unsigned char buffer[HEIGHT][WIDTH][3];

enum intersectObject {tri, sph, light, none};

double lightRadius = 0.25;

//http://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
double randDouble(double lowerBound, double upperBound)
{
	double temp = ((double)rand())/(RAND_MAX);
	return lowerBound + temp * (upperBound - lowerBound);
}

double radians(double degrees)
{
	double temp = 2 * PI / 360.0;
	return tan(temp * degrees);
}

struct Vertex
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double normal[3];
	double shininess;
};

typedef struct _Triangle
{
	struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
	double position[3];
	double color_diffuse[3];
	double color_specular[3];
	double shininess;
	double radius;
} Sphere;

typedef struct _Light
{
	double position[3];
	double color[3];

	_Light()
	{

	}

	_Light(double pos[3], double col[3])
	{
		for(int i = 0; i < 3; i++)
		{
			position[i] = pos[i];
			color[i] = col[i];
		}
	}
} Light;

struct vector 
{
	double x;
	double y;
	double z;

	vector()
	{
		x = y = z = 0;
	}

	vector(double a[3], double b[3])
	{
		x = a[0] - b[0];
		y = a[1] - b[1];
		z = a[2] - b[2];
	}

	vector(double x_, double y_, double z_)
	{
		x = x_;
		y = y_;
		z = z_;
	}

	void normalize()
	{
		double mag = x * x + y * y + z * z;
		mag = sqrt(mag);
		x /= mag;
		y /= mag;
		z /= mag;
	}

	double length()
	{
		double mag = x * x + y * y + z * z;
		mag = sqrt(mag);
		return mag;
	}

	void print()
	{
		std::cout<<"x: "<<this->x<<" y: "<<this->y<<" z: "<<this->z<<std::endl;
	}
};

double dotProduct(vector a, vector b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

vector cross(vector a, vector b)
{
	vector temp;
	temp.x = a.y * b.z - a.z * b.y;
	temp.y = -(a.x * b.z - a.z * b.x);
	temp.z = a.x * b.y - a.y * b.x;

	return temp;
}

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

double sphereIntersection(vector v, double start[3], int sphereIndex)
{
	v.normalize();
	Sphere temp = spheres[sphereIndex];
	double t0;
	double t1;
	double b;
	double c;

	//camera at (0, 0, 0) so no need for x0, y0, z0
	b = 2 * (v.x * (start[0] - temp.position[0]) + 
		v.y * (start[1] - temp.position[1]) + 
		v.z * (start[2] - temp.position[2]));
	c = pow(start[0] - temp.position[0], 2) + 
		pow(start[1] - temp.position[1], 2) + 
		pow(start[2] - temp.position[2], 2) - pow(temp.radius, 2);
	double l = pow(b, 2) - 4 * c;
	if(l < 0)
	{
		return -1;
	}
	t0 = -b - sqrt(pow(b, 2) - 4 * c);
	t0 /= 2;
	t1 = -b + sqrt(pow(b, 2) - 4 * c);
	t1 /= 2;

	if(min(t0, t1) > 0)
	{
		return min(t0, t1);
	}
	if(t0 <= 0 && t1 <= 0)
	{
		return -1;
	}
	if(t0 > 0)
	{
		return t0;
	}
	else
	{
		return t1;
	}
}

double triangleIntersectionRedone(vector v, double start[3], int triangleIndex, double &alpha, double &beta,
	double &gamma, vector &normal)
{
	//normalize v just in case i forget to do it outside of the function\

	//step 1: generate normal for the plane
	vector v1;
	vector v2;
	
	Vertex A = triangles[triangleIndex].v[0];
	Vertex B = triangles[triangleIndex].v[1];
	Vertex C = triangles[triangleIndex].v[2];


	//v1 is B-A
	v1.x = B.position[0] - A.position[0];
	v1.y = B.position[1] - A.position[1];
	v1.z = B.position[2] - A.position[2];

	//v2 is C-A
	v2.x = C.position[0] - A.position[0];
	v2.y = C.position[1] - A.position[1];
	v2.z = C.position[2] - A.position[2];

	normal = cross(v1, v2);

	normal.normalize();

	double d = normal.x * A.position[0] + normal.y * A.position[1] + normal.z * A.position[2];
	
	double t = d - (normal.x * start[0] + normal.y * start[1] + normal.z * start[2]);

	double denom = normal.x * v.x + normal.y * v.y + normal.z * v.z;
	if(denom == 0)
	{
		return -1;
	}

	t /= denom;

	if(t <= 0)
	{
		return -1;
	}

	double intersectionPoint[3] = {start[0] + v.x * t, start[1] + v.y * t, start[2] + v.z * t};

	//use winding method to determine inside/out
	vector Q_A(intersectionPoint[0] - A.position[0], intersectionPoint[1] - A.position[1],
		intersectionPoint[2] - A.position[2]);

	vector Q_B(intersectionPoint[0] - B.position[0], intersectionPoint[1] - B.position[1],
		intersectionPoint[2] - B.position[2]);

	vector Q_C(intersectionPoint[0] - C.position[0], intersectionPoint[1] - C.position[1],
		intersectionPoint[2] - C.position[2]);

	vector B_A(B.position[0] - A.position[0], B.position[1] - A.position[1], B.position[2] - A.position[2]);
	
	vector C_B(C.position[0] - B.position[0], C.position[1] - B.position[1], C.position[2] - B.position[2]);

	vector A_C(A.position[0] - C.position[0], A.position[1] - C.position[1], A.position[2] - C.position[2]);

	vector C_A(C.position[0] - A.position[0], C.position[1] - A.position[1], C.position[2] - A.position[2]);

	vector alphaComp = cross(C_B, Q_B);
	vector betaComp = cross(A_C, Q_C);
	vector gammaComp = cross(B_A, Q_A);

	if(dotProduct(alphaComp, normal) < 0)
	{
		return -1;
	}

	if(dotProduct(betaComp, normal) < 0)
	{
		return -1;
	}

	if(dotProduct(gammaComp, normal) < 0)
	{
		return -1;
	}

	//get alpha, beta, gamma values
	vector botVec = cross(B_A, C_A);

	alpha = dotProduct(alphaComp, normal)/dotProduct(botVec, normal);
	beta = dotProduct(betaComp, normal)/dotProduct(botVec, normal);
	gamma = dotProduct(gammaComp, normal)/dotProduct(botVec, normal);

	return t;
}

//same as other one, except used after we have determined that the particular sphere
//is the first object to intersect the ray from the camera
void IntersectionRecursiveRedone(vector v, double startPosition[3], double cIntensity[3], int objIndex_, intersectObject obj_, vector normal, int bounces = 0)
{
	if(bounces > maxBounces)
	{
		return;
	}

	double alpha, beta, gamma;
	intersectObject obj = none;
	int objIndex = 0;
	double minT = -1;
	std::vector<Light> lightList;
	vector n_;
	vector attenuation(0.1, 0.05, .005);

	//determines if there is a sphere that is in the way and closer than anything
	for(int i = 0; i < num_spheres; i++)
	{
		//do intersection code
		double a = sphereIntersection(v, startPosition, i);
		if(a != -1 && minT == -1 && !(objIndex_ == i && obj_ == sph))
		{
			minT = a;
			objIndex = i;
			obj = sph;
		}
		else if(a != -1 && minT != -1 && !(objIndex_ == i && obj_ == sph))
		{
			if(a < minT)
			{
				minT = a;
				objIndex = i;
				obj = sph;
			}
		}
	}

	//determines if there is a triangle that is in the way and closer than anything
	for(int i = 0; i < num_triangles; i++)
	{
		double alphaTemp, betaTemp, gammaTemp;
		vector nTemp;
		double a = triangleIntersectionRedone(v, startPosition, i, alphaTemp, betaTemp, gammaTemp, nTemp);
		if(a != -1 && minT == -1 && !(objIndex_ == i && obj_ == tri))
		{
			obj = tri;
			objIndex = i;
			minT = a;
			alpha = alphaTemp;
			beta = betaTemp;
			gamma = gammaTemp;
			n_.x = nTemp.x;
			n_.y = nTemp.y;
			n_.z = nTemp.z;
		}
		else if(a != -1 && minT != -1 && !(objIndex_ == i && obj_ == tri))
		{
			if(a < minT)
			{
				minT = a;
				objIndex = i;
				obj = tri;
				alpha = alphaTemp;
				beta = betaTemp;
				gamma = gammaTemp;
				n_.x = nTemp.x;
				n_.y = nTemp.y;
				n_.z = nTemp.z;
			}
		}
	}

	if(obj == tri)
	{
		double start[3] = {startPosition[0] + v.x * minT, startPosition[1] + v.y * minT, startPosition[2] + v.z * minT};

		vector v_(-v.x, -v.y, -v.z);
		v_.normalize();

		double dot = dotProduct(v_, n_);

		vector r;
		r.x = 2 * dot * n_.x - v_.x;
		r.y = 2 * dot * n_.y - v_.y;
		r.z = 2 * dot * n_.z - v_.z;

		r.normalize();

		IntersectionRecursiveRedone(r, start, cIntensity, objIndex, obj, n_, bounces + 1);

		//check to see which lights this object is exposed to
		for(int i = 0; i < num_lights; i++)
		{
			for(int k = 0; k < shadowRayCount; k++)
			{
				double distanceToLight;
				double distanceToObject;

				double q, f;
				q = 2 * PI * randDouble(0, 1);
				f = std::acos(2 * randDouble(0, 1) - 1);

				vector onSphere;

				onSphere.x = lightRadius * std::cos(q) * std::sin(f) + lights[i].position[0];
				onSphere.y = lightRadius * std::sin(q) * std::sin(f) + lights[i].position[1];
				onSphere.z = lightRadius * std::cos(f) + lights[i].position[2];

				vector lightV(onSphere.x - start[0], onSphere.y - start[1], onSphere.z - start[2]);
				distanceToLight = lightV.length();
				lightV.normalize();
				if(dotProduct(n_, lightV) < 0)
				{
					continue;
				}
			
				//check to see if there is a sphere in the way of the light
				bool sphereInFront = false;
				for(int j = 0; j < num_spheres; j++)
				{
					double t = sphereIntersection(lightV, start, j);
					if(t == -1)
					{
						continue;
					}
					vector temp(lightV.x * t, lightV.y * t, lightV.z * t);
					double distance = temp.length();
					if(distance < distanceToLight)
					{
						sphereInFront = true;
						break;
					}
				}
				if(sphereInFront)
				{
					continue;
				}
				//check to see if there is a triangle in the way of the light
				bool triangleInFront = false;
				for(int j = 0; j < num_triangles; j++)
				{
					if(j == objIndex)
						continue;
					double alphaTemp, betaTemp, gammaTemp;
					vector nTemp;
					double t = triangleIntersectionRedone(lightV, start, j, alphaTemp, betaTemp, gammaTemp, nTemp); 
					if(t == -1)
					{
						continue;
					}
					vector temp(lightV.x * t, lightV.y * t, lightV.z * t);
					double distance = temp.length();
					if(distance < distanceToLight)
					{
						triangleInFront = true;
						break;
					}
				}
				if(triangleInFront)
				{
					continue;
				}

				//if we've reached this point, that no triangle or sphere is in the way of the light, and its in view
				Light l;
				l.color[0] = lights[i].color[0];
				l.color[1] = lights[i].color[1];
				l.color[2] = lights[i].color[2];
				l.position[0] = onSphere.x;
				l.position[1] = onSphere.y;
				l.position[2] = onSphere.z;
				lightList.push_back(l);
			}
		}
		//n is normal at the particular point we are calculating things at
		//interpolate using barycentric coordinates
		vector n;

		n.x = alpha * triangles[objIndex].v[0].normal[0] + 
				beta * triangles[objIndex].v[1].normal[0] + 
				gamma * triangles[objIndex].v[2].normal[0];

		n.y = alpha * triangles[objIndex].v[0].normal[1] + 
				beta * triangles[objIndex].v[1].normal[1] + 
				gamma * triangles[objIndex].v[2].normal[1];

		n.z = alpha * triangles[objIndex].v[0].normal[2] + 
				beta * triangles[objIndex].v[1].normal[2] + 
				gamma * triangles[objIndex].v[2].normal[2];

		n.normalize();

		double colors[3] = {0, 0, 0};

		if(lightList.size() == 0)
		{ 
			cIntensity[0] = min(cIntensity[0] + ambient_light[0], 1);
			cIntensity[1] = min(cIntensity[1] + ambient_light[1], 1);
			cIntensity[2] = min(cIntensity[2] + ambient_light[2], 1);
			return;
		}

		std::vector<Light>::iterator it = lightList.begin();
		std::vector<Light>::iterator end = lightList.end();

		for( ; it != end; it++)
		{
			//vector to light
			vector l(it->position[0] - start[0], 
					 it->position[1] - start[1], 
					 it->position[2] - start[2]);

			double distanceFromLight = l.length();

			double att = ((double)1)/(attenuation.x + attenuation.y * distanceFromLight + attenuation.z * pow(distanceFromLight, 2));

			l.normalize();
			//v_ is vector TO the viewer, already flipped to him. dont do it again
			

			//now calc light reflection vector
			double L_N = dotProduct(l, n);
			vector r_;
			r_.x = 2 * L_N * n.x - l.x;
			r_.y = 2 * L_N * n.y - l.y;
			r_.z = 2 * L_N * n.z - l.z;

			r_.normalize();

			double r_v = dotProduct(v_, r_);

			colors[0] += att * (alpha * (triangles[objIndex].v[0].color_diffuse[0] * it->color[0] * max(L_N, 0) + 
					(triangles[objIndex].v[0].color_specular[0]) * it->color[0] * max(pow(r_v, triangles[objIndex].v[0].shininess), 0)) +
							     beta * (triangles[objIndex].v[1].color_diffuse[0] * it->color[0] * max(L_N, 0) + 
					(triangles[objIndex].v[1].color_specular[0]) * it->color[0] * max(pow(r_v, triangles[objIndex].v[1].shininess), 0)) + 
								 gamma * (triangles[objIndex].v[2].color_diffuse[0] * it->color[0] * max(L_N, 0) + 
					(triangles[objIndex].v[2].color_specular[0] * it->color[0]) * max(pow(r_v, triangles[objIndex].v[2].shininess), 0)));

			colors[1] += att * (alpha * (triangles[objIndex].v[0].color_diffuse[1] * it->color[1] * max(L_N, 0) + 
					(triangles[objIndex].v[0].color_specular[1]) * it->color[1] * max(pow(r_v, triangles[objIndex].v[0].shininess), 0)) +
							     beta * (triangles[objIndex].v[1].color_diffuse[1] * it->color[1] * max(L_N, 0) + 
					(triangles[objIndex].v[1].color_specular[1]) * it->color[1] * max(pow(r_v, triangles[objIndex].v[1].shininess), 0)) + 
								 gamma * (triangles[objIndex].v[2].color_diffuse[1] * it->color[1] * max(L_N, 0) + 
					(triangles[objIndex].v[2].color_specular[1]) * it->color[1] * max(pow(r_v, triangles[objIndex].v[2].shininess), 0)));

			colors[2] += att * (alpha * (triangles[objIndex].v[0].color_diffuse[2] * it->color[2] * max(L_N, 0) + 
					(triangles[objIndex].v[0].color_specular[2]) * it->color[2] * max(pow(r_v, triangles[objIndex].v[0].shininess), 0)) +
							     beta * (triangles[objIndex].v[1].color_diffuse[2] * it->color[2] * max(L_N, 0) + 
					(triangles[objIndex].v[1].color_specular[2]) * it->color[2] * max(pow(r_v, triangles[objIndex].v[1].shininess), 0)) + 
								 gamma * (triangles[objIndex].v[2].color_diffuse[2] * it->color[2] * max(L_N, 0) + 
					(triangles[objIndex].v[2].color_specular[2]) * it->color[2] * max(pow(r_v, triangles[objIndex].v[2].shininess), 0)));
		}

		for(int k = 0; k < 3; k++)
		{
			colors[k] /= shadowRayCount;
		}

		double finalColor[3] = {0, 0, 0};

		vector l(r.x, r.y, r.z);

		double distanceFromLight = l.length();

		l.normalize();
		//v_ is vector TO the viewer, already flipped to him. dont do it again
			

		//now calc light reflection vector
		double L_N = dotProduct(l, n);
		vector r_;
		r_.x = 2 * L_N * n.x - l.x;
		r_.y = 2 * L_N * n.y - l.y;
		r_.z = 2 * L_N * n.z - l.z;

		r_.normalize();

		double r_v = dotProduct(v_, r_);

		double rSpec = alpha * triangles[objIndex].v[0].color_specular[0] * max(pow(r_v, triangles[objIndex].v[0].shininess), 0) +
					   beta * triangles[objIndex].v[1].color_specular[0] * max(pow(r_v, triangles[objIndex].v[1].shininess), 0) +
					   gamma * triangles[objIndex].v[2].color_specular[0] * max(pow(r_v, triangles[objIndex].v[2].shininess), 0);

		double gSpec = alpha * triangles[objIndex].v[0].color_specular[1] * max(pow(r_v, triangles[objIndex].v[0].shininess), 0) +
					   beta * triangles[objIndex].v[1].color_specular[1] * max(pow(r_v, triangles[objIndex].v[1].shininess), 0) +
					   gamma * triangles[objIndex].v[2].color_specular[1] * max(pow(r_v, triangles[objIndex].v[2].shininess), 0);

		double bSpec = alpha * triangles[objIndex].v[0].color_specular[2] * max(pow(r_v, triangles[objIndex].v[0].shininess), 0) +
					   beta * triangles[objIndex].v[1].color_specular[2] * max(pow(r_v, triangles[objIndex].v[1].shininess), 0) +
				       gamma * triangles[objIndex].v[2].color_specular[2] * max(pow(r_v, triangles[objIndex].v[2].shininess), 0);

		//need to calculate l reflection vector for this, i can't just straight up add
		//the reflection light with a spec component, its over powering, plus it doesn't make sense to apply it
		//to the lights, but not the reflection, which is essentially a light of a very 
		//specific color.

		finalColor[0] = (1 - rSpec) * colors[0] + (rSpec) * cIntensity[0];
		finalColor[1] = (1 - gSpec) * colors[1] + (gSpec) * cIntensity[1];
		finalColor[2] = (1 - bSpec) * colors[2] + (bSpec) * cIntensity[2];

		cIntensity[0] = min(finalColor[0] + ambient_light[0], 1);
		cIntensity[1] = min(finalColor[1] + ambient_light[1], 1);
		cIntensity[2] = min(finalColor[2] + ambient_light[2], 1);
	}
	else if(obj == sph)
	{
		double start[3] = {startPosition[0] + v.x * minT, startPosition[1] + v.y * minT, startPosition[2] + v.z * minT};

		vector normal(start[0] - spheres[objIndex].position[0], 
					  start[1] - spheres[objIndex].position[1],
					  start[2] - spheres[objIndex].position[2]);

		normal.normalize();

		vector v_(-v.x, -v.y, -v.z);
		v_.normalize();
		
		double dot = dotProduct(v_, normal);

		vector r;
		r.x = 2 * dot * normal.x - v_.x;
		r.y = 2 * dot * normal.y - v_.y;
		r.z = 2 * dot * normal.z - v_.z;

		r.normalize();

		IntersectionRecursiveRedone(r, start, cIntensity, objIndex, obj, normal, bounces + 1);


		//check to see which lights this object is exposed to
		for(int i = 0; i < num_lights; i++)
		{
			for(int k = 0; k < shadowRayCount; k++)
			{
				double distanceToLight;
				double distanceToObject;
				double q, f;
				q = 2 * PI * randDouble(0, 1);
				f = std::acos(2 * randDouble(0, 1) - 1);

				vector onSphere;

				onSphere.x = lightRadius * std::cos(q) * std::sin(f) + lights[i].position[0];
				onSphere.y = lightRadius * std::sin(q) * std::sin(f) + lights[i].position[1];
				onSphere.z = lightRadius * std::cos(f) + lights[i].position[2];

				vector lightV(onSphere.x - start[0], onSphere.y - start[1], onSphere.z - start[2]);
				distanceToLight = lightV.length();
				lightV.normalize();

				if(dotProduct(lightV, normal) <= 0)
				{
					continue;
				}
			
				//check to see if there is a sphere in the way of the light
				bool sphereInFront = false;
				for(int j = 0; j < num_spheres; j++)
				{
					if(j == objIndex)
						continue;
					double t = sphereIntersection(lightV, start, j);
					if(t == -1)
					{
						continue;
					}
					//vector temp(start[0] + lightV.x * t, start[1] + lightV.y * t, start[2] + lightV.z * t);
					vector temp(lightV.x * t, lightV.y * t, lightV.z * t);
					double distance = temp.length();
					if(distance < distanceToLight)
					{
						sphereInFront = true;
						break;
					}
				}
				if(sphereInFront)
				{
					continue;
				}

				//check to see if there is a triangle in the way of the light
				bool triangleInFront = false;
				for(int j = 0; j < num_triangles; j++)
				{
					double alphaTemp, betaTemp, gammaTemp;
					vector nTemp;
					double t = triangleIntersectionRedone(lightV, start, j, alphaTemp, betaTemp, gammaTemp, nTemp); 
					if(t == -1)
					{
						continue;
					}
					//vector temp(start[0] + lightV.x * t, start[1] + lightV.y * t, start[2] + lightV.z * t);
					vector temp(lightV.x * t, lightV.y * t, lightV.z * t);
					double distance = temp.length();
					if(distance < distanceToLight)
					{
						triangleInFront = true;
						break;
					}
				}
				if(triangleInFront)
				{
					continue;
				}

				//if we've reached this point, that no triangle or sphere is in the way of the light, and its in view
				lightList.push_back(lights[i]);
			}
		}

		double colors[3] = {0, 0, 0};

		if(lightList.size() == 0)
		{
			cIntensity[0] = min(cIntensity[0] + ambient_light[0], 1);
			cIntensity[1] = min(cIntensity[1] + ambient_light[1], 1);
			cIntensity[2] = min(cIntensity[2] + ambient_light[2], 1);
			return;
		}

		std::vector<Light>::iterator it = lightList.begin();
		std::vector<Light>::iterator end = lightList.end();

		for( ; it != end; it++)
		{
			vector l(it->position[0] - start[0], 
					 it->position[1] - start[1], 
					 it->position[2] - start[2]);

			double distanceFromLight = l.length();

			double att = ((double)1)/(attenuation.x + attenuation.y * distanceFromLight + attenuation.z * pow(distanceFromLight, 2));

			l.normalize();
			//v_ is vector TO the viewer, already flipped to him. dont do it again
			

			//now calc light reflection vector
			double L_N = dotProduct(l, normal);
			vector r_;
			r_.x = 2 * L_N * normal.x - l.x;
			r_.y = 2 * L_N * normal.y - l.y;
			r_.z = 2 * L_N * normal.z - l.z;

			r_.normalize();

			double r_v = dotProduct(v_, r_);
			//now do lighting calcs with r_, v_, n, & l
			colors[0] += att * (spheres[objIndex].color_diffuse[0] * it->color[0] * max(L_N, 0) +
					(spheres[objIndex].color_specular[0]) * it->color[0] * max(pow(r_v, spheres[objIndex].shininess), 0));

			colors[1] += att * (spheres[objIndex].color_diffuse[1] * it->color[1] * max(L_N, 0) +
				(spheres[objIndex].color_specular[1]) * it->color[1] * max(pow(r_v, spheres[objIndex].shininess), 0));

			colors[2] += att * (spheres[objIndex].color_diffuse[2] * it->color[2] * max(L_N, 0) +
				(spheres[objIndex].color_specular[2]) * it->color[2] * max(pow(r_v, spheres[objIndex].shininess), 0));
		}

		for(int k = 0; k < 3; k++)
		{
			colors[k] /= shadowRayCount;
		}

		double finalColor[3] = {0, 0, 0};

		double rSpec = spheres[objIndex].color_specular[0];

		double gSpec = spheres[objIndex].color_specular[1];

		double bSpec = spheres[objIndex].color_specular[2];

		//need to calculate l reflection vector for this, i can't just straight up add
		//the reflection light with a spec component, its over powering, plus it doesn't make sense to apply it
		//to the lights, but not the reflection, which is essentially a light of a very 
		//specific color. 

		vector l(r.x, r.y, r.z);

		l.normalize();
		//v_ is vector TO the viewer, already flipped to him. dont do it again
			

		//now calc light reflection vector
		double L_N = dotProduct(l, normal);
		vector r_;
		r_.x = 2 * L_N * normal.x - l.x;
		r_.y = 2 * L_N * normal.y - l.y;
		r_.z = 2 * L_N * normal.z - l.z;

		r_.normalize();
		
		double r_v = dotProduct(v_, r_);

		finalColor[0] = (1 - rSpec) * colors[0] + (rSpec) * cIntensity[0] * max(pow(r_v, spheres[objIndex].shininess), 0);
		finalColor[1] = (1 - gSpec) * colors[1] + (gSpec) * cIntensity[1] * max(pow(r_v, spheres[objIndex].shininess), 0);
		finalColor[2] = (1 - bSpec) * colors[2] + (bSpec) * cIntensity[2] * max(pow(r_v, spheres[objIndex].shininess), 0);

		cIntensity[0] = min(finalColor[0] + ambient_light[0], 1);
		cIntensity[1] = min(finalColor[1] + ambient_light[1], 1);
		cIntensity[2] = min(finalColor[2] + ambient_light[2], 1);
	}
	else if(bounces == 0)
	{
		cIntensity[0] = 1;
		cIntensity[1] = 1;
		cIntensity[2] = 1;
	}
	return;
}

//MODIFY THIS FUNCTION
void draw_scene()
{
	unsigned int x,y;
	//simple output
	for(x=0;x < WIDTH;x++)
	{
		bool print = true;
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for(y=0;y < HEIGHT;y++)
		{
			double a = ((double)WIDTH/HEIGHT);
			double xLow = -a * std::tan(radians(fov)/2);
			double yHigh = -std::tan(radians(fov)/2);
			if(antiAliasing)
			{
				double a = ((double)WIDTH/HEIGHT);

				double xLow = -a * std::tan(radians(fov)/2);
				double yHigh = -std::tan(radians(fov)/2);

				double xLowerBound;
				double xUpperBound;
				double yLowerBound;
				double yUpperBound;

				//get x bounds
				if(x == 0)
				{
					double xCenter = xLow + 2 * ((double)x)/WIDTH * a * std::tan(radians(fov)/2);
					double xNextCenter = xLow + 2 * ((double)(x + 1))/WIDTH * a * std::tan(radians(fov)/2);
					double midDistance = (xNextCenter - xCenter)/2;
					xLowerBound = xCenter - midDistance;
					xUpperBound = xCenter + midDistance;
				}
				else if(x == WIDTH - 1)
				{
					double xCenter = xLow + 2 * ((double)x)/WIDTH * a * std::tan(radians(fov)/2);
					double xNextCenter = xLow + 2 * ((double)(x + 1))/WIDTH * a * std::tan(radians(fov)/2);
					double xPrevCenter = xLow + 2 * ((double)(x - 1))/WIDTH * a * std::tan(radians(fov)/2);
					double upperMidDistance = (xNextCenter - xCenter)/2;
					double lowerMidDistance = (xCenter - xPrevCenter)/2;
					xLowerBound = xCenter - lowerMidDistance;
					xUpperBound = xCenter + upperMidDistance;
				}
				else
				{
					double xCenter = xLow + 2 * ((double)x)/WIDTH * a * std::tan(radians(fov)/2);
					double xPrevCenter = xLow + 2 * ((double)(x - 1))/WIDTH * a * std::tan(radians(fov)/2);
					double midDistance = (xCenter - xPrevCenter)/2;
					xLowerBound = xCenter - midDistance;
					xUpperBound = xCenter + midDistance;
				}

				//get y bounds
				if(y == 0)
				{
					double yCenter = yHigh + 2 * ((double)y)/HEIGHT * std::tan(radians(fov)/2);
					double yNextCenter = yHigh + 2 * ((double)(y + 1))/HEIGHT * std::tan(radians(fov)/2);
					double midDistance = (yCenter - yNextCenter)/2;
					yUpperBound = yCenter - midDistance;
					yLowerBound = yCenter + midDistance;
				}
				else if(y == HEIGHT - 1)
				{
					double yCenter = yHigh + 2 * ((double)y)/HEIGHT * std::tan(radians(fov)/2);
					double yNextCenter = yHigh + 2 * ((double)(y + 1))/HEIGHT * std::tan(radians(fov)/2);
					double yPrevCenter = yHigh + 2 * ((double)(y - 1))/HEIGHT * std::tan(radians(fov)/2);
					double lowerMidDistance = (yCenter - yNextCenter)/2;
					double upperMidDistance = (yPrevCenter - yCenter)/2;
					yLowerBound = yCenter - lowerMidDistance;
					yUpperBound = yCenter + upperMidDistance;
				}
				else
				{
					double yCenter = yHigh + 2 * ((double)y)/HEIGHT * std::tan(radians(fov)/2);
					double yPrevCenter = yHigh + 2 * ((double)(y - 1))/HEIGHT * std::tan(radians(fov)/2);
					double midDistance = (yPrevCenter - yCenter)/2;
					yUpperBound = yCenter - midDistance;
					yLowerBound = yCenter + midDistance;
				}

				int xSubdivisions = 2;
				int ySubdivisions = 2;

				double xInc = (xUpperBound - xLowerBound)/xSubdivisions;
				double yInc = (yUpperBound - yLowerBound)/ySubdivisions;

				double color[3] = {0, 0, 0};

				for(int i = 0; i < xSubdivisions; i++)
				{
					for(int j = 0; j < ySubdivisions; j++)
					{
						intersectObject obj = none;
						int objIndex = 0;
						vector temp;

						temp.x = xLowerBound + i * xInc;
						temp.y = yUpperBound - j * yInc;
						temp.z = -1;
						temp.normalize();

						double start[3] = {0, 0, 0};
						double cIntensity[3] = {0, 0, 0};
						vector normTemp(0, 0, 0);
						double alpha = 0, beta = 0, gamma = 0;
						IntersectionRecursiveRedone(temp, start, cIntensity, 0, none, normTemp);
						//plot_pixel(x, y, floor(cIntensity[0] * 255), floor(cIntensity[1] * 255), floor(cIntensity[2] * 255));
						color[0] += cIntensity[0];
						color[1] += cIntensity[1];
						color[2] += cIntensity[2];
					}
				}

				color[0] /= (xSubdivisions * ySubdivisions);
				color[1] /= (xSubdivisions * ySubdivisions);
				color[2] /= (xSubdivisions * ySubdivisions);

				plot_pixel(x, y, floor(color[0] * 255), floor(color[1] * 255), floor(color[2] * 255));
			}
			else
			{
				vector temp;
				temp.x = xLow + 2 * ((double)x)/WIDTH * a * std::tan(radians(fov)/2);
				temp.y = yHigh + 2 * ((double)y)/HEIGHT * std::tan(radians(fov)/2);
				temp.z = -1;

				temp.normalize();

			

				double start[3] = {0, 0, 0};
				double cIntensity[3] = {0, 0, 0};
				vector normTemp(0, 0, 0);
				double alpha = 0, beta = 0, gamma = 0;
				IntersectionRecursiveRedone(temp, start, cIntensity, 0, none, normTemp);
				plot_pixel(x, y, floor(cIntensity[0] * 255), floor(cIntensity[1] * 255), floor(cIntensity[2] * 255));
			}
		}
		glEnd();
		glFlush();
	}
	printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
	glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
	glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
	buffer[HEIGHT-y-1][x][0]=r;
	buffer[HEIGHT-y-1][x][1]=g;
	buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
	plot_pixel_display(x,y,r,g,b);
	if(mode == MODE_JPEG)
		plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
	Pic *in = NULL;

	in = pic_alloc(640, 480, 3, NULL);
	printf("Saving JPEG file: %s\n", filename);

	memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
	if (jpeg_write(filename, in))
		printf("File saved Successfully\n");
	else
		printf("Error in Saving\n");

	pic_free(in);      

}

void parse_check(char *expected,char *found)
{
	if(stricmp(expected,found))
	{
		char error[100];
		printf("Expected '%s ' found '%s '\n",expected,found);
		printf("Parse error, abnormal abortion\n");
		exit(0);
	}

}

void parse_doubles(FILE*file, char *check, double p[3])
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check(check,str);
	fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
	printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
	char str[100];
	fscanf(file,"%s",str);
	parse_check("rad:",str);
	fscanf(file,"%lf",r);
	printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
	char s[100];
	fscanf(file,"%s",s);
	parse_check("shi:",s);
	fscanf(file,"%lf",shi);
	printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
	FILE *file = fopen(argv,"r");
	int number_of_objects;
	char type[50];
	int i;
	Triangle t;
	Sphere s;
	Light l;
	fscanf(file,"%i",&number_of_objects);

	printf("number of objects: %i\n",number_of_objects);
	char str[200];

	parse_doubles(file,"amb:",ambient_light);

	for(i=0;i < number_of_objects;i++)
	{
		fscanf(file,"%s\n",type);
		printf("%s\n",type);
		if(stricmp(type,"triangle")==0)
		{

			printf("found triangle\n");
			int j;

			for(j=0;j < 3;j++)
			{
				parse_doubles(file,"pos:",t.v[j].position);
				parse_doubles(file,"nor:",t.v[j].normal);
				parse_doubles(file,"dif:",t.v[j].color_diffuse);
				parse_doubles(file,"spe:",t.v[j].color_specular);
				parse_shi(file,&t.v[j].shininess);
			}

			if(num_triangles == MAX_TRIANGLES)
			{
				printf("too many triangles, you should increase MAX_TRIANGLES!\n");
				exit(0);
			}
			triangles[num_triangles++] = t;
		}
		else if(stricmp(type,"sphere")==0)
		{
			printf("found sphere\n");

			parse_doubles(file,"pos:",s.position);
			parse_rad(file,&s.radius);
			parse_doubles(file,"dif:",s.color_diffuse);
			parse_doubles(file,"spe:",s.color_specular);
			parse_shi(file,&s.shininess);

			if(num_spheres == MAX_SPHERES)
			{
				printf("too many spheres, you should increase MAX_SPHERES!\n");
				exit(0);
			}
			spheres[num_spheres++] = s;
		}
		else if(stricmp(type,"light")==0)
		{
			printf("found light\n");
			parse_doubles(file,"pos:",l.position);
			parse_doubles(file,"col:",l.color);

			if(num_lights == MAX_LIGHTS)
			{
				printf("too many lights, you should increase MAX_LIGHTS!\n");
				exit(0);
			}
			lights[num_lights++] = l;
		}
		else
		{
			printf("unknown type in scene description:\n%s\n",type);
			exit(0);
		}
	}
	return 0;
}

void display()
{

}

void init()
{
	glMatrixMode(GL_PROJECTION);
	glOrtho(0,WIDTH,0,HEIGHT,1,-1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
	//hack to make it only draw once
	static int once=0;
	if(!once)
	{
		draw_scene();
		if(mode == MODE_JPEG)
			save_jpg();
	}
	once=1;
}

int main (int argc, char ** argv)
{
	if (argc<2 || argc > 3)
	{  
		printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
		exit(0);
	}
	if(argc == 3)
	{
		mode = MODE_JPEG;
		filename = argv[2];
	}
	else if(argc == 2)
		mode = MODE_DISPLAY;

	/*
	int maxBounces = 5;
	bool antiAliasing = false;
	double shadowRayCount = 100;
	*/
	while(true)
	{
		int i;
		std::cout<<"Enter the number of bounces: ";
		std::cin>>i;
		if(!std::cin)
		{
			std::cin.clear();
			std::cin.ignore(INT_MAX, '\n');
			std::cout<<"Please use the correct data type\n";
			continue;
		}
		if(i < 0)
		{
			maxBounces = 0;
		}
		else if(i > 5)
		{
			maxBounces = 5;
		}
		else
		{
			maxBounces = i;
		}
		std::cin.ignore(INT_MAX, '\n');
		break;
	}

	while(true)
	{
		char i;
		std::cout<<"Do you want anti-aliasing (Y\\N): ";
		std::cin>>i;
		if(!std::cin)
		{
			std::cin.clear();
			std::cin.ignore(INT_MAX, '\n');
			std::cout<<"Please use the correct data type\n";
			continue;
		}
		if(i == 'n' || i == 'N')
		{
			antiAliasing = false;
		}
		else if(i == 'y' || i == 'Y')
		{
			antiAliasing = true;
		}
		else
		{
			std::cout<<"Since you didn't say yes or no, antialiasing is off\n";
		}
		std::cin.ignore(INT_MAX, '\n');
		break;
	}

	while(true)
	{
		int i;
		std::cout<<"Enter the number of shadow rays: ";
		std::cin>>i;
		if(!std::cin)
		{
			std::cin.clear();
			std::cin.ignore(INT_MAX, '\n');
			std::cout<<"Please use the correct data type\n";
			continue;
		}
		if(i < 1)
		{
			shadowRayCount = 1;
		}
		else if(i > 100)
		{
			shadowRayCount = 100;
		}
		else
		{
			shadowRayCount = i;
		}
		std::cin.ignore(INT_MAX, '\n');
		break;
	}

	std::srand(time(nullptr));

	glutInit(&argc,argv);
	loadScene(argv[1]);

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(0,0);
	glutInitWindowSize(WIDTH,HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	init();
	std::cout<<radians(fov)<<std::endl;
	glutMainLoop();
}
