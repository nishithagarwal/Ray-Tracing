/*
 CSCI 480
 Assignment 3 Raytracer
 
 Name: <Your name here>
 */

#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "pic.h"
#include <string.h>
#include <math.h>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

//pixel resolution settings
#define resolution_x 640
#define resolution_y 480

unsigned char buffer[HEIGHT][WIDTH][3];

float aspect = (float) WIDTH/HEIGHT;
double orig[3]= {0.0,0.0,0.0};

double r[3];

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

Vertex p1,p2,p3,p4;
Vertex vertices[resolution_y][resolution_x];

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
} Light;

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



void normalize(double *p)
{
    double length_n = fabs (sqrt((p[0]*p[0])+(p[1]*p[1])+(p[2]*p[2])));
    p[0]=(double) p[0]/length_n;
    p[1]=(double) p[1]/length_n;
    p[2]=(double) p[2]/length_n;
    
    //printf("\nv1 %f,%f,%f,%f ",p[0],p[1],p[2],length_n);
    
}

float dotProduct(double a[3],double b[3])
{
    float c = (a[0]*b[0])+(a[1]*b[1])+(a[2]*b[2]);
    return c;
}

void crossProduct(double a[3],double b[3], double *c)
{
    c[0] = a[1]*b[2]-b[1]*a[2];
    c[1] = -a[0]*b[2]+b[0]*a[2];
    c[2] = a[0]*b[1]-b[0]*a[1];
    
}


void setNormalSphere(double *normal,double p[3], float t,int k)
{
    normal[0] = (p[0]*t) - spheres[k].position[0];
    normal[1] = (p[1]*t) - spheres[k].position[1];
    normal[2] = (p[2]*t) - spheres[k].position[2];
    
    //printf("\nNormals1 %f,%f,%f ", normal[0], normal[1], normal[2]);
    normalize(normal);
    
    //printf("\nPoint %f,%f,%f ",d[0]*t,d[1]*t,d[2]*t);
    //printf("\nSphere %f,%f,%f ",spheres[0].position[0],spheres[0].position[1],spheres[0].position[2]);
    //printf("\nNormals %f,%f,%f ", normal[0], normal[1], normal[2]);
    
}

void setNormalTriangle(double *normal,double p[3], float t,int k)
{
    
    /*double e1[3]={
     triangles[k].v[1].position[0]-triangles[k].v[0].position[0],
     triangles[k].v[1].position[1]-triangles[k].v[0].position[1],
     triangles[k].v[1].position[2]-triangles[k].v[0].position[2]};
     
     double e2[3]={
     triangles[k].v[2].position[0]-triangles[k].v[0].position[0],
     triangles[k].v[2].position[1]-triangles[k].v[0].position[1],
     triangles[k].v[2].position[2]-triangles[k].v[0].position[2]};
     
     crossProduct(e2,e1,normal);
     
     normal[0] = fabs(normal[0]);
     normal[1] = fabs(normal[1]);
     normal[2] = fabs(normal[2]);
     
     normalize(normal);*/
    
    
    
    
}

float max(float a,float b)
{
    if(a>=b)
        return a;
    else
        return b;
}

float min(float a,float b)
{
    if(a>=b)
        return b;
    else
        return a;
}


void calcPixels()
{
    //Position of 4 points on the image plane
    p1.position[0]=(-aspect)*(tan(0.0174532925*fov/2));
    p1.position[1]=tan(0.0174532925*fov/2);
    p1.position[2]=-1.0;
    
    p2.position[0]=(aspect)*(tan(0.0174532925*fov/2));
    p2.position[1]=tan(0.0174532925*fov/2);
    p2.position[2]=-1.0;
    
    p3.position[0]=(aspect)*(tan(0.0174532925*fov/2));
    p3.position[1]=-tan(0.0174532925*fov/2);
    p3.position[2]=-1.0;
    
    p4.position[0]=(-aspect)*(tan(0.0174532925*fov/2));
    p4.position[1]=-tan(0.0174532925*fov/2);
    p4.position[2]=-1.0;
    
    printf("\n %f %f %f",p1.position[0],p1.position[1],p1.position[2]);
    printf("\n %f %f %f",p2.position[0],p2.position[1],p2.position[2]);
    printf("\n %f %f %f",p3.position[0],p3.position[1],p3.position[2]);
    printf("\n %f %f %f",p4.position[0],p4.position[1],p4.position[2]);
    
    for(int i=0;i<resolution_y;i++)
    {
        vertices[i][0].position[0]=p1.position[0];
        for(int j=0;j<resolution_x;j++)
        {
            vertices[i][j].position[2]=-1.0;
            if(j==0)
                vertices[i][j].position[0]=p4.position[0];
            else
                vertices[i][j].position[0]=vertices[i][j-1].position[0]+((p2.position[0]-p1.position[0])/(resolution_x-1));
            if(i==0)
                vertices[i][j].position[1]=p4.position[1];
            else
                vertices[i][j].position[1]=vertices[i-1][j].position[1]+((p1.position[1]-p4.position[1])/(resolution_y-1));
            
        }
    }
}



float intersectSphere(double orig[3], double d[3], Sphere s)
{
    double o[3] = {orig[0]-s.position[0],orig[1]-s.position[1],orig[2]-s.position[2]};
    
    float a = dotProduct(d, d);
    float b = 2 * dotProduct(d, o);
    float c = dotProduct(o, o) - (s.radius * s.radius);
    
    float disc = b * b - 4 * a * c;
    float distSqrt = sqrtf(disc);
    
    if (disc < 0)
        return 0;
    
    //printf("\na,b,c, t: %f,%f,%f,%f",a,b,c,disc);
    float q;
    if (b < 0)
        q = (-b - distSqrt)/2.0;
    else
        q = (-b + distSqrt)/2.0;
    
    float t0 = q / a;
    float t1 = c / q;
    
    if (t0 > t1)
    {
        float temp = t0;
        t0 = t1;
        t1 = temp;
    }
    
    if (t1 < 0)
        return 0;
    
    if (t0 < 0)
    {
        return t1;
    }
    else
    {
        return t0;
    }
    
    
}



float intersectTriangle(double orig[3],double d[3], Triangle *t)
{
    double e1[3]={t->v[1].position[0]-t->v[0].position[0],t->v[1].position[1]-t->v[0].position[1],t->v[1].position[2]-t->v[0].position[2]};
    double e2[3]={t->v[2].position[0]-t->v[0].position[0],t->v[2].position[1]-t->v[0].position[1],t->v[2].position[2]-t->v[0].position[2]};
    
    double e1e2[3];
    crossProduct(e1,e2,e1e2);
    
    
    double p[3];
    crossProduct(d,e2,p);
    
    normalize(e1e2);
    
    float a = dotProduct(e1, p);
    if (a > -0.00001 && a < 0.00001)
        return 0;
    
    float f =  1 / a;
    double s[3] ={orig[0]-t->v[0].position[0],orig[1]-t->v[0].position[1],orig[2]-t->v[0].position[2]};
    float u = f*(dotProduct(s, p));
    if(u < 0.0 || u > 1.0)
        return 0;
    
    double q[3];
    crossProduct(s, e1, q);
    float v = f * (dotProduct(d, q));
    if(v < 0.0 || u + v > 1.0)
        return 0;
    
    // at this stage we can compute t to find out where
	// the intersection point is on the line
	float t0 = f * dotProduct(e2,q);
    
	if (t0 > 0.00001) // ray intersection
		return t0;
    
	else // this means that there is a line intersection
        // but not a ray intersection
        return 0;
    
    
}

//SET COLOR FOR EACH SPHERE
void setColorSphere(Vertex *v,float t,int k)
{
    //SET DIFFUSE COLOR
    double l[3]={(-v->position[0]*t)+lights[0].position[0],(-v->position[1]*t)+lights[0].position[1],(-v->position[2]*t)+lights[0].position[2]};
    
    normalize(l);
    
    float ln = max(dotProduct(l, v->normal),0.0);
    
    v->color_diffuse[0] = spheres[k].color_diffuse[0]*ln;
    v->color_diffuse[1] = spheres[k].color_diffuse[1]*ln;
    v->color_diffuse[2] = spheres[k].color_diffuse[2]*ln;
    
    //SET SPECULAR COLOR
    double e[3] = {-v->position[0]*t,-v->position[1]*t,-v->position[2]*t};
    //double r[3] = {(2*ln*v->normal[0])-l[0],(2*ln*v->normal[1])-l[1],(2*ln*v->normal[2])-l[2]}; //reflection vector
    
    r[0] = (2*ln*v->normal[0])-l[0];
    r[1] = (2*ln*v->normal[1])-l[1];
    r[2] = (2*ln*v->normal[2])-l[2]; //reflection vector
    
    normalize(r);
    normalize(e);
    
    
    float re = max(dotProduct(r, e),0.0);
    
    v->color_specular[0] = spheres[k].color_specular[0]* pow(re,spheres[k].shininess);
    v->color_specular[1] = spheres[k].color_specular[1]* pow(re,spheres[k].shininess);
    v->color_specular[2] = spheres[k].color_specular[2]* pow(re,spheres[k].shininess);
    
    //ADD SHADOWS  - Set black if view direction intersects another object
    double p[3] = {v->position[0]*t,v->position[1]*t,v->position[2]*t};
    bool flag=false;
    for(int j=0;j<num_spheres;j++)
    {
        if(j!=k && intersectSphere(p, l, spheres[j]))
        {
            flag=true;
        }
    }
    
    
    
    for(int j=0;j<num_triangles;j++)
    {
        if(intersectTriangle(p, l, &triangles[k]) )
        {
            flag=true;
        }
    }
    
    
    if(flag==true){
        v->color_diffuse[0] = 0;
        v->color_diffuse[1] = 0;
        v->color_diffuse[2] = 0;
        
        v->color_specular[0] = 0;
        v->color_specular[1] = 0;
        v->color_specular[2] = 0;
        
        r[0] = 0;
        r[1] = 0;
        r[2] = 1;
    }
}


//SET COLOR FOR EACH TRAINGLE
void  setColorTraingle(Vertex *v,float t,int k,int s)
{
    
    float A = fabs(
                   triangles[k].v[0].position[0]*(triangles[k].v[1].position[1] - triangles[k].v[2].position[1])
                   +triangles[k].v[1].position[0]*(triangles[k].v[2].position[1] - triangles[k].v[0].position[1])
                   +triangles[k].v[2].position[0]*(triangles[k].v[0].position[1] - triangles[k].v[1].position[1])
                   )/2;
    
    float A1 = fabs(
                    v->position[0]*t*(triangles[k].v[1].position[1] - triangles[k].v[2].position[1])
                    +triangles[k].v[1].position[0]*(triangles[k].v[2].position[1] - v->position[1]*t)
                    +triangles[k].v[2].position[0]*(v->position[1]*t - triangles[k].v[1].position[1])
                    )/2;
    
    float A2 = fabs(
                    triangles[k].v[0].position[0]*(v->position[1]*t - triangles[k].v[2].position[1])
                    +v->position[0]*t*(triangles[k].v[2].position[1] - triangles[k].v[0].position[1])
                    +triangles[k].v[2].position[0]*(triangles[k].v[0].position[1] - v->position[1]*t)
                    )/2;
    
    float A3 = fabs(
                    triangles[k].v[0].position[0]*(triangles[k].v[1].position[1] - v->position[1]*t)
                    +triangles[k].v[1].position[0]*(v->position[1]*t - triangles[k].v[0].position[1])
                    +v->position[0]*t*(triangles[k].v[0].position[1] - triangles[k].v[1].position[1])
                    )/2;
    
    A=A1+A2+A3;
    
    v->normal[0] = (A1/A)*triangles[k].v[0].normal[0]+(A2/A)*triangles[k].v[1].normal[0]+(A3/A)*triangles[k].v[2].normal[0];
    v->normal[1] = (A1/A)*triangles[k].v[0].normal[1]+(A2/A)*triangles[k].v[1].normal[1]+(A3/A)*triangles[k].v[2].normal[1];
    v->normal[2] = (A1/A)*triangles[k].v[0].normal[2]+(A2/A)*triangles[k].v[1].normal[2]+(A3/A)*triangles[k].v[2].normal[2];
    
    normalize(v->normal);
    
    
    double l[3]={-v->position[0]*t+lights[s].position[0],-v->position[1]*t+lights[s].position[1],-v->position[2]*t+lights[s].position[2]}; //light vector
    normalize(l);
    
    float ln = max(dotProduct(l, v->normal),0.0);
    //ln =1.0;
    
    double p0_diffuse[3] = {
        triangles[k].v[0].color_diffuse[0]*ln,
        triangles[k].v[0].color_diffuse[1]*ln,
        triangles[k].v[0].color_diffuse[2]*ln}; //red
    
    double p1_diffuse[3] = {
        triangles[k].v[1].color_diffuse[0]*ln,
        triangles[k].v[1].color_diffuse[1]*ln,
        triangles[k].v[1].color_diffuse[2]*ln}; //green
    
    double p2_diffuse[3] = {
        triangles[k].v[2].color_diffuse[0]*ln,
        triangles[k].v[2].color_diffuse[1]*ln,
        triangles[k].v[2].color_diffuse[2]*ln}; //blue
    
    
    //printf("\nA A1 A2 A3 Sum %f %f %f %f %f",A,A1,A2,A3,A1+A2+A3);
    
    
    v->color_diffuse[0] = (A1/A)*p0_diffuse[0]+(A2/A)*p1_diffuse[0]+(A3/A)*p2_diffuse[0];
    v->color_diffuse[1] = (A1/A)*p0_diffuse[1]+(A2/A)*p1_diffuse[1]+(A3/A)*p2_diffuse[1];
    v->color_diffuse[2] = (A1/A)*p0_diffuse[2]+(A2/A)*p1_diffuse[2]+(A3/A)*p2_diffuse[2];
    
    
    //SET SPECULAR COLOR
    double e[3] = {-v->position[0]*t,-v->position[1]*t,-v->position[2]*t};
    
    r[0] = (2*ln*v->normal[0])-l[0];
    r[1] = (2*ln*v->normal[1])-l[1];
    r[2] = (2*ln*v->normal[2])-l[2]; //reflection vector
    
    
    normalize(r);
    normalize(e);
    
    
    float re = max(dotProduct(r, e),0.0);
    
    
    double p0_specular[3] = {
        triangles[k].v[0].color_specular[0]* pow(re,triangles[k].v[0].shininess),
        triangles[k].v[0].color_specular[1]* pow(re,triangles[k].v[0].shininess),
        triangles[k].v[0].color_specular[2]* pow(re,triangles[k].v[0].shininess)}; //red
    
    double p1_specular[3] = {
        triangles[k].v[1].color_specular[0]* pow(re,triangles[k].v[1].shininess),
        triangles[k].v[1].color_specular[1]* pow(re,triangles[k].v[1].shininess),
        triangles[k].v[1].color_specular[2]* pow(re,triangles[k].v[1].shininess)}; //green
    
    double p2_specular[3] = {
        triangles[k].v[2].color_specular[0]* pow(re,triangles[k].v[2].shininess),
        triangles[k].v[2].color_specular[1]* pow(re,triangles[k].v[2].shininess),
        triangles[k].v[2].color_specular[2]* pow(re,triangles[k].v[2].shininess)}; //blue
    
     v->color_specular[0] = (A1/A)*p0_specular[0]+(A2/A)*p1_specular[0]+(A3/A)*p2_specular[0];
     v->color_specular[1] = (A1/A)*p0_specular[1]+(A2/A)*p1_specular[1]+(A3/A)*p2_specular[1];
     v->color_specular[2] = (A1/A)*p0_specular[2]+(A2/A)*p1_specular[2]+(A3/A)*p2_specular[2];
    
    
    //ADD SHADOWS  - Set black if view direction intersects another object
    double p[3] = {v->position[0]*t,v->position[1]*t,v->position[2]*t};
    bool flag=false;
    for(int j=0;j<num_spheres;j++)
    {
        if(intersectSphere(p, l, spheres[j]))
        {
            flag=true;
        }
    }
    for(int j=0;j<num_triangles-1;j++)
    {
        if(j!=k && intersectTriangle(p, l, &triangles[j]))
        {
            flag=true;
        }
    }
    if(flag==true){
        v->color_diffuse[0] = 0;
        v->color_diffuse[1] = 0;
        v->color_diffuse[2] = 0;
        
        v->color_specular[0] = 0;
        v->color_specular[1] = 0;
        v->color_specular[2] = 0;
        
        r[0] = 0;
        r[1] = 0;
        r[2] = 1;
    }
    
    
    
}


//MODIFY THIS FUNCTION
void draw_scene()
{
    float min_t=-1;
    int min_k=-1;
    
    float min_t1=-1;
    int min_k1=-1;
    
    calcPixels();
    
    for(int i=0;i<resolution_y;i++)
    {
        glPointSize(2.0);
        glBegin(GL_POINTS);
        for(int j=0;j<resolution_x;j++)
        {
            bool flag=false;
            float min_dist=100000;
            
            bool flag1=false;
            float min_dist1=100000;
            
            //normalize the direction vector
            normalize(vertices[i][j].position);
            
            //Check if Ray intersects with any Traingle
            for(int k=0;k<num_triangles;k++)
            {
                float t = intersectTriangle(orig,vertices[i][j].position, &triangles[k]);
                if(t!=0)
                {
                    float dist  = sqrtf( pow(vertices[i][j].position[0]*t-vertices[i][j].position[0],2)
                                        +pow(vertices[i][j].position[1]*t-vertices[i][j].position[1],2)
                                        +pow(vertices[i][j].position[2]*t-vertices[i][j].position[2],2)
                                        );
                    flag=true;
                    
                    
                    if (dist<min_dist){
                        min_t=t;
                        min_k=k;
                        min_dist=dist;
                    }
                }
            }
            
            
            if(flag==true){
                //setNormalTriangle(vertices[i][j].normal,vertices[i][j].position,min_t,min_k);
                
                double finalColor[3]={0.0,0.0,0.0};
                
                for(int s=0;s<num_lights;s++){
                    setColorTraingle(&vertices[i][j],min_t,min_k,s);
                    
                    finalColor[0]+=lights[s].color[0]*(vertices[i][j].color_diffuse[0]+vertices[i][j].color_specular[0]);
                    finalColor[1]+=lights[s].color[1]*(vertices[i][j].color_diffuse[1]+vertices[i][j].color_specular[1]);
                    finalColor[2]+=lights[s].color[2]*(vertices[i][j].color_diffuse[2]+vertices[i][j].color_specular[2]);
                }
                
                finalColor[0]+=ambient_light[0];
                finalColor[1]+=ambient_light[1];
                finalColor[2]+=ambient_light[2];
                

                
                
                /*for(int k=0;k<num_triangles;k++)
                {
                    float t1 = intersectTriangle(orig,r, &triangles[k]);
                    if(t1!=0)
                    {
  
                        flag1=true;
                        
                        
                        //if (t1<min_t1){
                            min_t1=t1;
                            min_k1=k;
                        //}
                    }
                }
                
                if(flag1==true)
                {
                    printf("TRUE1");
                    flag1=false;
                    
                    finalColor[0]*= 1-triangles[min_k].v[0].color_specular[0];
                    finalColor[1]*= 1-triangles[min_k].v[0].color_specular[1];
                    finalColor[2]*= 1-triangles[min_k].v[0].color_specular[2];

                    
                    
                    finalColor[0]+= triangles[min_k].v[0].color_specular[0]*triangles[min_k1].v[0].color_diffuse[0];
                    finalColor[1]+= triangles[min_k].v[0].color_specular[1]*triangles[min_k1].v[0].color_diffuse[1];
                    finalColor[2]+= triangles[min_k].v[0].color_specular[2]*triangles[min_k1].v[0].color_diffuse[2];
                    
                }*/
                
                
                
                finalColor[0] = max(finalColor[0],0.0);
                finalColor[0] = min(finalColor[0],1.0);
                finalColor[1] = max(finalColor[1],0.0);
                finalColor[1] = min(finalColor[1],1.0);
                finalColor[2] = max(finalColor[2],0.0);
                finalColor[2] = min(finalColor[2],1.0);
                
                
                plot_pixel(j,i
                           ,finalColor[0]*255
                           ,finalColor[1]*255
                           ,finalColor[2]*255);
                
                //printf("\nFinal Color %f,%f,%f ",finalColor[0],finalColor[1],finalColor[2]);
                
            }
            
            
            
            //Check if Ray intersects with Sphere
            for(int k=0;k<num_spheres;k++)
            {
                float t=intersectSphere(orig,vertices[i][j].position,spheres[k]);
                if(t!=0)
                {
                    setNormalSphere(vertices[i][j].normal,vertices[i][j].position,t,k);
                    setColorSphere(&vertices[i][j],t,k);
                    
                    double finalColor[3]={
                        lights[0].color[0]*ambient_light[0]+vertices[i][j].color_diffuse[0]+vertices[i][j].color_specular[0],
                        lights[0].color[1]*ambient_light[1]+vertices[i][j].color_diffuse[1]+vertices[i][j].color_specular[1],
                        lights[0].color[2]*ambient_light[2]+vertices[i][j].color_diffuse[2]+vertices[i][j].color_specular[2]};
                    
                    //printf("\nDiffuse %f,%f,%f ", vertices[i][j].color_diffuse[0], vertices[i][j].color_diffuse[1], vertices[i][j].color_diffuse[2]);
                    //printf("\nSpecular %f,%f,%f ", vertices[i][j].color_specular[0], vertices[i][j].color_specular[1], vertices[i][j].color_specular[2]);
                    // printf("\nFinal Color %f,%f,%f ",finalColor[0],finalColor[1],finalColor[2]);
                    
      
                    
                    /*for(int k=0;k<num_triangles;k++)
                    {
                        float t1 = intersectTriangle(orig,r, &triangles[k]);
                        if(t1!=0)
                        {
  
                            flag1=true;
                            
                            
                            //if (t1<min_t1){
                                min_t1=t1;
                                min_k1=k;
                            //}
                        }
                    }
                    
                    if(flag1==true)
                    {
                        printf("TRUE2");
                        flag1=false;
                        
                        finalColor[0]*= 1-spheres[min_k].color_specular[0];
                        finalColor[1]*= 1-spheres[min_k].color_specular[1];
                        finalColor[2]*= 1-spheres[min_k].color_specular[2];
                        
                        
                        finalColor[0]+= spheres[min_k].color_specular[0]*triangles[min_k1].v[0].color_diffuse[0];
                        finalColor[1]+= spheres[min_k].color_specular[1]*triangles[min_k1].v[0].color_diffuse[1];
                        finalColor[2]+= spheres[min_k].color_specular[2]*triangles[min_k1].v[0].color_diffuse[2];
                        
                    }*/
                    
                  
                    
                    /*for(int k=0;k<num_spheres;k++)
                    {
                        float t1 = intersectSphere(orig,r, spheres[k]);
                        if(t1!=0)
                        {

                            flag1=true;
                            
                            
                            //if (t1<min_t1){
                            min_t1=t1;
                            min_k1=k;
                            //}
                        }
                    }
                    
                    if(flag1==true)
                    {
                        printf("TRUE3");
                        flag1=false;
                        finalColor[0]+= spheres[min_k].color_specular[0]*spheres[min_k1].color_diffuse[0];
                        finalColor[1]+= spheres[min_k].color_specular[1]*spheres[min_k1].color_diffuse[1];
                        finalColor[2]+= spheres[min_k].color_specular[2]*spheres[min_k1].color_diffuse[2];
                        
                    }*/
                    
                    finalColor[0] = max(finalColor[0],0.0);
                    finalColor[0] = min(finalColor[0],1.0);
                    finalColor[1] = max(finalColor[1],0.0);
                    finalColor[1] = min(finalColor[1],1.0);
                    finalColor[2] = max(finalColor[2],0.0);
                    finalColor[2] = min(finalColor[2],1.0);
                    
                    plot_pixel(j,i
                               ,finalColor[0]*255
                               ,finalColor[1]*255
                               ,finalColor[2]*255);
                }
                
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
    if(strcasecmp(expected,found))
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
        if(strcasecmp(type,"triangle")==0)
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
        else if(strcasecmp(type,"sphere")==0)
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
        else if(strcasecmp(type,"light")==0)
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

void handleKeypress(unsigned char key, int x, int y)
{
    switch(key) {
        case 27: exit(0);
    }
    
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
    
    glClearColor(1.0,1.0,1.0,0);
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
    
    glutInit(&argc,argv);
    loadScene(argv[1]);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    int window = glutCreateWindow("Ray Tracer");
    glutDisplayFunc(display);
    glutKeyboardFunc(handleKeypress);
    glutIdleFunc(idle);
    init();
    glutMainLoop();
}
