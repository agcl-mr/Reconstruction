#include <CGAL/enum.h>
#include <fstream>
#include <queue>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <stack>
#include <GL/glut.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangle_2.h>



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<int[10],K> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<int[3],K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Delaunay;

typedef Delaunay::Point   Point;
typedef Delaunay::Triangle   Triangle;
typedef std::vector<Point> point_vec;
point_vec ch,out;

#define PI 3.1415926535


const int N = 100000;
int num_edges=0;
std::vector<int> graph[N];
std::vector<int> cycles[N];
int no_of_cycles=0,no_ver_graph=0;
std::vector<Point> cycle_vertices(N);

//int y=0;
float rad=0;
Delaunay dt;
typedef std::vector<Point> vertex_arrray;
vertex_arrray vertices,boundaries;
class Edge
{
public:
  Point source;
  Point target;
 };

typedef std::vector<Edge> Edges;
Edges shape,bdry,marked_edges,skinny;



//variables for OpenGL
int minx=9999,miny=9999,maxx=0,maxy=0;
int window;
float diagonalDistance;
GLdouble width, height;
float minX = std::numeric_limits<int>::max(), maxX = std::numeric_limits<int>::min();
float minY = std::numeric_limits<int>::max(), maxY = std::numeric_limits<int>::min();



bool isFinite(Delaunay::All_faces_iterator afi) /*Checking whether the face is finite or not*/
{
  if(afi->vertex(0) != dt.infinite_vertex() && afi->vertex(1) != dt.infinite_vertex() && afi->vertex(2) != dt.infinite_vertex())
      return true;
  return false;
}


void computconvexhull()
{
   std::vector<Point>::iterator ptr;
   for (ptr = vertices.begin(); ptr != vertices.end(); ptr++)
   {
        ch.push_back(*ptr);
   }

   CGAL::convex_hull_points_2(ch.begin(),ch.end(),back_inserter(out));
   out.push_back(out[0]);

}

void init()
{

         for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); fit++)
         {
             for(int k=0;k<=2;k++)
             {
                 fit->info()[k]=0;
             }
         }
}

float distance(Point a, Point b) /*compute the distance between two points*/
{
  return (float)(sqrt(((a.x()-b.x())  *(a.x()-b.x()))+ ((a.y()-b.y())* (a.y()-b.y()))));
}


//Checking skinny triangles
void check_skinnytraingle(Delaunay::Face::Face_handle fh1,Point circum_center)
{

     Point a = fh1->vertex(0)->point();
     Point b = fh1->vertex(1)->point();
     Point c = fh1->vertex(2)->point();

     float s1 = distance(a,b);
     float s2 = distance(b,c);
     float s3 = distance(c,a);
//Calculating the coordinates of incentre
     float x_coordnate = ((s2 * a.x()) + (s3 * b.x()) + (s1 * c.x())) / (s1+s2+s3);
     float y_coordnate = ((s2 * a.y()) + (s3 * b.y()) + (s1 * c.y())) / (s1+s2+s3);
     float m = std::max({s1, s2, s3});
    
//Calculating the distance between Circumcenter and incenter
     float center_diff = distance(Point(x_coordnate,y_coordnate),circum_center);
     Edge e;

        if(s1 <= s2 && s1 <= s3)
        {
            if(s1 < center_diff)
            {
                fh1->info()[0]=1;
                fh1->info()[1]=1;
            e.source = fh1->vertex(1)->point(); // skinny
            e.target = fh1->vertex(2)->point();
            skinny.push_back(e);
	          e.source = fh1->vertex(2)->point(); // skinny
            e.target = fh1->vertex(0)->point();
            skinny.push_back(e);

                if(fh1->neighbor(0)->neighbor(0) == fh1)
                    fh1->neighbor(0)->info()[0] = 1;
                else if(fh1->neighbor(0)->neighbor(1) == fh1)
                    fh1->neighbor(0)->info()[1] = 1;
                else if(fh1->neighbor(0)->neighbor(2) == fh1)
                    fh1->neighbor(0)->info()[2] = 1;

                if(fh1->neighbor(1)->neighbor(0) == fh1)
                    fh1->neighbor(1)->info()[0] = 1;
                else if(fh1->neighbor(1)->neighbor(1) == fh1)
                    fh1->neighbor(1)->info()[1] = 1;
                else if(fh1->neighbor(1)->neighbor(2) == fh1)
                    fh1->neighbor(1)->info()[2] = 1;


            }
        }

       else if(s2 <= s1 && s2 <= s3)
        {
            if(s2 < center_diff)
            {
                fh1->info()[1]=1;
                fh1->info()[2]=1;
            e.source = fh1->vertex(2)->point(); // skinny
            e.target = fh1->vertex(0)->point();
            skinny.push_back(e);
	          e.source = fh1->vertex(0)->point(); // skinny
            e.target = fh1->vertex(1)->point();
            skinny.push_back(e);
                if(fh1->neighbor(1)->neighbor(0) == fh1)
                    fh1->neighbor(1)->info()[0] = 1;
                else if(fh1->neighbor(1)->neighbor(1) == fh1)
                    fh1->neighbor(1)->info()[1] = 1;
                else if(fh1->neighbor(1)->neighbor(2) == fh1)
                    fh1->neighbor(1)->info()[2] = 1;

                if(fh1->neighbor(2)->neighbor(0) == fh1)
                    fh1->neighbor(2)->info()[0] = 1;
                else if(fh1->neighbor(2)->neighbor(1) == fh1)
                    fh1->neighbor(2)->info()[1] = 1;
                else if(fh1->neighbor(2)->neighbor(2) == fh1)
                    fh1->neighbor(2)->info()[2] = 1;

             }
        }

       else if(s3 <= s1 && s3 <= s2)
        {
            if(s3 < center_diff)
            {
                fh1->info()[2]=1;
                fh1->info()[0]=1;
            e.source = fh1->vertex(0)->point(); // skinny
            e.target = fh1->vertex(1)->point();
            skinny.push_back(e);
	          e.source = fh1->vertex(1)->point(); // skinny
            e.target = fh1->vertex(2)->point();
            skinny.push_back(e);
                if(fh1->neighbor(2)->neighbor(0) == fh1)
                    fh1->neighbor(2)->info()[0] = 1;
                else if(fh1->neighbor(2)->neighbor(1) == fh1)
                    fh1->neighbor(2)->info()[1] = 1;
                else if(fh1->neighbor(2)->neighbor(2) == fh1)
                    fh1->neighbor(2)->info()[2] = 1;

                if(fh1->neighbor(0)->neighbor(0) == fh1)
                    fh1->neighbor(0)->info()[0] = 1;
                else if(fh1->neighbor(0)->neighbor(1) == fh1)
                    fh1->neighbor(0)->info()[1] = 1;
                else if(fh1->neighbor(0)->neighbor(2) == fh1)
                    fh1->neighbor(0)->info()[2] = 1;

            }
        }

}


void populatemarkedEdges(void)
{
    Edge e;
    Delaunay::Finite_faces_iterator it;

    for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
     {

        if(it->info()[0] == 1 )
        {
            e.source = it->vertex(1)->point();
            e.target = it->vertex(2)->point();
            marked_edges.push_back(e);

         }

        if(it->info()[1] == 1 )
        {
            e.source = it->vertex(2)->point();
            e.target = it->vertex(0)->point();
            marked_edges.push_back(e);

         }

        if(it->info()[2] == 1 )
        {
            e.source = it->vertex(0)->point();
            e.target = it->vertex(1)->point();
            marked_edges.push_back(e);

         }
     }
}

void insertToShape(Point a, Point b, std::vector<Edge>  & shape)/*Inserting a new  into ec-shape*/
{
  Edge e;
  e.source = a;
  e.target = b;
  shape.push_back(e);

}



void post_process(void)
{

    Delaunay::Finite_faces_iterator it;
    int edgesum = 0,e0 = 0,e1 = 0,e2 = 0;
    for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
     {
        edgesum = it->info()[0] +  it->info()[1] +  it->info()[2];
        if(edgesum > 0)
        {
            e0 = it->neighbor(0)->info()[0] + it->neighbor(0)->info()[1] + it->neighbor(0)->info()[2];
            e1 = it->neighbor(1)->info()[0] + it->neighbor(1)->info()[1] + it->neighbor(1)->info()[2];
            e2 = it->neighbor(2)->info()[0] + it->neighbor(2)->info()[1] + it->neighbor(2)->info()[2];

            if(isFinite(it->neighbor(0)) && isFinite(it->neighbor(1)) && isFinite(it->neighbor(2)))
            {
                if(edgesum == 1)
                {
                    if(it->info()[0] == 1)
                    {
                        if(e0 == 1 )
                        {
                            it->info()[0] = 0;
                            it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;
                        }

                    }
                    else if(it->info()[1] == 1)
                    {
                        if( e1 == 1 )
                        {
                            it->info()[1] = 0;
                            it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;
                        }
                    }
                    else if(it->info()[2] == 1)
                    {
                        if(e2 == 1 )
                        {
                            it->info()[2] = 0;
                            it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;
                        }
                    }
                }
                else if(edgesum == 2)
                {

                   if(it->info()[0] == 0)
                    {

                        if(e1 == 1 && e2 == 1 )
                        {
                            it->info()[1] = 0; it->info()[2] = 0;
                            it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;
                            it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;
                        }
                    }
                  else if(it->info()[1] == 0)
                    {

                       if(e0 == 1 && e2 == 1 )
                        {
                            it->info()[0] = 0; it->info()[2] = 0;
                            it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;
                            it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;
                        }

                    }
                   else if(it->info()[2] == 0)
                   {

                       if(e0 == 1 && e1 == 1 )
                       {
                           it->info()[0] = 0; it->info()[1] = 0;
                           it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;
                           it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;
                       }
                   }


                }
                else if(edgesum == 3)
                {

                    if(e0 == 1 && e1 == 1 && e2 == 1)
                    {
                        it->info()[0] = 0; it->info()[1] = 0, it->info()[2] = 0;
                        it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;
                        it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;
                        it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;
                    }

                }

            }
          else if(!isFinite(it->neighbor(0)) && isFinite(it->neighbor(1)) && isFinite(it->neighbor(2)))
            {
                if(edgesum == 1)
                {
                    if(it->info()[1] == 1)
                    {
                        if(e1 == 1 )
                        {
                            it->info()[1] = 0;
                            it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;
                        }
                    }
                    else if(it->info()[2] == 1)
                    {
                        if(e2 == 1 )
                        {
                            it->info()[2] = 0;
                            it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;
                        }
                    }
                }
                else if(edgesum == 2)
                {
                     if(e1 == 1 && e2 == 1 && it->info()[0] == 0)
                     {
                         it->info()[1] = 0; it->info()[2] = 0;
                         it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;
                         it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;
                     }
                }
            }
            else if(!isFinite(it->neighbor(1)) && isFinite(it->neighbor(0)) && isFinite(it->neighbor(2)))
            {
                if(edgesum == 1)
                {
                    if(it->info()[0] == 1)
                    {
                        if(e0 == 1 )
                        {
                            it->info()[0] = 0;
                            it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;
                        }
                    }
                    else if(it->info()[2] == 1)
                    {
                        if(e2 == 1 )
                        {
                            it->info()[2] = 0;
                            it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;
                        }
                    }
                }
                else if(edgesum == 2)
                {
                     if(e0 == 1 && e2 == 1 && it->info()[1] == 0)
                     {
                         it->info()[0] = 0; it->info()[2] = 0;
                         it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;
                         it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;
                     }
                }
            }
            else if(!isFinite(it->neighbor(2)) && isFinite(it->neighbor(0)) && isFinite(it->neighbor(1)))
            {
                if(edgesum == 1)
                {
                    if(it->info()[0] == 1)
                    {
                        if(e0 == 1 )
                        {
                            it->info()[0] = 0;
                            it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;
                        }
                    }
                    else if(it->info()[1] == 1)
                    {
                        if(e1 == 1 )
                        {
                            it->info()[1] = 0;
                            it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;
                        }
                    }
                }
                else if(edgesum == 2)
                {
                     if(e0 == 1 && e1 == 1 && it->info()[2] == 0)
                     {
                         it->info()[0] = 0; it->info()[1] = 0;
                         it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;
                         it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;
                     }
                }
            }
            else if(!isFinite(it->neighbor(0)) && !isFinite(it->neighbor(1)) && isFinite(it->neighbor(2)))
            {
                if(edgesum == 1)
                {
                    if(e2 == 1 && it->info()[0] == 0 && it->info()[1] == 0)
                    {
                        it->info()[2] = 0;
                        it->neighbor(2)->info()[0] = 0; it->neighbor(2)->info()[1] = 0; it->neighbor(2)->info()[2] = 0;

                    }
                }
            }
            else if(!isFinite(it->neighbor(0)) && isFinite(it->neighbor(1)) && !isFinite(it->neighbor(2)))
            {
                if(edgesum == 1)
                {
                    if(e1 == 1 && it->info()[0] == 0 && it->info()[2] == 0)
                    {
                        it->info()[1] = 0;
                        it->neighbor(1)->info()[0] = 0; it->neighbor(1)->info()[1] = 0; it->neighbor(1)->info()[2] = 0;

                    }
                }
            }
            else if(isFinite(it->neighbor(0)) && !isFinite(it->neighbor(1)) && !isFinite(it->neighbor(2)))
            {
                if(edgesum == 1)
                {
                    if(e0 == 1 && it->info()[1] == 0 && it->info()[2] == 0)
                    {
                        it->info()[0] = 0;
                        it->neighbor(0)->info()[0] = 0; it->neighbor(0)->info()[1] = 0; it->neighbor(0)->info()[2] = 0;

                    }
                }
            }
        }
    }
}


bool check_ccentre_hull(Delaunay::Face::Face_handle it)
{
    Point cc = CGAL::circumcenter(Triangle(it->vertex(0)->point(),it->vertex(1)->point(),it->vertex(2)->point()));
    std::vector<Point>::iterator ptr;
    int flg = 0;
    for (ptr = out.begin(); ptr < out.end()-1; ptr++)
    {
           if(CGAL::orientation(*ptr,*(ptr+1),cc) == CGAL::COUNTERCLOCKWISE)
           {
               flg++;
           }
           else
           {
               flg=0;
               break;
           }
      }
    if(flg == 0)
      return 0;
     else
        return 1;

}

bool check_ccentre_in_or_out(Delaunay::Face::Face_handle it)
{

    Point cc = CGAL::circumcenter(Triangle(it->vertex(0)->point(),it->vertex(1)->point(),it->vertex(2)->point()));
    if( ((CGAL::orientation( it->vertex(0)->point(),it->vertex(1)->point(),cc) == CGAL::COUNTERCLOCKWISE )) || ((CGAL::orientation( it->vertex(0)->point(),it->vertex(1)->point(),cc) == CGAL::COLLINEAR)))
        if(((CGAL::orientation( it->vertex(1)->point(),it->vertex(2)->point(),cc) == CGAL::COUNTERCLOCKWISE )) || ((CGAL::orientation( it->vertex(1)->point(),it->vertex(2)->point(),cc) == CGAL::COLLINEAR)))
           if(((CGAL::orientation( it->vertex(2)->point(),it->vertex(0)->point(),cc) == CGAL::COUNTERCLOCKWISE )) || ((CGAL::orientation( it->vertex(2)->point(),it->vertex(0)->point(),cc) == CGAL::COLLINEAR)))
              return 1;
   return 0;


}

void collectboundaryedges()
{

    // collecting the boundary edges
      for (Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); fit++)
      {
          int edgesum = fit->info()[0] + fit->info()[1] + fit->info()[2];
          if(edgesum == 0)
          {
               if(!isFinite(fit->neighbor(0)))
               {
                   Edge e;
                   e.source = fit->vertex(1)->point();
                   e.target = fit->vertex(2)->point();
                   bdry.push_back(e);
               }
               if(!isFinite(fit->neighbor(1)))
               {
                   Edge e;
                   e.source = fit->vertex(2)->point();
                   e.target = fit->vertex(0)->point();
                   bdry.push_back(e);
               }
               if(!isFinite(fit->neighbor(2)))
               {
                   Edge e;
                   e.source = fit->vertex(0)->point();
                   e.target = fit->vertex(1)->point();
                   bdry.push_back(e);
               }
          }
          else if(edgesum == 1)
          {
              if(fit->info()[0] ==1 )
              {
                  Edge e,e1;
                  e.source = fit->vertex(0)->point();
                  e.target = fit->vertex(1)->point();
                  bdry.push_back(e);
                  e1.source = fit->vertex(2)->point();
                  e1.target = fit->vertex(0)->point();
                  bdry.push_back(e1);

              }
              else if(fit->info()[1] ==1 )
              {
                  Edge e,e1;
                  e.source = fit->vertex(1)->point();
                  e.target = fit->vertex(2)->point();
                  bdry.push_back(e);
                  e1.source = fit->vertex(0)->point();
                  e1.target = fit->vertex(1)->point();
                  bdry.push_back(e1);

              }
              else if(fit->info()[2] ==1 )
              {
                  Edge e,e1;
                  e.source = fit->vertex(1)->point();
                  e.target = fit->vertex(2)->point();
                  bdry.push_back(e);
                  e1.source = fit->vertex(2)->point();
                  e1.target = fit->vertex(0)->point();
                  bdry.push_back(e1);

              }
          }
         else if(edgesum ==2)
          {
              if(fit->info()[0] !=1 )
              {
                    Edge e;
                   e.source = fit->vertex(1)->point();
                   e.target = fit->vertex(2)->point();
                    bdry.push_back(e);
              }
              else if(fit->info()[1] !=1 )
              {
                    Edge e;
                   e.source = fit->vertex(2)->point();
                   e.target = fit->vertex(0)->point();
                    bdry.push_back(e);
              }
              else if(fit->info()[2] !=1 )
              {
                    Edge e;
                   e.source = fit->vertex(0)->point();
                   e.target = fit->vertex(1)->point();
                    bdry.push_back(e);
              }

          }

      }

}


//############### Drawing functions##############
/* For displaying points as filled circles */
void drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius)
{
	int triangleAmount = 100;
	GLfloat twicePi = 2.0f * 3.14159265;
	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(x, y);
	for(int i = 0; i <= triangleAmount; i++)
	   glVertex2f(x + (radius * cos(i *  twicePi / triangleAmount)), y + (radius * sin(i * twicePi / triangleAmount)));
	glEnd();
}

void pointset(void)
{
  /***************Open GL commands to smoothen the line segments **************/
  glEnable (GL_LINE_SMOOTH);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  glPointSize( 1.0 );


  float r2=rad;
  for(int i = 0; i < vertices.size(); i++) /*For displaying point set*/
   {
        glColor3f(0.0, 0.0, 1.0);
        drawFilledCircle(vertices[i].x(), vertices[i].y(), r2);
   }

      glLineWidth(3.0);
      glColor3f(1.0, 0.0, 0.0);
      std::vector<int> ls;
      for(int k=0;k<no_of_cycles;k++)
        {
           ls = cycles[k];
           if(ls.size()>10)
            {
              int i=0;
              for(; i < ls.size()-1; i++)
                {
                   glBegin(GL_LINES);
                   glVertex2f(cycle_vertices[ls[i]].x(),cycle_vertices[ls[i]].y());
                   glVertex2f(cycle_vertices[ls[i+1]].x(),cycle_vertices[ls[i+1]].y());
                   glEnd();
                }
               glBegin(GL_LINES);
               glVertex2f(cycle_vertices[ls[i]].x(),cycle_vertices[ls[i]].y());
               glVertex2f(cycle_vertices[ls[0]].x(),cycle_vertices[ls[0]].y());
               glEnd();
            }
         }
}


void reshape(int w, int h)
{
  width = (GLdouble) w;
    height = (GLdouble) h;
    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glRotatef(-180,180,0,1);
    glOrtho(minx-20.0,maxx+20.0,miny-20.0,maxy+20.0, -1.f, 1.f);
    return;
}
void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);
  pointset();
  glFlush();
  return;
}
// CYCLE REMOVAL CODE/////////////////////////////////////
void remove(std::vector<int> &v) //remove the duplicate
{
        auto end = v.end();
        for (auto it = v.begin(); it != end; ++it) {
                end = std::remove(it + 1, end, *it);
        }

        v.erase(end, v.end());
}


int point_already_exist(Point p,int indx)
{
    for(int k=0;k<indx;k++)
    {
        if(cycle_vertices[k]==p)
        {
            return k;
        }
    }
    return -1;
}
// add the edges to the graph
void addEdge(int u, int v)
{
        graph[u].push_back(v);
        graph[v].push_back(u);
}


void dfs_cycle(int u, int p, int color[], int par[])
  {

         // already (completely) visited vertex.
          if (color[u] == 2) {
                  return;
          }

          // seen vertex, but was not completely visited -> cycle detected.
          // backtrack based on parents to find the complete cycle.
        //  std::cout << "color of U "<< color[u]<< std::endl;
          if (color[u] == 1)
          {

                  int cur = p;
                  cycles[no_of_cycles].push_back(cur);

                  while (cur != u)
                  {
                          cur = par[cur];
                          cycles[no_of_cycles].push_back(cur);
                  }
                  no_of_cycles++;
                  return;
           }

          par[u] = p;

          // partially visited.
          color[u] = 1;
          // simple dfs on graph
          for (int v : graph[u]) {

                  // if it has not been visited previously
                  if (v == par[u])
                  {
                          continue;
                  }
                 dfs_cycle(v, u, color, par);
          }

          // completely visited.
          color[u] = 2;

  }
void create_graph()
{
    int loc1=-1,loc2=-1,indx =0;

    for(int i=0;i<bdry.size();i++)
    {
        loc1 = point_already_exist(bdry[i].source,indx);
        loc2 = point_already_exist(bdry[i].target,indx);

        if(loc1 == -1)
        {
            cycle_vertices[indx] = bdry[i].source;
            indx++;

        }
        if(loc2 == -1)
        {
            cycle_vertices[indx] = bdry[i].target;
            indx++;

        }

        if(loc1 == -1 && loc2 == -1)
        {
             addEdge((indx-2), (indx-1));
             num_edges++;

        }
        else if(loc1 == -1 && loc2 != -1)
        {
            addEdge((indx-1), loc2);
            num_edges++;
        }
        else if(loc1 != -1 && loc2 == -1)
        {
            addEdge(loc1,(indx-1));
            num_edges++;
        }
        else if(loc1 != -1 && loc2 != -1)
        {
                 addEdge(loc1,loc2);
                 num_edges++;
        }
    }
    no_ver_graph = indx;

}

void apply_degree_constraint()
{
    float smallest = 10000, small = 10000,temp = 0;
    int smallest_indx = 0, small_indx = 0,temp_indx = 0;

    create_graph(); // Creating the graph with the boundary vertices

    // remove the duplicate edges and remove the degree one vertices
    for(int i=0;i<no_ver_graph;i++)
    {
        std::vector<int>::iterator position;
        remove(graph[i]);
        if(graph[i].size() == 1)
        {
            position = std::find(graph[graph[i][0]].begin(), graph[graph[i][0]].end(), i);
            if (position != graph[graph[i][0]].end())
                graph[graph[i][0]].erase(position);

            graph[i].clear();
        }

    }
//    if(choice != 2)
  //  {
        for(int i=0;i<no_ver_graph;i++)
         {

            if(graph[i].size()>2)
             {
                 smallest = distance(cycle_vertices[i],cycle_vertices[graph[i][0]]);
                 smallest_indx = graph[i][0];
                 small = distance(cycle_vertices[i],cycle_vertices[graph[i][1]]);
                 small_indx =graph[i][1];

                 for(int j=2;j<graph[i].size();j++)
                 {

                     if(small < smallest)
                     {
                        temp = small; temp_indx = small_indx;
                        small = smallest; small_indx = smallest_indx;
                        smallest = temp;  smallest_indx = temp_indx;
                     }
                    temp = distance(cycle_vertices[i],cycle_vertices[graph[i][j]]);
                    temp_indx =  graph[i][j];

                   if(temp < smallest)
                    {
                        small = smallest; small_indx = smallest_indx;
                        smallest = temp;  smallest_indx = temp_indx;
                    }
                    else if(temp < small)
                    {
                        small = temp; small_indx = temp_indx;
                    }

                }
                graph[i].clear();
                graph[i].push_back(smallest_indx);
                graph[i].push_back(small_indx);
             }

          }
  //    }

}

void remove_cycle(void)
{
    int color[N];
    int par[N];

    for(int i =0;i<10000;i++)
    {
        color[i]=0;
        par[i]=0;
    }
    for (int k=0;k<no_ver_graph;k++)
     {
         if (color[k]==2)
            continue;
         dfs_cycle(k, 1, color, par);
     }

}

void ct_shape()
{
    Delaunay::Finite_faces_iterator it;

    Point circ_cent_1,circ_cent_2,edgend_1,edgend_2;
    int  flg_1=-1,flg_2=-1;

    for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
     {

        Delaunay::Face& fce = *it;

        circ_cent_1 = CGAL::circumcenter(Triangle(it->vertex(0)->point(),it->vertex(1)->point(),it->vertex(2)->point()));

        for(int k=0;k<=2;k++)
                {

                    if( isFinite( it->neighbor(k) ) )
                    {

                             circ_cent_2 = CGAL::circumcenter(Triangle(it->neighbor(k)->vertex(0)->point(),it->neighbor(k)->vertex(1)->point(),it->neighbor(k)->vertex(2)->point()));

                             edgend_1= Point(fce.vertex(fce.cw(k))->point().x(),fce.vertex(fce.cw(k))->point().y());
                             edgend_2= Point(fce.vertex(fce.ccw(k))->point().x(),fce.vertex(fce.ccw(k))->point().y());

                             flg_1=-1;
                             flg_2=-1;

                             if((CGAL::orientation(edgend_1,edgend_2,circ_cent_1)) == CGAL::COUNTERCLOCKWISE)
                                flg_1 = 1;
                              else if((CGAL::orientation(edgend_1,edgend_2,circ_cent_1)) == CGAL::CLOCKWISE)
                                flg_1 = 0;

                             if((CGAL::orientation(edgend_1,edgend_2,circ_cent_2)) == CGAL::COUNTERCLOCKWISE)
                                flg_2 = 1;
                             else if((CGAL::orientation(edgend_1,edgend_2,circ_cent_2)) == CGAL::CLOCKWISE)
                                flg_2 = 0;

                             if((flg_1 == 1 && flg_2 == 1)||(flg_1 == 0 && flg_2 == 0))
                             {
                                 it->info()[k] = 1;

                             }

                    }
                    else
                    {
                        if(!check_ccentre_hull(it))
                            it->info()[k] = 1;
                    }
                }


       if(check_ccentre_in_or_out(it))
       {
          check_skinnytraingle(it,circ_cent_1);

       }

    }
}

//################################################
int main(int argc, char **argv)
{
    std::string input = "";
    if (argc < 2)
    {
        std::cout << "Refer the ReadMe file for furthur details." << std::endl;
        exit(0);
    }
    else
    {
        input = std::string(argv[1]);
    }

    std::ifstream inputFile(input.c_str());
    std::istream_iterator<Point> begin(inputFile);
    std::istream_iterator<Point> end;

    dt.insert(begin, end);

    init();
    std::cout <<"Number of faces ="<< dt.number_of_faces() << std::endl;
    Delaunay::Vertex_iterator vi=dt.vertices_begin();
        do{
            if(vi->point().x()>maxx)/*Finding the maximum and minimum x,y coordinates*/
                maxx=vi->point().x();//////////////////////////////////////////////////
            if(vi->point().y()>maxy)///////////////////////////////////////////////////
                maxy=vi->point().y();//////////////////////////////////////////////////
            if(vi->point().x()<minx)///////////////////////////////////////////////////
                minx=vi->point().x();//////////////////////////////////////////////////
            if(vi->point().y()<miny)///////////////////////////////////////////////////
                miny=vi->point().y();//////////////////////////////////////////////////
           vertices.push_back(vi->point());

           if(vertices.size()>500000)
            {
                printf("Input size is too high, please increase the size of input_points and shape array\n");
                exit(1);
            }
            vi++;

        }while(vi!=dt.vertices_end());


  diagonalDistance = distance(Point(minX, minY, 1), Point(maxX, maxY, 1));
  std::cout <<"Number of point ="<< vertices.size() << std::endl;

  computconvexhull(); //computing convex hull


  Delaunay::Edge_iterator ei=dt.edges_begin();
  int i;
    do{
        Delaunay::Face& f = *(ei->first);
        i = ei->second;
        insertToShape(f.vertex(f.cw(i))->point(),f.vertex(f.ccw(i))->point(),shape);
        ei++;
    }while(ei!=dt.edges_end());





    ct_shape();
    post_process();
    populatemarkedEdges();
    collectboundaryedges();
    apply_degree_constraint();
    remove_cycle();

     if(argc > 2)
        rad = std::atof(argv[2]);
     else
         rad =1;

     glutInit(&argc, argv);
     glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
     glutInitWindowSize(650, 650);
     float w = glutGet(GLUT_WINDOW_WIDTH);
     float h = glutGet(GLUT_WINDOW_HEIGHT);
     window = glutCreateWindow("CT Shape");
     glutReshapeFunc(reshape);
     glClearColor(1.0, 1.0, 1.0, 0.0);
     glutDisplayFunc(display);
     glutMainLoop();
     return 0;
}
